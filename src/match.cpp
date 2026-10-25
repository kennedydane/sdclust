/**************************************
<match.cpp: A class for storing and doing sd_heuristic matching>
    Copyright (C) <2008>  <Dane Kennedy>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

**************************************/


#include "match.h"

//using namespace std;

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cassert>
#include <sys/types.h>
//#include "clusters.h"
//#include "fasta.h"

#define MAX(i,j) (i>j?i:j)
#define MIN(i,j) (i<j?i:j)
#define SQR(i) (i)*(i)

//globals global;


/************************ MATCH ************************/


bool match::thresholds_set = false;
uint8_t match::matches_threshold = 0;
float match::var_threshold = 0;

match::match() {
  this->positions=NULL;
  this->matches_count=0;
  this->match_found = false;
}

void match::initialise(uint8_t m_thresh, float sd_thresh, position *positions)
{
  if (!this->thresholds_set) {
    this->thresholds_set = true;
    this->matches_threshold = m_thresh;
    this->var_threshold = sd_thresh*sd_thresh;
  }
  this->positions = positions;
  this->matches_count=0;
  this->match_found = false;
}

void match::reset() {
  this->matches_count=0;
  this->match_found = false;
}

bool match::add_match(position where) {
  matches_count++;

  if (this->matches_count == 1) {
    //this->sd = 0;
    this->var=0;
    this->mean = where.source - where.target;
    //this->min_pos = where;
    //this->max_pos = where;
    this->positions[0] = where;

    return true;
  }
  else if (this->matches_count <= this->matches_threshold) {
    this->positions[this->matches_count-1] = where;
    float offset = where.source - where.target;
    //float var = sd*sd;
    float total = this->mean * (this->matches_count-1);
    float sum_sqr = this->var*(this->matches_count - 2) + total * this->mean;
    this->mean = (total + offset) / this->matches_count;
    sum_sqr += offset*offset;
    total += offset;
    this->var = (sum_sqr - total*this->mean)/(this->matches_count-1);
    //this->sd = sqrt(var);
    
    //if (where.source < min_pos.source)      min_pos.source = where.source;
    //else if (where.source > max_pos.source) max_pos.source = where.source;
    //if (where.target < min_pos.target)      min_pos.target = where.target;
    //else if (where.target > max_pos.target) max_pos.target = where.target;

    if ((this->matches_count == matches_threshold) && /*(sd <= sd_threshold)*/ (this->var <= var_threshold)) {
      this->match_found = true;
    }
    return false;
  }
  else {
    position old = this->positions[(this->matches_count-1) % this->matches_threshold];
    this->positions[(this->matches_count-1) % this->matches_threshold] = where;

    //if (where.source > max_pos.source)      max_pos.source = where.source;
    //else if (where.source < min_pos.source) min_pos.source = where.source; // Shouldn't really happen...
    //if (where.target > max_pos.target)      max_pos.target = where.target;
    //else if (where.target < min_pos.target) min_pos.target = where.target;

    if (!this->match_found) {
      float old_offset = old.source - old.target;
      float offset = where.source - where.target;
      //float var = sd*sd;
      float total = this->mean * (this->matches_threshold);
      float sum_sqr = this->var*(this->matches_threshold-1) + total * this->mean;
      this->mean = (total + offset - old_offset) / this->matches_threshold;  
      sum_sqr += offset*offset - old_offset*old_offset;
      total += offset - old_offset;
      if (this->matches_threshold > 1)
        this->var = (sum_sqr - total * this->mean) / (this->matches_threshold - 1);
      else
        this->var=0;
      if (this->var < 0)
        this->var = -this->var;
      //this->sd = sqrt(var);

      if (/*sd <= sd_threshold*/this->var <= var_threshold) {
        this->match_found = true;
      }
    }
    return false;
  }
}

/*void match::show_match(uint32_t source, uint32_t target) {
  if ((matches_count >= matches_threshold) && (sd <= sd_threshold)) {
    cout << source << ", " << target << ": count = " << (int)matches_count << " mean = " << mean << " sd = " << sd << "\n";
  } 
};*/

void match::cluster(clusters *the_clusters, uint32_t source, sequence* s_seq, uint32_t target, sequence* t_seq, bool rc, d1 *d1_calc, d2 *d2_calc, match_stats *stats, ostream &out) {
  if (this->match_found && (the_clusters->find_parent(source) != the_clusters->find_parent(target))) {

    stats->heuristic_success++;
    uint16_t s_min, s_max, t_min, t_max;

    //sequence *s_seq = the_sequences->get_sequence(source);
    //sequence *t_seq = the_sequences->get_sequence(target);
    
    /*if (rc) {
      s_seq=s_seq->rc(); 
    }*/

    /*s_min = MAX(0, MIN(min_pos.source-globals::skip, max_pos.source - globals::d2_window_length));
    t_min = MAX(0, MIN(min_pos.target-globals::skip, max_pos.target - globals::d2_window_length));
    s_max = MIN(MAX(min_pos.source + globals::d2_window_length, max_pos.source+globals::mer_size-1), s_seq->get_len()-globals::d2_word_length);
    t_max = MIN(MAX(min_pos.target + globals::d2_window_length, max_pos.target+globals::mer_size-1), t_seq->get_len()-globals::d2_word_length);*/

    //uint16_t s_min, s_max, t_min, t_max;
    uint32_t d1_score = d1_calc->calc_score(t_seq, &s_min,&s_max, &t_min, &t_max);

    //uint32_t c_words = d1_calc->calc_score(t_seq, &t_min, &t_max);

    int32_t d2_score;
    
    if (d1_score>=globals::d1_threshold) {
      stats->d1_success++;
      if (!globals::heuristic_only) {
        //d2_score = d2_calc->calc_score(s_seq, t_seq, 0, s_seq->get_len()-globals::d2_word_length, MAX(0,(int16_t)t_min-20), MIN(t_seq->get_len()-globals::d2_word_length,t_max+20));
        //d2_score = d2_calc->calc_score(s_seq, t_seq, 0, s_seq->get_len()-globals::d2_word_length, dummy1, dummy2);
        d2_score = d2_calc->calc_score(s_seq, t_seq, s_min, s_max, t_min, t_max);

      }
      else {
        d2_score = 0;
      }

    }
    else d2_score=globals::d2_threshold*10;

    /*if (rc) {
      delete s_seq;
    }*/
    if (d2_score <= globals::d2_threshold) {
      if (!globals::heuristic_only)
        stats->d2_success++;
      if (!globals::complete) //Only put them in same cluster if doing single linkage, otherwise we're only interested in the extended table anyway
        the_clusters->make_union(source, target);
      if (globals::show_extended) {
        out << source << "\t" << target << "\t" << rc << "\t" << d2_score << "\n";
      }
    }
  };
}

/************************ MATCH LIST ************************/

match_list::match_list(uint32_t list_size, uint8_t m_thresh, float sd_thresh, clusters *the_clusters) {
  this->list_size = list_size;
  this->the_list = new match[list_size];
  this->position_pool = new position[(size_t)list_size * m_thresh];
  for (uint32_t i=0; i < list_size; i++)
    this->the_list[i].initialise(m_thresh, sd_thresh, &this->position_pool[(size_t)i * m_thresh]);
  this->matches=new vector<uint32_t>;
  this->the_clusters = the_clusters;
  this->d1_calc = new d1();
  this->d2_calc = new d2();
  this->stats.heuristic_success = 0;
  this->stats.d1_success = 0;
  this->stats.d2_success = 0;
}

match_list::~match_list() {
  delete[] this->the_list;
  delete[] this->position_pool;
  delete this->matches;
  delete this->d1_calc;
  delete this->d2_calc;
}

vector<uint32_t>* match_list::get_matches() {
  return this->matches;
}

void match_list::reset_list() {
  for (vector<uint32_t>::iterator i = this->matches->begin(); i != this->matches->end(); ++i) {
    this->the_list[*i].reset();
  };
  this->matches->clear();
  
}

void match_list::add_match(uint32_t who_am_i, uint32_t which_seq, position where) {
  if ((the_clusters->find_parent(who_am_i) != the_clusters->find_parent(which_seq)) && (this->the_list[which_seq].add_match(where))) {
    this->matches->push_back(which_seq);
  }
  
}

void match_list::find_clusters(uint32_t who_am_i, sequence* the_sequence, sequence_list *the_sequences, bool rc, ostream &out) {
  this->d1_calc->set_num_words1(the_sequence);
  for (vector<uint32_t>::iterator i = this->matches->begin(); i != this->matches->end(); ++i) {
    this->the_list[*i].cluster(this->the_clusters, who_am_i, the_sequence, *i, the_sequences->get_sequence(*i), rc, this->d1_calc, this->d2_calc, &this->stats, out);
  }
}
