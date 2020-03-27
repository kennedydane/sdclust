/**************************************
<d1and2.cpp: classes implementing the d1 and d2 distance measures>

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

#include "d1and2.h"
#include "globals.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <sys/types.h>

#define MAX(i,j) (i>j?i:j)
#define MIN(i,j) (i<j?i:j)

using namespace std;

/****************** d1 ******************/

/*void d1::set_num_words1(sequence *src_seq) {
  this->src_seq=src_seq;
  memset(this->num_words1, 0, sizeof(int16_t)*65536);
  //uint8_t old_mer_size = src_seq->get_mer_size();
  //src_seq->set_mer_size(globals::d2_word_length);
  for (uint16_t i=0; i <= src_seq->get_len()-globals::d1_word_length; i++) {
    seq_word the_word = src_seq->get_word(i, globals::d1_word_length);
    this->num_words1[the_word]++;

  };
  //src_seq->set_mer_size(old_mer_size);
  //cout << "--" << wc << "--";
};*/

void d1::set_num_words1(sequence *src_seq) {
  this->src_seq=src_seq;
  //this->src_seq_len=src_seq->get_len();
  memset(this->num_words1,0,sizeof(d1_match)*this->array_size);
  for (uint16_t i=0; i <= src_seq->get_len()-globals::d1_word_length; i++) {
    seq_word the_word = src_seq->get_word(i, globals::d1_word_length);
    if (this->num_words1[the_word].count++==0) {
      this->num_words1[the_word].min=i;
      this->num_words1[the_word].max=i;
    }
    else {
      this->num_words1[the_word].max=i;
    }
    //cerr << num_words1[the_word].count << " " << num_words1[the_word].min << " " << num_words1[the_word].max << "\n";

  };
}

uint32_t d1::calc_score(sequence *target_seq,uint16_t *s_min_pos, uint16_t *s_max_pos, uint16_t *t_min_pos, uint16_t *t_max_pos) {
  if ((this->src_seq->get_len() < globals::d1_window_length) || (target_seq->get_len() < globals::d1_window_length)) {
    return 0;
  }
  int32_t common=0;
  int32_t common_max=0;
  uint16_t s_min, s_max, t_min, t_max;
  t_min=target_seq->get_len()-globals::d1_word_length;
  t_max=0;
  
  s_min=this->src_seq->get_len()-globals::d1_word_length;
  s_max=0;

  memset(num_words2, 0, sizeof(uint16_t)*this->array_size);
  for (uint16_t i=0; i <= target_seq->get_len()-globals::d1_word_length; i++) {
    seq_word the_word = target_seq->get_word(i, globals::d1_word_length);
    
    if (num_words1[the_word].count) {
      s_min=MIN(s_min, num_words1[the_word].min);
      s_max=MAX(s_max, num_words1[the_word].max);

      //  *min_pos = i;
      //*min_pos = MIN(*min_pos, i);
      if (common+1>=globals::d1_threshold) {
        t_max = i;
        t_min=MIN(t_min,MAX(0,i-globals::d1_window_length));
      }
      //*max_pos = MAX(MIN(*min_pos+globals::d1_window_length-globals::d1_word_length,target_seq->get_len()-globals::d1_word_length),i);
    }

    common -= MIN(num_words1[the_word].count,num_words2[the_word]);
    num_words2[the_word]++;
    common += MIN(num_words1[the_word].count,num_words2[the_word]);

    if (i>globals::d2_window_length) {
      seq_word old_word = target_seq->get_word(i-globals::d1_window_length, globals::d1_word_length);
      
      common-=MIN(num_words1[old_word].count,num_words2[old_word]);
      num_words2[old_word]--;
      common+=MIN(num_words1[old_word].count,num_words2[old_word]);
    }
    common_max = MAX(common, common_max);
    //if (common_max >= globals::d1_threshold)
    //  return common_max;

  }
  *s_max_pos=MAX(MIN(src_seq->get_len()-globals::d1_word_length,s_min+globals::d1_window_length-globals::d1_word_length),MIN(s_max+globals::d1_word_length,src_seq->get_len()-globals::d1_word_length));
  *s_min_pos=MIN(MAX(0,s_max-globals::d1_window_length+globals::d1_word_length),MAX(s_min-globals::d1_word_length, 0));
  
  *t_max_pos=MAX(MIN(target_seq->get_len()-globals::d1_word_length,t_min+globals::d1_window_length-globals::d1_word_length),MIN(t_max+globals::d1_word_length,target_seq->get_len()-globals::d1_word_length));
  *t_min_pos=MIN(MAX(0,t_max-globals::d1_window_length+globals::d1_word_length),MAX(t_min-globals::d1_word_length, 0));
  return common_max;
  //return (globals::d1_window_length-globals::d1_word_length+1)+(this->src_seq_len-globals::d1_word_length+1)-(2*common_max);
}


/****************** d2 ******************/

inline void d2::add_word(int16_t delta1[], int16_t delta2[], int32_t *score, seq_word the_word) {
/*  *score = *score - SQR(delta1[the_word]-delta2[the_word]);
  delta1[the_word]++;
  *score = *score + SQR(delta1[the_word]-delta2[the_word]);*/
  *score +=  2 * ((delta1[the_word]++)-delta2[the_word]) + 1;
}

inline void d2::del_word(int16_t *delta1, int16_t *delta2, int32_t *score, seq_word the_word) {
/*  *score = *score - SQR(delta1[the_word]-delta2[the_word]);
  delta1[the_word]--;
  *score = *score + SQR(delta1[the_word]-delta2[the_word]);*/
  *score -= 2 * ((delta1[the_word]--)-delta2[the_word]) - 1;
}

int32_t d2::calc_score(sequence* s_seq, sequence*  t_seq, uint16_t s_min, uint16_t s_max, uint16_t t_min, uint16_t t_max) {
  int num_words_win = globals::d2_window_length - globals::d2_word_length+1;

  if ((s_max - s_min +1 < num_words_win) || (t_max - t_min + 1 < num_words_win)) {
    cerr << ":-O " << s_max << " " << s_min << " " << t_max << " " << t_min << " :-O\n" <<flush;
    return globals::d2_threshold*10;
  }
  memset(delta1,0,sizeof(int16_t)*this->array_size); //Zero the two tables.
  memset(delta2,0,sizeof(int16_t)*this->array_size);
  int32_t score=0;
  int32_t min_score;
  
  for (uint16_t i=0; i < num_words_win; ++i) {
    seq_word s_word = s_seq->get_word(i+s_min, globals::d2_word_length);
    seq_word t_word = t_seq->get_word(i+t_min, globals::d2_word_length);

    add_word(delta1, delta2, &score,  s_word);
    add_word(delta2, delta1, &score,  t_word);
  }

  if (score <= globals::d2_threshold) {
    return score; //UNCOMMENT FOR SPEED!
  }
  min_score = score;
  
  this->src_pos = s_min;
  this->tgt_pos = t_min;
  
  while (this->src_pos <= s_max-num_words_win+1)  {
  
  //for (int16_t i=s_min/*+1*/; i <= s_max-num_words_win+1; ++i) { 

    if (this->src_pos != s_min) {
      del_word(delta1, delta2, &score, s_seq->get_word(this->src_pos-1, globals::d2_word_length));
      add_word(delta1, delta2, &score, s_seq->get_word(this->src_pos+num_words_win-1, globals::d2_word_length));
    }
    min_score = MIN(min_score, score);

    if (min_score <= globals::d2_threshold) return min_score; //UNCOMMENT FOR SPEED!

    while (++this->tgt_pos <= t_max-num_words_win+1) {
    //for (int16_t j=t_min+1; j <= t_max-num_words_win+1; ++j) {
      del_word(delta2, delta1, &score, t_seq->get_word(this->tgt_pos-1, globals::d2_word_length));
      add_word(delta2, delta1, &score, t_seq->get_word(this->tgt_pos+num_words_win-1, globals::d2_word_length));
      min_score = MIN(min_score, score);
      if (min_score <= globals::d2_threshold) return min_score; //UNCOMMENT FOR SPEED!
    }

    ++this->src_pos;

    if (this->src_pos <= s_max-num_words_win+1) {
      del_word(delta1, delta2, &score, s_seq->get_word(this->src_pos-1, globals::d2_word_length));
      add_word(delta1, delta2, &score, s_seq->get_word(this->src_pos+num_words_win-1, globals::d2_word_length));

      min_score = MIN(min_score, score);
      if (min_score <= globals::d2_threshold) return min_score; //UNCOMMENT FOR SPEED!

      //for (int16_t j=t_max-(num_words_win)+1; j>t_min; --j) {
      while (--this->tgt_pos > t_min) {
        del_word(delta2, delta1, &score, t_seq->get_word(this->tgt_pos+num_words_win-1, globals::d2_word_length));
        add_word(delta2, delta1, &score, t_seq->get_word(this->tgt_pos-1, globals::d2_word_length));

        min_score = MIN(min_score, score);
        if (min_score <= globals::d2_threshold) return min_score; //UNCOMMENT FOR SPEED!
      }
    }
    ++this->src_pos;
  }
  return min_score;
}
