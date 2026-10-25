/**************************************
<match.h: A class for storing and doing sd_heuristic matching>
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


#ifndef MATCH_H
#define MATCH_H

#include <vector>
#include <sys/types.h>
#include <fstream>
#include "clusters.h"
#include "fasta.h"
#include "globals.h"
#include "d1and2.h"

class sequence;

using namespace std;

/**
	@author Dane Kennedy,,, <kennedy.dane@gmail.com>
*/

class sequence_list;
class d1;
class d2;

class position {
  public:
    uint16_t source, target;
};

// Counters that each thread keeps for itself and adds to globals at the end.
struct match_stats {
  uint64_t heuristic_success;
  uint64_t d1_success;
  uint64_t d2_success;
};

class match{
  public:
    match();
    void initialise(uint8_t m_thresh, float sd_thresh, position *positions);
    void reset();
    bool add_match(position where);
    void cluster(clusters* the_clusters, uint32_t source, sequence* src_seq, uint32_t target, sequence* t_seq, bool rc, d1 *d1_calc, d2 *d2_calc, match_stats *stats, ostream &out);
  private:

    static bool thresholds_set;
    static uint8_t matches_threshold;
    static float var_threshold;

    uint16_t matches_count;
    bool match_found;
    position* positions; // matches_threshold slots in the owning match_list's pool
    float mean, var;
};



// The matching state for one thread: a match per possible target sequence,
// plus its own d1 and d2 calculators. The cluster table is shared.
class match_list{
  public:
    match_list(uint32_t list_size, uint8_t m_thresh, float sd_thresh, clusters *the_clusters);
    ~match_list();
    vector<uint32_t>* get_matches();
    void reset_list();
    void add_match(uint32_t who_am_i, uint32_t which_seq, position where);
    void find_clusters(uint32_t who_am_i, sequence* the_sequence, sequence_list *the_sequences, bool rc, ostream &out);
    match_stats* get_stats() {return &this->stats;};

  private:
    match* the_list;
    position* position_pool;
    uint32_t list_size;
    vector<uint32_t> *matches;
    clusters* the_clusters;
    d1 *d1_calc;
    d2 *d2_calc;
    match_stats stats;
};

#endif
