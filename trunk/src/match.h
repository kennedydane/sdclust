/**************************************
match.h: Wadda wadda

**************************************/
#ifndef MATCH_H
#define MATCH_H

//#include "common.h"

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

class match{
  public:
    match();
    match(uint8_t m_thresh, float sd_thresh);
    ~match();
    void reset();
    bool add_match(position where);
    //void show_match(uint32_t source, uint32_t target);
    void cluster(clusters* the_clusters, uint32_t source, sequence* src_seq, uint32_t target, sequence* t_seq, bool rc, d1 *d1_calc, d2 *d2_calc);
    //uint32_t common_words(sequence *s_seq, sequence *t_seq);
    
    //void set_num_words1(sequence *s_seq);
    //int32_t d2_score(sequence *s_seq, sequence *t_seq, uint16_t s_min, uint16_t s_max, uint16_t t_min, uint16_t t_max);
  private:

    static bool thresholds_set;
    static uint8_t matches_threshold;
    //static float sd_threshold;
    static float var_threshold;

    uint16_t matches_count;
    bool match_found;
    position* positions;
    float /*sd,*/ mean, var;
    //float min_sd;
    //position min_pos, max_pos;
    static int16_t delta1[65536];
    static int16_t delta2[65536];
    static int16_t num_words1[65536];
    static int16_t num_words2[65536];
//    float min_sd;
};




class match_list{
  public:
    match_list(uint32_t list_size, uint8_t m_thresh, float sd_thresh);
    ~match_list();
    vector<uint32_t>* get_matches();
    void reset_list();
    void add_match(uint32_t who_am_i, uint32_t which_seq, position where);
    //void show_matches(uint32_t who_am_i);
    void find_clusters(uint32_t who_am_i, sequence* the_sequence, sequence_list *the_sequences, bool rc, d1* d1_calc, d2* d2_calc);
    void show_clusters();
    

  private:
    match** the_list;
    uint32_t list_size;
    vector<uint32_t> *matches;
    clusters* the_clusters;
};

#endif
