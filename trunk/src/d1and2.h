/**************************************
<d1and2.h: classes implementing the d1 and d2 distance measures>

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

#ifndef D1AND2_H
#define D1AND2_H

#include <sys/types.h>
#include "globals.h"
#include "fasta.h"

/**
	@author Dane Kennedy,,, <kennedy.dane@gmail.com>
 */

using namespace std;

class sequence;

struct d1_match {
  uint16_t count;
  int16_t min;
  int16_t max;
};

class d1 {
  public:
    d1(){array_size=1<<(globals::d1_word_length<<1); num_words1=new d1_match[array_size]; num_words2 = new uint16_t[array_size];};
    ~d1(){delete[] num_words1; delete[] num_words2;};
    void set_num_words1(sequence* src_seq);
    uint32_t calc_score(sequence* target_seq, uint16_t *s_min_pos, uint16_t *s_max_pos, uint16_t *min_pos, uint16_t *max_pos);
  private:
    uint32_t array_size;
    sequence* src_seq;
    //uint16_t src_seq_len;
    //int16_t num_words1[65536];
    d1_match *num_words1;
    uint16_t *num_words2;
};

class d2 {
  public:
    d2(){array_size=1<<(globals::d2_word_length<<1);delta1=new int16_t[array_size];delta2=new int16_t[array_size];};
    ~d2(){delete[] delta1; delete[] delta2;};
    int32_t calc_score(sequence *src_seq, sequence* target_seq, uint16_t s_min, uint16_t s_max, uint16_t t_min, uint16_t t_max);
  private:
    void add_word(int16_t delta1[], int16_t delta2[], int32_t *score, seq_word the_word);
    void del_word(int16_t *delta1, int16_t *delta2, int32_t *score, seq_word the_word);
    int16_t *delta1;//[65536];
    int16_t *delta2;//[65536];
    uint32_t array_size;
    uint16_t src_pos, tgt_pos;
};


#endif
