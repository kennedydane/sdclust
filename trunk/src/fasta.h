/**************************************
<fasta.h: class for accessing fasta files. It also has stuff for storing them and
then doing the clustering. This should be separated in the future.>

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

#ifndef FASTA_H
#define FASTA_H

/**
	@author Dane Kennedy,,, <kennedy.dane@gmail.com>
 */

#include <string>
#include <iostream>
#include <fstream>
#include <sys/types.h>

//#include "common.h"
#include "match.h"
#include "globals.h"
#include "d1and2.h"

using namespace std;

// forward declarations because of cyclic dependencies

class merray; // forward declaration because of heinous class interdependencies.
class sequence;
class sequence_list;
class match_list;

//typedef uint64_t seq_word;


class fasta{
  public:
    fasta(char *source_file_name);
    uint32_t count_sequences();
    sequence *get_sequence(uint8_t mer_size);
    void reset();

    ~fasta();
    
  private:
    ifstream source_file;
    char *buffer;
    uint32_t size;
    uint32_t where;
};

struct boundary{
  uint16_t begin, end;
};

const boundary zero_boundary = {0,0};

class boundary_list{
  public:
    boundary_list(boundary the_boundary) {this->next = NULL; this->the_boundary = the_boundary;};
    ~boundary_list();
    void set_next(boundary_list *who) {this->next = who;};
    boundary_list* get_next() {return this->next;};
    uint16_t count();
    boundary get_boundary() {return this->the_boundary;};
  private:
    boundary_list *next;
    boundary the_boundary;
};

class boundary_array{
  public:
    boundary_array(boundary_list *the_list);
    ~boundary_array();
    uint16_t len() {return this->how_many;};
    boundary get_boundary(uint16_t which);
  private:
    boundary *the_array;
    uint16_t how_many;
};


class sequence{
  public:
    sequence() {};
    sequence(string sequence_string, uint8_t mer_size);
    ~sequence() {delete[] this->the_seq; delete this->boundaries;};
    sequence* rc();
    
    //static void set_mer_size(uint8_t mer_size) {sequence::mer_size = mer_size;};
    //static uint8_t get_mer_size() {return sequence::mer_size;};

    uint16_t get_len(){return this->len;};
    void print_boundaries();
    uint16_t count_mers(uint8_t mer_size);
    void update_mer_count_array(uint32_t *the_array, uint8_t array_mer_size, uint8_t mer_size);
    void add_mers(uint32_t my_index, merray* the_array, uint32_t *mer_count_array, uint8_t array_mer_size, uint8_t mer_size);
    seq_word get_word(uint16_t which_word, uint8_t mer_size);
    string get_word_s(uint16_t which_word, uint8_t mer_size);
    string get_sequence();
    bool find_matches(uint32_t my_index, match_list *matches, sequence_list *seq_list, merray *the_array, uint8_t mer_size); //, bool update_search = false);
  private:
    seq_word *the_seq; //The sequence stored in compact form
    uint16_t len;
    boundary_array *boundaries;
    //static uint8_t mer_size;
};

class sequence_list{
  public:
    sequence_list(fasta* the_source, uint8_t mer_size, uint8_t matches_threshold, float sd_threshold);
    void add_sequence(sequence *the_seq);
    sequence* get_sequence(uint32_t which) {return this->sequences[which];};
    uint32_t count_mers(uint8_t mer_size);
    
    void merge_sort_mers(uint8_t array_mer_size, uint8_t mer_size);
    void print_mers(uint8_t mer_size);
    void create_unique_mer_array(uint8_t mer_size); // Creates a list of indexes into a sorted mer_array of all the unique words
    void print_unique_mers(uint8_t mer_size);
    uint32_t find_mer(seq_word the_mer, uint8_t mer_size); //, bool update_search = false); //returns index into mer_array
    seq_word get_mer(uint32_t merr_index, uint8_t mer_size);
    //uint8_t get_mer_size() {return this->mer_size;};
    uint8_t get_array_mer_size() {return this->array_mer_size;};
    string get_mer_s(uint32_t merr_index, uint8_t mer_size);
    uint32_t* get_mer_count_array() {return this->mer_count_array;};
    uint32_t len() {return this->how_many;};
    void find_matches(uint32_t start, uint32_t end, uint8_t bucket_size, uint8_t mer_size);
    ~sequence_list();
    //int32_t d2_score(sequence *s_seq, sequence *t_seq, uint16_t s_min, uint16_t s_max, uint16_t t_min, uint16_t t_max);
  private:
    void create_merray(uint8_t array_mer_size, uint8_t mer_size);
    void create_mer_count_array(uint8_t mer_size);
    void update_unique_mer_array(uint32_t min_index);
    uint8_t array_mer_size;
    uint32_t unique_mers;
    uint32_t* mer_count_array;
    uint32_t *unique_mer_array;
    uint32_t how_many;
    merray *mer_array;
    //uint8_t mer_size;
    sequence **sequences;
    match_list *matches;
    d1 *d1_calc;
    d2 *d2_calc;

};

class merray{
  public:
    merray(uint32_t how_many);
    ~merray();
    void add_mer(uint32_t mer_index, seq_word the_mer, uint32_t *mer_count_array, uint8_t mer_size, uint8_t array_mer_size);
    void show_mers();
    void show_mers(sequence_list *the_sequences, uint8_t mer_size);
    //uint64_t get_mer(sequence_list *the_sequences, uint32_t mer_array_index, uint8_t which_mer);
    uint32_t get_mer_index(uint32_t mer_array_index) {return this->the_array[mer_array_index];};
    uint32_t len() {return this->array_len;};
    void merge_sort_mers(sequence_list* the_sequences, uint8_t mer_size);
    uint32_t count_unique_mers(sequence_list* the_sequences, uint8_t mer_size);
    uint32_t* create_unique_mer_list(sequence_list* the_sequences, uint32_t unique_mer_count, uint8_t mer_size);
    void print_mer(sequence_list* the_sequences, uint32_t mer_array_index, uint8_t mer_size);
  private:
    void insertion_sort_mers(sequence_list *the_sequences, uint32_t start, uint32_t finish, uint8_t mer_size);
    void insert(sequence_list *the_sequences, uint32_t start, uint32_t finish, uint8_t mer_size);
    void merge_sort_mers(sequence_list* the_sequences, uint32_t start, uint32_t finish, uint8_t mer_size);
    void merge(sequence_list* the_sequences, uint32_t left, uint32_t right, uint8_t mer_size);
    
    uint32_t *the_array;
    uint32_t array_len;
    //uint32_t *mer_counts;

};



#endif
