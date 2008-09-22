/**************************************
<fasta.cpp: class for accessing fasta files. It also has stuff for storing them and
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

#include "fasta.h"
#include "match.h"

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cassert>
#include <sys/types.h>
#include "globals.h"

//globals global_vars;


using namespace std;

#define MAX(i,j) (i>j?i:j)
#define MIN(i,j) (i<j?i:j)
#define SQR(i) (i)*(i)


/************************ BOUNDARY LIST ************************/

boundary_list::~boundary_list() {if (this->next !=0) delete this->next;}

uint16_t boundary_list::count() {
  if (this->next == NULL)
    return 1;
  else
    return 1 + this->next->count();
}

/************************ BOUNDARY ARRAY ************************/

boundary_array::boundary_array(boundary_list *the_list) {
  this->how_many = the_list->count();
  this->the_array = new boundary[this->how_many];
  boundary_list *temp_list = the_list;
  if (how_many > 0) {
    for (int32_t i = this->how_many - 1; i >= 0; i--) {
      the_array[i] = temp_list->get_boundary();
      temp_list = temp_list->get_next();
    }
  }
}

boundary_array::~boundary_array() {
  delete[] the_array;
}

boundary boundary_array::get_boundary(uint16_t which) {
  if (which < this->how_many)
    return this->the_array[which];
  return zero_boundary;
}

/************************ SEQUENCE ************************/

//uint8_t sequence::mer_size = 12;

sequence::sequence(string sequence_string, uint8_t mer_size){
  //boundaries are created based on the mer size. There is no point in storing boundaries
  //that are smaller than mer size...
  this->len = sequence_string.length();
  //this->mer_size = mer_size;
  uint16_t how_many;
  if (this->len % 32 == 0)
    how_many = (this->len / 32);
  else
    how_many = (this->len / 32) + 1;
  this->the_seq = new seq_word[how_many];
  for (uint16_t i = 0; i < how_many; i++) {
    this->the_seq[i] = 0;
  }
  bool valid_sequence = false;
  uint16_t start=0;
  boundary_list *b_list = NULL;
  for (uint16_t i = 0; i < this->len; i++) {

    char val = 0;
    switch (sequence_string[i]) {
      case 'a' :
      case 'A' : val = 0; ++globals::A_count;
      if (!valid_sequence) {
        valid_sequence = true;
        start = i;
      }
      break;
      
      case 'c' :
      case 'C' : val = 1; ++globals::C_count;
      if (!valid_sequence) {
        valid_sequence = true;
        start = i;
      }
      break;
      
      case 'g' :
      case 'G' : val = 2; ++globals::G_count;
      if (!valid_sequence) {
        valid_sequence = true;
        start = i;
      }
      break;
      
      case 't' :
      case 'T' : val = 3; ++globals::T_count;
      if (!valid_sequence) {
        valid_sequence = true;
        start = i;
      }
      break;
      
      default : val = (random()>>5)&3; ++globals::N_count;
      if (valid_sequence) {
        valid_sequence = false;
        if ((i-start) >= mer_size) {
          boundary the_boundary = {start, i};
          if (b_list == NULL) {
            b_list = new boundary_list(the_boundary);
          }
          else {
            boundary_list *temp_list = new boundary_list(the_boundary);
            temp_list->set_next(b_list);
            b_list = temp_list;
          }
        }
      }
    }
    this->the_seq[i / 32] = (this->the_seq[i / 32]<<2) | val;
  };
  for (uint16_t i = this->len; i % 32 != 0; i++) { // makes sure the last word is in the right place.
    this->the_seq[i / 32] = this->the_seq[i/32]<<2;
  }
  if (valid_sequence) {
    if ((this->len - start >= mer_size)) {
      boundary the_boundary = {start, len};
      if (b_list == 0) {
        b_list = new boundary_list(the_boundary);
      }
      else {
        boundary_list *temp_list = new boundary_list(the_boundary);
        temp_list->set_next(b_list);
        b_list = temp_list;
      }
    }
  }
  if (b_list == 0) {
    b_list = new boundary_list(zero_boundary);
  }
  this->boundaries = new boundary_array(b_list);
  delete b_list;
}

sequence* sequence::rc() {
  sequence *temp = new sequence();
  uint16_t how_many;
  uint8_t offset=0;
  if (this->len % 32 == 0)
    how_many = (this->len / 32);
  else {
    how_many = (this->len / 32) + 1;
    offset = 32 - (this->len % 32);
  }
  temp->len = len;
  temp->the_seq = new seq_word[how_many];

  boundary the_boundary, t_boundary=this->boundaries->get_boundary(this->boundaries->len()-1);
  the_boundary.begin = this->len - t_boundary.end;
  the_boundary.end = this->len - t_boundary.begin;
  boundary_list *b_list = new boundary_list(the_boundary);

  for (int32_t i=this->boundaries->len()-2; i >= 0; i--) {
    t_boundary = this->boundaries->get_boundary(i);
    the_boundary.begin = this->len - t_boundary.end;
    the_boundary.end = this->len - t_boundary.begin;
    boundary_list *temp_list = new boundary_list(the_boundary);
    temp_list->set_next(b_list);
    b_list = temp_list;
  };
  temp->boundaries = new boundary_array(b_list);
  delete b_list;

  for (uint16_t j = 0; j < how_many; j++) {
    seq_word word = ~this->the_seq[j];
    uint16_t which = how_many-1-j;
    temp->the_seq[which] = 0;
    for (uint8_t k = 0; k< 32; k++) {
      temp->the_seq[which] = (temp->the_seq[which] << 2) | (word & 3);
      word = word >> 2;
    }
  }
  if (offset != 0) {
    for (uint16_t j=0; j < how_many-1; j++) {
      temp->the_seq[j] = (temp->the_seq[j] << (offset*2)) | (temp->the_seq[j+1] >> (64-(offset*2)));
    }
    temp->the_seq[how_many-1] = temp->the_seq[how_many-1] << (offset*2);
  }
  return temp;

}


void sequence::print_boundaries() {
  if (this->boundaries->len() != 0) {
    for (int32_t i = this->boundaries->len()-1; i >= 0; i--) {
      boundary temp_boundary = this->boundaries->get_boundary(i);
      cout << "(" << temp_boundary.begin << ", " << temp_boundary.end << ")";
    }
  }
}

uint16_t sequence::count_mers(uint8_t mer_size) {
  uint16_t count = 0;
  for (uint16_t i = 0; i < this->boundaries->len(); i++) {
    boundary the_boundary = this->boundaries->get_boundary(i);
    if (the_boundary.end - the_boundary.begin >= mer_size)
      count += the_boundary.end - the_boundary.begin - mer_size + 1;
  }
  return count;
}

void sequence::update_mer_count_array(uint32_t *the_array, uint8_t array_mer_size, uint8_t mer_size) {
  for (uint16_t i = 0; i < this->boundaries->len(); i++) {
    boundary the_boundary = this->boundaries->get_boundary(i);
    for (uint16_t j = the_boundary.begin; j < the_boundary.end - mer_size +1; j++) {
      the_array[this->get_word(j, mer_size) >> ((mer_size - array_mer_size) * 2)]++;
    }
  }
}

void sequence::add_mers(uint32_t my_index, merray* the_array, uint32_t *mer_count_array, uint8_t array_mer_size, uint8_t mer_size) {
  for (uint16_t i = 0; i < this->boundaries->len(); i++) {
    boundary the_boundary = this->boundaries->get_boundary(i);
    if (the_boundary.end - the_boundary.begin >= mer_size) {
      for (uint16_t j=the_boundary.begin; j < the_boundary.end - mer_size + 1; j++) {
        if ((my_index < 1048575) && (j<4095)) { //just changed my_index < 1000000
          uint32_t index = (my_index << 12) + j; //20 bits for my_index = 10^6 possible ESTs, 12 bits for position = 4000 max EST len
          the_array->add_mer(index, this->get_word(j, mer_size), mer_count_array, mer_size, array_mer_size);
        }
      }
    }
  }
}

seq_word sequence::get_word(uint16_t which_word, uint8_t mer_size) {
  //unsigned long long int is 8 bytes long. i.e. we are limited to mer sizes of 32. Which should be plenty...
//  if (which_word + mer_size <= this->len) {
  uint16_t which_byte = which_word >> 5;// / 32;
  //uint8_t which_bit = (which_word % 32) * 2;
  uint8_t which_bit = (which_word & 31) << 1;
  seq_word the_val;
  if (which_bit != 0) {
    uint8_t shift_it = 64-which_bit;
    if ((shift_it >> 1) < mer_size) {
      the_val = (the_seq[which_byte] << which_bit) | this->the_seq[which_byte+1] >> (shift_it);
    }
    else {
      the_val=the_seq[which_byte] << which_bit;
    }
  }
  else {
    the_val = this->the_seq[which_byte];
  }
  return (the_val >> (64 - (mer_size << 1)));
}


string sequence::get_word_s(uint16_t which_word, uint8_t mer_size) {
  seq_word the_word = this->get_word(which_word, mer_size);
  string the_string = "";
  the_word = the_word << (64 - (2 * mer_size));
  
  for (uint16_t i = 0; i < mer_size; i++) {
    switch (the_word >> 62) {
      case 0 : the_string.append("A"); break;
      case 1 : the_string.append("C"); break;
      case 2 : the_string.append("G"); break;
      case 3 : the_string.append("T"); break;
    };
    the_word = the_word << 2;
  }
  return the_string;
}

string sequence::get_sequence() {
  //uint8_t temp_mer_size = this->mer_size;
  //this->mer_size = 1;
  string the_string = "";
  for (uint16_t i = 0; i < this->boundaries->len(); i++) {
    boundary the_boundary = this->boundaries->get_boundary(i);
    for (uint16_t j = the_boundary.begin; j < the_boundary.end; j++) {
      seq_word the_word = this->get_word(j,1);
      switch (the_word) {
        case 0 : the_string.append("A"); break;
        case 1 : the_string.append("C"); break;
        case 2 : the_string.append("G"); break;
        case 3 : the_string.append("T"); break;
      }
    }
    the_string.append("x");
  }
  //this->mer_size = temp_mer_size;
  return the_string.substr(0,the_string.length()-1); // Remove the last x.
}

bool sequence::find_matches(uint32_t my_index, match_list *matches, sequence_list *seq_list, merray *the_array, uint8_t mer_size) { //, bool update_search) {
  uint32_t mer_array_index;
  uint32_t mer_index;
  uint32_t which_sequence;
  uint16_t which_word;
  uint32_t array_len = the_array->len();
  for (uint16_t i = 0; i < this->boundaries->len(); ++i) {
    boundary the_boundary = this->boundaries->get_boundary(i);
    for (int32_t j = the_boundary.begin; j < the_boundary.end - mer_size + 1; j += globals::skip) {
      seq_word the_word = this->get_word(j, mer_size);
      mer_array_index = seq_list->find_mer(the_word, mer_size); //, update_search);
      mer_index = the_array->get_mer_index(mer_array_index);
      while ((seq_list->get_mer(mer_index, mer_size) == the_word)/* && (++mer_array_index < array_len)*/) {
        which_sequence = mer_index >> 12;
        if ((my_index < which_sequence)) {
          which_word = mer_index & 4095;
          position pos = {j, which_word};
          matches->add_match(my_index, which_sequence, pos);
        }
        //mer_array_index++;
        if (++mer_array_index < array_len)
          mer_index = the_array->get_mer_index(mer_array_index);
        else
          break;
      }
    }
  }
  return true;
}


/************************ SEQUENCE LIST ************************/

sequence_list::sequence_list(fasta *the_source, uint8_t mer_size, uint8_t matches_threshold, float sd_threshold) {
  this->how_many = the_source->count_sequences();
  //this->mer_size = mer_size;
  this->mer_count_array = NULL;
  this->mer_array = NULL;
  this->unique_mer_array = NULL;
  the_source->reset();
  this->sequences = new sequence*[this->how_many];
  for (uint32_t i = 0; i < this->how_many; i++) {
    this->sequences[i] = the_source->get_sequence(mer_size);
  }
  this->matches = new match_list(this->how_many, matches_threshold, sd_threshold);
  this->d1_calc= new d1();
  this->d2_calc = new d2();

}

sequence_list::~ sequence_list() {
  for (uint32_t i = 0; i < this->how_many; i++) {
    delete this->sequences[i];
  }
  delete[] this->sequences;
  delete this->matches;
  if (this->mer_count_array != NULL)
    delete[] this->mer_count_array;
  if (this->mer_array != NULL)
    delete this->mer_array;
  if (this->unique_mer_array != NULL)
    delete[] this->unique_mer_array;
  delete this->d1_calc;
  delete this->d2_calc;
}

uint32_t sequence_list::count_mers(uint8_t mer_size) {
  uint32_t total = 0;
  for (uint32_t i = 0; i < this->how_many; i++)
    total += this->sequences[i]->count_mers(mer_size);
  return total;
}

void sequence_list::create_mer_count_array(uint8_t mer_size) {
  //This makes mer_count_array contain the starting index for each mer...
  if (this->mer_count_array != NULL) {
    delete [] this->mer_count_array;
    this->mer_count_array = NULL;
  }
  uint32_t array_size = 1 << (this->array_mer_size * 2);
  this->mer_count_array = new uint32_t [array_size];
  for (uint32_t i = 0; i < array_size; i++) {
    this->mer_count_array[i] = 0;
  }
  for (uint32_t i = 0; i < this->how_many; i++) {
    this->sequences[i]->update_mer_count_array(this->mer_count_array, array_mer_size, mer_size);
  };
  for (uint32_t i = array_size-1; i > 0; i--) {
    mer_count_array[i] = mer_count_array[i-1];
  }
  mer_count_array[0]=0;
  for (uint32_t i = 1; i < array_size; i++) {
    mer_count_array[i] += mer_count_array[i-1];
  }
}

void sequence_list::create_merray(uint8_t array_mer_size, uint8_t mer_size) {
  this->array_mer_size = array_mer_size;

  uint32_t total = this->count_mers(mer_size);
  this->create_mer_count_array(mer_size);
  if (this->mer_array != NULL) {
    delete this->mer_array;
    this->mer_array = NULL;
  }
  this->mer_array = new merray(total);
  for (uint32_t i = 0; i < this->how_many; i++) {
    this->sequences[i]->add_mers(i, this->mer_array, this->mer_count_array, array_mer_size, mer_size);
  }
}

void sequence_list::merge_sort_mers(uint8_t array_mer_size, uint8_t mer_size) {
  if (array_mer_size > mer_size)
    array_mer_size = mer_size;
  //this->create_mer_count_array();
  this->create_merray(array_mer_size, mer_size);
  if (array_mer_size != mer_size)
    this->mer_array->merge_sort_mers(this, mer_size );
  delete[] this->mer_count_array;
  this->mer_count_array = NULL;
}

void sequence_list::print_mers(uint8_t mer_size) {
  this->mer_array->show_mers(this, mer_size);
}

void sequence_list::create_unique_mer_array(uint8_t mer_size) {
  this->unique_mers = this->mer_array->count_unique_mers(this, mer_size);
  if (this->unique_mer_array != NULL)
    delete[] this->unique_mer_array;
  this->unique_mer_array = this->mer_array->create_unique_mer_list(this, this->unique_mers, mer_size);
}

void sequence_list::print_unique_mers(uint8_t mer_size) {
  for (uint32_t i = 0; i < this->unique_mers; i++) {
    this->mer_array->print_mer(this, unique_mer_array[i], mer_size);

  }
}



uint32_t sequence_list::find_mer(seq_word the_mer, uint8_t mer_size) { //, bool update_search) {
  //binary search that gives the index into mer_array of the nearest matching mer. Should also try implementing
  //with an interpolated search...
  uint32_t left   = 0;
  uint32_t right  = this->unique_mers-1;
  uint32_t middle = 0;
  if (this->get_mer(this->mer_array->get_mer_index(unique_mer_array[left]), mer_size) == the_mer) {
//    if ((update_search) && (unique_mer_array[left]+1 < this->mer_array->len())) {
//      unique_mer_array[left]++;
//    }
    return unique_mer_array[left];
  }
  else if (this->get_mer(this->mer_array->get_mer_index(unique_mer_array[right]), mer_size) == the_mer) {
//    if ((update_search) && (unique_mer_array[right]+1 < this->mer_array->len())) {
//      unique_mer_array[right]++;
//    }
    return unique_mer_array[right];
  }
  
  while (left+1 < right) {
    middle = left + ((right - left) / 2);
    seq_word found_mer = this->get_mer(this->mer_array->get_mer_index(unique_mer_array[middle]), mer_size);
    if (found_mer < the_mer) {
      left = middle;
    }
    else if (found_mer > the_mer) {
      right = middle;
    }
    else {
//      if ((update_search) && (unique_mer_array[middle]+1 < this->mer_array->len())) {
//        ++unique_mer_array[middle];
//      }
      return unique_mer_array[middle];
    }
  }
  return 0;
}


seq_word sequence_list::get_mer(uint32_t merr_index, uint8_t mer_size) {
//  uint32_t which_sequence = merr_index >> 12;
//  uint16_t which_word = merr_index & 4095;
//  return this->sequences[which_sequence]->get_word(which_word);
  return this->sequences[merr_index >> 12]->get_word(merr_index & 4095, mer_size);
}

string sequence_list::get_mer_s(uint32_t merr_index, uint8_t mer_size) {
  uint32_t which_sequence = merr_index >> 12;
  uint16_t which_word = merr_index & 4095;
  return this->sequences[which_sequence]->get_word_s(which_word, mer_size);
}

void sequence_list::update_unique_mer_array(uint32_t min_index) {
  // This is very clever, if I do say so myself.

  uint32_t mer_array_len = this->mer_array->len()-1;
  while (((this->mer_array->get_mer_index(this->unique_mer_array[0]) >> 12) < min_index) && (this->unique_mer_array[0] < mer_array_len)) {
    ++this->unique_mer_array[0];
  }
  uint32_t u_mers = 1;

  for (uint32_t i=1; i < this->unique_mers; ++i) {
    if (this->unique_mer_array[i] <= this->unique_mer_array[i-1]) {
      this->unique_mer_array[i] = this->unique_mer_array[i-1];
    }
    else {
      while (((this->mer_array->get_mer_index(this->unique_mer_array[i]) >> 12) < min_index) && (this->unique_mer_array[i] < mer_array_len)) {
        ++this->unique_mer_array[i];
      }
      //++u_mers;
      this->unique_mer_array[u_mers++] = this->unique_mer_array[i];

    }
  }
  this->unique_mers = u_mers;
}



void sequence_list::find_matches(uint32_t start, uint32_t end, uint8_t bucket_size, uint8_t mer_size) {
//  this->merge_sort_mers(bucket_size, 0);
  if (start >= how_many) {
    start = 0;
  }
  if (end >= how_many) {
    end = how_many-1;
  }
  uint32_t step = how_many/100 + 1;

  if (globals::progress)
    cerr << "Creating the dictionary." << flush;
  this->merge_sort_mers(bucket_size, mer_size);
  if (globals::progress)
    cerr << "." << flush;

  this->create_unique_mer_array(mer_size);
  if (globals::progress)
    cerr << "." << flush;
  
  //this->mer_array->show_mers(this);

  //cerr << this->unique_mers << "----" << flush;

  for (uint32_t i = start; i <= end; i++) {
    if (i % step == 0) {
      if (i % (step*2)==0) {
        if (globals::progress)
          cerr << "Updating dictionary." << flush;
        update_unique_mer_array(i);
        if (globals::progress)
          cerr << ".. " << flush;
      };
      if (globals::progress)
        cerr << (i*100)/how_many << "%... " << flush;
    };

    this->sequences[i]->find_matches(i, this->matches, this, this->mer_array, mer_size); //, true);

    //sequence::set_mer_size(8);
    //match* a_match;
    //a_match->set_num_words1(this->sequences[i]);
    //sequence::set_mer_size(globals::mer_size);

    this->matches->find_clusters(i, this->sequences[i], this, false, this->d1_calc, this->d2_calc);
    this->matches->reset_list();

    if (globals::do_rc) {
      sequence *temp = sequences[i]->rc();
      temp->find_matches(i, this->matches, this, this->mer_array, mer_size);
      //sequence::set_mer_size(8);
      //a_match->set_num_words1(temp);
      
      //sequence::set_mer_size(globals::mer_size);

      this->matches->find_clusters(i, temp, this, true, this->d1_calc, this->d2_calc);
      this->matches->reset_list();
      
      delete temp;
    };
  }
  if (globals::show_compact)
    this->matches->show_clusters();
}

/************************ FASTA ************************/

fasta::fasta(char *source_file_name)
{ this->source_file.open(source_file_name, ios::in);
  source_file.seekg(0,ios::end);
  this->size = source_file.tellg();
  source_file.seekg(0,ios::beg);
  this->buffer = new char[size];
  source_file.read(this->buffer, this->size);
  source_file.close();
  this->where = 0;
}

uint32_t fasta::count_sequences() {
  //NOT BULLETPROOF IN THE SLIGHTEST
  uint32_t counter = 0;
  for (uint32_t i = 0; i < this->size; i++)
    if (this->buffer[i] == '>')
      counter++;
  return counter;
  
}

sequence* fasta::get_sequence(uint8_t mer_size) {
  // also not bullet proof. But I'll assume that people are using decent fasta files.
  while ((buffer[where] == '>') || (buffer[where] == ';')) {
    where ++;
    while ((this->where < this->size) && (buffer[where-1] != '\n'))
      where++;
  }
  uint32_t start = where;
  string the_sequence="";
  while ((this->where < this->size) && (buffer[where] != '>')) {
    if (buffer[where] == '\n' || buffer[where] == ' ' || buffer[where] == '\t' || buffer[where] == '\r') {
      if (where-start >= 1)
        the_sequence.append(&buffer[start], where-start);
      start = where+1;
    }
    where++;
  }
  sequence* new_sequence = new sequence(the_sequence, mer_size);
  return new_sequence;
}

void fasta::reset() {
  this->where=0;
}

fasta::~fasta()
{
  delete[] this->buffer;
}

/************************ MERRAY ************************/

merray::merray(uint32_t how_many)
{
  this->array_len = how_many;
  this->the_array = new unsigned int[how_many];
}

merray::~merray()
{
  delete[] this->the_array;
}

void merray::add_mer(uint32_t mer_index, seq_word the_mer, uint32_t *mer_count_array, uint8_t mer_size, uint8_t array_mer_size) {
  this->the_array[mer_count_array[the_mer >> (2*(mer_size - array_mer_size))]++] = mer_index;

}

void merray::show_mers() {
  for (uint32_t i = 0; i < this->array_len; i++) {
    cout << this->the_array[i] << " ";
  };
}

void merray::show_mers(sequence_list *the_sequences, uint8_t mer_size) {
  for (uint32_t i = 0; i < this->array_len; i++) {
    cout << (the_array[i]>>12) << " " << (the_array[i] & 4095) << "___";
    cout << the_sequences->get_mer_s(this->the_array[i], mer_size) << " " << flush;
  }
}

/*
uint64_t merray::get_mer(sequence_list *the_sequences, uint32_t mer_array_index, uint8_t which_mer) {
  //cout << mer_index << "--" << (int)which_mer << "..." << flush;
  return the_sequences->get_mer(this->the_array[mer_array_index], which_mer);
}
*/


/*uint32_t merray::get_mer_index(uint32_t mer_array_index) {
//   if (mer_array_index >= this->array_len) {
//     cerr << "what the crap -- " << mer_array_index << " " << this->array_len << " " << flush;
//   }
  //assert(mer_array_index < this->array_len);
   return this->the_array[mer_array_index];
}*/



void merray::print_mer(sequence_list* the_sequences, uint32_t mer_array_index, uint8_t mer_size) {
  cout << the_sequences->get_mer_s(this->the_array[mer_array_index], mer_size) << "\n";
}

void merray::insert(sequence_list *the_sequences, uint32_t start, uint32_t finish, uint8_t mer_size) {
  uint32_t the_val = this->the_array[finish];
  int64_t i = finish-1;
  while (i>=start && (the_sequences->get_mer(this->the_array[i], mer_size) > the_sequences->get_mer(the_val, mer_size))) {
    this->the_array[i+1] = this->the_array[i];
    i--;
  }
  this->the_array[i+1] = the_val;
}

void merray::insertion_sort_mers(sequence_list *the_sequences, uint32_t start, uint32_t finish, uint8_t mer_size) {
  for (uint32_t i=start+1; i<=finish; ++i) {
    this->insert(the_sequences, start, i, mer_size);
  }
}

void merray::merge(sequence_list *the_sequences, uint32_t left, uint32_t right, uint8_t mer_size) {
  //uint8_t which_mer = the_sequences->get_mer_size();
  uint32_t middle = ((left+right)/2) + 1;
  uint32_t t_left = left, t_middle = middle;
  uint32_t *temp_list = new uint32_t[right-left+1];
  for (uint32_t i = 0; i<=right-left; ++i) {
    if ((t_left < middle) && ((t_middle > right) || (the_sequences->get_mer(this->the_array[t_left], mer_size) <= the_sequences->get_mer(this->the_array[t_middle], mer_size)))) {
      temp_list[i] = the_array[t_left++];
    }
    else {
      temp_list[i]=the_array[t_middle++];
    }
  }
  for (uint32_t j = 0; j <= right-left; ++j) {
    this->the_array[j+left] = temp_list[j];
  };
  delete [] temp_list;
}

void merray::merge_sort_mers(sequence_list *the_sequences, uint32_t start, uint32_t finish, uint8_t mer_size) {
  //uint8_t which_mer = the_sequences->get_mer_size();
  if (finish > start) {
    uint32_t range = finish - start;
    if (range == 1) { // move this to insertion_sort_mers????
      if (the_sequences->get_mer(this->the_array[start], mer_size) > the_sequences->get_mer(this->the_array[finish], mer_size)) {
        uint32_t temp = this->the_array[start];
        this->the_array[start] = this->the_array[finish];
        this->the_array[finish] = temp;
      }
    }
    else if (range < 10) {
      this->insertion_sort_mers(the_sequences, start, finish, mer_size);
    }

    else {
      uint32_t middle = (finish + start)/2;
      this->merge_sort_mers(the_sequences, start, middle, mer_size);
      this->merge_sort_mers(the_sequences, middle+1, finish, mer_size);
      this->merge(the_sequences, start, finish, mer_size);
    }
  };
}

void merray::merge_sort_mers(sequence_list *the_sequences, uint8_t mer_size) {
  uint32_t *mer_count_array = the_sequences->get_mer_count_array();
  uint8_t array_mer_size = the_sequences->get_array_mer_size();
  uint32_t array_size = 1 << (array_mer_size * 2);

  if (mer_count_array[0] != 0)
    this->merge_sort_mers(the_sequences, 0, mer_count_array[0]-1, mer_size);

  for (uint32_t i=1; i < array_size; ++i) {

    if (mer_count_array[i-1] != mer_count_array[i])
       this->merge_sort_mers(the_sequences, mer_count_array[i-1], mer_count_array[i]-1, mer_size);
  }

}

uint32_t merray::count_unique_mers(sequence_list *the_sequences, uint8_t mer_size) {
  uint32_t temp=1;
  seq_word the_mer = the_sequences->get_mer(this->the_array[0], mer_size);
  for (uint32_t i = 1; i < this->array_len; ++i) {
    if (the_sequences->get_mer(this->the_array[i], mer_size) != the_mer) {
      the_mer = the_sequences->get_mer(this->the_array[i], mer_size);
      temp++;
    }
  }
  return temp;
}

uint32_t* merray::create_unique_mer_list(sequence_list* the_sequences, uint32_t unique_mer_count, uint8_t mer_size) {
  uint32_t *temp = new uint32_t[unique_mer_count];
  seq_word the_mer = the_sequences->get_mer(this->the_array[0], mer_size);
  temp[0]=0;
  uint32_t counter=1;
  for (uint32_t i = 1; i < this->array_len; ++i) {
    if (the_sequences->get_mer(this->the_array[i], mer_size) != the_mer) {
      the_mer = the_sequences->get_mer(this->the_array[i], mer_size);
      temp[counter++]=i;
    }
  }
  return temp;
}



