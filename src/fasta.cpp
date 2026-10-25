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
#include <cstdlib>
#include <cctype>
#include <sys/types.h>
#include "globals.h"
#include "clusters.h"

#include <algorithm>
#include <mutex>
#include <sstream>
#include <thread>
#include <vector>

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

// Only k-mers that fit in an index entry are indexed: 20 bits for the
// sequence (about 10^6 ESTs) and 12 bits for the position (4095 bases).
#define INDEXABLE(my_index, j) (((my_index) < 1048575) && ((j) < 4095))

void sequence::count_index_mers(uint32_t my_index, merray* the_array, uint8_t mer_size, uint32_t *partition_counts) {
  for (uint16_t i = 0; i < this->boundaries->len(); i++) {
    boundary the_boundary = this->boundaries->get_boundary(i);
    if (the_boundary.end - the_boundary.begin >= mer_size) {
      for (uint16_t j = the_boundary.begin; j < the_boundary.end - mer_size + 1; j++) {
        if (INDEXABLE(my_index, j))
          the_array->count_mer(this->get_word(j, mer_size), partition_counts);
      }
    }
  }
}

void sequence::add_mers(uint32_t my_index, merray* the_array, uint8_t mer_size, uint32_t *partition_cursors) {
  for (uint16_t i = 0; i < this->boundaries->len(); i++) {
    boundary the_boundary = this->boundaries->get_boundary(i);
    if (the_boundary.end - the_boundary.begin >= mer_size) {
      for (uint16_t j=the_boundary.begin; j < the_boundary.end - mer_size + 1; j++) {
        if (INDEXABLE(my_index, j)) {
          uint32_t index = (my_index << 12) + j; //20 bits for my_index = 10^6 possible ESTs, 12 bits for position = 4000 max EST len
          the_array->add_mer(index, this->get_word(j, mer_size), partition_cursors);
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

bool sequence::find_matches(uint32_t my_index, match_list *matches, sequence_list *seq_list, merray *the_array, uint8_t mer_size) {
  uint32_t from, to;
  bool exact = the_array->exact();
  for (uint16_t i = 0; i < this->boundaries->len(); ++i) {
    boundary the_boundary = this->boundaries->get_boundary(i);
    for (int32_t j = the_boundary.begin; j < the_boundary.end - mer_size + 1; j += globals::skip) {
      seq_word the_word = this->get_word(j, mer_size);
      // Every occurrence of the_word in a later sequence, in sequence order.
      the_array->find_mers(the_word, my_index, &from, &to);
      for (uint32_t k = from; k < to; ++k) {
        uint32_t mer_index = the_array->get_mer_index(k);
        if (!exact && seq_list->get_mer(mer_index, mer_size) != the_word)
          continue;
        position pos = {(uint16_t)j, (uint16_t)(mer_index & 4095)};
        matches->add_match(my_index, mer_index >> 12, pos);
      }
    }
  }
  return true;
}


/************************ SEQUENCE LIST ************************/

sequence_list::sequence_list(fasta *the_source, uint8_t mer_size, uint8_t matches_threshold, float sd_threshold) {
  this->how_many = the_source->count_sequences();
  this->mer_array = NULL;
  this->the_clusters = NULL;
  this->matches_threshold = matches_threshold;
  this->sd_threshold = sd_threshold;
  the_source->reset();
  this->sequences = new sequence*[this->how_many];
  for (uint32_t i = 0; i < this->how_many; i++) {
    this->sequences[i] = the_source->get_sequence(mer_size);
  }
}

sequence_list::~ sequence_list() {
  for (uint32_t i = 0; i < this->how_many; i++) {
    delete this->sequences[i];
  }
  delete[] this->sequences;
  if (this->mer_array != NULL)
    delete this->mer_array;
  if (this->the_clusters != NULL)
    delete this->the_clusters;
}

seq_word sequence_list::get_mer(uint32_t merr_index, uint8_t mer_size) {
  return this->sequences[merr_index >> 12]->get_word(merr_index & 4095, mer_size);
}

string sequence_list::get_mer_s(uint32_t merr_index, uint8_t mer_size) {
  uint32_t which_sequence = merr_index >> 12;
  uint16_t which_word = merr_index & 4095;
  return this->sequences[which_sequence]->get_word_s(which_word, mer_size);
}

// Runs work(thread number) on each of threads threads and waits for them all.
template <class F> static void run_threads(uint32_t threads, F work) {
  vector<thread> workers;
  for (uint32_t t = 1; t < threads; t++)
    workers.push_back(thread(work, t));
  work(0);
  for (uint32_t t = 0; t < workers.size(); t++)
    workers[t].join();
}

static uint32_t thread_count() {
  if (globals::threads > 0)
    return globals::threads;
  uint32_t cores = thread::hardware_concurrency();
  return cores > 0 ? cores : 1;
}

void sequence_list::find_matches(uint32_t start, uint32_t end, uint8_t bucket_size, uint8_t mer_size) {
  if (start >= how_many) {
    start = 0;
  }
  if (end >= how_many) {
    end = how_many-1;
  }
  uint32_t threads = thread_count();
  uint32_t step = how_many/100 + 1;

  if (globals::progress)
    cerr << "Creating the dictionary (" << threads << " threads)..." << flush;
  this->mer_array = new merray(this, bucket_size, mer_size, threads);
  if (globals::progress)
    cerr << ". " << flush;

  this->the_clusters = new clusters(how_many);
  vector<match_list*> lists;
  for (uint32_t t = 0; t < threads; t++)
    lists.push_back(new match_list(how_many, this->matches_threshold, this->sd_threshold, this->the_clusters));

  // Threads take ESTs a chunk at a time, roughly in order. The extended table
  // (-C) is buffered per chunk and printed in chunk order, so its lines come
  // out in the same order as a single-threaded run.
  const uint32_t chunk = 16;
  uint32_t chunks = (end - start) / chunk + 1;
  atomic<uint32_t> next_chunk(0);
  mutex output_lock;
  vector<string> done(chunks);
  vector<bool> is_done(chunks, false);
  uint32_t next_print = 0;

  run_threads(threads, [&](uint32_t t) {
    match_list *matches = lists[t];
    uint32_t c;
    while ((c = next_chunk.fetch_add(1)) < chunks) {
      ostringstream out;
      uint32_t first = start + c * chunk;
      uint32_t last = first + chunk - 1 < end ? first + chunk - 1 : end;
      for (uint32_t i = first; i <= last; i++) {
        if (globals::progress && i % step == 0) {
          lock_guard<mutex> guard(output_lock);
          cerr << ((uint64_t)i*100)/how_many << "%... " << flush;
        }

        this->sequences[i]->find_matches(i, matches, this, this->mer_array, mer_size);
        matches->find_clusters(i, this->sequences[i], this, false, out);
        matches->reset_list();

        if (globals::do_rc) {
          sequence *temp = sequences[i]->rc();
          temp->find_matches(i, matches, this, this->mer_array, mer_size);
          matches->find_clusters(i, temp, this, true, out);
          matches->reset_list();
          delete temp;
        };
      }
      if (globals::show_extended) {
        lock_guard<mutex> guard(output_lock);
        done[c] = out.str();
        is_done[c] = true;
        while (next_print < chunks && is_done[next_print]) {
          cout << done[next_print];
          string().swap(done[next_print]);
          next_print++;
        }
      }
    }
  });

  for (uint32_t t = 0; t < threads; t++) {
    match_stats *stats = lists[t]->get_stats();
    globals::heuristic_success += stats->heuristic_success;
    globals::d1_success += stats->d1_success;
    globals::d2_success += stats->d2_success;
    delete lists[t];
  }

  if (globals::show_compact)
    this->the_clusters->show_clusters();
}

/************************ FASTA ************************/

fasta::fasta(const char *source_file_name)
{
  this->buffer = NULL;
  this->size = 0;
  this->where = 0;
  this->record = 0;
  this->source_file.open(source_file_name, ios::in | ios::binary);
  if (!this->source_file.is_open()) {
    cerr << "sd_clust: cannot open '" << source_file_name << "'\n";
    exit(EXIT_FAILURE);
  }
  source_file.seekg(0,ios::end);
  streamoff file_size = source_file.tellg();
  if (file_size < 0) {
    cerr << "sd_clust: cannot read '" << source_file_name << "'\n";
    exit(EXIT_FAILURE);
  }
  if (file_size == 0) {
    cerr << "sd_clust: '" << source_file_name << "' is empty\n";
    exit(EXIT_FAILURE);
  }
  if ((uint64_t)file_size > UINT32_MAX) {
    cerr << "sd_clust: '" << source_file_name << "' is too large (limit is 4 GB)\n";
    exit(EXIT_FAILURE);
  }
  this->size = (uint32_t)file_size;
  source_file.seekg(0,ios::beg);
  this->buffer = new char[size];
  source_file.read(this->buffer, this->size);
  if ((uint32_t)source_file.gcount() != this->size) {
    cerr << "sd_clust: failed to read '" << source_file_name << "'\n";
    exit(EXIT_FAILURE);
  }
  source_file.close();

  // The first thing in the file (after blank lines and ';' comments) must be a '>' header,
  // otherwise the sequence data before it would be silently merged into the first record.
  uint32_t i = 0;
  while (i < this->size) {
    if (buffer[i] == ' ' || buffer[i] == '\t' || buffer[i] == '\r' || buffer[i] == '\n') {
      i++;
    }
    else if (buffer[i] == ';' && (i == 0 || buffer[i-1] == '\n')) {
      while (i < this->size && buffer[i] != '\n')
        i++;
    }
    else
      break;
  }
  if (i == this->size) {
    cerr << "sd_clust: '" << source_file_name << "' contains no sequences\n";
    exit(EXIT_FAILURE);
  }
  if (buffer[i] != '>') {
    cerr << "sd_clust: '" << source_file_name << "' is not a FASTA file (expected a header line starting with '>')\n";
    exit(EXIT_FAILURE);
  }
}

uint32_t fasta::count_sequences() {
  uint32_t counter = 0;
  for (uint32_t i = 0; i < this->size; i++)
    if (this->buffer[i] == '>') {
      if (i==0) {
        counter++;
      }
      else if (buffer[i-1] == '\n') {
        counter++;

      }
    }
  return counter;
  
}

sequence* fasta::get_sequence(uint8_t mer_size) {
  string header = "";
  while ((this->where < this->size) && ((buffer[where] == '>') || (buffer[where] == ';') || isspace((unsigned char)buffer[where]))) {
    if (isspace((unsigned char)buffer[where])) { // blank lines before or between headers
      where++;
      continue;
    }
    uint32_t line_start = where;
    where ++;
    while ((this->where < this->size) && (buffer[where-1] != '\n'))
      where++;
    if (buffer[line_start] == '>') {
      header.assign(&buffer[line_start], where-line_start);
      while (!header.empty() && (header[header.length()-1] == '\n' || header[header.length()-1] == '\r'))
        header.erase(header.length()-1);
    }
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
  if (where-start >= 1) // last line without a trailing newline
    the_sequence.append(&buffer[start], where-start);

  if (the_sequence.empty()) {
    cerr << "sd_clust: sequence " << this->record << " (" << header << ") is empty\n";
    exit(EXIT_FAILURE);
  }
  if (the_sequence.length() > MAX_SEQUENCE_LEN) {
    cerr << "sd_clust: sequence " << this->record << " (" << header << ") is " << the_sequence.length()
         << " bases long; the maximum supported length is " << MAX_SEQUENCE_LEN << "\n";
    exit(EXIT_FAILURE);
  }
  this->record++;
  sequence* new_sequence = new sequence(the_sequence, mer_size);
  return new_sequence;
}

void fasta::reset() {
  this->where=0;
  this->record=0;
}

fasta::~fasta()
{
  delete[] this->buffer;
}

/************************ MERRAY ************************/

merray::merray(sequence_list *the_sequences, uint8_t bucket_size, uint8_t mer_size, uint32_t threads)
{
  if (bucket_size > mer_size)
    bucket_size = mer_size;
  this->suffix_bits = 2 * (mer_size - bucket_size);
  this->drop_bits = this->suffix_bits > 32 ? this->suffix_bits - 32 : 0;
  this->bucket_count = 1u << (2 * bucket_size);
  const uint8_t partition_bits = 2 * bucket_size < 10 ? 2 * bucket_size : 10;
  this->partition_shift = 2 * bucket_size - partition_bits;
  const uint32_t partitions = 1u << partition_bits;
  uint32_t how_many = the_sequences->len();

  // Each thread takes a fixed, contiguous share of the sequences, so it can
  // be given its own place to write in every partition.
  vector<uint32_t> first_seq(threads + 1);
  for (uint32_t t = 0; t <= threads; t++)
    first_seq[t] = (uint32_t)(((uint64_t)how_many * t) / threads);
  vector<uint32_t> per_thread((size_t)threads * partitions, 0);

  // Count the k-mers in each bucket and, per thread, in each partition...
  this->bucket_fill = new atomic<uint32_t>[this->bucket_count];
  for (uint32_t b = 0; b < this->bucket_count; b++)
    this->bucket_fill[b].store(0, memory_order_relaxed);
  run_threads(threads, [&](uint32_t t) {
    for (uint32_t j = first_seq[t]; j < first_seq[t+1]; j++)
      the_sequences->get_sequence(j)->count_index_mers(j, this, mer_size, &per_thread[(size_t)t * partitions]);
  });

  // ...turn the counts into where each bucket starts, and where each thread
  // starts writing in each partition...
  this->bucket_start = new uint32_t[this->bucket_count + 1];
  uint32_t total = 0;
  for (uint32_t b = 0; b < this->bucket_count; b++) {
    this->bucket_start[b] = total;
    total += this->bucket_fill[b].load(memory_order_relaxed);
  }
  this->bucket_start[this->bucket_count] = total;
  delete[] this->bucket_fill;
  this->bucket_fill = NULL;
  for (uint32_t p = 0; p < partitions; p++) {
    uint32_t where = this->bucket_start[p << this->partition_shift];
    for (uint32_t t = 0; t < threads; t++) {
      uint32_t count = per_thread[(size_t)t * partitions + p];
      per_thread[(size_t)t * partitions + p] = where;
      where += count;
    }
  }
  this->array_len = total;
  this->the_array = new uint64_t[total];

  // ...write every entry into its partition...
  run_threads(threads, [&](uint32_t t) {
    for (uint32_t j = first_seq[t]; j < first_seq[t+1]; j++)
      the_sequences->get_sequence(j)->add_mers(j, this, mer_size, &per_thread[(size_t)t * partitions]);
  });

  // ...then, one partition at a time, move each entry into its bucket and
  // sort the buckets. Entries are unique, so the result doesn't depend on the
  // order they went in.
  atomic<uint32_t> next_partition(0);
  run_threads(threads, [&](uint32_t) {
    vector<uint64_t> scratch;
    vector<uint32_t> cursor(1u << this->partition_shift);
    uint32_t p;
    while ((p = next_partition.fetch_add(1)) < partitions) {
      uint32_t first_bucket = p << this->partition_shift;
      uint32_t last_bucket = first_bucket + (1u << this->partition_shift);
      uint64_t *begin = this->the_array + this->bucket_start[first_bucket];
      uint64_t *end = this->the_array + this->bucket_start[last_bucket];
      if (this->partition_shift > 0 && end - begin > 1) {
        scratch.assign(begin, end);
        for (uint32_t b = first_bucket; b < last_bucket; b++)
          cursor[b - first_bucket] = this->bucket_start[b];
        for (size_t k = 0; k < scratch.size(); k++) {
          seq_word the_mer = the_sequences->get_mer((uint32_t)scratch[k], mer_size);
          this->the_array[cursor[(the_mer >> this->suffix_bits) - first_bucket]++] = scratch[k];
        }
      }
      for (uint32_t b = first_bucket; b < last_bucket; b++)
        if (this->bucket_start[b+1] - this->bucket_start[b] > 1)
          sort(this->the_array + this->bucket_start[b], this->the_array + this->bucket_start[b+1]);
    }
  });
}

merray::~merray()
{
  delete[] this->the_array;
  delete[] this->bucket_start;
}

void merray::count_mer(seq_word the_mer, uint32_t *partition_counts) {
  uint32_t b = the_mer >> this->suffix_bits;
  this->bucket_fill[b].fetch_add(1, memory_order_relaxed);
  partition_counts[b >> this->partition_shift]++;
}

void merray::add_mer(uint32_t mer_index, seq_word the_mer, uint32_t *partition_cursors) {
  uint64_t suffix = (the_mer & ((((seq_word)1) << this->suffix_bits) - 1)) >> this->drop_bits;
  uint32_t b = the_mer >> this->suffix_bits;
  this->the_array[partition_cursors[b >> this->partition_shift]++] = (suffix << 32) | mer_index;
}

void merray::find_mers(seq_word the_mer, uint32_t after_seq, uint32_t *from, uint32_t *to) {
  uint64_t suffix = (the_mer & ((((seq_word)1) << this->suffix_bits) - 1)) >> this->drop_bits;
  uint32_t b = the_mer >> this->suffix_bits;
  uint64_t *bucket_begin = this->the_array + this->bucket_start[b];
  uint64_t *bucket_end = this->the_array + this->bucket_start[b+1];
  // First entry for this k-mer from a sequence after after_seq. If after_seq+1
  // doesn't fit in 20 bits the sum carries into the k-mer bits, past them all.
  uint64_t *lo = lower_bound(bucket_begin, bucket_end, (suffix << 32) + (((uint64_t)after_seq + 1) << 12));
  uint64_t *hi = lo;
  while (hi < bucket_end && (*hi >> 32) == suffix)
    ++hi;
  *from = lo - this->the_array;
  *to = hi - this->the_array;
}
