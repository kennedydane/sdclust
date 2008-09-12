/**************************************
sd_clust.cpp: Wadda wadda

**************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


/**
	@author Dane Kennedy,,, <kennedy.dane@gmail.com>
 */

#include <cstdlib>
#include <string>
#include "fasta.h"
#include "match.h"
#include "globals.h"
#include <fstream>
#include <iomanip>

#include <sys/types.h>

#include <getopt.h>

using namespace std;


char* source_file = "";

//time_t ts, te;
uint64_t ts, te; //timing variables (time start, time end)

void show_version(); //forward declerations
void show_help();

bool fasta_stats=false;
void show_fasta_stats(sequence_list *the_sequences);

static char *allowed_options = "cCeEfhHnpqQvs:w:RS:t:T:W:L:1:2:";

bool compare_two = false;
bool compare_all = false;
uint32_t first=0, second=1;

static struct option long_options[] = {
  {"progress",0,0, 'p'},
  {"sd",1,0,'s'},
  {"show_compact", 0, 0, 'c'},
  {"show_extended",0,0,'C'},
  {"single",0,0,'q'},
  {"complete",0,0,'Q'},
  {"compare",0,0,'e'},
  {"compare_all",0,0,'E'},
  {"fasta",0,0,'f'},
  {"no_rc",0,0,'n'},
  {"help",0,0,'h'},
  {"heur_only",0,0,'H'},
  {"version",0,0,'v'},
  {"wc_thresh",1,0,'t'},
  {"word_len",1,0,'w'},
  {"skip",1,0,'S'},
  {"d2_w",1,0,'W'},
  {"d2_t",1,0,'T'},
  {"d2_l",1,0,'L'},
  {"range",0,0,'R'},
  {"1",1,0,'1'},
  {"2",1,0,'2'},
  {0, 0, 0, 0}
};

void parseCommandLine(int argc, char* argv[]) {
  int opt = getopt_long(argc, argv, allowed_options, long_options, 0);
  while (opt != -1) {
    switch (opt) {
      case 'c' : //show compact clusters
        globals::show_compact = true;
        break;
      case 'C': //show extended cluster table
        globals::show_extended = true;
        break;
      case 'h': //show_help
        show_help();
        exit(0);
      case 'H': //only do heuristics (no d_2)
        globals::heuristic_only = true;
        break;
      case 'R': //cluster range
        globals::range = true;
        break;
      case 'e' : //compare two sequences
        compare_two=true;
        break;
      case 'E' : //compare every single sequence against every other sequence
        compare_all=true;
        break;
      case 'f': // give some stats about the sequences
        fasta_stats=true;
        break;
      case '1': //index of first sequence
        first = atoi(optarg);
        break;
      case '2': //index of second sequence
        second = atoi(optarg);
        break;
      case 'p' : //show _progress
        globals::progress=true;
        break;
      case 'q': //do single linkage
        globals::complete=false;
        break;
      case 'Q': //do complete linkage
        globals::show_compact=false; //It's going to output garbage anyway
        globals::complete=true;
        break;
      case 'n': // reverse compliment
        globals::do_rc = false;
        break;
      case 'v': //print version;
        show_version();
        exit(0);
      case 's' : // sd threshold
        globals::sd_threshold = atof(optarg);
        break;
      case 't' : // match count threshold
        globals::matches_threshold = atoi(optarg);
        break;
      case 'w' : //word_len
        globals::mer_size = atoi(optarg);
        break;
      case 'S': //skip size
        globals::skip = atoi(optarg);
        break;
      case 'T': //d2_threshold
        globals::d2_threshold = atoi(optarg);
        globals::d1_threshold = (2*(globals::d2_window_length-globals::d2_word_length+1)-globals::d2_threshold)/2;
        break;
      case 'W': //d2_word length
        globals::d2_word_length = atoi(optarg);
        globals::d1_word_length = globals::d2_word_length;
        globals::d1_threshold = (2*(globals::d2_window_length-globals::d2_word_length+1)-globals::d2_threshold)/2;
        break;
      case 'L': //d2 window length
        globals::d2_window_length = atoi(optarg);
        globals::d1_window_length = globals::d2_window_length;
        globals::d1_threshold = (2*(globals::d2_window_length-globals::d2_word_length+1)-globals::d2_threshold)/2;
      break;
    }
    opt = getopt_long(argc, argv, allowed_options, long_options, 0);
  }
  if (optind<argc) {
    source_file = argv[optind];
  }
  else {
    show_help();
    exit(1);
  }
}

void show_version() {
  cout << "sd_clust version " << globals::version << "\n" << flush;
}

void show_help() {
  show_version();
  cout << "\nUsage: sd_clust [options] <filename>\n\n";
  cout << "  [options]\n";
  cout << "    -v, --version        : show version number (stdout)\n";
  cout << "    -h, --help           : show this help screen (stdout)\n";
  cout << "    -p, --progress       : show percentage progress (stderr)\n";
  cout << "    -f, --fasta          : show summary about the fasta sequences (stderr)\n"; 
  cout << "\n";
  cout << "    -c, --show_compact   : show the compact cluster table (stdout)\n";
  cout << "    -C, --show_extended  : show the extended cluster table (stdout)\n";
  cout << "    -q, --single         : do single linkage clustering (default)\n";
  cout << "    -Q, --complete       : find all overlapping ESTs\n";
  cout << "    -n, --no_rc          : don't evaluate reverse compliment\n";
  cout << "    -H, --heur_only      : cluster based on heuristics only\n";
  cout << "\n";
  cout << "    -S, --skip <int>     : set number of words to skip (default = " << (int)globals::skip << ")\n";
  cout << "    -s, --sd <float>     : set standard deviation threshold (default = " << globals::sd_threshold << ")\n";
  cout << "    -t, --wc_thresh <int>: set word count threshold (default = " << (int)globals::matches_threshold << ")\n";
  cout << "    -w, --word_len <int> : set word length (default = " << (int)globals::mer_size << ")\n";
  cout << "\n";
  cout << "    -L, --d2_l <int>     : set d2 window length (default = " << (int)globals::d2_window_length << ")\n";
  cout << "    -T, --d2_t <int>     : set d2 threshold (default = " << (int)globals::d2_threshold << ")\n";
  cout << "    -W, --d2_w <int>     : set d2 word length (default = " << (int)globals::d2_word_length << ")\n";
  cout << "\n";
  cout << "    -e, --compare        : Find d2 score of two sequences (stdout)\n";
  cout << "    -E, --compare_all    : Find d2 score of every pair of sequences and\n";
  cout << "                           outputs successfull matches in -C format (stdout)\n";
  cout << "    -1, --1 <int>        : first sequence to be compared (default = 0)\n";
  cout << "    -2, --2 <int>        : second sequence to be compared (default = 1)\n";
  cout << "\n";
  cout << "  <filename> is a fasta formated file of sequences\n";
  cout << "\n";
  cout << "Notes:\n";
  cout << "  1. The format of the compact cluster table (-c) is simply that each line\n";
  cout << "     contains all the sequences in the cluster. Sequences are numbered from 0\n";
  cout << "     and each line is terminated with a \".\".\n";
  cout << "  2. The format of the extended cluster table (-C) is each line is a tab\n";
  cout << "     separated list of:\n";
  cout << "       Sequence1\tSequence2\tReverseComplement\td2Score\n";
  cout << "     where Sequence1 and Sequence2 are the numbers of the sequences in the\n";
  cout << "     file (starting at 0), ReverseComplement indicates whether the match was\n";
  cout << "     with the reverse complement (0=no, 1=yes) and d2Score indicates the d2\n";
  cout << "     score between the sequences.\n";
  cout << "  3. Each time sd_clust finishes an analysis it appends the log file\n";
  cout << "     (sd_clust.log). The line it appends the following info:\n";
  cout << "       f_name w_thresh w_size sd_thresh skip_val secs heur_pass d1_pass d2_pass\n";
  cout << "     where f_name is the fasta file's name, w_thresh is the heuristic's\n";
  cout << "     word count threshold, w_size is the heuristic's word size, sd_thresh\n";
  cout << "     is the standard deviation threshold, skip_val specifies every\n";
  cout << "     skip_valth word to check, secs is the time to complete, heur_pass\n";
  cout << "     is the number of potential matches the heuristic found, d1 pass is\n";
  cout << "     is the number of times the d1 heuristic passes and d2_pass is the\n";
  cout << "     number of matches that were found.\n";
  cout << "\n";
  cout << "Examples:\n";
  cout << "  sd_clust -c my_seqences.fasta\n";
  cout << "    This will use the default d2 and sd options and output the compact cluster\n";
  cout << "    table to standard output.\n";
  cout << "\n";
  cout << "  sd_clust -C -Q mysequences.fasta\n";
  cout << "    This will output the extended cluster table showing all overlaps between\n";
  cout << "    ESTs.\n";
  cout << "\n";
  cout << "  sd_clust -c -t 20 -w 16 -s 2.0 -S 3 -n my_sequences.fasta\n";
  cout << "    This changes the heuristic to use a word length of 16, a word count\n";
  cout << "    threshold of 20, a standard deviation threshold of 2 and  it only checks\n";
  cout << "    every 3rd word, but doesn't check the reverse complement. Compact cluster\n";
  cout << "    table is output.\n";
  cout << "\n";
  cout << "  sd_clust -e -1 34 -2 83 my_sequences.fasta\n";
  cout << "    This will calculate the d2 score between sequence #34 and sequence #83.\n";
  cout << "    Note that sequences are numbered from 0.\n" << flush;
}

void echoOptions() {
  cerr << "\nRunning sd_clust with the following options:\n";
  cerr << "  Source file      =\t" << source_file << "\n";
  cerr << "  Match Threshhold =\t" << (int)globals::matches_threshold << "\n";
  cerr << "  Word Length      =\t" << (int)globals::mer_size << "\n";
  cerr << "  SD Threshold     =\t" << globals::sd_threshold << "\n";
  cerr << "  Skip value       =\t" << (int)globals::skip << "\n";
  cerr << "  RC               =\t" << (int)globals::do_rc << "\n";
  cerr << "  d2 threshold     =\t" << (int)globals::d2_threshold << "\n";
  cerr << "  d2 word length   =\t" << (int)globals::d2_word_length << "\n";
  cerr << "  d2 window length =\t" << (int)globals::d2_window_length << "\n";
  cerr << "  d1 threshold     =\t" << (int)globals::d1_threshold << "\n";
  cerr << "  d1 word length   =\t" << (int)globals::d1_word_length << "\n";
  cerr << "  d1 window length =\t" << (int)globals::d1_window_length << "\n";
  if (!globals::complete) 
    cerr << "  Single Linkage\n\n" << flush;
  else
    cerr << "  Complete Linkage\n\n" << flush;
  
}

void show_fasta_stats(sequence_list *the_sequences) {
  uint32_t no_seqs = the_sequences->len();
  uint64_t total_len = 0;
  
  sequence *temp_seq;
  uint16_t shortest = the_sequences->get_sequence(0)->get_len();
  uint16_t longest = shortest;
  
  for (uint32_t i=0; i < no_seqs; ++i) {
    temp_seq = the_sequences->get_sequence(i);
    uint16_t seq_len = temp_seq->get_len();
    total_len += seq_len;
    if (seq_len < shortest) shortest = seq_len;
    if (seq_len > longest) longest = seq_len;
  }

  cerr << "Summary of                 :\t" << source_file << "\n";
  cerr << "  Number of sequences      :\t" << no_seqs << "\n";
  cerr << "  Number of A's            :\t" << globals::A_count << "\n";
  cerr << "  Number of C's            :\t" << globals::C_count << "\n";
  cerr << "  Number of G's            :\t" << globals::G_count << "\n";
  cerr << "  Number of T's            :\t" << globals::T_count << "\n";
  cerr << "  Number of N's            :\t" << globals::N_count << "\n";
  cerr << "  Total number of bases    :\t" << total_len << "\n";
  cerr << "  Average sequence length  :\t" << (int)((float)total_len)/((float)no_seqs) << "\n";
  cerr << "  Longest sequence length  :\t" << longest << "\n";
  cerr << "  Shortest sequence length :\t" << shortest << "\n\n" << flush;
}

void do_log_file() {
  ofstream logfile;
  logfile.open("sd_clust.log",ios::app);
  if (!logfile.is_open()) {
    cerr << "LOG FILE NOT OPEN...\n";
  }
  else {
    logfile << setiosflags(ios::fixed) << setprecision(2);
    logfile << source_file << "\t" << (int)globals::matches_threshold << "\t" << (int)globals::mer_size << "\t" << globals::sd_threshold << "\t" << (int)globals::skip << "\t" << ((double)(te-ts))/CLOCKS_PER_SEC << "\t" << globals::heuristic_success << "\t" << globals::d1_success << "\t" << globals::d2_success << "\n";

    logfile.close();
  }
}

int main(int argc, char *argv[]) { 
  srandom(563573);
  parseCommandLine(argc,argv);

  //ts = time((time_t *)NULL);
  ts = clock();

  fasta *fasta_file = new fasta(source_file);
  sequence_list *the_sequences = new sequence_list(fasta_file, globals::mer_size, globals::matches_threshold, globals::sd_threshold);

  delete fasta_file;

  if (fasta_stats) {
    show_fasta_stats(the_sequences);
    exit(0);
  }

  if (globals::progress) {
    echoOptions();
    show_fasta_stats(the_sequences);
  }

  if (compare_two) {
    //d1* d1_calc = new d1();
    d2* d2_calc = new d2();
    //sequence::set_mer_size(globals::d2_word_length);
    //match *a_match = new match();
    sequence *one = the_sequences->get_sequence(first);
    cout << "Sequence #" << first << " [" <<  one->get_len() << "] =\n" << one->get_sequence() << "\n\n" << flush;
    sequence *two = the_sequences->get_sequence(second);
    cout << "Sequence #" << second << " [" <<  two->get_len() << "] =\n" << two->get_sequence() << "\n\n" << flush;
    int32_t forward_d2_score = d2_calc->calc_score(one, two, 0, one->get_len()-globals::d2_word_length, 0, two->get_len()-globals::d2_word_length);
    sequence* two_rc = two->rc();
    cout << "RC of Sequence #" << second << " [" <<  two_rc->get_len() << "] =\n" << two_rc->get_sequence() << "\n\n" << flush;
    int32_t rc_d2_score = d2_calc->calc_score(one, two_rc, 0, one->get_len()-globals::d2_word_length, 0, two_rc->get_len()-globals::d2_word_length);
    cout << "Forward d2 score = " << forward_d2_score << ", Reverse d2 score = " << rc_d2_score << "\n" <<flush;
    delete two_rc;
    //delete a_match;
    delete the_sequences;
    exit(0);
  }
  else if (compare_all) {
    d1 *d1_calc = new d1();
    d2 *d2_calc = new d2();
    uint16_t dummy1, dummy2, dummy3, dummy4;
    /*this->set_num_words1(s_seq);
    uint32_t in_common = this->common_words(s_seq, t_seq);
    if (in_common>=65) {cout << source << " " << target << " " << in_common << "\n";};*/
    
    //sequence::set_mer_size(globals::d2_word_length);
    //match *a_match = new match();
    for (uint32_t i=0; i < the_sequences->len()-1; i++) {
      if (globals::progress)
          cerr << i << " " << flush;
      sequence *one = the_sequences->get_sequence(i);
      d1_calc->set_num_words1(one);
      for (uint32_t j=i+1; j < the_sequences->len(); j++) {
        sequence *two = the_sequences->get_sequence(j);
        sequence *two_rc = two->rc();
        if (d1_calc->calc_score(two, &dummy1, &dummy2, &dummy3, &dummy4)>=globals::d1_threshold) {
          globals::d1_success++;
          int32_t forward_score = d2_calc->calc_score(one, two, 0, one->get_len()-globals::d2_word_length, 0, two->get_len()-globals::d2_word_length);
          if (forward_score <= globals::d2_threshold) {
            globals::d2_success++;
            cout << i << "\t" << j << "\t0\t" << forward_score << "\n" << flush;
          }
        };
        if (d1_calc->calc_score(two_rc, &dummy1, &dummy2, &dummy3, &dummy4)>=globals::d1_threshold) {
          globals::d1_success++;
          int32_t reverse_score = d2_calc->calc_score(one, two_rc, 0, one->get_len()-globals::d2_word_length, 0, two_rc->get_len()-globals::d2_word_length);
        
          if (reverse_score <= globals::d2_threshold) {
            globals::d2_success++;
            cout << i << "\t" << j << "\t1\t" << reverse_score << "\n" << flush;
          }
        }
        delete(two_rc);

      }
    }
  }
  else {
    if (!globals::range) {
      first = 0;
      second = the_sequences->len();
    }
      the_sequences->find_matches(first, second, globals::bucket_size, globals::mer_size);
  }

  //te = time((time_t *)NULL);
  te = clock();
  do_log_file();

  delete the_sequences;

  return EXIT_SUCCESS;
}
