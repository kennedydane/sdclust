#ifndef GLOBALS_H
#define GLOBALS_H

#include <sys/types.h>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

typedef uint64_t seq_word;


class globals {
public:
    static uint8_t d1_word_length;
    static uint8_t d1_threshold;
    static uint8_t d1_window_length;

    static uint8_t  d2_threshold;
    static uint8_t  d2_word_length;
    static uint8_t  d2_window_length;
    static uint8_t  mer_size;
    static uint8_t  bucket_size;
    static uint8_t  matches_threshold;
    static float    sd_threshold;

    static uint64_t heuristic_success;
    static uint64_t d1_success;
    static uint64_t d2_success;
    static ofstream scores_comparison;
    static uint8_t  skip;
    static bool     progress;
    static bool     show_compact;
    static bool     show_extended;
    static bool     complete;
    static bool     range;
    static bool     heuristic_only;
    static bool     do_rc;
    static string   version;

    static uint32_t A_count;
    static uint32_t C_count;
    static uint32_t G_count;
    static uint32_t T_count;
    static uint32_t N_count;
};
#endif
