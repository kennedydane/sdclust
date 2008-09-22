/**************************************
<globals.h: A class for storing some global settings.>
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
