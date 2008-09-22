/**************************************
<globals.cpp: A class for storing some global settings.>
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


#include "globals.h"
#include <fstream>
//#include <sys/types.h>

using namespace std;

uint8_t  globals::bucket_size       = 12; // Recommended max of 13. // Could be done by check free memory though.

uint8_t globals::d1_word_length     = 6;
uint8_t globals::d1_threshold       = 75;
uint8_t globals::d1_window_length   = 100;

uint8_t  globals::d2_threshold      = 40;
uint8_t  globals::d2_word_length    = 6;
uint8_t  globals::d2_window_length  = 100;

//uint8_t  globals::d1_threshold = (2*(globals::d2_window_length-globals::d2_word_length+1)-globals::d2_threshold)/2;

uint8_t  globals::mer_size          = 14;
uint8_t  globals::matches_threshold = 3;
float    globals::sd_threshold      = 4;
uint8_t  globals::skip              = 10;

uint64_t globals::heuristic_success = 0;
uint64_t globals::d1_success        = 0;
uint64_t globals::d2_success        = 0;
bool     globals::progress          = false;
bool     globals::show_compact      = false;
bool     globals::show_extended     = false;
bool     globals::complete          = false;
bool     globals::range             = false;
bool     globals::heuristic_only    = false;
bool     globals::do_rc             = true;
ofstream globals::scores_comparison;
string   globals::version           = "0.1.0 pre-alpha";

uint32_t globals::A_count=0;
uint32_t globals::C_count=0;
uint32_t globals::G_count=0;
uint32_t globals::T_count=0;
uint32_t globals::N_count=0;


