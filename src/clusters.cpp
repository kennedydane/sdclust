/**************************************
<clusters.cpp: class that uses union-find data structures for managing clusters>

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
#include "clusters.h"
#include "globals.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <sys/types.h>


using namespace std;

/************************clusters************************/

clusters::clusters(){
  this->parent = 0;
  this->size=0;
  this->total=0;
}

clusters::clusters(int32_t entries) {
  this->parent = new atomic<int32_t>[entries];
  for (int32_t i=0; i < entries; ++i) {
    parent[i].store(i, memory_order_relaxed);
  }
  this->size=entries;
  this->total=0;
}

clusters::~clusters() {
  if (this->parent != 0) {
    delete[] this->parent;
  }
}

int32_t clusters::find_parent(int32_t i) {
  int32_t p = parent[i].load(memory_order_acquire);
  while (p != i) {
    int32_t gp = parent[p].load(memory_order_acquire);
    if (gp != p) // path halving: skip a level if nobody else has changed it
      parent[i].compare_exchange_weak(p, gp, memory_order_release, memory_order_relaxed);
    i = gp;
    p = parent[i].load(memory_order_acquire);
  }
  return i;
}

void clusters::make_union(int32_t i, int32_t j) {
  this->total.fetch_add(1, memory_order_relaxed);
  while (true) {
    i = find_parent(i);
    j = find_parent(j);
    if (i == j)
      return;
    // Always hang the larger root under the smaller one. A fixed order means
    // two threads can never link roots into a cycle.
    if (i < j) {
      int32_t t = i; i = j; j = t;
    }
    int32_t expected = i;
    if (parent[i].compare_exchange_strong(expected, j, memory_order_acq_rel, memory_order_relaxed))
      return;
    // i stopped being a root while we were looking: try again from the top.
  }
}

void clusters::show_clusters() {
  // Clusters are printed in a fixed order whatever order the unions happened
  // in: members in ascending order on each line, and lines ordered by their
  // smallest member.
  vector<int32_t> first(this->size, -1), last(this->size, -1), next(this->size, -1);
  for (int32_t i=0; i < this->size; i++) {
    int32_t r = find_parent(i);
    if (first[r] < 0)
      first[r] = i;
    else
      next[last[r]] = i;
    last[r] = i;
  }
  for (int32_t i=0; i < this->size; i++) {
    int32_t r = find_parent(i);
    if (first[r] == i) {
      cout << i;
      for (int32_t j=next[i]; j >= 0; j = next[j]) {
        cout << " " << j;
      }
      cout << ".\n";
    }
  }
  if (globals::progress)
    cerr << "Total CLUSTER operations = " << this->total << "\n";
}
