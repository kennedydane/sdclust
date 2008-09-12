/**************************************
clusters.cpp: Wadda wadda

**************************************/
#include "clusters.h"
#include "globals.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <sys/types.h>


using namespace std;

void entry::initialise(int32_t i) {
  this->cluster=i;
  this->last=i;
  this->match=i;
  this->next=-1;
  this->rank=0;
}

/************************clusters************************/

clusters::clusters(){
  this->the_entries = 0;
  this->size=0;
}

clusters::clusters(int32_t entries) {
  this->the_entries = new entry[entries];
  for (int32_t i=0; i < entries; ++i) {
    the_entries[i].initialise(i);
  }
  this->size=entries;
  this->total=0;
}

clusters::~clusters() {
  if (this->the_entries != 0) {
    delete[] this->the_entries;
  }
}

int32_t clusters::find_parent(int32_t i) {
  int32_t k, r;
  r = i;
  while (r != the_entries[r].cluster)
    r = the_entries[r].cluster;
  k = i;
  while (k != the_entries[k].cluster) {
    i = k;
    k = the_entries[k].cluster;
    the_entries[i].cluster = r;
  }

  return r;
}

void clusters::make_union(int32_t i, int32_t j) {
  i = find_parent(i);
  j = find_parent(j);
  this->total++;
  if (i != j) {
    int32_t q;
    if (the_entries[i].rank > the_entries[j].rank) {
      q = the_entries[i].last;
      the_entries[q].next = j;
      the_entries[i].last = the_entries[j].last;
      the_entries[j].cluster = i;
    }
    else {
      q = the_entries[j].last;
      the_entries[q].next = i;
      the_entries[j].last = the_entries[i].last;
      the_entries[i].cluster = j;

      if (the_entries[i].rank == the_entries[j].rank)
        the_entries[j].rank++;

    }
  }
}

void clusters::show_clusters() {
  int32_t i,j;
  for (i=0; i< this->size; i++) {
    if (the_entries[i].cluster == i) {
      cout << i;
      for (j=the_entries[i].next; j >=0; j = the_entries[j].next) {
        cout << " " << j;
      }
      cout << ".\n";
    }
  }
  if (globals::progress)
    cerr << "Total CLUSTER operations = " << this->total << "\n";
}


