/**************************************
clusters.h: Just some union find data structures

**************************************/
#ifndef CLUSTERS_H
#define CLUSTERS_H

#include <sys/types.h>
//#include "common.h"

/**
	@author Dane Kennedy,,, <kennedy.dane@gmail.com>
*/

using namespace std;

class entry {
  public:
    void initialise(int32_t i);
  protected:
    int32_t cluster; // refers to index of root
    int32_t rank; // for find-union
    int32_t next;
    int32_t last;
    int32_t match;
    friend class clusters;
};



class clusters{
  public:
    clusters();
    clusters(int32_t entries);
    ~clusters();
    int32_t find_parent(int32_t i);
    void make_union(int32_t i, int32_t j);
    void show_clusters();
  private:
    entry *the_entries;
    int32_t size;
    int32_t total;
};

#endif
