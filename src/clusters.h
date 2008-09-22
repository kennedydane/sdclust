/**************************************
<clusters.h: class that uses union-find data structures for managing clusters>

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
