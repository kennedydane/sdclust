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
#include <cstdint>
#include <atomic>
//#include "common.h"

/**
	@author Dane Kennedy,,, <kennedy.dane@gmail.com>
*/

using namespace std;

// A union-find table that many threads can share without locks. Each entry
// points at its parent; a root points at itself. Unions only ever change a
// root (with a compare-and-swap), and path halving only ever points an entry
// at one of its ancestors, so concurrent finds and unions are safe.
class clusters{
  public:
    clusters();
    clusters(int32_t entries);
    ~clusters();
    int32_t find_parent(int32_t i);
    void make_union(int32_t i, int32_t j);
    void show_clusters();
  private:
    atomic<int32_t> *parent;
    int32_t size;
    atomic<int32_t> total;
};

#endif
