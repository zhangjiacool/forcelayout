/* Copyright (C) 2013-2014 Kari Pahula

   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:
   
   The above copyright notice and this permission notice (including the next
   paragraph) shall be included in all copies or substantial portions of the
   Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/

#ifndef _WORLD_H
#define _WORLD_H

#define RELAX_EXTRA 1

struct pair {
  double x, y;
};

struct vertex {
  struct pair pos;
  float radius;
  float weight;
};

struct edge {
  float weight;
};
 
struct world_work {
  int start, end;
  double energy;
  void *extra;
  struct pair data[];
};

struct world {
  struct thread_control *pool;
  double allforces;
  struct edge *edges;
  struct vertex *vertices;
  double *dist;
  int *mapping;		//index: internal id > 0
  int nitems;
  int maxid;
  int *r_mapping;	//index: id
  double energy;
  double maxmove;
  double repulsioncap;
  double world_weight_inv;
  struct world_work **world_work;
  struct options *options;
};

void load_world_positions(struct world *, struct vertex *, const char *);

#endif
