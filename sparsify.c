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

#include <config.h>

#include <math.h>
#include "world.h"
#include "worker.h"

static double resolve_overlap(struct world *world, struct pair *newpos, int i) {
  double overlap = 0;
  struct vertex *v1 = &world->vertices[i];
  if (v1->weight <= 0)
    return 0;
  struct pair force = {0};
  for (int j = 0; j < world->nitems; ++j) {
    struct vertex *v2 = &world->vertices[j];
    if (i == j || v2->weight <= 0)
      continue;
    
    float relax = v1->radius+v2->radius+RELAX_EXTRA/2;
    // hypot is expensive.
    double taxidist = fabs(v1->pos.x-v2->pos.x)+fabs(v1->pos.y-v2->pos.y);
    if (relax*2 < taxidist)
      continue;
    double dist = hypot(v1->pos.x-v2->pos.x, v1->pos.y-v2->pos.y);
    if (dist < relax) {
      overlap += relax-dist;
      double nudge = -(relax+RELAX_EXTRA-dist)/2;
      if (v1->weight > v2->weight) {
	nudge *= v2->weight/v1->weight;
      }
      double normX = (v2->pos.x-v1->pos.x)/dist, normY = (v2->pos.y-v1->pos.y)/dist;
      force.x += nudge*normX;
      force.y += nudge*normY;
    }
  }

  newpos->x = v1->pos.x+force.x;
  newpos->y = v1->pos.y+force.y;

  return overlap;
}

static double count_overlap(double dist, struct vertex *v1, struct vertex *v2)
{
  double overlap;
  overlap = dist-v1->radius;
  if (overlap < 0) {
    return 3;
  } else {
    overlap = (v2->radius+RELAX_EXTRA)/overlap;
    return overlap > 3 ? 3 : overlap;
  }
}

static void sparsify_work(void *cfg, void *data)
{
  struct world *world = cfg;
  struct world_work *work = data;
  work->energy = 0;
  struct pair *newpos = work->data;
  for (int i = work->start; i < work->end; ++i, ++newpos) {
    work->energy += resolve_overlap(world, newpos, i);
  }
}

// Make the map sparser to resolve overlaps
void sparsify_world(struct world *world)
{
  // First, shift everything by a constant factor
  double total_overlap = 0;
  int noverlap = 0;
  for (int i = 0; i < world->nitems; ++i) {
    struct vertex *v1 = &world->vertices[i], *v2;
    if (v1->weight <= 0)
      continue;
    for (int j = i+1; j < world->nitems; ++j) {
      v2 = &world->vertices[j];
      if (v2->weight <= 0)
	continue;
      double dist = hypot(v1->pos.x-v2->pos.x, v1->pos.y-v2->pos.y);
      double relax = v1->radius + v2->radius + RELAX_EXTRA;
      if (dist < relax) {
	noverlap += v1->weight+v2->weight;
	total_overlap += v2->weight*count_overlap(dist, v1, v2);
	total_overlap += v1->weight*count_overlap(dist, v2, v1);
      }
    }
  }
  total_overlap /= noverlap;

  for (int i = 0; i < world->nitems; ++i) {
    world->vertices[i].pos.x *= total_overlap;
    world->vertices[i].pos.y *= total_overlap;
  }
}

  // Then, bump vertices around until overlaps are resolved
double sparsify_step(struct world *world)
{
  struct work_phase work_ops = {
    .work = &sparsify_work
  };
  double energy = 0;
  give_work(world->pool, &work_ops, world, world->world_work);
  struct world_work **workptr = world->world_work;
  do {
    struct world_work *work = *workptr;
    struct pair *newpos = work->data;
    for (int i = work->start; i < work->end; ++i, ++newpos) {
      world->vertices[i].pos = *newpos;
    }
    energy += work->energy;
  } while (*(++workptr));
      
  return energy;
}
