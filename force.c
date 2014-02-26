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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/mman.h>
#include "world.h"
#include "worker.h"

#define COOLING 0.995
#define REPULSION_CAP_CHANGE 1.15

struct barycenter {
  long double x, y;
};

static double count_energy(struct world *, struct pair *, struct pair *, int);

void init_force(struct world *world)
{
  int start = 0;
  long pagesize = sysconf(_SC_PAGESIZE);
  pagesize=512;
  int items_per_unit =(pagesize-sizeof(struct world_work))/sizeof(struct pair);
  int nbufs = world->nitems/items_per_unit+1;
  struct world_work **work = malloc((nbufs+1)*sizeof(struct world_work **));
  work[nbufs] = NULL;
  for (struct world_work **workptr = work; start < world->nitems; ++workptr) {
    struct world_work *buf = mmap(NULL, pagesize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    buf->extra = malloc(sizeof(struct barycenter));
    buf->start = start;
    start += items_per_unit;
    buf->end = start >= world->nitems ? world->nitems : start;
    *workptr = buf;
  }
  world->world_work = work;
}

void work_map(void *cfg, void *data)
{
  struct world_work *work = data;
  struct world *world = cfg;
  struct barycenter *barycenter = work->extra;
  work->energy = 0;
  barycenter->x = barycenter->y = 0;
  struct pair *newpos = work->data;
  for (int i = work->start; i < work->end; ++i, ++newpos) {
    if (world->vertices[i].weight <= 0)
      continue;
    work->energy += count_energy(world, &world->vertices[i].pos, newpos, i);
    float weight = world->vertices[i].weight;
    barycenter->x += newpos->x*weight;
    barycenter->y += newpos->y*weight;
  }
}

struct copy_data {
  struct barycenter barycenter;
  struct world *world;
};

void work_copy(void *cfg, void *data)
{
  struct world_work *work = data;
  struct copy_data *copy_data = cfg;
  struct pair *newpos = work->data;
  struct vertex *tgt = copy_data->world->vertices;
  for (int i = work->start; i < work->end; ++i, ++newpos) {
    struct pair pos = *newpos;
    pos.x -= copy_data->barycenter.x;
    pos.y -= copy_data->barycenter.y;
    tgt[i].pos = pos;
  }
}

double world_step(struct world *world)
{
  double energy = 0;
  struct work_phase work_ops = {
    .work = &work_map
  };
  give_work(world->pool, &work_ops, world, world->world_work);
  struct world_work **workptr = world->world_work;
  struct barycenter barycenter = {0, 0};
  do {
    struct world_work *work = *workptr;
    barycenter.x += ((struct barycenter *)work->extra)->x;
    barycenter.y += ((struct barycenter *)work->extra)->y;
    energy += work->energy;
  } while (*(++workptr));
  barycenter.x *= world->world_weight_inv;
  barycenter.y *= world->world_weight_inv;
#if 1
  workptr = world->world_work;
  do {
    struct world_work *work = *workptr;
    struct pair *newpos = work->data;
    for (int i = work->start; i < work->end; ++i, ++newpos) {
      struct pair pos = *newpos;
      pos.x -= barycenter.x;
      pos.y -= barycenter.y;
      world->vertices[i].pos = pos;
    }
  } while (*(++workptr));
#else
  work_ops.work = &work_copy;
  struct copy_data *copy_data = malloc(sizeof(struct copy_data));
  copy_data->barycenter = barycenter;
  copy_data->world = world;
  give_work(world->pool, &work_ops, copy_data, world->world_work);
  free(copy_data);
#endif
  world->maxmove *= COOLING;
  world->repulsioncap *= REPULSION_CAP_CHANGE;
  return energy;
}

static double count_energy(struct world *world, struct pair *pos, struct pair *newpos, int i) {
  double maxmove = world->maxmove;
  double totalenergy = 0;
  struct vertex *v1 = &world->vertices[i], *v2 = world->vertices;
  struct edge *edge = world->edges+i*world->nitems;
  struct pair force = {0};
  for (int j = 0; j < world->nitems; ++j, ++v2, ++edge) {
    double energy = 0;
    if (i == j || world->vertices[j].weight < 0)
      continue;

    double relax = v1->radius+v2->radius+RELAX_EXTRA;
    double dist = hypot(pos->x-v2->pos.x, pos->y-v2->pos.y);

    double repulsionenergy;
    if (edge->weight > 0) {
      energy = edge->weight / v2->weight * pow(dist - relax, 2) / (v2->weight + relax);
      if (dist < relax)
	energy = -energy;
    }

    repulsionenergy = pow(v2->weight+relax, 2)/dist;
    double cap = world->repulsioncap*v2->weight;
    if (repulsionenergy > cap)
      repulsionenergy = cap;
    repulsionenergy -= 0.01;
    energy -= repulsionenergy;

    double normX = (v2->pos.x-pos->x)/dist, normY = (v2->pos.y-pos->y)/dist;
    energy /= v1->weight;
    force.x += energy*normX;
    force.y += energy*normY;
  }

  double energy = hypot(force.x, force.y);
  if (energy > world->maxmove) {
    double scale = world->maxmove/energy;
    force.x *= scale;
    force.y *= scale;
  }

  newpos->x = pos->x+force.x;
  newpos->y = pos->y+force.y;

  return energy;
}
