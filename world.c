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

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jansson.h>
#include <assert.h>
#include <sys/mman.h>
#include <math.h>
#include <pthread.h>
#include <dirent.h>

#define ITERATIONS 1000
#define MAXTHREADS 16

#include "world.h"
#include "force.h"
#include "adjust.h"
#include "worker.h"
#include "sparsify.h"

struct options {
  int threads;
  int verbose;
  const char *output;
  const char *initial_positions;
  int iterations;
  const char *rotate_to;
};

static void inc_weight(struct world *world, int i, int j)
{
  struct edge *edge = world->edges + i + j*world->nitems;
  ++edge->weight;
  edge = world->edges + j + i*world->nitems;
  ++edge->weight;
}

static void count_edge_closure(int nitems, struct edge *edges, char *closure_map, int i)
{
  struct edge *edge = edges+i*nitems;
  closure_map[i] = 1;
  for (int j = 0; j < nitems; ++j, ++edge) {
    if (i == j)
      continue;
    if (edge->weight > 0 && !closure_map[j]) {
      count_edge_closure(nitems, edges, closure_map, j);
    }
  }
}

static int key_comparator(const void *k1, const void *k2)
{
  return strcmp(* (char * const *) k1, * (char * const *) k2);
}

void init_world(struct world *world, json_t *json) {
  json_t *items = json_object_get(json, "items"), *picks = json_object_get(json, "picks"), *c;
  const char *key;
  struct vertex *ptr;
  int i = 0;
  int mult = 0;
  int heaviestitem, maxweight = 0;
  int nthreads;
  const char **keys;

  assert(items);
  if (world->options->threads <= 0)
    nthreads = sysconf(_SC_NPROCESSORS_ONLN);
  else if (world->options->threads > MAXTHREADS)
    nthreads = MAXTHREADS;
  else
    nthreads = world->options->threads;

  world->pool = init_workers(nthreads);
  world->maxmove = 30;
  world->repulsioncap = 10;
  world->nitems = json_object_size(items);
  world->edges = calloc(world->nitems*world->nitems, sizeof(struct edge));
  world->mapping = malloc((1+world->nitems)*sizeof(int));
  ptr = world->vertices = malloc(world->nitems*sizeof(struct vertex));
  world->maxid = 0;
  keys = malloc((world->nitems+1)*sizeof(char *));
  keys[world->nitems] = NULL;
  json_object_foreach (items, key, c) {
    keys[i++] = key;
    int id = atoi(key);
    world->maxid = id > world->maxid ? id : world->maxid;
  }
  world->r_mapping = malloc((1+world->maxid)*sizeof(int));
  qsort(&keys[0], world->nitems, sizeof(char *), key_comparator);

  for (i = 0; i < world->nitems; ++i) {
    json_t *val = json_object_get(items, keys[i]);
    int id = atoi(keys[i]);
    ptr->weight = 1+json_integer_value(json_object_get(val, "weight"));
    ptr->radius = sqrtf(ptr->weight)/M_PI;
    if (ptr->weight > maxweight) {
      heaviestitem = id;
      maxweight = ptr->weight;
    }

    // Initial placement: concentric rings around 0,0
    ptr->pos.x = (mult+8)*10*sin((double)i/world->nitems*2*M_PI);
    ptr->pos.y = (mult+8)*10*cos((double)i/world->nitems*2*M_PI);
    ++mult;
    mult %= 16;
    world->r_mapping[id] = i+1;
    world->mapping[i+1] = id;
    ++ptr;
  }
  if (world->options->initial_positions) {
    load_world_positions(world, world->vertices, world->options->initial_positions);
  }

  json_object_foreach (picks, key, c) {
    for (i = 0; i < json_array_size(c); ++i) {
      int ref1 = world->r_mapping[json_integer_value(json_array_get(c, i))];
      if (ref1) {
	for (int j = i+1; j < json_array_size(c); ++j) {
	  int ref2 = world->r_mapping[json_integer_value(json_array_get(c, j))];
	  if (ref2) {
	    inc_weight(world, ref1-1, ref2-1);
	  }
	}
      }
    }
  }

  // Sanity check: only pick items which are connected to the heaviest item
  char *closure_map = calloc(world->nitems+1, sizeof(char));
  count_edge_closure(world->nitems, world->edges, closure_map, world->r_mapping[heaviestitem]-1);
  for (int i = 0; i < world->nitems; ++i) {
    if (!closure_map[i]) {
      world->vertices[i].weight = -INFINITY;
    }
  }
  free(closure_map);

  world->world_weight_inv = 0;
  for (i = 0; i < world->nitems; ++i) {
    if (isfinite(world->vertices[i].weight)) {
      world->world_weight_inv += world->vertices[i].weight;
    }
  }
  world->world_weight_inv = 1/world->world_weight_inv;
  init_force(world);
}

static json_t * world_to_json(struct world *world) {
  json_t *res = json_object();
  for (int i = 0; i < world->nitems; ++i) {
    char id[8];
    struct vertex *par = &world->vertices[i];
    if (par->weight <= 0)
      continue;
    json_t *comic = json_object();
    json_object_set_new(comic, "x", json_real(par->pos.x));
    json_object_set_new(comic, "y", json_real(par->pos.y));
    json_object_set_new(comic, "radius", json_real(par->radius));
    json_object_set_new(comic, "weight", json_integer(par->weight));
    snprintf(id, 8, "%i", world->mapping[i+1]);
    json_object_set_new(res, id, comic);
  }
  return res;
}

void load_world_positions(struct world *world, struct vertex *vertices, const char *path) {
  json_error_t error;
  json_t *pos, *positions = json_load_file(path, 0, &error);
  const char *key;
  if (!positions) {
    fprintf(stderr, "forcelayout: %s %s", error.text, error.source);
    exit(1);
  }
  json_object_foreach(positions, key, pos) {
    int keyval = atoi(key);
    if (keyval > world->maxid)
      continue;
    int rid = world->r_mapping[atoi(key)];
    if (rid != 0) {
      struct vertex *par = &vertices[rid-1];
      par->pos.x = json_real_value(json_object_get(pos, "x"));
      par->pos.y = json_real_value(json_object_get(pos, "y"));
      par->weight = json_integer_value(json_object_get(pos, "weight"));
      par->radius = sqrtf(par->weight)/M_PI;
    }
  }
  json_decref(positions);
}


static void usage() {
  fprintf(stderr, "usage: forcelayout [-j threads] [-i iterations] [-r reference] [-q] input.json output.json\n");
  exit(1);
}

int main(int argc, char *argv[]) {
  double energy;
  size_t len;
  json_t *json;
  struct world world;
  struct options options = {
    .threads = 0,
    .verbose = 1,
    .output = NULL,
    .initial_positions = NULL,
    .iterations = 0,
    .rotate_to = NULL
  };
  int opt;
  while ((opt = getopt(argc, argv, "j:p:i:qr:")) != -1) {
    switch (opt) {
    case 'j':
      options.threads = atoi(optarg);
      break;
    case 'q':
      options.verbose = 0;
      break;
    case 'p':
      options.initial_positions = optarg;
      break;
    case 'r':
      options.rotate_to = optarg;
      break;
    case 'i':
      options.iterations = atoi(optarg);
      break;
    default:
      usage();
    }
  }

  if (optind+2 != argc)
    usage();
  if (options.iterations <= 0)
    options.iterations = ITERATIONS;
  json = json_load_file(argv[optind], 0, NULL);
  options.output = argv[optind+1];
  assert(json_is_object(json));

  world.options = &options;
  init_world(&world, json);
  pthread_t rotate_loader_thread;
  struct compare_init compare_init;
  if (options.rotate_to) {
    compare_init.world = &world;
    compare_init.filepath = options.rotate_to;
    pthread_create(&rotate_loader_thread, NULL, compare_initer, &compare_init);
  }

  // Main force-directed graph algorithm
  for (int i = 0; i < options.iterations; ++i) {
#ifdef DEBUG
    char tmpname[100];
    snprintf(tmpname, 100, "/tmp/world%i.json", i);
    json_dump_file(world_to_json(&world), tmpname, JSON_INDENT(2));
#endif
    energy = world_step(&world);
    if (options.verbose)
      fprintf(stderr, "%i forces %f\n", i, energy);
  }

  sparsify_world(&world);
  do {
    energy = sparsify_step(&world);
    if (options.verbose)
      fprintf(stderr, "overlap %f\n", energy);
  } while (energy > 0);

  if (options.rotate_to) {
    void *retval;
    struct compare_data *compare_data;
    pthread_join(rotate_loader_thread, &retval);
    compare_data = retval;
    compare_world(compare_data);
  }

  json_dump_file(world_to_json(&world), options.output, JSON_INDENT(2));
}
