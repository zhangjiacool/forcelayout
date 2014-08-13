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

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "adjust.h"
#include "worker.h"

#define COMPARE_STEPS 4000

struct compare {
  double rotate;
  int mirror, i;
  double badness;
};

struct compare_data {
  struct world *world;
  struct vertex *compare_vertices;
  struct compare *array, **work;
};

/*
  The force-directed layout algorithm is fully deterministic and
  similar inputs would produce similar layouts, but the result might
  be mirrored and rotated compared to a run with slightly different
  input.  This step rotates and mirrors a result to find a most
  similar configuration to a given result.
*/

void *compare_initer(void *ptr)
{
  struct compare_data *compare_data = malloc(sizeof(struct compare_data));
  compare_data->array = malloc(COMPARE_STEPS*sizeof(struct compare));
  compare_data->work = malloc((COMPARE_STEPS+1)*sizeof(struct compare *));
  compare_data->work[COMPARE_STEPS] = NULL;
  for (int i = 0; i < COMPARE_STEPS; ++i) {
    compare_data->work[i] = &compare_data->array[i];
    if (i < COMPARE_STEPS/2) {
      struct compare values = {
	.rotate = i*M_PI*2/(COMPARE_STEPS/2),
	.mirror = 0,
	.i = i
      };
      compare_data->array[i] = values;
    } else {
      compare_data->array[i] = compare_data->array[i-COMPARE_STEPS/2];
      compare_data->array[i].mirror = 1;
    }
  }
  struct compare_init *init = ptr;
  struct world *world = init->world;
  compare_data->compare_vertices = calloc(world->nitems, sizeof(struct vertex));
  load_world_positions(world, compare_data->compare_vertices, init->filepath);
  compare_data->world = world;
  return compare_data;
}

struct translate_matrix {
  double a, b, c, d;
};

static void translate_world(struct world *, const struct compare *);

static void make_translate(struct translate_matrix *m, const struct compare *compare) {
  double sin_val, cos_val;
#ifdef _GNU_SOURCE
  sincos(compare->rotate, &sin_val, &cos_val);
#else
  sin_val = sin(compare->rotate);
  cos_val = cos(compare->rotate);
#endif
  if (!compare->mirror) {
    m->a = cos_val;
    m->b = -sin_val;
    m->c = sin_val;
    m->d = cos_val;
  } else {
    m->a = -sin_val;
    m->b = cos_val;
    m->c = cos_val;
    m->d = sin_val;
  }
}


static void work_adjust(void *cfg, void *data)
{
  struct compare_data *compare_data = cfg;
  struct compare *compare = data;
  struct world *world = compare_data->world;
  // Start by filling out the matrix
  struct translate_matrix m;
  make_translate(&m, compare);
  double badness = 0;
  struct vertex *v1 = world->vertices, *v2 = compare_data->compare_vertices;
  for (int i = 0; i < world->nitems; ++i, ++v1, ++v2) {
    if (v2->weight == 0 || isinf(v1->weight) || isinf(v2->weight))
      continue;
    double x = v1->pos.x*m.a+v1->pos.y*m.b;
    double y = v1->pos.x*m.c+v1->pos.y*m.d;
    badness += hypot(x-v2->pos.x,y-v2->pos.y)*v1->weight;
  }
  compare->badness = badness;
}

void compare_world(struct compare_data *compare)
{
  // First, try out a full circle and a mirrored setup.
  struct work_phase compare_ops = {
    .work = &work_adjust
  };
  give_work(compare->world->pool, &compare_ops, compare, compare->work);
  struct compare best_compare;
  best_compare.badness = DBL_MAX;
  for (int i = 0; i < COMPARE_STEPS; ++i) {
    if (compare->array[i].badness < best_compare.badness) {
      best_compare = compare->array[i];
    }
  }
  
  // Then, fine tune by scanning around the best compare
  double tune_base = best_compare.i*M_PI*2/(COMPARE_STEPS/2);
  int mirror = best_compare.mirror;
  for (int i = -COMPARE_STEPS/2; i < COMPARE_STEPS/2; ++i) {
    struct compare value = {
      .rotate = tune_base+i*M_PI*2/(COMPARE_STEPS*COMPARE_STEPS/4),
      .mirror = mirror
    };
    compare->array[i+COMPARE_STEPS/2] = value;
  }
  give_work(compare->world->pool, &compare_ops, compare, compare->work);
  for (int i = 0; i < COMPARE_STEPS; ++i) {
    if (compare->array[i].badness < best_compare.badness) {
      best_compare = compare->array[i];
    }
  }

  translate_world(compare->world, &best_compare);
}

static void translate_world(struct world *world, const struct compare *compare) {
  struct translate_matrix m;
  make_translate(&m, compare);
  struct vertex *v = world->vertices;
  for (int i = 0; i < world->nitems; ++i, ++v) {
    double x = v->pos.x*m.a+v->pos.y*m.b;
    double y = v->pos.x*m.c+v->pos.y*m.d;
    v->pos.x = x;
    v->pos.y = y;
  }
}
