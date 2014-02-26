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

#include "worker.h"
#include <stdlib.h>
#include <pthread.h>
#include <sched.h>
#include <stdio.h>

struct thread_control {
  int nthreads, nthreads_working;
  pthread_mutex_t mutex;
  pthread_cond_t work_available, work_taken, work_done;
  pthread_t *threads;
  void *worker_arg, *job_item;
  struct work_phase work;
};

struct job_data {
  int dummy;
};

struct worker_data {
  int i;
  struct thread_control *control;
};

static void *worker(void *ptr)
{
  struct worker_data *data = ptr;
  struct thread_control *control = data->control;
  void *job_item;
  pthread_mutex_lock(&control->mutex);
  for (;;) {
    while (!control->job_item)
      pthread_cond_wait(&control->work_available, &control->mutex);
    job_item = control->job_item;
    control->job_item = NULL;
    if (control->work.init_phase)
      control->work.init_phase(control->worker_arg, job_item);
    ++control->nthreads_working;
    pthread_cond_signal(&control->work_taken);
    pthread_mutex_unlock(&control->mutex);
    control->work.work(control->worker_arg, job_item);
    //printf("thread %i done work\n", data->i);
    pthread_mutex_lock(&control->mutex);
    if (control->work.end_phase)
      control->work.end_phase(control->worker_arg, job_item);
    --control->nthreads_working;
    pthread_cond_signal(&control->work_done);
  }
}

struct thread_control *init_workers(int nthreads)
{
  struct thread_control *control = malloc(sizeof(struct thread_control));
  control->job_item = NULL;
  control->nthreads = nthreads;
  control->nthreads_working = 0;
  control->threads = malloc(nthreads*sizeof(pthread_t *));
  pthread_cond_init(&control->work_available, NULL);
  pthread_cond_init(&control->work_taken, NULL);
  pthread_cond_init(&control->work_done, NULL);
  pthread_mutex_init(&control->mutex, NULL);
  pthread_attr_t init_attr;
  pthread_attr_init(&init_attr);
#ifdef SCHED_BATCH
  pthread_attr_setschedpolicy(&init_attr, SCHED_BATCH);
#endif
  for (int t = 0; t < nthreads; ++t) {
    struct worker_data *worker_data = malloc(sizeof(struct worker_data));
    worker_data->i = t;
    worker_data->control = control;
    pthread_create(&control->threads[t], &init_attr, worker, worker_data);
  }
  return control;
}

void give_work(struct thread_control *control, struct work_phase *phase, void *arg, void *work)
{
  pthread_mutex_lock(&control->mutex);
  control->work = *phase;
  control->worker_arg = arg;
  for (struct job_data **ptr = work; *ptr; ++ptr) {
    while (control->job_item != NULL)
      pthread_cond_wait(&control->work_taken, &control->mutex);
    control->job_item = *ptr;
    pthread_cond_signal(&control->work_available);
  }
  while (control->nthreads_working > 0)
    pthread_cond_wait(&control->work_done, &control->mutex);
  pthread_mutex_unlock(&control->mutex);
}
