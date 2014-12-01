/*
  Copyright (c) 2014 Dorle Osterode, Stefan Dang, Lukas Götz
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <stdio.h>
#include "core/types_api.h"
#include "core/error.h"
#include "core/str_api.h"
#include "genometools.h"

#ifndef GT_SCAFFOLDER_GRAPH_H
#define GT_SCAFFOLDER_GRAPH_H

typedef struct GtScaffolderGraph GtScaffolderGraph;

/* basic scaffold graph functions */
GtScaffolderGraph *gt_scaffolder_graph_new(GtUword nofvertices,
                                         GtUword nofedges);
void gt_scaffolder_graph_delete(GtScaffolderGraph *graph);

/* for test purposes: low level graph functions */
void gt_scaffolder_graph_add_vertex(GtScaffolderGraph *graph,
                                    const GtStr *header_seq,
                                    GtUword seqlen,
                                    float astat,
                                    float copynum);
void gt_scaffolder_graph_add_edge(GtScaffolderGraph *graph,
                                  GtUword vstartID,
                                  GtUword vendID,
                                  GtWord dist,
                                  float stddev,
                                  GtUword numpairs,
                                  bool dir,
                                  bool comp);
int gt_scaffolder_test_graph(GtUword max_nof_vertices,
                             GtUword max_nof_edges,
                             bool init_vertices,
                             GtUword nof_vertices,
                             bool init_edges,
                             GtUword nof_edges,
                             bool print_graph);


/* graph display functions
SK: sga Format unterstützen (asqg) */
int gt_scaffolder_graph_print(const GtScaffolderGraph *g,
			      const char *filename,
                              GtError *err);
void gt_scaffolder_graph_print_generic(const GtScaffolderGraph *g,
				       GtFile *f);

/* extended scaffold graph functions */
GtScaffolderGraph *gt_scaffolder_graph_new_from_file(const char *ctg_filename,
                                                   GtUword min_ctg_len,
                                                   const char *dist_filename,
                                                   GtError *err);
int gt_scaffolder_graph_mark_repeats(const char *filename,
                                     GtScaffolderGraph *graph,
                                     float copy_num_cutoff,
                                     float astat_cutoff,
                                     GtError *err);
int gt_scaffolder_graph_filter(GtScaffolderGraph *graph,
                               float pcutoff,
                               float cncutoff,
                               GtUword ocutoff);
/* SK: nicht-öffentlich? */
void gt_scaffolder_makescaffold(GtScaffolderGraph *graph);

#endif
