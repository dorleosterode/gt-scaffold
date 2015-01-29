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

#include "core/array_api.h"

#include "gt_scaffolder_graph.h"
#include "extended/assembly_stats_calculator.h"

#ifndef GT_SCAFFOLDER_ALGORITHMS_H
#define GT_SCAFFOLDER_ALGORITHMS_H

/* write scaffold into file */
int gt_scaffolder_graph_write_scaffold(GtArray *records,
				       const char *file_name,
				       GtError *err);

/* iterate graph to construct scaffold-records */
GtArray *gt_scaffolder_graph_iterate_scaffolds(const GtScaffolderGraph *graph,
                                        GtAssemblyStatsCalculator *scaf_stats);

/* load astatics and copy number of every contig and mark repeated contigs */
int gt_scaffolder_graph_mark_repeats(const char *filename,
                                     GtScaffolderGraph *graph,
                                     float copy_num_cutoff,
                                     float astat_cutoff,
                                     GtError *err);

/* mark polymorphic edges/vertices and inconsistent edges in scaffold graph */
void gt_scaffolder_graph_filter(GtScaffolderGraph *graph,
                                float pcutoff,
                                float cncutoff,
                                GtWord ocutoff);

/* check if vertex holds just sense or antisense edges */
bool
gt_scaffolder_graph_isterminal(const GtScaffolderGraphVertex *vertex);

/* create new walk */
GtScaffolderGraphWalk *gt_scaffolder_walk_new(void);

/* remove walk <*walk> */
void gt_scaffolder_walk_delete(GtScaffolderGraphWalk *walk);

/* add edge <*edge> to walk <*walk> */
void gt_scaffolder_walk_addegde(GtScaffolderGraphWalk *walk,
                                       GtScaffolderGraphEdge *edge);

GtScaffolderGraphRecord *
gt_scaffolder_graph_record_new(GtScaffolderGraphVertex *root);

void gt_scaffolder_graph_record_add_edge(GtScaffolderGraphRecord *rec,
					 GtScaffolderGraphEdge *edge);

void gt_scaffolder_graph_record_delete(GtScaffolderGraphRecord *rec);

void gt_scaffolder_calc_cc_and_terminals(const GtScaffolderGraph *graph,
                                         GtArray *ccs);

GtScaffolderGraphWalk *gt_scaffolder_create_walk(GtScaffolderGraph *graph,
                 GtScaffolderGraphVertex *start);

/* Konstruktion des Scaffolds mit groesster Contig-Gesamtlaenge */
void gt_scaffolder_makescaffold(GtScaffolderGraph *graph);

void gt_scaffolder_removecycles(GtScaffolderGraph *graph);

/* generates the fasta-string for the given scaffold-record */
GtStr *gt_scaffolder_graph_generate_string(GtScaffolderGraphRecord *rec,
                                           GtStr *ids);

#endif
