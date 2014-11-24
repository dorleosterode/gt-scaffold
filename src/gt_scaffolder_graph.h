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

#ifndef GT_SCAFFOLDER_GRAPH_H
#define GT_SCAFFOLDER_GRAPH_H

typedef struct GtScaffoldGraph GtScaffoldGraph;

/* basic scaffold graph functions */
GtScaffoldGraph *gt_scaffolder_graph_new(GtUword nofvertices, GtUword nofedges);
void graph_add_vertex(GtScaffoldGraph *graph, GtUword seqlen, float astat,
  float copynum);
void graph_add_edge(GtScaffoldGraph *graph, GtUword vstartID, GtUword vendID,
  GtWord dist, float stddev, GtUword numpairs, bool dir, bool comp);
void gt_scaffolder_graph_delete(GtScaffoldGraph *graph);

/* Darstellung des Graphen
SK: gt_file benutzen, print_graph_generic / print_graph
SK: sga Format unterstützen asgq (?) */
int write_graph(const GtScaffoldGraph *g, const char *filename);
void print_graph(const GtScaffoldGraph *g, FILE *f);

/* extended scaffold graph functions
SK: Fasta-Iterator (core/fastareaderseqiterator) */
GtScaffoldGraph *gt_scaffolder_graph_new_from_file(const char *ctgfilename,
                 GtUword minctglen, const char *distfilename, GtError *err);
int gt_scaffolder_graph_filtering(GtScaffoldGraph *graph, float pcutoff,
    float cncutoff, GtUword ocutoff); /* SK: GtError Objekt; nicht benötigt */
void gt_scaffolder_makescaffold(GtScaffoldGraph *graph); /* SK: nicht-öffentlich? */


#endif
