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



#include "core/error.h"
#include "core/file_api.h"
#include "core/str_api.h"
#include "core/types_api.h"

#ifndef GT_SCAFFOLDER_GRAPH_H
#define GT_SCAFFOLDER_GRAPH_H

typedef enum { GIS_UNVISITED, GIS_POLYMORPHIC, GIS_INCONSISTENT,
               GIS_REPEAT, GIS_VISITED, GIS_PROCESSED,
               GIS_SCAFFOLD} GraphItemState;

/* vertex of scaffold graph (describes one contig) */
typedef struct GtScaffolderGraphVertex {
  /* vertex index in array */
  GtUword index;
  /* header sequence of corresponding contig */
  GtStr *header_seq;
  /* sequence length of corresponding contig */
  GtUword seq_len;
  /* a-statistics value for classifying contig as repeat or unique contig */
  float astat;
  /* estimated copy number of corresponding contig */
  float copy_num;
  GtUword nof_edges;
  struct GtScaffolderGraphEdge **edges;
  /* vertex state (vertex can adapt every state except GIS_INCONSISTENT) */
  GraphItemState state;
} GtScaffolderGraphVertex;

/* edge of scaffold graph (describes orientation of two contigs) */
typedef struct GtScaffolderGraphEdge {
  /* pointer to end vertex of edge */
  GtScaffolderGraphVertex *end;
  /* pointer to start vertex of edge */
  GtScaffolderGraphVertex *start;
  /* estimated distance between contigs of start and end vertex */
  GtWord dist;
  /* standard deviation of estimated distance */
  float std_dev;
  /* number of read pairs resulting that distance */
  GtUword num_pairs;
  /* edge state */
  GraphItemState state;
  /* describes direction of corresponding contigs
     sense = true & same = true: ctg1 & ctg2 in sense direction
     sense = true & same = false: ctg1 in sense & ctg2 in antisense direction
     sense = false & same = true: ctg1 & ctg2 in antisense direction
     sense = false & same = false: ctg1 in antisense & ctg2 in sense direction*/
  bool sense;
  bool same;
} GtScaffolderGraphEdge;

/* scaffold graph */
typedef struct GtScaffolderGraph {
  GtScaffolderGraphVertex *vertices;
  GtUword nof_vertices;
  GtUword max_nof_vertices;
  GtScaffolderGraphEdge *edges;
  GtUword nof_edges;
  GtUword max_nof_edges;
} GtScaffolderGraph;

/* linear scaffold */
typedef struct GtScaffolderGraphWalk {
  GtUword nof_edges;
  GtUword size;
  GtUword total_contig_len;
  GtScaffolderGraphEdge **edges;
}GtScaffolderGraphWalk;

/* Construct graph data structure <*GtScaffolderGraph>. Wrap around two
   seperate constructor functions, which allocate memory for <max_nof_edges>
   edges and <max_nof_vertices> vertices. */
GtScaffolderGraph *gt_scaffolder_graph_new(GtUword max_nof_vertices,
                                           GtUword max_nof_edges);

/* Free all memory allocated for <*graph> including vertices and edges */
void gt_scaffolder_graph_delete(GtScaffolderGraph *graph);

/* Initialize a new vertex in <*graph>. Each vertex represents a contig and
   contains information about the sequence header <*header_seq>, sequence
   length <seq_len>, A-statistics <astat> and estimated copy number <copy_num>*/
void gt_scaffolder_graph_add_vertex(GtScaffolderGraph *graph,
                                    const GtStr *header_seq,
                                    GtUword seq_len,
                                    float astat,
                                    float copy_num);

/* Initialize a new edge in <*graph>. Each edge between two contig
   vertices <vstartID> and <vendID> contains information about the distance
   <dist>, standard deviation <std_dev>, number of pairs <num_pairs> and the
   direction of <vstartID> <dir> and corresponding <vendID> <same> */
void gt_scaffolder_graph_add_edge(GtScaffolderGraph *graph,
                                  GtUword vstartID,
                                  GtUword vendID,
                                  GtWord dist,
                                  float std_dev,
                                  GtUword num_pairs,
                                  bool dir,
                                  bool same);

GtScaffolderGraphEdge *gt_scaffolder_graph_find_edge(
                                    const GtScaffolderGraph *graph,
                                    GtUword vertexid_1,
                                    GtUword vertexid_2);

/* determines corresponding vertex id to contig header */
int gt_scaffolder_graph_get_vertex_id(const GtScaffolderGraph *graph,
                                      GtUword *vertex_id,
                                      const GtStr *header_seq);

/* assign edge <*edge> new attributes */
void gt_scaffolder_graph_alter_edge(GtScaffolderGraphEdge *edge,
                                    GtWord dist,
                                    float std_dev,
                                    GtUword num_pairs,
                                    bool sense,
                                    bool same);

/* print graphrepresentation in dot-format into file filename */
int gt_scaffolder_graph_print(const GtScaffolderGraph *g,
                              const char *filename,
                              GtError *err);

/* print graphrepresentation in dot-format into gt-filestream f */
void gt_scaffolder_graph_print_generic(const GtScaffolderGraph *g,
                                       GtFile *f);

/* create scaffold graph from file */
/* TODO: include a-statistics, copy number */
GtScaffolderGraph *gt_scaffolder_graph_new_from_file(const char *ctg_filename,
                                                     GtUword min_ctg_len,
                                                     const char *dist_filename,
                                                     GtError *err);

/* Function to test basic graph functionality on different scenarios:
- Create graph and allocate space for <max_nof_vertices> vertices and
  <max_nof_edges> edges.
- Initialize vertices if <init_vertices> is true and create <nof_vertices>
  vertices.
- Initialize edges if <init_edges> is true and create <nof_edges> edges.
- Delete graph. */
int gt_scaffolder_graph_test(GtUword max_nof_vertices,
                             GtUword max_nof_edges,
                             bool init_vertices,
                             GtUword nof_vertices,
                             bool init_edges,
                             GtUword nof_edges,
                             bool print_graph);

#endif
