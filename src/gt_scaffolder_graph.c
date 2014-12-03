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
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <core/assert_api.h>
#include "core/queue_api.h"
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/error.h"
#include "gt_scaffolder_graph.h"
#include "gt_scaffolder_parser.h"

/* Initialize vertex portion inside <*graph>. Allocate memory for
   <max_nof_vertices> vertices. */
void gt_scaffolder_graph_init_vertices(GtScaffolderGraph *graph,
                                       GtUword max_nof_vertices)
{
  gt_assert(graph != NULL);
  gt_assert(graph->vertices == NULL);
  gt_assert(max_nof_vertices > 0);
  graph->vertices = gt_malloc(sizeof (*graph->vertices) * max_nof_vertices);
  graph->nof_vertices = 0;
  graph->max_nof_vertices = max_nof_vertices;
}

/* Initialize edge portion inside <*graph>. Allocate memory for
   <max_nof_edges> edges. */
void gt_scaffolder_graph_init_edges(GtScaffolderGraph *graph,
                                    GtUword max_nof_edges)
{
  gt_assert(graph != NULL);
  gt_assert(graph->edges == NULL);
  gt_assert(max_nof_edges > 0);
  graph->edges = gt_malloc(sizeof (*graph->edges) * max_nof_edges);
  graph->nof_edges = 0;
  graph->max_nof_edges = max_nof_edges;
}

/* Construct graph data structure <*GtScaffolderGraph>. Wrap around two
   seperate constructor functions, which allocate memory for <max_nof_edges>
   edges and <max_nof_vertices> vertices. */
GtScaffolderGraph *gt_scaffolder_graph_new(GtUword max_nof_vertices,
                                           GtUword max_nof_edges)
{
  GtScaffolderGraph *graph;

  graph = gt_malloc(sizeof (*graph));
  graph->vertices = NULL;
  graph->edges = NULL;
  gt_scaffolder_graph_init_vertices(graph, max_nof_vertices);
  gt_scaffolder_graph_init_edges(graph, max_nof_edges);

  return graph;
}

/* Free all memory allocated for <*graph> including vertices and edges */
void gt_scaffolder_graph_delete(GtScaffolderGraph *graph)
{
  GtScaffolderGraphVertex *vertex;
  /*LG: NULL pointer has to be freed !?
    gt_assert(graph != NULL);*/

  if (graph != NULL) {
    if (graph->vertices != NULL) {
      /* If existent, free header_seq and pointer to outgoing edges first */
      for ( vertex = graph->vertices;
            vertex < (graph->vertices + graph->nof_vertices);
            vertex++
          )
      {
        gt_str_delete(vertex->header_seq);
        gt_free(vertex->edges);
      }
    }

    gt_free(graph->vertices);
    gt_free(graph->edges);
  }

  gt_free(graph);
}

/* Initialize a new vertex in <*graph>. Each vertex represents a contig and
   contains information about the sequence header <*header_seq>, sequence
   length <seq_len>, A-statistics <astat> and estimated copy number <copy_num>*/
void gt_scaffolder_graph_add_vertex(GtScaffolderGraph *graph,
                                    const GtStr *header_seq,
                                    GtUword seq_len,
                                    float astat,
                                    float copy_num)
{
  GtUword nextfree;

  gt_assert(graph != NULL);
  gt_assert(graph->nof_vertices < graph->max_nof_vertices);

  nextfree = graph->nof_vertices;

  /* Initialize vertex */
  graph->vertices[nextfree].index = nextfree; /* SD: remove without breaking */
  graph->vertices[nextfree].seq_len = seq_len;
  graph->vertices[nextfree].astat = astat;
  graph->vertices[nextfree].copy_num = copy_num;
  graph->vertices[nextfree].nof_edges = 0;
  if (header_seq != NULL) {
    graph->vertices[nextfree].header_seq = gt_str_clone(header_seq);
  }
  graph->vertices[nextfree].state = GIS_UNVISITED;

  /* Allocate initial space for pointer to outgoing edges */
  graph->vertices[nextfree].edges = gt_malloc(sizeof (*graph->vertices->edges));

  graph->nof_vertices++;
}

void gt_scaffolder_graph_add_edge_ptr_to_vertex(GtScaffolderGraph *graph,
                                                GtUword edgeID,
                                                GtUword vertexID)
{
  /* Allocate new space for pointer to this edge */
  if (graph->vertices[vertexID].nof_edges > 0) {
    graph->vertices[vertexID].edges =
      /* SK: realloc zu teuer? Besser: DistEst parsen und gezielt allokieren */
      gt_realloc( graph->vertices[vertexID].edges, sizeof (*graph->edges) *
                  (graph->vertices[vertexID].nof_edges + 1) );
  }
  /* Assign adress of this edge to the pointer */
  graph->vertices[vertexID].edges[graph->vertices[vertexID].nof_edges] =
    &graph->edges[edgeID];

  graph->vertices[vertexID].nof_edges++;
}

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
                                  bool same)
{

  gt_assert(graph != NULL);
  gt_assert(graph->nof_edges < graph->max_nof_edges);

  GtUword nextfree = graph->nof_edges;

  /* Inititalize edge */
  graph->edges[nextfree].start = graph->vertices + vstartID;
  graph->edges[nextfree].end = graph->vertices + vendID;
  graph->edges[nextfree].dist = dist;
  graph->edges[nextfree].std_dev = std_dev;
  graph->edges[nextfree].num_pairs = num_pairs;
  graph->edges[nextfree].sense = dir;
  graph->edges[nextfree].same = same;
  graph->edges[nextfree].state = GIS_UNVISITED;

  /* Add ptr to edge to start and end vertex */
  gt_scaffolder_graph_add_edge_ptr_to_vertex(graph, nextfree, vstartID);
  gt_scaffolder_graph_add_edge_ptr_to_vertex(graph, nextfree, vendID);

  graph->nof_edges++;
}

GtScaffolderGraphEdge
*gt_scaffolder_graph_find_edge(const GtScaffolderGraph *graph,
                               GtUword vertexid_1,
                               GtUword vertexid_2)
{
  GtScaffolderGraphEdge **edge;
  GtScaffolderGraphVertex *v1 = graph->vertices + vertexid_1;

  for (edge = v1->edges; edge < (v1->edges + v1->nof_edges); edge++) {
    if ((*edge)->end->index == vertexid_2)
      return *edge;
  }
  return NULL;
}

/* determines corresponding vertex id to contig header */
/* SD: Binaersuche separat testen */
int gt_scaffolder_graph_get_vertex_id(const GtScaffolderGraph *graph,
                                      GtUword *vertex_id,
                                      const GtStr *header_seq)
{
  GtScaffolderGraphVertex *min_vertex, *max_vertex, *mid_vertex;
  int had_err, cmp;
  bool found;

  had_err = 0;
  found = false;

  /* binary search */
  min_vertex = graph->vertices;
  max_vertex = graph->vertices + graph->nof_vertices - 1;
  while (max_vertex >= min_vertex)
    {
      /* calculate midpoint */
      mid_vertex = min_vertex + ((max_vertex - min_vertex) / 2);
      cmp = gt_str_cmp(mid_vertex->header_seq, header_seq);
      if (cmp == 0)
      {
        found = true;
        *vertex_id = mid_vertex->index;
        break;
      }
      else if (cmp < 0)
        min_vertex = mid_vertex + 1;
      else
        max_vertex = mid_vertex - 1;
    }

  /* contig header was not found */
  if (found == false)
    had_err = -1;

  return had_err;
}

/* assign edge <*edge> new attributes */
void gt_scaffolder_graph_alter_edge(GtScaffolderGraphEdge *edge,
                                    GtWord dist,
                                    float std_dev,
                                    GtUword num_pairs,
                                    bool sense,
                                    bool same)
{
  /* check if edge exists */
  gt_assert(edge != NULL);

  /* assign edge new attributes */
  edge->dist = dist;
  edge->std_dev = std_dev;
  edge->num_pairs = num_pairs;
  edge->sense = sense;
  edge->same = same;
}

/* print graphrepresentation in dot-format into file filename */
int gt_scaffolder_graph_print(const GtScaffolderGraph *g,
                              const char *filename,
                              GtError *err)
{
  gt_assert(g != NULL);

  int had_err = 0;

  GtFile *f = gt_file_new(filename, "w", err);
  if (f == NULL)
    had_err = -1;

  if (!had_err) {
    gt_scaffolder_graph_print_generic(g, f);
    gt_file_delete(f);
  }

  return had_err;
}

/* print graphrepresentation in dot-format into gt-filestream f */
void gt_scaffolder_graph_print_generic(const GtScaffolderGraph *g,
                                       GtFile *f)
{
  gt_assert(g != NULL);

  GtScaffolderGraphVertex *v;
  GtScaffolderGraphEdge *e;
  /* 0: GIS_UNVISITED, 1: GIS_POLYMORPHIC, 2: GIS_INCONSISTENT,
     3: GIS_VISITED, 4: GIS_PROCESSED */
  const char *color_array[] = {"black", "gray", "gray", "red", "green"};

  /* print first line into f */
  gt_file_xprintf(f, "graph {\n");

  /* iterate over all vertices and print them. add attribute color according
     to the current state */
  for (v = g->vertices; v < (g->vertices + g->nof_vertices); v++) {
    gt_file_xprintf(f, GT_WU " [color=\"%s\" label=\"%s\"];\n", v->index,
                    color_array[v->state], gt_str_get(v->header_seq));
  }

  /* iterate over all edges and print them. add attribute color according to
     the current state and label the edge with the distance*/
  for (e = g->edges; e < (g->edges + g->nof_edges); e++) {
    gt_file_xprintf(f,
                    GT_WU " -- " GT_WU " [color=\"%s\" label=\"" GT_WD "\"];\n",
                    e->start->index, e->end->index,
                    color_array[e->state], e->dist);
  }

  /* print the last line into f */
  gt_file_xprintf(f, "}\n");
}

/* create scaffold graph from file */
/* TODO: include a-statistics, copy number */
GtScaffolderGraph *gt_scaffolder_graph_new_from_file(const char *ctg_filename,
                                                     GtUword min_ctg_len,
                                                     const char *dist_filename,
                                                     GtError *err)
{
  GtScaffolderGraph *graph;
  int had_err;
  GtUword nof_distances, nof_contigs;

  graph = NULL;
  had_err = 0;
  nof_contigs = 0;
  nof_distances = 0;

  /* count contigs */
  had_err = gt_scaffolder_parser_count_contigs(ctg_filename, min_ctg_len,
            &nof_contigs, err);

  if (had_err == 0)
  {
    /* allocate memory for vertices of scaffolder graph */
    graph = gt_malloc(sizeof (*graph));
    graph->vertices = NULL;
    graph->edges = NULL;
    gt_scaffolder_graph_init_vertices(graph, nof_contigs);

    /* parse contigs in FASTA-format and save them as vertices of
       scaffold graph */
    gt_scaffolder_parser_read_contigs(graph, ctg_filename, min_ctg_len, err);

    if (had_err == 0)
    {
      /* count distance information */
      nof_distances = 0;
      had_err = gt_scaffolder_parser_count_distances(graph, dist_filename,
              &nof_distances, err);

      if (had_err == 0)
      {
        /* allocate memory for edges of scaffolder graph */
        gt_scaffolder_graph_init_edges(graph, nof_distances);
        /* parse distance information of contigs in abyss-dist-format and
           save them as edges of scaffold graph */
        had_err =
          gt_scaffolder_parser_read_distances(dist_filename, graph, false, err);
      }
    }
  }
  /* SK: loeschen: gt_error_check(err);*/
  /* SK: graph / callback loeschen und auf NULL setzen */

  if (had_err != 0)
  {
    gt_scaffolder_graph_delete(graph);
    graph = NULL;
  }

  return graph;
}

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
                             bool print_graph)
{
  int had_err = 0;
  GtScaffolderGraph *graph;

  /* Construct graph and init vertices / edges */
  if (!init_vertices && !init_edges )
    graph = gt_scaffolder_graph_new(max_nof_vertices, max_nof_edges);
  /* Construct graph, don't init as this will be done later */
  else {
    graph = gt_malloc(sizeof (*graph));
    graph->vertices = NULL;
    graph->edges = NULL;
  }

  if (graph == NULL)
    had_err = -1;

  /* Init vertex portion of graph. <nof_vertices> Create vertices. */
  if (init_vertices) {
    gt_scaffolder_graph_init_vertices(graph, max_nof_vertices);
    if (graph->vertices == NULL)
      had_err = -1;

    unsigned i;
    for (i = 0; i < nof_vertices; i++) {
      gt_scaffolder_graph_add_vertex(graph, gt_str_new_cstr("foobar"), 100, 20, 40);
    }
  }

  /* Init edge portion of graph. Connect every vertex with another vertex until
  <nof_edges> is reached. */
  if (init_edges) {
    unsigned vertex1 = 0, vertex2 = 0;

    gt_scaffolder_graph_init_edges(graph, max_nof_edges);

    if (graph->edges == NULL)
      had_err = -1;

    /* Connect 1st vertex with every other vertex, then 2nd one, etc */
    unsigned i;
    for (i = 0; i < nof_edges; i++) {
      if (vertex2 < nof_vertices - 1)
        vertex2++;
      else if (vertex1 < nof_vertices - 2) {
        vertex1++;
        vertex2 = vertex1 + 1;
      }
      gt_scaffolder_graph_add_edge(graph, vertex1, vertex2, 2, 1.5, 4, true,
                                   true);
    }
  }

  /* Print the graph for diff comparison */
  /* SD: Ask Dorle about error object and (!ma) assertion */
  if (print_graph) {
    GtError *err;
    char outfile[] = "gt_scaffolder_graph_test.dot";
    err = gt_error_new();
    gt_scaffolder_graph_print(graph, outfile, err);
    gt_error_delete(err);
  }

  gt_scaffolder_graph_delete(graph);
  if (graph != NULL)
    had_err = -1;

  return had_err;
}
