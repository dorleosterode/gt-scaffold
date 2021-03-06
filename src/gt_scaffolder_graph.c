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

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/file_api.h"
#include "core/ma_api.h"
#include "core/str_api.h"

#include "gt_scaffolder_graph.h"
#include "gt_scaffolder_parser.h"

/* Initialize vertex portion inside <*graph>. Allocate memory for
   <max_nof_vertices> vertices. */
static void gt_scaffolder_graph_init_vertices(GtScaffolderGraph *graph,
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
static void gt_scaffolder_graph_init_edges(GtScaffolderGraph *graph,
                                           GtUword max_nof_edges)
{
  gt_assert(graph != NULL);
  gt_assert(graph->edges == NULL);
  gt_assert(max_nof_edges > 0);
  graph->edges = gt_malloc(sizeof (*graph->edges) * max_nof_edges);
  graph->nof_edges = 0;
  graph->max_nof_edges = max_nof_edges;
}

/* Construct graph data structure <*GtScaffolderGraph>. Initialize edges and
   allocate memory for <max_nof_edges> if != 0. Initialize vertices and allocate
   memory for <max_nof_vertices> vertices if != 0. Edges and/or vertices will
   have to be initialized seperately otherwise. */
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

  if (graph != NULL) {

    /* Iterate over vertices and free header_seq and pointer to outgoing edges
       first */
    if (graph->vertices != NULL) {
      for ( vertex = graph->vertices;
            vertex < (graph->vertices + graph->nof_vertices);
            vertex++
          )
      {
        gt_str_delete(vertex->header_seq);
        gt_free(vertex->edges);
      }
    }

    /* Now delete vertices and edges*/
    gt_free(graph->vertices);
    gt_free(graph->edges);
  }

  gt_free(graph);
}

/* Initialize a new vertex in <*graph>. Each vertex represents a contig and
   contains information about the sequence header <*header_seq>, sequence
   length <seq_len>, A-statistics <astat> and estimated copy number <copy_num>*/
void gt_scaffolder_graph_add_vertex(GtScaffolderGraph *graph,
                                    GtStr *header_seq,
                                    GtUword seq_len,
                                    float astat,
                                    float copy_num)
{
  GtUword nextfree;

  gt_assert(graph != NULL);
  gt_assert(graph->vertices != NULL);
  gt_assert(graph->nof_vertices < graph->max_nof_vertices);

  nextfree = graph->nof_vertices;

  /* Initialize vertex */
  graph->vertices[nextfree].seq_len = seq_len;
  graph->vertices[nextfree].astat = astat;
  graph->vertices[nextfree].copy_num = copy_num;
  graph->vertices[nextfree].nof_edges = 0;
  if (header_seq != NULL) {
    graph->vertices[nextfree].header_seq = header_seq;
  }
  graph->vertices[nextfree].state = GIS_UNVISITED;
  graph->vertices[nextfree].edges = NULL;

  graph->nof_vertices++;
}

/* Initialize a new edge in <*graph>. Each edge between two contig
   vertices <vstartID> and <vendID> contains information about the distance
   <dist>, standard deviation <std_dev>, number of pairs <num_pairs> and the
   direction of <vstartID> <dir> and corresponding <vendID> <same> */
void gt_scaffolder_graph_add_edge(GtScaffolderGraph *graph,
                                  GtScaffolderGraphVertex *vstart,
                                  GtScaffolderGraphVertex *vend,
                                  GtWord dist,
                                  float std_dev,
                                  GtUword num_pairs,
                                  bool dir,
                                  bool same)
{
  GtUword nextfree = graph->nof_edges;

  gt_assert(graph != NULL);
  gt_assert(graph->vertices != NULL);
  gt_assert(vstart != NULL);
  gt_assert(graph->edges != NULL);
  gt_assert(graph->nof_edges < graph->max_nof_edges);

  /* Inititalize edge */
  graph->edges[nextfree].start = vstart;
  graph->edges[nextfree].end = vend;
  graph->edges[nextfree].dist = dist;
  graph->edges[nextfree].std_dev = std_dev;
  graph->edges[nextfree].num_pairs = num_pairs;
  graph->edges[nextfree].sense = dir;
  graph->edges[nextfree].same = same;
  graph->edges[nextfree].state = GIS_UNVISITED;

  /* Add ptr to edge to start vertex */
  gt_assert(vstart != NULL);
  vstart->edges[vstart->nof_edges] = graph->edges + nextfree;
  vstart->nof_edges++;

  graph->nof_edges++;
}

/* Returns pointer to edge between <*vertex1> and <*vertex2> */
GtScaffolderGraphEdge
*gt_scaffolder_graph_find_edge(const GtScaffolderGraphVertex *vertex_1,
                               const GtScaffolderGraphVertex *vertex_2)
{
  GtUword eid;

  for (eid = 0; eid < vertex_1->nof_edges; eid++) {
    if (vertex_1->edges[eid]->end == vertex_2)
      return vertex_1->edges[eid];
  }
  return NULL;
}

/* determines corresponding vertex to contig header */
bool gt_scaffolder_graph_get_vertex(const GtScaffolderGraph *graph,
                                    GtScaffolderGraphVertex **vertex,
                                    const GtStr *header_seq)
{
  GtScaffolderGraphVertex *min_vertex, *max_vertex, *mid_vertex;
  int cmp;
  bool found;

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
        *vertex = mid_vertex;
        break;
      }
      else if (cmp < 0)
        min_vertex = mid_vertex + 1;
      else
        max_vertex = mid_vertex - 1;
    }
  return found;
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

/* determine vertex id*/
GtUword gt_scaffolder_graph_get_vertex_id(const GtScaffolderGraph *graph,
                                       const GtScaffolderGraphVertex *vertex)
{
  gt_assert(graph != NULL);
  /* LG: gt_assert(graph->vertices != NULL); */
  return (vertex - graph->vertices);
}

/* print graphrepresentation in dot-format into file filename */
int gt_scaffolder_graph_print(const GtScaffolderGraph *g,
                              const char *filename,
                              GtError *err)
{
  int had_err = 0;
  GtFile *f;

  gt_assert(g != NULL);

  f = gt_file_new(filename, "w", err);
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
  GtScaffolderGraphVertex *v;
  GtScaffolderGraphEdge *e;
  /* 0: GIS_UNVISITED, 1: GIS_POLYMORPHIC, 2: GIS_INCONSISTENT,
     3: GIS_REPEAT, 4: GIS_VISITED, 5: GIS_PROCESSED, 6: GIS_SCAFFOLD
     7: GIS_CYCLIC */
  const char *color_array[] = {"black", "gray80", "gainsboro", "ivory3", "red",
                               "green", "magenta", "blue"};

  gt_assert(g != NULL);

  /* print first line into f */
  gt_file_xprintf(f, "digraph {\n");

  /* iterate over all vertices and print them. add attribute color according
     to the current state */
  for (v = g->vertices; v < (g->vertices + g->nof_vertices); v++) {
    gt_file_xprintf(f, GT_WU " [color=\"%s\" label=\"%s\"];\n",
                    gt_scaffolder_graph_get_vertex_id(g, v),
                    color_array[v->state], gt_str_get(v->header_seq));
  }

  /* iterate over all edges and print them. add attribute color according to
     the current state and label the edge with the distance*/
  for (e = g->edges; e < (g->edges + g->nof_edges); e++) {
    gt_assert(e != NULL);
    gt_file_xprintf(f,
                    GT_WU " -> " GT_WU " [color=\"%s\" label="
                    "\"" GT_WD "\" arrowhead=\"%s\"];\n",
                    gt_scaffolder_graph_get_vertex_id(g, e->start),
                    gt_scaffolder_graph_get_vertex_id(g, e->end),
                    color_array[e->state], e->dist, e->sense?"normal":"inv");
  }

  /* print the last line into f */
  gt_file_xprintf(f, "}\n");
}

/* print graphrepresentation in dot-format into gt-filestream f */
void gt_scaffolder_graph_print_scaffold(const GtScaffolderGraph *g,
                                        GtFile *f)
{
  GtScaffolderGraphVertex *v;
  GtScaffolderGraphEdge *e;

  gt_assert(g != NULL);

  /* print first line into f */
  gt_file_xprintf(f, "digraph {\n");

  /* iterate over all vertices and print just the scaffold vertices */
  for (v = g->vertices; v < (g->vertices + g->nof_vertices); v++) {
    if (v->state == GIS_SCAFFOLD)
      gt_file_xprintf(f, GT_WU " [label=\"%s\"];\n",
                      gt_scaffolder_graph_get_vertex_id(g, v),
                      gt_str_get(v->header_seq));
  }

  /* iterate over all edges and print just the scaffold edges */
  for (e = g->edges; e < (g->edges + g->nof_edges); e++) {
    gt_assert(e != NULL);
    if (e->state == GIS_SCAFFOLD)
      gt_file_xprintf(f,
                      GT_WU " -> " GT_WU " [label="
                      "\"" GT_WD "\" arrowhead=\"%s\"];\n",
                      gt_scaffolder_graph_get_vertex_id(g, e->start),
                      gt_scaffolder_graph_get_vertex_id(g, e->end),
                      e->dist, e->sense?"normal":"inv");
  }

  /* print the last line into f */
  gt_file_xprintf(f, "}\n");
}

/* create scaffold graph from file */
int gt_scaffolder_graph_new_from_file(GtScaffolderGraph **graph_par,
                                      const char *ctg_filename,
                                      GtUword min_ctg_len,
                                      const char *dist_filename,
                                      bool astat_is_annotated,
                                      GtError *err)
{
  GtScaffolderGraph *graph;
  int had_err;
  GtUword nof_distances, nof_contigs;

  graph = NULL;

  /* count contigs */
  had_err = gt_scaffolder_parser_count_contigs(ctg_filename, min_ctg_len,
            &nof_contigs, err);

  if (had_err == 0)
  {
    /* allocate memory for vertices of scaffolder graph */

    /*graph = gt_scaffolder_graph_new(0, 0);*/
    graph = gt_malloc(sizeof (*graph));
    graph->vertices = NULL;
    graph->edges = NULL;

    gt_scaffolder_graph_init_vertices(graph, nof_contigs);
    /* parse contigs in FASTA-format and save them as vertices of
     scaffold graph */

    had_err = gt_scaffolder_parser_read_contigs(graph, ctg_filename,
              min_ctg_len, astat_is_annotated, err);
  }

  if (had_err == 0)
  {
    /* count distance information */
    had_err = gt_scaffolder_parser_count_distances(graph, dist_filename,
            &nof_distances, err);
  }

  if (had_err == 0)
  {
    if (nof_contigs == 1 && nof_distances == 0) {
      fprintf(stderr, "Graph only contains 1 vertex and no edges: "
                      "Did not perform scaffolding!\n");
      exit(0);
    }

    /* allocate memory for edges of scaffolder graph */
    gt_scaffolder_graph_init_edges(graph, nof_distances);
    /* parse distance information of contigs in abyss-dist-format and
       save them as edges of scaffold graph */
    had_err = gt_scaffolder_parser_read_distances(dist_filename,
              graph, false, err);
  }

  if (had_err != 0)
  {
    gt_scaffolder_graph_delete(graph);
    graph = NULL;
  }

  *graph_par = graph;

  return had_err;
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
                             bool print_graph,
                             GtError *err)
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
    unsigned i;

    gt_scaffolder_graph_init_vertices(graph, max_nof_vertices);
    if (graph->vertices == NULL)
      had_err = -1;

    for (i = 0; i < nof_vertices; i++) {
      gt_scaffolder_graph_add_vertex(graph, gt_str_new_cstr("foobar"),
                                     100, 20, 40);
      /* Simply allocate maximum amount of memory for pointer to potentially
         outgoing edges */
      if (nof_edges > 0)
        graph->vertices[i].edges
          = gt_malloc(sizeof (*graph->vertices->edges) * nof_edges);
    }
  }

  /* Init edge portion of graph. Connect every vertex with another vertex until
  <nof_edges> is reached. */
  if (init_edges) {
    GtScaffolderGraphVertex *vertex1, *vertex2;
    unsigned i;
    vertex1 = graph->vertices;
    vertex2 = graph->vertices;

    gt_scaffolder_graph_init_edges(graph, max_nof_edges);

    if (graph->edges == NULL)
      had_err = -1;

    /* Connect 1st vertex with every other vertex, then 2nd one, etc */
    for (i = 0; i < nof_edges; i++) {
      if (vertex2-graph->vertices < nof_vertices - 1)
        vertex2++;
      else if (vertex1-graph->vertices < nof_vertices - 2) {
        vertex1++;
        vertex2 = vertex1 + 1;
      }
      gt_scaffolder_graph_add_edge(graph, vertex1, vertex2, 2, 1.5, 4, true,
                                   true);
    }
  }

  /* Print the graph for diff comparison */
  if (print_graph) {
    char outfile[] = "gt_scaffolder_graph_test.dot";
    gt_scaffolder_graph_print(graph, outfile, err);
  }

  gt_scaffolder_graph_delete(graph);

  return had_err;
}
