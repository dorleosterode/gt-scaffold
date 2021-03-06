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

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/array_api.h"
#include "core/fasta_reader_rec.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/queue_api.h"

#include "extended/assembly_stats_calculator.h"

#include "gt_scaffolder_graph.h"

/* increment size for realloc of walk */
#define INCREMENT_SIZE 32

/* Check if vertex already has been filtered out of the graph */
static bool
vertex_is_marked(const GtScaffolderGraphVertex *vertex) {
  gt_assert(vertex != NULL);

  if (vertex->state == GIS_POLYMORPHIC || vertex->state == GIS_REPEAT
     || vertex->state == GIS_CYCLIC)
    return true;
  else
    return false;
}

/* Check if edge already has been filtered out of the graph */
static bool edge_is_marked(const GtScaffolderGraphEdge *edge) {
  gt_assert(edge != NULL);

  if (edge->state == GIS_INCONSISTENT || edge->state == GIS_POLYMORPHIC
     || edge->state == GIS_CYCLIC || edge->state == GIS_REPEAT)
    return true;
  else
    return false;
}

/* Marks edge and its twin as GIS_STATE */
static void mark_edge(GtScaffolderGraphEdge *edge, GraphItemState state) {
  GtUword eid;

  gt_assert(edge != NULL);

  /* mark edge as state*/
  edge->state = state;
  /* search twin-edge and mark twin-edge as state */
  for (eid = 0; eid < edge->end->nof_edges; eid++) {
    if (edge->end->edges[eid]->end == edge->start)
      edge->end->edges[eid]->state = state;
  }
}

/* Marks vertex and all its edges and twins as GIS_STATE */
static void mark_vertex(GtScaffolderGraphVertex *vertex, GraphItemState state) {
  GtUword eid;

  gt_assert(vertex != NULL);

  /* mark vertex as state */
  vertex->state = state;

  /* mark all edges and their twins as state */
  for (eid = 0; eid < vertex->nof_edges; eid++)
    mark_edge(vertex->edges[eid], state);
}

/* load a and copy number of every contig and mark repeated contigs */
int gt_scaffolder_graph_mark_repeats(const char *filename,
                                     GtScaffolderGraph *graph,
                                     float copy_num_cutoff,
                                     float astat_cutoff,
                                     GtError *err)
{
  const GtUword BUFSIZE_2 = 1024;
  FILE *file;
  char line[BUFSIZE_2+1], ctg_header[BUFSIZE_2+1];
  GtWord num1, num2, num3;
  float astat, copy_num;
  bool valid_contig;
  int had_err;
  GtStr *gt_str_field;
  GtScaffolderGraphVertex *ctg;

  had_err = 0;

  if (strlen(filename) != 0) {
    gt_str_field = gt_str_new();
    file = fopen(filename, "rb");
    if (file == NULL) {
      had_err = -1;
      gt_error_set(err, "can not read A-statistic file %s", filename);
    }
    if (had_err != -1)
    {
      /* iterate over each line of file until eof (contig record) */
      while (fgets(line, BUFSIZE_2, file) != NULL) {

        /* remove '\n' from end of line */
        line[strlen(line)-1] = '\0';

        num1 = 0;
        num2 = 0;
        num3 = 0;
        copy_num = 0.0;
        astat = 0.0;

        /*parse record consisting of ctg_header, a-statistics and copy number*/
        if (sscanf(line,"%s\t" GT_WD "\t" GT_WD "\t" GT_WD "\t%f\t%f",
          ctg_header, &num1, &num2, &num3, &copy_num, &astat) == 6)
        {
          /* get vertex id corresponding to root contig header */
          gt_str_set(gt_str_field, ctg_header);
          /* gt_scaffolder_graph_get_vertex_id in if Statement schieben */
          valid_contig = gt_scaffolder_graph_get_vertex(graph, &ctg,
                     gt_str_field);

          if (valid_contig) {
            ctg->astat = astat;
            ctg->copy_num = copy_num;
          }
        }
        else {
          had_err = -1;
          gt_error_set(err, "Invalid record in A-statistic file %s", filename);
          break;
        }
      }
      fclose(file);
    }
    gt_str_delete(gt_str_field);
  }

  if (had_err != -1)
  {
    GtScaffolderGraphVertex *vertex;

    /*iterate over all vertices and mark them as repeats if cutoff is exceeded*/
    for (vertex = graph->vertices;
      vertex < (graph->vertices + graph->nof_vertices); vertex++)
    {
       if (vertex->astat <= astat_cutoff ||
          (strlen(filename) != 0 && vertex->copy_num < copy_num_cutoff))
        mark_vertex(vertex, GIS_REPEAT);
    }
  }

  return had_err;
}

/* check if unique order of edges <*edge1>, <*edge2> with probability
   <cutoff> exists */
static bool
gt_scaffolder_graph_ambiguousorder(const GtScaffolderGraphEdge *edge1,
                                   const GtScaffolderGraphEdge *edge2,
                                   float cutoff)
{
  float expval, variance, interval, prob12, prob21, p_wrong;

  gt_assert(edge1 != NULL);
  gt_assert(edge2 != NULL);

  expval = edge1->dist - edge2->dist;
  variance = 2 * ((edge1->std_dev * edge1->std_dev) +
                  (edge2->std_dev * edge2->std_dev));
  interval = (0 - expval) / sqrt(variance);
  prob12 = 0.5 * (1 + erf(interval) );
  prob21 = 1.0 - prob12;

  p_wrong = 1.0 - MAX(prob12, prob21);
  return (p_wrong > cutoff) ? true : false;
}

/* calculate overlap of edge1->end and edge2->end with respect to the
   distance to edge1->start and edge2->start, which must be equal. */
static GtWord
gt_scaffolder_calculate_overlap(const GtScaffolderGraphEdge *edge1,
                                const GtScaffolderGraphEdge *edge2)
{
  GtWord overlap = 0;
  GtWord start1, start2, end1, end2;

  gt_assert(edge1 != NULL);
  gt_assert(edge2 != NULL);

  start1 = edge1->dist;
  start2 = edge2->dist;
  end1 = edge1->dist + edge1->end->seq_len - 1;
  end2 = edge2->dist + edge2->end->seq_len - 1;

  if (start2 <= end1 &&
      start1 <= end2)
  {
    GtWord intersect_start = MAX(start1, start2);
    GtWord intersect_end = MIN(end1, end2);
    overlap = intersect_end - intersect_start + 1;
  }
  return overlap;
}

/* checks if edge1->end and edge2->end are polymorphic */
static void
gt_scaffolder_graph_check_mark_polymorphic(const GtScaffolderGraphEdge *edge1,
                                           const GtScaffolderGraphEdge *edge2,
                                           float pcutoff,
                                           float cncutoff)
{
  GtScaffolderGraphVertex *poly_vertex;

  gt_assert(edge1 != NULL);
  gt_assert(edge2 != NULL);

  if (gt_scaffolder_graph_ambiguousorder(edge1, edge2, pcutoff) &&
      (edge1->end->copy_num + edge2->end->copy_num) < cncutoff) {
    /* mark vertex with lower copy number as polymorphic */
    if (edge1->end->copy_num < edge2->end->copy_num)
      poly_vertex = edge1->end;
    else
      poly_vertex = edge2->end;
    /* mark all edges of the polymorphic vertex as polymorphic */
    if (!vertex_is_marked(poly_vertex)) {
      mark_vertex(poly_vertex, GIS_POLYMORPHIC);
    }
  }
}

/* marks all edges of vertex in direction <sense> as inconsistent */
static void mark_edges_in_twin_dir(GtScaffolderGraphVertex *vertex,
                                   bool sense)
{
  GtUword eid;

  for (eid = 0; eid < vertex->nof_edges; eid++) {
    if (vertex->edges[eid]->sense == sense)
      vertex->edges[eid]->state = GIS_INCONSISTENT;
  }
}

/* mark polymorphic edges/vertices and inconsistent edges in scaffold graph */
void gt_scaffolder_graph_filter(GtScaffolderGraph *graph,
                                float pcutoff,
                                float cncutoff,
                                GtWord ocutoff)
{
  GtScaffolderGraphVertex *vertex;
  GtScaffolderGraphEdge *edge1, *edge2;
  GtUword eid1, eid2;
  GtWord sense_maxoverlap, antisense_maxoverlap, overlap;
  bool twin_dir;

  gt_assert(graph != NULL);

  /* iterate over all vertices */
  for (vertex = graph->vertices;
       vertex < (graph->vertices + graph->nof_vertices); vertex++) {

    /* ignore repeat vertices */
    if (vertex_is_marked(vertex))
      continue;

    /* iterate over all pairs of edges */
    for (eid1 = 0; eid1 < vertex->nof_edges; eid1++) {
      for (eid2 = eid1 + 1; eid2 < vertex->nof_edges; eid2++) {
        edge1 = vertex->edges[eid1];
        edge2 = vertex->edges[eid2];

        if (edge1->sense == edge2->sense) {
          /* check if edge1->end and edge2->end are polymorphic */
          gt_scaffolder_graph_check_mark_polymorphic(edge1, edge2,
                                                     pcutoff, cncutoff);
          /* SD: Nur das erste Paar polymoprh markieren? */
        }
      }
    }

    /* no need to check inconsistent edges for polymorphic vertices */
    if (vertex_is_marked(vertex))
      continue;

    sense_maxoverlap = 0;
    antisense_maxoverlap = 0;
    /* iterate over all pairs of edges, that are not polymorphic */
    for (eid1 = 0; eid1 < vertex->nof_edges; eid1++) {
      for (eid2 = eid1 + 1; eid2 < vertex->nof_edges; eid2++) {
        edge1 = vertex->edges[eid1];
        edge2 = vertex->edges[eid2];
        if ((edge1->sense == edge2->sense) &&
            (!edge_is_marked(edge1) && !edge_is_marked(edge2))) {
          overlap = gt_scaffolder_calculate_overlap(edge1, edge2);

          /* differentiate between maximal overlap of sense edge
             pairs and antisense edge pairs */
          if (edge1->sense && overlap > sense_maxoverlap)
            sense_maxoverlap = overlap;
          if (!edge1->sense && overlap > antisense_maxoverlap)
            antisense_maxoverlap = overlap;
        }
      }
    }

    /* check if maxoverlap is larger than ocutoff and mark edges
       as inconsistent */
    if (sense_maxoverlap > ocutoff || antisense_maxoverlap > ocutoff) {

      for (eid1 = 0; eid1 < vertex->nof_edges; eid1++) {
        if (sense_maxoverlap > ocutoff && vertex->edges[eid1]->sense) {
          vertex->edges[eid1]->state = GIS_INCONSISTENT;
          /* adapted from SGA, is it necessary to mark all edges
             in twin dir? */
          twin_dir = !vertex->edges[eid1]->same;
          mark_edges_in_twin_dir(vertex->edges[eid1]->end, twin_dir);
        }
        if (antisense_maxoverlap > ocutoff && !vertex->edges[eid1]->sense) {
          vertex->edges[eid1]->state = GIS_INCONSISTENT;
          twin_dir = vertex->edges[eid1]->same;
          mark_edges_in_twin_dir(vertex->edges[eid1]->end, twin_dir);
        }
      }

    }
  }
}

/* check if vertex holds just sense or antisense edges */
bool gt_scaffolder_graph_isterminal(const GtScaffolderGraphVertex *vertex)
{
  bool dir;
  bool set_dir = false;
  GtUword eid;

  gt_assert(vertex != NULL);

  if (vertex->nof_edges == 0)
    return true;

  for (eid = 0; eid < vertex->nof_edges; eid++) {
    if (set_dir) {
      if (vertex->edges[eid]->sense != dir &&
          !edge_is_marked(vertex->edges[eid]))
        return false;
    }
    else {
      if (!edge_is_marked(vertex->edges[eid])) {
        dir = vertex->edges[eid]->sense;
        set_dir = true;
      }
    }

  }

  return true;
}

/* traverse the graph and calculate all connected components. all
   terminal vertices for each connected component are stored in an
   GtArray. all GtArray for the terminial vertices are stored in
   ccs. */
void gt_scaffolder_calc_cc_and_terminals(const GtScaffolderGraph *graph,
                                         GtArray *ccs)
{
  GtArray *terminal_vertices = NULL;
  GtQueue *vqueue = NULL;
  GtScaffolderGraphVertex *vertex, *currentvertex, *nextvertex;
  GtUword eid;

  gt_assert(graph != NULL);
  gt_assert(ccs != NULL);

  for (vertex = graph->vertices; vertex <
      (graph->vertices + graph->nof_vertices); vertex++) {
    if (!vertex_is_marked(vertex))
      vertex->state = GIS_UNVISITED;
  }

  vqueue = gt_queue_new();

  for (vertex = graph->vertices; vertex <
      (graph->vertices + graph->nof_vertices); vertex++) {
    if (vertex_is_marked(vertex) ||
        vertex->state == GIS_VISITED)
      continue;

    vertex->state = GIS_PROCESSED;
    gt_queue_add(vqueue, vertex);
    /* create a new gt_array-object to store the terminal vertices for
       the next cc */
    terminal_vertices = gt_array_new(sizeof (vertex));

    while (gt_queue_size(vqueue) != 0) {
      currentvertex = (GtScaffolderGraphVertex*)gt_queue_get(vqueue);

      /* store all terminal vertices */
      if (gt_scaffolder_graph_isterminal(currentvertex))
        gt_array_add(terminal_vertices, currentvertex);

      currentvertex->state = GIS_VISITED;
      for (eid = 0; eid < currentvertex->nof_edges; eid++) {
        if (!edge_is_marked(currentvertex->edges[eid])) {
          nextvertex = currentvertex->edges[eid]->end;
          /* just take vertices, that are consistent */
          if (vertex_is_marked(nextvertex))
            continue;
          if (nextvertex->state == GIS_UNVISITED) {
            nextvertex->state = GIS_PROCESSED;
            gt_queue_add(vqueue, nextvertex);
          }
        }
      }
    }
    /* save the terminal_vertices for this cc in ccs */
    gt_array_add(ccs, terminal_vertices);
  }

  gt_queue_delete(vqueue);
}

/* checks, if e1->start = e2->end and e1->end = e2->start */
static bool
is_twin(const GtScaffolderGraphEdge *e1, const GtScaffolderGraphEdge *e2) {
  if (e1->start == e2->end &&
      e1->end == e2->start)
    return true;
  return false;
}

/* DFS to detect Cycles given a starting vertex */
GtScaffolderGraphEdge *
gt_scaffolder_detect_cycle_recursive(GtScaffolderGraphVertex *v,
                                     GtScaffolderGraphVertex *p,
                                     bool dir,
                                     GtArray *visited)
{
  GtUword eid;
  GtScaffolderGraphEdge *back;
  bool next_dir;

  gt_assert(v != NULL);

  gt_array_add(visited, v);
  v->state = GIS_VISITED;

  for (eid = 0; eid < v->nof_edges; eid++) {
    back = v->edges[eid];
    if (back->sense == dir &&
        !edge_is_marked(back) &&
        back->end != p) {

      if (!vertex_is_marked(back->end)) {
        if (back->end->state == GIS_VISITED)
          return back;

        if (back->end->state == GIS_UNVISITED) {
          /* SGA: set cur_dir to !back->twin->dir */
          if (back->same)
            next_dir = back->sense ? true : false;
          else
            next_dir = back->sense ? false : true;

          back = gt_scaffolder_detect_cycle_recursive(back->end, v, next_dir,
                                                      visited);

          if (back != NULL)
            return back;
        }
      }
    }
  }

  v->state = GIS_PROCESSED;
  return NULL;
}

/*  remove cycles */
void gt_scaffolder_removecycles(GtScaffolderGraph *graph) {
  bool found_cycle = true;
  GtUword i, j, k;
  GtArray *ccs, *terminal_vertices, *visited;
  GtScaffolderGraphEdge *back_edge;
  GtScaffolderGraphVertex *start, *v;

  gt_assert(graph != NULL);

  ccs = gt_array_new(sizeof (GtArray *));
  visited = gt_array_new(sizeof (GtScaffolderGraphVertex *));

  while (found_cycle) {
    found_cycle = false;

    gt_scaffolder_calc_cc_and_terminals(graph, ccs);

    /* initialize all vertices as not visited */
    for (v = graph->vertices; v < (graph->vertices + graph->nof_vertices);
         v++) {
      if (!vertex_is_marked(v))
        v->state = GIS_UNVISITED;
    }

    /* iterate over all ccs */
    for (i = 0; i < gt_array_size(ccs); i++) {
      terminal_vertices = *(GtArray **) gt_array_get(ccs, i);

      /* iterate over all terminal vertices of this cc */
      for (j = 0; j < gt_array_size(terminal_vertices); j++) {
        start = *(GtScaffolderGraphVertex **)
                gt_array_get(terminal_vertices, j);
        /* search for a cycle, if terminal vertex has edges */
        if (start->nof_edges > 0) {
          GtUword eid;
          bool dir = true;
          bool set_dir = false;

          for (eid = 0; eid < start->nof_edges; eid++) {
            if (!edge_is_marked(start->edges[eid])) {
              dir = start->edges[eid]->sense;
              set_dir = true;
            }
          }

          if (!set_dir)
            continue;

          if (vertex_is_marked(start))
            continue;

          back_edge = gt_scaffolder_detect_cycle_recursive(start, NULL,
                                                           dir, visited);

          /* mark all visited vertices as unvisited for the next search */
          for (k = 0; k < gt_array_size(visited); k++) {
            v = *(GtScaffolderGraphVertex **) gt_array_get(visited, k);
            v->state = GIS_UNVISITED;
          }

          gt_array_reset(visited);

          if (back_edge != NULL) {
            found_cycle = true;
            mark_vertex(back_edge->start, GIS_CYCLIC);
            mark_vertex(back_edge->end, GIS_CYCLIC);
          }
        }
      }
    }

    /* clean up of the ccs */
    for (i = 0; i < gt_array_size(ccs); i++)
      gt_array_delete(*(GtArray **) gt_array_get(ccs, i));

    gt_array_reset(ccs);
  }

  for (i = 0; i < gt_array_size(ccs); i++)
    gt_array_delete(*(GtArray **) gt_array_get(ccs, i));

  gt_array_delete(ccs);
  gt_array_delete(visited);
}

/* create new walk */
GtScaffolderGraphWalk *gt_scaffolder_walk_new(void)
{
  GtScaffolderGraphWalk *walk;

  walk = gt_malloc(sizeof (*walk));
  walk->total_contig_len = 0;
  walk->size = 0;
  walk->nof_edges = 0;
  walk->edges = NULL;
  return walk;
}

/* remove walk <*walk> */
void gt_scaffolder_walk_delete(GtScaffolderGraphWalk *walk)
{
  if (walk != NULL)
    gt_free(walk->edges);
  gt_free(walk);
}

/* add edge <*edge> to walk <*walk> */
void gt_scaffolder_walk_addegde(GtScaffolderGraphWalk *walk,
                                GtScaffolderGraphEdge *edge)
{
  gt_assert(walk != NULL);
  gt_assert(edge != NULL);

  if (walk->size == walk->nof_edges) {
    walk->size += INCREMENT_SIZE;
    walk->edges = gt_realloc(walk->edges, walk->size*sizeof (*walk->edges));
  }
  walk->edges[walk->nof_edges] = edge;
  walk->total_contig_len += edge->end->seq_len;
  walk->nof_edges++;
}

/* creates all minimal walks with respect to edge->dist from start to
   every other terminal vertex in the current cc. the walk with the
   greatest total contig length is returned */
GtScaffolderGraphWalk
*gt_scaffolder_create_walk(GtScaffolderGraph *graph,
                           GtScaffolderGraphVertex *start)
{
  GtQueue *wqueue;
  GtArray *terminal_vertices;
  GtScaffolderGraphEdge *edge, *reverseedge, *nextedge, **edgemap;
  GtScaffolderGraphVertex *endvertex, *currentvertex, *nextendvertex;
  GtUword lengthbestwalk, lengthcwalk, eid, i;
  GtScaffolderGraphWalk *bestwalk, *currentwalk;
  float distance, *distancemap;
  GtScaffolderGraphNode *node;
  bool dir;

  gt_assert(graph != NULL);
  gt_assert(start != NULL);

  lengthbestwalk = 0;

  /* do we need this initialization? isn't bestwalk
     just a pointer to the best currentwalk? */
  bestwalk = gt_scaffolder_walk_new();

  wqueue = gt_queue_new();
  terminal_vertices = gt_array_new(sizeof (start));

  /* initialize distancemap with GT_WORD_MAX, we want to minimize
     over distancemap */
  distancemap = gt_malloc(sizeof (*distancemap) * graph->nof_vertices);
  for (i = 0; i < graph->nof_vertices; i++)
    distancemap[i] = GT_WORD_MAX;

  edgemap = gt_malloc(sizeof (*edgemap)*graph->nof_vertices);

  /* check if node has edges */
  if (start->nof_edges == 0) {
    gt_free(edgemap);
    gt_free(distancemap);
    return NULL;
  }

  for (eid = 0; eid < start->nof_edges; eid++) {
    edge = start->edges[eid];
    if (!edge_is_marked(edge) &&
        !vertex_is_marked(edge->end))
    {
      GtScaffolderGraphNode *initial_node = gt_malloc(sizeof (*initial_node));
      endvertex = edge->end;

      GtUword endvertex_index = gt_scaffolder_graph_get_vertex_id(graph,
                                                                  endvertex);
      distancemap[endvertex_index] = edge->dist;
      edgemap[endvertex_index] = edge;

      initial_node->edge = edge;
      initial_node->dist = edge->dist;

      gt_queue_add(wqueue, initial_node);
    }
  }

  while (gt_queue_size(wqueue) != 0) {
    node = (GtScaffolderGraphNode*)gt_queue_get(wqueue);
    edge = node->edge;
    endvertex = edge->end;

    /* determine opposite direction of twin of edge egde */
    /* according to SGA: EdgeDir yDir = !pXY->getTwin()->getDir(); */
    if (edge->same)
      dir = edge->sense;
    else
      dir = !edge->sense;

    /* store all terminal vertices */
    if (gt_scaffolder_graph_isterminal(endvertex))
      gt_array_add(terminal_vertices, endvertex);

    for (eid = 0; eid < endvertex->nof_edges; eid++) {
      nextedge = endvertex->edges[eid];
      if (nextedge->sense == dir) {
        if (!edge_is_marked(nextedge) &&
            !vertex_is_marked(nextedge->end) &&
            !is_twin(edge, nextedge))
        {
          nextendvertex = nextedge->end;

          distance = node->dist + nextedge->dist;

          /* GT_WORD_MAX is the initial value */
          GtUword next_endvertex_index =
            gt_scaffolder_graph_get_vertex_id(graph, nextendvertex);
          if (distancemap[next_endvertex_index] == GT_WORD_MAX ||
              distancemap[next_endvertex_index] > distance)
            {
              GtScaffolderGraphNode *current_node =
              gt_malloc(sizeof (*current_node));
              distancemap[next_endvertex_index] = distance;
              edgemap[next_endvertex_index] = nextedge;

              current_node->edge = nextedge;
              current_node->dist = distance;

              gt_queue_add(wqueue, current_node);
            }
        }
      }
    }
    gt_free(node);
  }

  /* create walk for each found terminal vertex and choose the walk
     with the greatest total contig length */
  while (gt_array_size(terminal_vertices) != 0) {
    currentvertex =
      *(GtScaffolderGraphVertex **) gt_array_pop(terminal_vertices);
    gt_assert(currentvertex < graph->vertices + graph->nof_vertices);
    gt_assert(currentvertex >= graph->vertices);

    currentwalk = gt_scaffolder_walk_new();
    while (currentvertex != start) {
      reverseedge = edgemap[gt_scaffolder_graph_get_vertex_id(graph,
                            currentvertex)];
      currentvertex = reverseedge->start;
      gt_scaffolder_walk_addegde(currentwalk, reverseedge);
    }

    currentwalk->total_contig_len += start->seq_len;
    lengthcwalk = currentwalk->total_contig_len;

    if (lengthcwalk > lengthbestwalk) {
      gt_scaffolder_walk_delete(bestwalk);
      bestwalk = currentwalk;
      lengthbestwalk = lengthcwalk;
    }
    else
      gt_scaffolder_walk_delete(currentwalk);
  }
  gt_free(distancemap);
  gt_free(edgemap);
  gt_array_delete(terminal_vertices);
  gt_queue_delete(wqueue);

  return bestwalk;
}

/* constructs the scaffolds for every cc. all vertices and edges in a
   scaffold are marked as GIS_SCAFFOLD. */
void gt_scaffolder_makescaffold(GtScaffolderGraph *graph)
{
  GtUword max_num_bases, i, j;
  GtScaffolderGraphWalk *walk, *bestwalk;
  GtScaffolderGraphVertex *start;
  GtArray *terminal_vertices, *cc_walks, *ccs;

  gt_assert(graph != NULL);

  /* remove cycles */
  gt_scaffolder_removecycles(graph);

  cc_walks = gt_array_new(sizeof (walk));

  /* create GtArray to store all ccs in it */
  ccs = gt_array_new(sizeof (GtArray *));

  gt_scaffolder_calc_cc_and_terminals(graph, ccs);

  for (i = 0; i < gt_array_size(ccs); i++) {
    terminal_vertices = *(GtArray **) gt_array_get(ccs, i);

    /* mark all lonesome vertices as scaffold */
    if (gt_array_size(terminal_vertices) == 1) {
      GtScaffolderGraphVertex *v;
      v = *(GtScaffolderGraphVertex **) gt_array_get(terminal_vertices, 0);
      if (v->nof_edges > 0) {
        bool lonesome = true;
        GtUword eid;
        for (eid = 0; eid < v->nof_edges; eid++) {
          if (!edge_is_marked(v->edges[eid])) {
            lonesome = false;
            break;
          }
        }
        if (lonesome)
          v->state = GIS_SCAFFOLD;
      }
      else
        v->state = GIS_SCAFFOLD;
    }

    if (gt_array_size(terminal_vertices) > 1) {
      /* calculate all paths between terminal vertices in this cc */
      for (j = 0; j < gt_array_size(terminal_vertices); j++) {
        start = *(GtScaffolderGraphVertex **)
                gt_array_get(terminal_vertices, j);
        gt_assert(start >= graph->vertices);
        gt_assert(start < graph->vertices + graph->nof_vertices);
        walk = gt_scaffolder_create_walk(graph, start);
        if (walk != NULL) {
          gt_array_add(cc_walks, walk);
        }
      }
    }

    /* the best walk in this cc is chosen */
    max_num_bases = 0;
    bestwalk = NULL;
    for (j = 0; j < gt_array_size(cc_walks); j++) {
      walk = *(GtScaffolderGraphWalk **) gt_array_get(cc_walks, j);
      if (walk->total_contig_len > max_num_bases) {
        bestwalk = walk;
        max_num_bases = walk->total_contig_len;
      }
    }

    /* mark all nodes and edges in the best walk as GIS_SCAFFOLD */
    if (bestwalk != NULL) {
      GtWord id;
      GtUword ed;
      bestwalk->edges[bestwalk->nof_edges - 1]->start->state = GIS_SCAFFOLD;
      for (id = (bestwalk->nof_edges - 1); id >= 0; id--) {
        bestwalk->edges[id]->state = GIS_SCAFFOLD;
        /* mark also the twin edges! */
        for (ed = 0; ed < bestwalk->edges[id]->end->nof_edges; ed++) {
          if (is_twin(bestwalk->edges[id], bestwalk->edges[id]->end->edges[ed]))
              bestwalk->edges[id]->end->edges[ed]->state = GIS_SCAFFOLD;
        }
        bestwalk->edges[id]->end->state = GIS_SCAFFOLD;
      }
    }

    /* free all walks in cc_walks */
    for (j = 0; j < gt_array_size(cc_walks); j++)
      gt_scaffolder_walk_delete(*(GtScaffolderGraphWalk **)
                                gt_array_get(cc_walks, j));

    gt_array_reset(cc_walks);
  }

  for (i = 0; i < gt_array_size(ccs); i++)
    gt_array_delete(*(GtArray **) gt_array_get(ccs, i));

  /* free all walks in cc_walks */
  for (i = 0; i < gt_array_size(cc_walks); i++)
    gt_scaffolder_walk_delete(*(GtScaffolderGraphWalk **)
                  gt_array_get(cc_walks, i));

  gt_array_delete(ccs);
  gt_array_delete(cc_walks);
}

/* functions to use GtScaffolderGraphRecords */
GtScaffolderGraphRecord *
gt_scaffolder_graph_record_new(GtScaffolderGraphVertex *root) {
  GtScaffolderGraphRecord *rec;

  gt_assert(root != NULL);

  rec = gt_malloc(sizeof (*rec));
  rec->root = root;
  rec->edges = gt_array_new(sizeof (GtScaffolderGraphEdge *));
  return rec;
}

void gt_scaffolder_graph_record_add_edge(GtScaffolderGraphRecord *rec,
                                         GtScaffolderGraphEdge *edge)
{
  gt_assert(rec != NULL);
  gt_assert(edge != NULL);

  gt_array_add(rec->edges, edge);
}

void gt_scaffolder_graph_record_delete(GtScaffolderGraphRecord *rec) {
  gt_assert(rec != NULL);

  gt_array_delete(rec->edges);

  gt_free(rec);
}

/* iterate over graph and return each scaffold in a scaffold record */
GtArray *gt_scaffolder_graph_iterate_scaffolds(const GtScaffolderGraph *graph,
                                          GtAssemblyStatsCalculator *scaf_stats)
{
  GtScaffolderGraphVertex *vertex, *next_edge_end;
  GtScaffolderGraphEdge *next_edge, *edge, *unmarked_edge;
  GtScaffolderGraphRecord *rec;
  GtArray *records;
  bool dir;
  GtUword eid, nof_edges_in_dir, scaf_seqlen, nof_scaffold_edges;

  next_edge = NULL;
  edge = NULL;
  unmarked_edge = NULL;

  records = gt_array_new(sizeof (rec));

  /* initialize all vertices as not visited */
  for (vertex = graph->vertices;
       vertex < (graph->vertices + graph->nof_vertices); vertex++) {
    if (!vertex_is_marked(vertex) && vertex->state != GIS_SCAFFOLD)
      vertex->state = GIS_UNVISITED;
  }

  /* iterate over all vertices */
  for (vertex = graph->vertices;
       vertex < (graph->vertices + graph->nof_vertices); vertex++) {

    if (vertex->state == GIS_VISITED || vertex_is_marked(vertex))
      continue;

    /* count unmarked edges and save one of them */
    nof_scaffold_edges = 0;
    for (eid = 0; eid < vertex->nof_edges; eid++) {
      edge = vertex->edges[eid];
      if (!edge_is_marked(edge)) {
        if (edge->state == GIS_SCAFFOLD) {
          nof_scaffold_edges++;
          unmarked_edge = edge;
        }
      }
    }

    if (nof_scaffold_edges <= 1) {
      /* found new scaffold */
      rec = gt_scaffolder_graph_record_new(vertex);
      scaf_seqlen = vertex->seq_len;

      vertex->state = GIS_VISITED;

      if (nof_scaffold_edges == 1) {

        next_edge = unmarked_edge;

        while (1) {
          /* store edge in scaffold-record */
          gt_scaffolder_graph_record_add_edge(rec, next_edge);
          scaf_seqlen += next_edge->dist;

          next_edge_end = next_edge->end;
          scaf_seqlen += next_edge_end->seq_len;

          if (next_edge_end->state == GIS_VISITED)
            break;

          next_edge_end->state = GIS_VISITED;

          /* according to SGA: EdgeDir nextDir = !pXY->getTwin()->getDir(); */
          if (next_edge->same)
            dir = next_edge->sense;
          else
            dir = !next_edge->sense;

          /* count valid edges (unmarked, no twin) in direction dir
             and save one of them */
          nof_edges_in_dir = 0;
          for (eid = 0; eid < next_edge_end->nof_edges; eid++) {
            edge = next_edge_end->edges[eid];
            if (edge->sense == dir && !edge_is_marked(edge)
                && !is_twin(next_edge, edge) && edge->state == GIS_SCAFFOLD) {
              nof_edges_in_dir++;
              unmarked_edge = edge;
            }
          }

          if (nof_edges_in_dir == 1)
            next_edge = unmarked_edge;
          else
            break;
        }
      }
      /* store scaffold-record rec */
      gt_array_add(records, rec);
      gt_assembly_stats_calculator_add(scaf_stats, scaf_seqlen);
    }
  }
  return records;
}

/* write scaffold into file */
int gt_scaffolder_graph_write_scaffold(GtArray *records,
                                       const char *file_name,
                                       GtError *err)
{
  GtFile *file;
  GtScaffolderGraphRecord *rec;
  GtScaffolderGraphEdge *e;
  GtUword i, j;
  int had_err = 0;

  /* create file */
  file = gt_file_new(file_name, "w", err);
  if (file == NULL) {
    had_err = -1;
    gt_error_set(err,"can not create file %s", file_name);
  }

  if (had_err == 0) {
    for (i = 0; i < gt_array_size(records); i++) {
      rec = *(GtScaffolderGraphRecord **) gt_array_get(records, i);

      gt_file_xprintf(file, "%s", gt_str_get(rec->root->header_seq));

      for (j = 0; j < gt_array_size(rec->edges); j++) {
        e = *(GtScaffolderGraphEdge **) gt_array_get(rec->edges, j);

        gt_file_xprintf(file, "\t%s," GT_WD ",%f,%d,%d,",
                        gt_str_get(e->end->header_seq),
                        e->dist,
                        e->std_dev,
                        e->sense,
                        e->same);

      }

      gt_file_xprintf(file, "\n");
    }

    gt_file_delete(file);
  }

  return had_err;
}
