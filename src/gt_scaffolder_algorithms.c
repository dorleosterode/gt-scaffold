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

#include "gt_scaffolder_graph.h"

/* increment size for realloc of walk */
/* SK: Increment evaluieren, zB mit Faktor 1.2 multiplizieren */
#define INCREMENT_SIZE 32

/* SK: BUFSIZE_2 lokal in Funktion deklarieren */
const GtUword BUFSIZE_2 = 1024;

/* Check if vertex already has been filtered out of the graph */
static bool
vertex_is_marked(const GtScaffolderGraphVertex *vertex) {
  if (vertex->state == GIS_POLYMORPHIC || vertex->state == GIS_REPEAT
     || vertex->state == GIS_CYCLIC)
    return true;
  else
    return false;
}

/* Check if edge already has been filtered out of the graph */
static bool edge_is_marked(const GtScaffolderGraphEdge *edge) {
  if (edge->state == GIS_INCONSISTENT || edge->state == GIS_POLYMORPHIC
     || edge->state == GIS_CYCLIC)
    return true;
  else
    return false;
}

/* load a and copy number of every contig and mark repeated contigs */
int gt_scaffolder_graph_mark_repeats(const char *filename,
                                     GtScaffolderGraph *graph,
                                     float copy_num_cutoff,
                                     float astat_cutoff,
                                     GtError *err)
{
  FILE *file;
  char line[BUFSIZE_2+1], ctg_header[BUFSIZE_2+1];
  GtWord num1, num2, num3;
  float astat, copy_num;
  bool valid_contig;
  int had_err;
  GtStr *gt_str_field;
  GtScaffolderGraphVertex *vertex, *ctg;

  had_err = 0;
  valid_contig = false;
  gt_str_field = gt_str_new();
  file = fopen(filename, "rb");
  if (file == NULL) {
    had_err = -1;
    gt_error_set(err, "can not read file %s", filename);
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

      /* parse record consisting of ctg_header, a-statistics and copy number */
      /* SD: %[^>,] failed, parsed the whole line instead */
      if (sscanf(line,"%s\t" GT_WD "\t" GT_WD "\t" GT_WD "\t%f\t%f",
          ctg_header, &num1, &num2, &num3, &copy_num, &astat) == 6)
      {
        /* get vertex id corresponding to root contig header */
        gt_str_set(gt_str_field, ctg_header);
        /* gt_scaffolder_graph_get_vertex_id in if Statement schieben */
        valid_contig = gt_scaffolder_graph_get_vertex(graph, &ctg,
                     gt_str_field);

        if (valid_contig) {
          /* SK: Evaluieren, ob Knoten hier als Repeat gesetzt werden können */
          ctg->astat = astat;
          ctg->copy_num = copy_num;
        }
      }
      else
        had_err = -1;
        /* SK: Fehlermeldung setzen */
    }
  }
  fclose(file);
  gt_str_delete(gt_str_field);

  if (had_err != -1)
  {
    /* iterate over all vertices */
    /* SK: Vertex lokal deklarieren */
    for (vertex = graph->vertices;
      vertex < (graph->vertices + graph->nof_vertices); vertex++)
    {
      /* SK: Cutoffs und Bedingungen nochmal evaluieren */
      if (vertex->astat <= astat_cutoff || vertex->copy_num < copy_num_cutoff)
        vertex->state = GIS_REPEAT;
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
  gt_assert(edge1 != NULL);
  gt_assert(edge2 != NULL);

  float expval, variance, interval, prob12, prob21;

  expval = edge1->dist - edge2->dist;
  variance = 2 * (edge1->std_dev * edge1->std_dev) -
             (edge2->std_dev * edge2->std_dev);
  interval = -expval / sqrt(variance);
  prob12 = 0.5 * (1 + erf(interval) );
  prob21 = 1.0 - prob12;

  return (prob12 <= cutoff && prob21 <= cutoff) ? true : false;
}

/* LG:overlap can be negative!? */
static GtWord gt_scaffolder_calculate_overlap(GtScaffolderGraphEdge *edge1,
                                              GtScaffolderGraphEdge *edge2)
{
  gt_assert(edge1 != NULL);
  gt_assert(edge2 != NULL);

  GtWord overlap = 0;

  if (edge2->dist >= edge1->end->seq_len ||
      edge1->dist >= edge2->end->seq_len)
  {
    GtWord intersect_start = MAX(edge1->dist, edge2->dist);
    GtWord intersect_end = MAX(edge1->end->seq_len, edge2->end->seq_len);
    /* SK: gt_assert(intersect_end >= intersect_start); */
    overlap = intersect_end - intersect_start;
  }
  return overlap;
}

static void
gt_scaffolder_graph_check_mark_polymorphic(GtScaffolderGraphEdge *edge1,
                                           GtScaffolderGraphEdge *edge2,
                                           float pcutoff,
                                           float cncutoff)
{
  gt_assert(edge1 != NULL);
  gt_assert(edge2 != NULL);

  GtScaffolderGraphVertex *poly_vertex;

  if (gt_scaffolder_graph_ambiguousorder(edge1, edge2, pcutoff) &&
      (edge1->end->copy_num + edge2->end->copy_num) < cncutoff) {
    /* mark vertex with lower copy number as polymorphic */
    if (edge1->end->copy_num < edge2->end->copy_num)
      poly_vertex = edge1->end;
    else
      poly_vertex = edge2->end;
    /* mark all edges of the polymorphic vertex as polymorphic */
    if (!vertex_is_marked(poly_vertex)) {
      GtUword eid;
      for (eid = 0; eid < poly_vertex->nof_edges; eid++)
        poly_vertex->edges[eid]->state = GIS_POLYMORPHIC;
      poly_vertex->state = GIS_POLYMORPHIC;
    }
  }
}

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
  gt_assert(graph != NULL);

  GtScaffolderGraphVertex *vertex;
  GtScaffolderGraphEdge *edge1, *edge2;
  GtUword eid1, eid2;
  GtWord sense_maxoverlap, antisense_maxoverlap, overlap;
  bool twin_dir;

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
          /* LG: adapted from SGA, is it necessary to mark all edges
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
  GtUword eid;

  gt_assert(vertex != NULL);

  if (vertex->nof_edges == 0)
    return true;

  dir = vertex->edges[0]->sense;
  for (eid = 1; eid < vertex->nof_edges; eid++) {
    if (vertex->edges[eid]->sense != dir)
      return false;
  }

  return true;
}

void gt_scaffolder_calc_cc_and_terminals(const GtScaffolderGraph *graph,
                                         GtArray *ccs)
{
  gt_assert(graph != NULL);
  gt_assert(ccs != NULL);

  GtArray *terminal_vertices = NULL;
  GtQueue *vqueue = NULL;
  GtScaffolderGraphVertex *vertex, *currentvertex, *nextvertex;
  GtUword eid;

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

static bool
is_twin(const GtScaffolderGraphEdge *e1, const GtScaffolderGraphEdge *e2) {
  if (e1->start == e2->end &&
      e1->end == e2->start)
    return true;
  return false;
}

/* DFS to detect Cycles given a starting vertex */
GtScaffolderGraphEdge
*gt_scaffolder_detect_cycle(GtScaffolderGraphVertex *v,
                            bool dir,
                            GtArray *visited) {
  GtUword eid, beid;
  GtScaffolderGraphEdge *back;
  GtArray *stack;

  gt_assert(v != NULL);

  stack = gt_array_new(sizeof (GtScaffolderGraphEdge *));

  gt_array_add(visited, v);

  for (eid = 0; eid < v->nof_edges; eid++) {
    if (v->edges[eid]->sense == dir &&
        !edge_is_marked(v->edges[eid])) {
      gt_array_add(stack, v->edges[eid]);

      while (gt_array_size(stack) != 0) {
        back = *(GtScaffolderGraphEdge **) gt_array_pop(stack);
        if (!vertex_is_marked(back->end)) {
          if (back->end->state == GIS_VISITED)
            return back; /* SK: Stack aufräumen */
          back->end->state = GIS_VISITED;
          for (beid = 0; beid < back->end->nof_edges; beid++) {
            if (back->end->edges[beid]->sense == dir &&
                !edge_is_marked(back->end->edges[beid]) &&
                !is_twin(back, back->end->edges[beid]))
              gt_array_add(stack, back->end->edges[beid]);
          }
        }
      }
    }
  }

  gt_array_delete(stack);

  return NULL;

}

/*  remove cycles */
void gt_scaffolder_removecycles(GtScaffolderGraph *graph) {
  bool done = false, found_cycle;
  GtUword i, j, k;
  GtArray *ccs, *terminal_vertices, *visited;
  GtScaffolderGraphEdge *back_edge;
  GtScaffolderGraphVertex *start, *v;

  ccs = gt_array_new(sizeof (GtArray *));
  visited = gt_array_new(sizeof (GtScaffolderGraphVertex *));

  while (!done) {
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
          back_edge = gt_scaffolder_detect_cycle(start,
                      start->edges[0]->sense, visited);

          /* mark all visited vertices as unvisited for the next search */
          for (k = 0; k < gt_array_size(visited); k++) {
            v = *(GtScaffolderGraphVertex **) gt_array_get(visited, k);
            v->state = GIS_UNVISITED;
          }

          gt_array_reset(visited);

          if (back_edge != NULL) {
            found_cycle = true;
            back_edge->state = GIS_CYCLIC;
            back_edge->start->state = GIS_CYCLIC;
            back_edge->end->state = GIS_CYCLIC;
          }
        }
      }
    }

    for (i = 0; i < gt_array_size(ccs); i++)
      gt_array_delete(*(GtArray **) gt_array_get(ccs, i));

    gt_array_reset(ccs);

    /* if a cycle was found, we have to search again */
    /* SK: done = !found_cycle, eine Variable verwenden */
    done = found_cycle ? false : true;
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

GtScaffolderGraphWalk
*gt_scaffolder_create_walk(GtScaffolderGraph *graph,
                           GtScaffolderGraphVertex *start)
{
  gt_assert(graph != NULL);
  gt_assert(start != NULL);

  /* BFS-Traversierung innerhalb aktueller Zusammenhangskomponente
     ausgehend von terminalen Knoten zu terminalen Knoten */
  GtQueue *wqueue;
  GtArray *terminal_vertices;
  GtScaffolderGraphEdge *edge, *reverseedge, *nextedge, **edgemap;
  GtScaffolderGraphVertex *endvertex, *currentvertex, *nextendvertex;
  GtUword lengthbestwalk, lengthcwalk, eid, i;
  GtScaffolderGraphWalk *bestwalk, *currentwalk;
  float distance, *distancemap;
  bool dir;
  lengthbestwalk = 0;
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
    return NULL;
  }

  dir = start->edges[0]->sense;
  for (eid = 0; eid < start->nof_edges; eid++) {
    edge = start->edges[eid];
    if (!edge_is_marked(edge) &&
        !vertex_is_marked(edge->end))
    {
      endvertex = edge->end;

      /* SK: genometools hashes verwenden, Dichte evaluieren
         SK: DistEst beim Einlesen prüfen, Basisadresse verwenden fuer Index */
      if (!vertex_is_marked(endvertex)) {
        GtUword endvertex_index = gt_scaffolder_graph_get_vertex_id(graph,
                                  endvertex);
        distancemap[endvertex_index] = edge->dist;
        edgemap[endvertex_index] = edge;

        gt_queue_add(wqueue, edge);
      }
    }
  }

  while (gt_queue_size(wqueue) != 0) {
    edge = (GtScaffolderGraphEdge*)gt_queue_get(wqueue);
    endvertex = edge->end;

    /* store all terminal vertices */
    if (gt_scaffolder_graph_isterminal(endvertex))
      gt_array_add(terminal_vertices, endvertex);

    for (eid = 0; eid < endvertex->nof_edges; eid++) {
      nextedge = endvertex->edges[eid];
      if (nextedge->sense == dir) {
        if (!edge_is_marked(nextedge) &&
            !vertex_is_marked(nextedge->end))
        {
          nextendvertex = nextedge->end;

          if (!vertex_is_marked(nextendvertex)) {

            distance = edge->dist + nextedge->dist;

            /* GT_WORD_MAX is the initial value */
            GtUword next_endvertex_index =
                    gt_scaffolder_graph_get_vertex_id(graph, nextendvertex);
            if (distancemap[next_endvertex_index] == GT_WORD_MAX ||
                distancemap[next_endvertex_index] > distance)
              {
                distancemap[next_endvertex_index] = distance;
                edgemap[next_endvertex_index] = nextedge;
                gt_queue_add(wqueue, nextedge);
              }
          }
        }
      }
    }
  }

  /* Ruecktraversierung durch EdgeMap für alle terminalen Knoten
     Konstruktion des Walks  */
  while (gt_array_size(terminal_vertices) != 0) {
    currentvertex =
      *(GtScaffolderGraphVertex **) gt_array_pop(terminal_vertices);
    gt_assert(currentvertex < graph->vertices + graph->nof_vertices);
    gt_assert(currentvertex >= graph->vertices);

    currentwalk = gt_scaffolder_walk_new();
    while (currentvertex != start) {
      reverseedge = edgemap[gt_scaffolder_graph_get_vertex_id(graph,
                            currentvertex)];
      /* Start NICHT end */
      currentvertex = reverseedge->start;
      /* Speicherung des aktuellen Walks */
      gt_scaffolder_walk_addegde(currentwalk, reverseedge);
    }

    /* Ermittelung Contig-Gesamtlaenge des aktuellen Walks  */
    lengthcwalk = currentwalk->total_contig_len;

    if (lengthcwalk > lengthbestwalk) {
      gt_scaffolder_walk_delete(bestwalk);
      bestwalk = currentwalk;
      lengthbestwalk = lengthcwalk;
    }
    else
      gt_scaffolder_walk_delete(currentwalk);
  }

  gt_array_delete(terminal_vertices);
  gt_queue_delete(wqueue);

  return bestwalk;
}

/* Konstruktion des Scaffolds mit groesster Contig-Gesamtlaenge */
void gt_scaffolder_makescaffold(GtScaffolderGraph *graph)
{
  gt_assert(graph != NULL);

  GtUword eid, max_num_bases, i, j;
  GtScaffolderGraphWalk *walk, *bestwalk;
  GtScaffolderGraphVertex *start;
  GtArray *terminal_vertices, *cc_walks, *ccs;

  /* Entfernung von Zyklen
     gt_scaffolder_removecycles(graph); */

  cc_walks = gt_array_new(sizeof (walk));

  /* create GtArray to store all ccs in it */
  ccs = gt_array_new(sizeof (GtArray *));

  gt_scaffolder_calc_cc_and_terminals(graph, ccs);

  for (i = 0; i < gt_array_size(ccs); i++) {
    terminal_vertices = *(GtArray **) gt_array_get(ccs, i);

    /* calculate all paths between terminal vertices in this cc */
    for (j = 0; j < gt_array_size(terminal_vertices); j++) {
      start = *(GtScaffolderGraphVertex **) gt_array_get(terminal_vertices, j);
      gt_assert(start >= graph->vertices);
      gt_assert(start < graph->vertices + graph->nof_vertices);
      walk = gt_scaffolder_create_walk(graph, start);
      if (walk != NULL) {
        gt_array_add(cc_walks, walk);
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
      bestwalk->edges[0]->start->state = GIS_SCAFFOLD;
      for (eid = 0; eid < bestwalk->nof_edges; eid++) {
        bestwalk->edges[eid]->state = GIS_SCAFFOLD;
        bestwalk->edges[eid]->end->state = GIS_SCAFFOLD;
      }
    }
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
