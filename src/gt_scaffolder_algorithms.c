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
#define INCREMENT_SIZE 32

/* load a and copy number of every contig and mark repeated contigs */
int gt_scaffolder_graph_mark_repeats(const char *filename,
                                     GtScaffolderGraph *graph,
                                     float copy_num_cutoff,
                                     float astat_cutoff,
                                     GtError *err)
{
  FILE *file;
  char line[1024], *field;
  GtUword root_ctg_id, field_counter;
  float astat, copy_num;
  bool astat_found, copy_num_found;
  int had_err, had_err_2;
  GtStr *gt_str_field;
  GtScaffolderGraphVertex *vertex;

  had_err = 0;
  had_err_2 = 0;
  file = fopen(filename, "rb");
  if (file == NULL) {
    had_err = -1;
    gt_error_set(err, " can not read file %s ", filename);
  }

  if (had_err != -1)
  {
    /* iterate over each line of file until eof (contig record) */
    while (fgets(line, 1024, file) != NULL) {

      /* remove '\n' from end of line */
      line[strlen(line)-1] = '\0';

      /* split line by first tab delimiter */
      field = strtok(line,"\t");

      /* get vertex id corresponding to root contig header */
      root_ctg_id = 0;
      gt_str_field = gt_str_new_cstr(field);
      had_err_2 = gt_scaffolder_graph_get_vertex_id(graph, &root_ctg_id,
                gt_str_field);
      gt_str_delete(gt_str_field);

      field_counter = 0;
      copy_num = 0.0;
      astat = 0.0;
      copy_num_found = false;
      astat_found = false;

      if (had_err_2 == 0) {
        /* iterate over tab delimited records */
        while (field != NULL)
        {

          /* parse record consisting of a-statistics and copy number */
          if (field_counter == 4 && sscanf(field,"%f", &copy_num) == 1)
            copy_num_found = true;
          if (field_counter == 5 && sscanf(field,"%f", &astat) == 1)
            astat_found = true;

          /* split line by next tab delimiter */
          field = strtok(NULL,"\t");
          field_counter++;
        }

        /* save a-statistics and copy number */
        if (copy_num_found && astat_found)
        {
          graph->vertices[root_ctg_id].astat = astat;
          graph->vertices[root_ctg_id].copy_num = copy_num;
        }
        else
          had_err = -1;
      }
    }
  }
  fclose(file);

  if (had_err != -1)
  {
    /* iterate over all vertices */
    for (vertex = graph->vertices;
      vertex < (graph->vertices + graph->nof_vertices); vertex++) {
      if (vertex->astat <= astat_cutoff || vertex->copy_num < copy_num_cutoff)
        vertex->state = GIS_REPEAT;
    }
  }

  return had_err;
}

/* check if unique order of edges <*edge1>, <*edge2> with probability
   <cutoff> exists */
bool
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

GtUword gt_scaffolder_calculate_overlap(GtScaffolderGraphEdge *edge1,
                                               GtScaffolderGraphEdge *edge2)
{
  gt_assert(edge1 != NULL);
  gt_assert(edge2 != NULL);

  GtUword overlap = 0;

  if (edge2->dist > (edge1->end->seq_len - 1) ||
      edge1->dist > (edge2->end->seq_len - 1)) {
    GtWord intersect_start = MAX(edge1->dist, edge2->dist);
    GtWord intersect_end = MAX(edge1->end->seq_len, edge2->end->seq_len);
    overlap = intersect_end - intersect_start;
  }
  return overlap;
}

void
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
    if (poly_vertex->state != GIS_POLYMORPHIC) {
      /* SK: Ueber korrekten Pointer iterieren */
      GtUword eid;
      for (eid = 0; eid < poly_vertex->nof_edges; eid++)
        poly_vertex->edges[eid]->state = GIS_POLYMORPHIC;
      poly_vertex->state = GIS_POLYMORPHIC;
    }
  }
}

/* mark polymorphic edges/vertices and inconsistent edges in scaffold graph */
int gt_scaffolder_graph_filter(GtScaffolderGraph *graph,
                               float pcutoff,
                               float cncutoff,
                               GtUword ocutoff)
{
  gt_assert(graph != NULL);

  GtScaffolderGraphVertex *vertex;
  GtScaffolderGraphEdge *edge1, *edge2;
  GtUword overlap, eid1, eid2;
  GtUword maxoverlap = 0;
  unsigned int dir; /* int statt bool, weil Iteration bislang nicht möglich */
  int had_err = 0;

  /* iterate over all vertices */
  for (vertex = graph->vertices;
       vertex < (graph->vertices + graph->nof_vertices); vertex++) {
    /* iterate over directions (sense/antisense) */
    for (dir = 0; dir < 2; dir++) {
      /* iterate over all pairs of edges */
      for (eid1 = 0; eid1 < vertex->nof_edges; eid1++) {
        for (eid2 = eid1 + 1; eid2 < vertex->nof_edges; eid2++) {
          edge1 = vertex->edges[eid1];
          edge2 = vertex->edges[eid2];
          /* SK: edge->sense == edge->sense pruefen? */
          if (edge1->sense == dir && edge2->sense == dir) {
            /* check if edge1->end and edge2->end are polymorphic */
            gt_scaffolder_graph_check_mark_polymorphic(edge1, edge2,
                                                       pcutoff, cncutoff);
            /* SD: Nur das erste Paar polymoprh markieren? */
          }
        }
      }

      /* no need to check inconsistent edges for polymorphic vertices */
      if (vertex->state == GIS_POLYMORPHIC)
        break;
      /* iterate over all pairs of edges, that are not polymorphic */
      for (eid1 = 0; eid1 < vertex->nof_edges; eid1++) {
        for (eid2 = eid1 + 1; eid2 < vertex->nof_edges; eid2++) {
          edge1 = vertex->edges[eid1];
          edge2 = vertex->edges[eid2];
          if (edge1->sense == dir && edge2->sense == dir &&
              edge1->state != GIS_POLYMORPHIC &&
              edge2->state != GIS_POLYMORPHIC) {
            overlap = gt_scaffolder_calculate_overlap(edge1, edge2);
            if (overlap > maxoverlap)
              maxoverlap = overlap;
          }
        }
      }

     /* check if maxoverlap is larger than ocutoff and mark edges
        as inconsistent */
      if (maxoverlap > ocutoff) {
        for (eid1 = 0; eid1 < vertex->nof_edges; eid1++)
          vertex->edges[eid1]->state = GIS_INCONSISTENT;
      }
    }
  }
  return had_err;
}

/* check if vertex holds just sense or antisense edges */
bool
gt_scaffolder_graph_isterminal(const GtScaffolderGraphVertex *vertex)
{
  gt_assert(vertex != NULL);

  bool dir;
  GtUword eid;

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
                                         GtArray *ccs) {
  GtArray *terminal_vertices = NULL;
  GtQueue *vqueue = NULL;
  GtScaffolderGraphVertex *vertex, *currentvertex, *nextvertex;
  GtUword eid;

  vqueue = gt_queue_new();

  for (vertex = graph->vertices; vertex <
	 (graph->vertices + graph->nof_vertices); vertex++) {
    if (vertex->state != GIS_POLYMORPHIC && vertex->state != GIS_REPEAT)
      vertex->state = GIS_UNVISITED;
  }

  for (vertex = graph->vertices; vertex <
	 (graph->vertices + graph->nof_vertices); vertex++) {
    if (vertex->state == GIS_POLYMORPHIC || vertex->state == GIS_VISITED
	|| vertex->state == GIS_REPEAT)
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
	if (currentvertex->edges[eid]->state != GIS_INCONSISTENT) {
	  nextvertex = currentvertex->edges[eid]->end;
	  /* just take vertices, that are consistent */
	  if (nextvertex->state == GIS_POLYMORPHIC
	      || nextvertex->state == GIS_REPEAT)
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

/*  remove cycles */
/* void gt_scaffolder_removecycles(GtScaffolderGraph *graph) { */

/*     /\* found one cc first remove all cycles starting at terminal */
/*        vertices. keep the array with terminal vertices clean. *\/ */
/*     while (gt_array_size(terminal_vertices) != 0) { */
/*       /\* search for cycles starting at every terminal vertex. if a */
/*       backedge is found, mark both vertices connected with this */
/*       backedge. *\/ */
/*     } */

/* } */

/* DFS to detect Cycles given a starting vertex */
/* TODO: remember all visited vertices to change the state to
   GIS_UNVISITED before search with the next starting vertex. */

/* SD: Commented because compiler complaines about it being unused function */

/*GtScaffolderGraphEdge
*gt_scaffolder_detect_cycle(GtScaffolderGraphVertex *v,
                            bool dir) {
  GtUword eid;
  GtScaffolderGraphVertex *end;
  GtScaffolderGraphEdge *back;

  gt_assert(v != NULL);

  v->state = GIS_VISITED;
  for (eid = 0; eid < v->nof_edges; eid++) {
    if (v->edges[eid]->sense == dir) {*/
      /* maybe we want just to mark the corresponding vertices at this
         point and return a boolean or something like that */ /*
      end = v->edges[eid]->end;
      if (end->state == GIS_VISITED)
        return v->edges[eid];
      if (end->state == GIS_UNVISITED) {
        back = gt_scaffolder_detect_cycle(end, dir);
        if (back != NULL)
          return back;
      }
    }
  }

  end->state = GIS_PROCESSED;
  return NULL;
}*/

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

GtScaffolderGraphWalk *gt_scaffolder_create_walk(GtScaffolderGraph *graph,
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
  GtUword lengthbestwalk, lengthcwalk, eid;
  GtScaffolderGraphWalk *bestwalk, *currentwalk;
  float distance, *distancemap;
  bool dir;
  lengthbestwalk = 0;
  bestwalk = gt_scaffolder_walk_new();

  wqueue = gt_queue_new();
  terminal_vertices = gt_array_new(sizeof (start));

  /* SK: Mit GtWord_Min / erwarteter Genomlaenge statt 0 initialisieren */
  distancemap = calloc(graph->nof_vertices, sizeof (*distancemap));
  edgemap = gt_malloc(sizeof (*edgemap)*graph->nof_vertices);

  /* check if node has edges */
  if (start->nof_edges == 0) {
    return NULL;
  }

  dir = start->edges[0]->sense;
  for (eid = 0; eid < start->nof_edges; eid++) {
    edge = start->edges[eid];
    endvertex = edge->end;

    /* SK: genometools hashes verwenden, Dichte evaluieren
       SK: DistEst beim Einlesen prüfen */
    distancemap[endvertex->index] = edge->dist;
    edgemap[endvertex->index] = edge;

    gt_queue_add(wqueue, edge);
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
        nextendvertex = nextedge->end;
        distance = edge->dist + nextedge->dist;

        /* SK: 0 steht fuer unitialisiert*/
        if (distancemap[nextendvertex->index] == 0 ||
        distancemap[nextendvertex->index] > distance) {
          distancemap[nextendvertex->index] = distance;
          edgemap[nextendvertex->index] = nextedge;
          gt_queue_add(wqueue, nextedge);
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
    while (currentvertex->index != start->index) {
      reverseedge = edgemap[currentvertex->index];
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

  GtScaffolderGraphVertex *vertex, *currentvertex, *nextvertex, *start;
  GtUword eid, max_num_bases;
  GtScaffolderGraphWalk *walk, *bestwalk;
  GtUword ccnumber;
  GtQueue *vqueue;
  GtArray *terminal_vertices, *cc_walks;

  /* Entfernung von Zyklen
     gt_scaffolder_removecycles(graph); */

  /* Iteration ueber alle Knoten, Makierung aller Knoten als unbesucht */
  /* SK: Schleife evtl nich mehr benoetigt, da vor-initialisiert */
  for (vertex = graph->vertices;
       vertex < (graph->vertices + graph->nof_vertices); vertex++) {
    if (vertex->state != GIS_POLYMORPHIC)
      vertex->state = GIS_UNVISITED;
  }

 /* BFS-Traversierung durch Zusammenhangskomponenten des Graphen,
    siehe GraphSearchTree.h */
  ccnumber = 0;
  vqueue = gt_queue_new();
  terminal_vertices = gt_array_new(sizeof (vertex));
  cc_walks = gt_array_new(sizeof (walk));

  for (vertex = graph->vertices; vertex <
       (graph->vertices + graph->nof_vertices); vertex++) {
    if (vertex->state == GIS_POLYMORPHIC || vertex->state == GIS_VISITED)
      continue;
    ccnumber += 1;
    vertex->state = GIS_PROCESSED;
    gt_queue_add(vqueue, vertex);
    gt_array_reset(terminal_vertices);
    gt_array_reset(cc_walks);

    while (gt_queue_size(vqueue) != 0) {
      currentvertex = (GtScaffolderGraphVertex*)gt_queue_get(vqueue);
      /*currentvertex.cc = ccnumber;*/

      /* store all terminal vertices to calculate all paths between them */
      if (gt_scaffolder_graph_isterminal(currentvertex)) {
        gt_assert(currentvertex >= graph->vertices);
        gt_assert(currentvertex < graph->vertices + graph->nof_vertices);
        gt_array_add(terminal_vertices, currentvertex);
        gt_assert(
          ( *(GtScaffolderGraphVertex **) gt_array_get_last(terminal_vertices) )
            == currentvertex );
      }

      currentvertex->state = GIS_VISITED;
      for (eid = 0; eid < currentvertex->nof_edges; eid++) {
        nextvertex = currentvertex->edges[eid]->end;
        /* why vertex->state? */
        if (nextvertex->state == GIS_POLYMORPHIC)
          continue;
        if (nextvertex->state == GIS_UNVISITED) {
          nextvertex->state = GIS_PROCESSED;
          gt_queue_add(vqueue, nextvertex);
        }
      }
    }

    /* calculate all paths between terminal vertices in this cc */
    while (gt_array_size(terminal_vertices) != 0) {
      start = *(GtScaffolderGraphVertex **) gt_array_pop(terminal_vertices);
      gt_assert(start >= graph->vertices);
      gt_assert(start < graph->vertices + graph->nof_vertices);
      walk = gt_scaffolder_create_walk(graph, start);
      gt_array_add(cc_walks, walk);
    }

    /* the best walk in this cc is chosen */
    max_num_bases = 0;
    bestwalk = NULL;
    while (gt_array_size(cc_walks) != 0) {
      walk = *(GtScaffolderGraphWalk **) gt_array_pop(cc_walks);
      if (walk->total_contig_len > max_num_bases) {
	bestwalk = walk;
	max_num_bases = walk->total_contig_len;
      }
    }

    /* mark all nodes and edges in the best walk as GIS_SCAFFOLD */
    bestwalk->edges[0]->start->state = GIS_SCAFFOLD;
    for (eid = 0; eid < bestwalk->nof_edges; eid++) {
      bestwalk->edges[0]->state = GIS_SCAFFOLD;
      bestwalk->edges[0]->end->state = GIS_SCAFFOLD;
    }
  }

  gt_array_delete(terminal_vertices);
  gt_array_delete(cc_walks);
  gt_queue_delete(vqueue);
}
