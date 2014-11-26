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
#include "gt_scaffolder_graph.h"
#include <core/assert_api.h>
#include "core/queue_api.h"
#include "core/types_api.h"
#include "core/minmax.h"
#include "core/fasta_reader_rec.h"
#include "core/error.h"

/* SK: Gt-Namenskonvention für Zustände einhalten (docs/ oder manuals/developermanual)
       Automatische Prüfung durch scripts/src_check */
typedef enum { GIS_UNVISITED, GIS_POLYMORPHIC, GIS_INCONSISTENT,
               GIS_VISITED, GIS_PROCESSED } GraphItemState;

/* vertex of scaffold graph (describes one contig) */
typedef struct GtScaffoldGraphVertex {
  /* unique vertex ID */
  GtUword id;
  /* header sequence of corresponding contig */
  char *headerseq;
  /* sequence length of corresponding contig */
  GtUword seqlen;
  /* a-statistics value for classifying contig as repeat or unique contig */
  float astat;
  /* estimated copy number of corresponding contig */
  float copynum;
  GtUword nofedges;
  struct GtScaffoldGraphEdge **edges;
  /* vertex state (vertex can adapt every state except GIS_INCONSISTENT) */
  GraphItemState state;
} GtScaffoldGraphVertex;

/* edge of scaffold graph (describes orientation of two contigs) */
typedef struct GtScaffoldGraphEdge {
  /* pointer to end vertex of edge */
  GtScaffoldGraphVertex *end;
  /* pointer to start vertex of edge */
  GtScaffoldGraphVertex *start;
  /* estimated distance between contigs of start and end vertex */
  GtWord dist;
  /* standard deviation of estimated distance */
  float stddev;
  /* number of read pairs resulting that distance */
  GtUword numpairs;
  /* edge state */
  GraphItemState state;
  /* describes direction of corresponding contigs
     sense = true & same = true: ctg1 & ctg2 in sense direction
     sense = true & same = false: ctg1 in sense & ctg2 in antisense direction
     sense = false & same = true: ctg1 & ctg2 in antisense direction
     sense = false & same = false: ctg1 in antisense & ctg2 in sense direction */
  bool sense;
  bool same;
} GtScaffoldGraphEdge;

/* scaffold graph */
struct GtScaffoldGraph {
  GtScaffoldGraphVertex *vertices;
  GtUword nofvertices;
  GtUword maxnofvertices;
  GtScaffoldGraphEdge *edges;
  GtUword nofedges;
  GtUword maxnofedges;
};

/* linear scaffold */
typedef struct GtScaffoldGraphWalk {
  GtUword nofedges;
  GtUword size;
  GtUword totalcontiglen;
  GtScaffoldGraphEdge **edges;
}GtScaffoldGraphWalk;

/* for parsing valid contigs */
typedef struct GtScaffoldValidCtg {
  GtUword nof;
  GtUword minctglen;
  GtUword headerlen;
  const char *headerseq;
  GtScaffoldGraph *graph;
}GtScaffoldValidCtg;

/* DistanceMap */
/* EdgeMap */

/* TODO: Funktionen gt_scaffolder_graph_create_vertices,
         gt_scaffolder_graph_create_edges erstellen */
/* Return pointer to data structure <*GtScaffoldGraph>. Allocate space for
   <totalnofvertices> vertices and <totalnofedges> edges. Initialize current
   number of vertices and edges with 0. */
GtScaffoldGraph *gt_scaffolder_graph_new(GtUword totalnofvertices, GtUword totalnofedges) {
  GtScaffoldGraph *graph;

  gt_assert(totalnofvertices > 0);
  gt_assert(totalnofedges > 0);

  graph = gt_malloc(sizeof(*graph));
  graph->vertices = gt_malloc(sizeof(*graph->vertices) * totalnofvertices);
  graph->edges = gt_malloc(sizeof(*graph->edges) * totalnofedges);
  graph->nofvertices = 0;
  graph->maxnofvertices = totalnofvertices;
  graph->nofedges = 0;
  graph->maxnofedges = totalnofedges;

  return graph;
}

void gt_scaffolder_graph_delete(GtScaffoldGraph *graph){
  gt_assert(graph != NULL);

  if (graph->vertices != NULL){
    gt_free(graph->vertices->headerseq);
    gt_free(graph->vertices->edges[0]);
    gt_free(graph->vertices->edges);
  }
  if (graph->edges != NULL){
    gt_free(graph->edges->end);
    gt_free(graph->edges->start);
  }
  gt_free(graph->edges);
  gt_free(graph->vertices);
  gt_free(graph);
}

/* Initialize a new vertex in <*graph>. Each vertex represents a contig and
   contains information about the sequence length <seqlen>, A-statistics
   <astat> and estimated copy number <copynum> */
void graph_add_vertex(GtScaffoldGraph *graph, GtUword seqlen, float astat,
  float copynum) {
  gt_assert(graph != NULL);
  gt_assert(graph->nofvertices < graph->maxnofvertices);

  GtUword nextfree = graph->nofvertices;

  /* Initialize vertex */
  graph->vertices[nextfree].id = nextfree; /* SD: redundant, needed for compiling */
  graph->vertices[nextfree].seqlen = seqlen;
  graph->vertices[nextfree].astat = astat;
  graph->vertices[nextfree].copynum = copynum;
  graph->vertices[nextfree].nofedges = 0;
  /* SD: Initialize state? graph->vertices[nextfree].state = GIS_UNVISITED; */

  /* Allocate initial space for pointer to outgoing edges */
  graph->vertices[nextfree].edges = gt_malloc(sizeof(*graph->edges));

  graph->nofvertices++;
}

/* Initialize a new, directed edge in <*graph>. Each edge represents a contig
   and contains information about the sequence length <seqlen>, A-statistics
   <astat> and estimated copy number <copynum> */
void graph_add_edge(GtScaffoldGraph *graph, GtUword vstartID, GtUword vendID,
  GtWord dist, float stddev, GtUword numpairs, bool dir, bool comp) {

  gt_assert(graph != NULL);
  gt_assert(graph->nofedges < graph->maxnofedges);

  GtUword nextfree = graph->nofedges;

  /* Inititalize edge */
  graph->edges[nextfree].start = graph->vertices + vstartID;
  graph->edges[nextfree].end = graph->vertices + vendID;
  graph->edges[nextfree].dist = dist;
  graph->edges[nextfree].stddev = stddev;
  graph->edges[nextfree].numpairs = numpairs;
  graph->edges[nextfree].sense = dir;
  graph->edges[nextfree].same = comp;

  /* Assign edge to start vertice */
  gt_assert(vstartID < graph->nofvertices && vendID < graph->nofvertices);
  /* Allocate new space for pointer to this edge */
  if(graph->vertices[vstartID].nofedges > 0) {
    graph->vertices[vstartID].edges =
      /* SK: realloc zu teuer? Besser: DistEst parsen und gezielt allokieren */
      gt_realloc( graph->vertices[vstartID].edges, sizeof(*graph->edges) *
                  (graph->vertices[vstartID].nofedges + 1) );
  }
  /* Assign adress of this edge to the pointer */
  graph->vertices[vstartID].edges[graph->vertices[vstartID].nofedges] =
    &graph->edges[nextfree];

  graph->vertices[vstartID].nofedges++;

  graph->nofedges++;
}

static GtScaffoldGraphEdge *graph_find_edge(GtScaffoldGraph *graph,
  GtUword vertexid_1, GtUword vertexid_2)
{
  GtScaffoldGraphEdge *edge;

  for (edge = (graph->vertices + vertexid_1)->edges[0];
       edge < ((graph->vertices + vertexid_1)->edges[0] +
               (graph->vertices + vertexid_1)->nofedges); edge++) {
    if (edge->end->id == vertexid_2)
      return edge;
  }
  return NULL;
}

static GtUword graph_get_vertex_id(GtScaffoldGraph *graph, const char* headerseq){
  GtUword vid;

  for (vid = 0; vid < graph->nofvertices; vid++){
    /*TODO: buchstabeweisen Vergleich einfuegen!!*/
    if (graph->vertices[vid].headerseq == headerseq)
      break;
  }

  return vid;

}

/* print graphrepresentation in dot-format into file filename */
int gt_scaffolder_graph_print(const GtScaffoldGraph *g,
                              const char *filename, GtError *err) {
  int had_err = 0;

  GtFile *f = gt_file_new(filename, "w", err);
  if (f == NULL)
    had_err = 1;

  if (!had_err) {
    gt_scaffolder_graph_print_generic(g, f);
    gt_file_delete(f);
  }

  return had_err;
}

/* print graphrepresentation in dot-format into gt-filestream f */
void gt_scaffolder_graph_print_generic(const GtScaffoldGraph *g,
                                       GtFile *f) {
  GtScaffoldGraphVertex *v;
  GtScaffoldGraphEdge *e;
  /* 0: GIS_UNVISITED, 1: GIS_POLYMORPHIC, 2: GIS_INCONSISTENT,
     3: GIS_VISITED, 4: GIS_PROCESSED */
  const char *color_array[5] = {"black", "gray", "gray", "red", "green"};

  /* print first line into f */
  gt_file_xprintf(f, "graph {\n");

  /* iterate over all vertices and print them. add attribute color according
     to the current state */
  for (v = g->vertices; v < (g->vertices + g->nofvertices); v++) {
    gt_file_xprintf(f, "%lu [color=\"%s\"];\n", v->id, color_array[v->state]);
  }

  /* iterate over all edges and print them. add attribute color according to
     the current state and label the edge with the distance*/
  for (e = g->edges; e < (g->edges + g->nofedges); e++) {
    gt_file_xprintf(f, "%lu -- %lu [color=\"%s\" label=\"%ld\"];\n",
                    e->start->id, e->end->id,
                    color_array[e->state], e->dist);
  }

  /* print the last line into f */
  gt_file_xprintf(f, "}\n");
}

/* parse distance information of contigs in abyss-dist-format and
   save them as edges of scaffold graph */
/* check for "mate-flag", composition is missing! */
static int gt_scaffolder_graph_read_distances(const char *filename,
  GtScaffoldGraph *graph, bool ismatepair, GtError *err)
{
  FILE *infile;
  char *buffer, *c, *field;
  GtUword pos, bufferlen, result, numpairs, num5, fieldsize;
  GtUword rootctgid, ctgid;
  GtScaffoldGraphEdge *edge = NULL;
  GtWord dist;
  float stddev;
  bool firstfield, nextfirstfield, curdir;
  int had_err;

  had_err = 0;
  infile = fopen(filename, "rb");
  if (infile == NULL)
  {
    gt_error_set (err , " invalid file %s ", filename);
    had_err = -1;
  }

  /* determine size of file */
  fseek(infile , 0 , SEEK_END);
  bufferlen = ftell (infile);
  rewind (infile);
  /* save content of file */
  buffer = gt_malloc(sizeof(*buffer) * bufferlen);
  result = fread(buffer, 1, bufferlen, infile);
  fieldsize = 64;
  field = gt_malloc(sizeof(*field)*fieldsize);

  if (result != bufferlen)
  {
    gt_error_set (err , " incomplete read of file %s ", filename);
    had_err = -1;
  }

  pos = 0;
  firstfield = true;
  nextfirstfield = true;
  /* sense direction as default */
  curdir = true;

  for (c = buffer; c != buffer + bufferlen; c++)
  {
    field[pos] = *c;
    pos++;

    if (fieldsize == pos)
    {
      fieldsize *= 2;
      field = gt_realloc(field, sizeof(*field)*fieldsize);
    }

    if (*c == ' ' || *c == '\n')
    {
      field[pos-1] = '\0';
      pos = 0;

      /* parse root contig ID */
      if (firstfield)
      {
        firstfield = false;
        rootctgid = graph_get_vertex_id(graph, field);
        /* Debbuging:
          printf("rootctgid: %s\n",field);*/
      }
      else
      {
        /* parse contig-record with attributes distance, number of pairs,
           standard deviation of distance */
        if (sscanf(field,"%ld,%lu,%f,%lu",&dist,&numpairs,&stddev,&num5) == 4)
        {
          /* check if edge between vertices already exixts */
          edge = graph_find_edge(graph, rootctgid, ctgid);
          if (edge != NULL)
          {
            /*  LG: laut SGA edge->stddev < stddev,  korrekt? */
            if (ismatepair == false && edge->stddev < stddev)
            {
              /* TODO alter_edge einfuegen*/
              /* LG: Ueberpruefung Kantenrichtung notwendig? */
              graph_add_edge(graph, rootctgid, ctgid, dist, stddev, numpairs,
                                     curdir, true);
            }
            else
            {
            /*Conflicting-Flag?*/
            }
          }
          else
            graph_add_edge(graph, rootctgid, ctgid, dist, stddev, numpairs,
                                   curdir, true);
          /* Debbuging:
             printf("dist: %ld\n numpairs: %lu\n stddev:"
                 "%f\n num5: %lu\n curdir: %d\n\n",dist, numpairs, stddev,
                 num5,curdir);*/

        }
        /* switch direction if semicolon occurs */
        else if (*field == ';')
          curdir = !curdir;
        nextfirstfield = true;
      }
    }
    /* parse contig ID */
    else if (*c == ',' && nextfirstfield)
    {
      nextfirstfield = false;
      field[pos-1] = '\0';
      pos = 0;
      /* TODO parsing composition */
      ctgid = graph_get_vertex_id(graph, field);
      /* Debbuging:
         printf("ctgid: %s\n",field);*/
    }
    if (*c == '\n')
    {
      firstfield = true;
      /* sense direction as default */
      curdir = true;
    }
  }
  return had_err;
}



static int gt_scaffolder_graph_count_ctg(GtUword length, void *data, GtError* err)
{
  int had_err;
  GtScaffoldValidCtg *validctg = (GtScaffoldValidCtg*) data;

  had_err = 0;
  if (length >= validctg->minctglen)
    validctg->nof++;
  if (length == 0)
  {
    gt_error_set (err , " invalid sequence length ");
    had_err = -1;
  }
  return had_err;
}

static int gt_scaffolder_graph_save_header(const char *description, GtUword length,
 void *data, GtError *err)
{
  int had_err;
  GtScaffoldValidCtg *validctg = (GtScaffoldValidCtg*) data;

  had_err = 0;
  validctg->headerseq = description;
  validctg->headerlen = length;
  if (length == 0)
  {
    gt_error_set (err , " invalid header length ");
    had_err = -1;
  }
  return had_err;
}

static int gt_scaffolder_graph_save_ctg(GtUword length, void *data, GtError* err)
{
  int had_err;
  GtScaffoldValidCtg *validctg = (GtScaffoldValidCtg*) data;

  had_err = 0;
  if (length > validctg->minctglen)
  {
    validctg->graph->vertices[validctg->graph->nofvertices].headerseq =
    gt_malloc(sizeof(char) * validctg->headerlen);
    strncpy(validctg->graph->vertices[validctg->graph->nofvertices].headerseq,
    validctg->headerseq, validctg->headerlen);
    graph_add_vertex(validctg->graph, length,0.0,0.0);
    /*validctg->graph->vertices[validctg->graph->nofvertices].seqlen = length;*/
    validctg->graph->nofvertices++;
  }
  if (length == 0)
  {
    gt_error_set (err , " invalid sequence length ");
    had_err = -1;
  }
  return had_err;
}

/* creation of scaffold graph */
GtScaffoldGraph *gt_scaffolder_graph_new_from_file(const char *ctgfilename,
                 GtUword minctglen, const char *distfilename, GtError *err)
{
  GtScaffoldGraph *graph;
  GtFastaReader* reader;
  GtStr *str_filename;
  int had_err;
  GtScaffoldValidCtg *validctg;

  had_err = 0;
  str_filename = gt_str_new();
  gt_str_set(str_filename, ctgfilename);
  /*TODO andere Name fuer validctg? callback_data*/
  validctg = gt_malloc(sizeof(*validctg));
  validctg->nof = 0;
  validctg->minctglen = minctglen;
  /* parse contigs in FASTA-format and save them as vertices of
     scaffold graph */
  reader = gt_fasta_reader_rec_new(str_filename);
  had_err = gt_fasta_reader_run(reader, NULL, NULL,
            gt_scaffolder_graph_count_ctg, validctg, err);
  /*TODO graph_new Funktion einfuegen */
  graph = gt_malloc(sizeof(*graph));
  graph->nofvertices = 0;
  graph->maxnofvertices = validctg->nof;
  graph->vertices = gt_malloc(sizeof(*graph->vertices) * validctg->nof);
  validctg->graph = graph;
  had_err = gt_fasta_reader_run(reader, gt_scaffolder_graph_save_header, NULL,
            gt_scaffolder_graph_save_ctg, validctg, err);
  gt_fasta_reader_delete(reader);
  gt_str_delete(str_filename);

  /* parse distance information of contigs in abyss-dist-format and
     save them as edges of scaffold graph */
  had_err = gt_scaffolder_graph_read_distances(distfilename, graph, false, err);

  if (had_err != 0)
    printf("error");

  return graph;
}



/* Pruefung auf eindeutige Ordnung der Kanten edge1, edge 2 mit Wahrscheinlichkeit
   cutoff */
static bool gt_scaffolder_graph_ambiguousorder(const GtScaffoldGraphEdge *edge1,
      const GtScaffoldGraphEdge *edge2, float cutoff) {

  float expval, variance, interval, prob12, prob21;

  expval = edge1->dist - edge2->dist;
  /* sichere Multiplikation, Division */
  /* SK: pow nicht verwenden, stattdessen mit sich selbst multiplizieren */
  variance = 2 * pow(edge1->stddev,2) - pow(edge2->stddev,2);
  interval = -expval / sqrt(variance);
  prob12 = 0.5 * (1 + erf(interval) );
  prob21 = 1.0 - prob12;

  return (prob12 <= cutoff && prob21 <= cutoff);
}

static GtUword calculate_overlap (GtScaffoldGraphEdge *edge1,
                                  GtScaffoldGraphEdge *edge2) {
  GtUword overlap = 0;
  GtWord intersect_start, intersect_end;
  if (edge2->dist > (edge1->end->seqlen - 1) ||
      edge1->dist > (edge2->end->seqlen - 1)) {
    intersect_start = MAX(edge1->dist, edge2->dist);
    intersect_end = MAX(edge1->end->seqlen - 1, edge2->end->seqlen - 1);
    overlap = intersect_end - intersect_start + 1;
  }
  return overlap;
}

/* Makierung polymorpher Kanten/Knoten und inkonsistenter Kanten im
   Scaffold Graphen */
int gt_scaffolder_graph_filtering(GtScaffoldGraph *graph, float pcutoff,
    float cncutoff, GtUword ocutoff) {
  GtScaffoldGraphVertex *vertex, *polymorphic_vertex;
  GtScaffoldGraphEdge *edge1, *edge2;
  GtUword eid_3, overlap;
  GtUword maxoverlap = 0;
  float sum_copynum;
  unsigned int dir; /* int statt bool, weil Iteration bislang nicht möglich */
  int had_err = 0;

  /* iterate over all vertices */
  for (vertex = graph->vertices;
       vertex < (graph->vertices + graph->nofvertices); vertex++) {
    /* iterate over directions (sense/antisense) */
    for (dir = 0; dir < 2; dir++) {
      /* iterate over all pairs of edges */
      for (edge1 = vertex->edges[0];
           edge1 < (vertex->edges[0] + vertex->nofedges); edge1++) {
        for (edge2 = edge1 + 1; edge2 < (vertex->edges[0] + vertex->nofedges);
             edge2++) {
          if (edge1->sense == dir && edge2->sense == dir) {
            /* check if edge1->end and edge2->end are polymorphic */
            sum_copynum = edge1->end->copynum + edge2->end->copynum;
            if (gt_scaffolder_graph_ambiguousorder(edge1, edge2, pcutoff) &&
                sum_copynum < cncutoff) {
              /* mark vertex with lower copynumber as polymorphic */
              if (edge1->end->copynum < edge2->end->copynum)
                polymorphic_vertex = edge1->end;
              else
                polymorphic_vertex = edge2->end;
              /* mark all edges of the polymorphic vertex as polymorphic */
              if (polymorphic_vertex->state != GIS_POLYMORPHIC) {
                for (eid_3 = 0; eid_3 < polymorphic_vertex->nofedges; eid_3++)
                  polymorphic_vertex->edges[eid_3]->state = GIS_POLYMORPHIC;
                polymorphic_vertex->state = GIS_POLYMORPHIC;
              }
            }
            /* SD: Nur das erste Paar polymoprh markieren? */
          }
        }
      }

      /* no need to check inconsistent edges for polymorphic vertices */
      if (vertex->state == GIS_POLYMORPHIC)
        break;
      /* iterate over all pairs of edges, that are not polymorphic */
      for (edge1 = vertex->edges[0];
           edge1 < (vertex->edges[0] + vertex->nofedges); edge1++) {
        for (edge2 = edge1 + 1; edge2 < (vertex->edges[0] + vertex->nofedges);
             edge2++) {
          if (edge1->sense == dir && edge2->sense == dir &&
              edge1->state != GIS_POLYMORPHIC && edge2->state != GIS_POLYMORPHIC) {

            overlap = calculate_overlap(edge1, edge2);
            if (overlap > maxoverlap)
              maxoverlap = overlap;
          }
        }
      }

     /* check if maxoverlap is larger than ocutoff and mark edges
        as inconsistent */
      if (maxoverlap > ocutoff) {
        for (edge1 = vertex->edges[0];
             edge1 < (vertex->edges[0] + vertex->nofedges); edge1++)
          edge1->state = GIS_INCONSISTENT;
      }
    }
  }
  /* SD: Knoten & Kanten werden nicht gelöscht, sondern der Einfachheit halber später geprüft*/
  return had_err;
}

/* check if vertex holds just sense or antisense edges */
static bool gt_scaffolder_graph_isterminal(const GtScaffoldGraphVertex *vertex) {
  GtScaffoldGraphEdge *edge;
  bool dir;

  dir = vertex->edges[0]->sense;
  for (edge = (vertex->edges[0] + 1); edge < (vertex->edges[0] + vertex->nofedges);
       edge++) {
    if (edge->sense != dir)
      return false;
  }

  return true;
}

/* Entfernung von Zyklen
static void gt_scaffolder_removecycles(GtScaffoldGraph *graph) {
}*/

/* Erstellung eines neuen Walks */
static GtScaffoldGraphWalk *gt_scaffolder_walk_new(void) {
  GtScaffoldGraphWalk *walk;

  walk = gt_malloc(sizeof(*walk));
  walk->totalcontiglen = 0;
  walk->size = 0;
  walk->nofedges = 0;
  return walk;
}

/* Loeschen eines Walks */
static void gt_scaffolder_walk_delete(GtScaffoldGraphWalk *walk) {
  if (walk != NULL)
    gt_free(walk->edges);
  gt_free(walk);
}

/* Ausgabe der Contig-Gesamtlaenge eines Walks */
static GtUword gt_scaffolder_walk_getlength(GtScaffoldGraphWalk *walk) {
  return walk->totalcontiglen;
}

/* Hinzufuegen einer Kante zum Walk */
static void gt_scaffolder_walk_addegde(GtScaffoldGraphWalk *walk, GtScaffoldGraphEdge *edge) {
  if (walk->size == walk->nofedges) {
    /* SK: 10 als Konstante definieren, const unsigned long increment_size 2er Potenz */
    walk->size += 10;
    walk->edges = gt_realloc(walk->edges, walk->size*sizeof(*walk->edges));
  }
  walk->edges[walk->nofedges] = edge;
  walk->totalcontiglen += edge->end->seqlen;
  walk->nofedges++;
}

/* Konstruktion des Scaffolds mit groesster Contig-Gesamtlaenge */
void gt_scaffolder_makescaffold(GtScaffoldGraph *graph) {
  GtScaffoldGraphVertex *vertex, *currentvertex, *nextvertex, *nextendvertex,
                        *endvertex;
  GtScaffoldGraphEdge *edge, *nextedge, *reverseedge, **edgemap;
  GtUword ccnumber, lengthcwalk, lengthbestwalk;
  GtQueue *vqueue, *wqueue;
  float distance, *distancemap;
  bool dir;
  GtScaffoldGraphWalk *bestwalk, *currentwalk;

  /* Entfernung von Zyklen
     gt_scaffolder_removecycles(graph); */

  /* Iteration ueber alle Knoten, Makierung aller Knoten als unbesucht */
  for (vertex = graph->vertices;
       vertex < (graph->vertices + graph->nofvertices); vertex++) {
    if (vertex->state != GIS_POLYMORPHIC)
      vertex->state = GIS_UNVISITED;
  }

 /* BFS-Traversierung durch Zusammenhangskomponenten des Graphen,
    siehe GraphSearchTree.h */
  ccnumber = 0;
  vqueue = gt_queue_new();
  /* SK: Mit GtWord_Min / erwarteter Genomlaenge statt 0 initialisieren */
  distancemap = calloc(graph->nofvertices, sizeof(*distancemap));
  edgemap = gt_malloc(sizeof(*edgemap)*graph->nofvertices);

  for (vertex = graph->vertices; vertex < (graph->vertices + graph->nofvertices); vertex++) {
    if (vertex->state == GIS_POLYMORPHIC || vertex->state == GIS_VISITED)
      continue;
    ccnumber += 1;
    vertex->state = GIS_PROCESSED;
    gt_queue_add(vqueue, vertex); // TODO: hier moechte evtl nicht die adresse von dem pointer uebergeben werden

    while (gt_queue_size(vqueue) != 0) {
      currentvertex = (GtScaffoldGraphVertex*)gt_queue_get(vqueue);
      //currentvertex.cc = ccnumber;

      /* BFS-Traversierung innerhalb aktueller Zusammenhangskomponente
         ausgehend von terminalen Knoten zu terminalen Knoten */
      lengthbestwalk = 0;
      bestwalk = gt_scaffolder_walk_new();
      /* SK: In weitere Funktion mit eigener Queue auslagern */
      wqueue = gt_queue_new();

      if (gt_scaffolder_graph_isterminal(vertex)) {
        dir = vertex->edges[0]->sense;
        for (edge = vertex->edges[0];
             edge < (vertex->edges[0] + vertex->nofedges); edge++) {
          endvertex = edge->end;

          /* SK: genometools hashes verwenden, Dichte evaluieren
             SK: DistEst beim Einlesen prüfen */
          distancemap[endvertex->id] = edge->dist;
          edgemap[endvertex->id] = edge;

          gt_queue_add(wqueue, edge);
        }
        while(gt_queue_size(wqueue) != 0) {
          edge = (GtScaffoldGraphEdge*)gt_queue_get(wqueue);
          endvertex = edge->end;

          /* Ruecktraversierung durch EdgeMap wenn terminaler Knoten erreicht,
             Konstruktion des Walks  */
          if (gt_scaffolder_graph_isterminal(endvertex)) {

            currentwalk = gt_scaffolder_walk_new();
            currentvertex = endvertex;
            while (currentvertex->id != vertex->id) {
              reverseedge = edgemap[currentvertex->id];
             /* Start NICHT end */
              currentvertex = reverseedge->end;
             /* Speicherung des aktuellen Walks */
             /* Kante vorher duplizieren */
              gt_scaffolder_walk_addegde(currentwalk, reverseedge);
            }

            /* Ermittelung Contig-Gesamtlaenge des aktuellen Walks  */
            lengthcwalk = gt_scaffolder_walk_getlength(currentwalk);
            if ( lengthcwalk > lengthbestwalk) {
              gt_scaffolder_walk_delete(bestwalk);
              bestwalk = currentwalk;
              lengthbestwalk = lengthcwalk;
            }
            else
              gt_scaffolder_walk_delete(currentwalk);
          }

          /* SD: Terminal Set implementieren, bestWalk über Rücktraversierung */
          for (nextedge = endvertex->edges[0];
               nextedge < (endvertex->edges[0] + endvertex->nofedges);
               nextedge++) {
            if (nextedge->sense == dir) {
              nextendvertex = nextedge->end;
              distance = edge->dist + nextedge->dist;

              if (distancemap[nextendvertex->id] == 0 ||
                  distancemap[nextendvertex->id] > distance) {
                distancemap[nextendvertex->id] = distance;
                edgemap[nextendvertex->id] = nextedge;
                gt_queue_add(wqueue, nextedge);
              }
            }
          }
        }
      }


      currentvertex->state = GIS_VISITED;
      for (edge = currentvertex->edges[0];
           edge < (currentvertex->edges[0] + currentvertex->nofedges);
           edge++) {
        nextvertex = edge->end;
        if (vertex->state == GIS_POLYMORPHIC)
          continue;
        if (nextvertex->state == GIS_UNVISITED) {
          nextvertex->state = GIS_PROCESSED;
          gt_queue_add(vqueue, nextvertex);
        }
      }
    }
  }
}



