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

const GtUword INCREMENT_SIZE = 32;

typedef enum { GIS_UNVISITED, GIS_POLYMORPHIC, GIS_INCONSISTENT,
               GIS_VISITED, GIS_PROCESSED } GraphItemState;

/* vertex of scaffold graph (describes one contig) */
typedef struct GtScaffoldGraphVertex {
  /* unique vertex ID */
  GtUword id;
  /* header sequence of corresponding contig */
  char *header_seq;
  /* sequence length of corresponding contig */
  GtUword seq_len;
  /* a-statistics value for classifying contig as repeat or unique contig */
  float astat;
  /* estimated copy number of corresponding contig */
  float copy_num;
  GtUword nof_edges;
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
} GtScaffoldGraphEdge;

/* scaffold graph */
struct GtScaffoldGraph {
  GtScaffoldGraphVertex *vertices;
  GtUword nof_vertices;
  GtUword max_nof_vertices;
  GtScaffoldGraphEdge *edges;
  GtUword nof_edges;
  GtUword max_nof_edges;
};

/* linear scaffold */
typedef struct GtScaffoldGraphWalk {
  GtUword nof_edges;
  GtUword size;
  GtUword total_contig_len;
  GtScaffoldGraphEdge **edges;
}GtScaffoldGraphWalk;

/* for parsing valid contigs,
   e.g. contigs with minimum length <min_ctg_len> */
typedef struct GtScaffoldGraphCallbackData {
  GtUword nof_valid_ctg;
  GtUword min_ctg_len;
  GtUword header_len;
  const char *header_seq;
  GtScaffoldGraph *graph;
} GtScaffoldGraphCallbackData;

/* DistanceMap */
/* EdgeMap */

/* Initialize vertex portion inside <*graph>. Allocate memory for
   <max_nof_vertices> vertices. */
void gt_scaffolder_graph_create_vertices(GtScaffoldGraph *graph,
                                         GtUword max_nof_vertices)
{
  gt_assert(graph != 0);
  gt_assert(graph->vertices == NULL);
  gt_assert(max_nof_vertices > 0);
  graph->vertices = gt_malloc(sizeof(*graph->vertices) * max_nof_vertices);
  graph->nof_vertices = 0;
  graph->max_nof_vertices = max_nof_vertices;
}

/* Initialize edge portion inside <*graph>. Allocate memory for
   <max_nof_edges> edges. */
void gt_scaffolder_graph_create_edges(GtScaffoldGraph *graph,
                                      GtUword max_nof_edges)
{
  gt_assert(graph != 0);
  gt_assert(graph->edges == NULL);
  gt_assert(max_nof_edges > 0);
  graph->edges = gt_malloc(sizeof(*graph->edges) * max_nof_edges);
  graph->nof_edges = 0;
  graph->max_nof_edges = max_nof_edges;
}

/* Construct graph data structure <*GtScaffoldGraph>. Wraps around two seperate
   constructor functions, which allocate memory for <max_nof_edges> edges and
   <maxnoefvertices> vertices. */
GtScaffoldGraph *gt_scaffolder_graph_new(GtUword max_nof_vertices,
                                         GtUword max_nof_edges)
{
  GtScaffoldGraph *graph;

  graph = gt_malloc(sizeof(*graph));
  gt_scaffolder_graph_create_vertices(graph, max_nof_vertices);
  gt_scaffolder_graph_create_edges(graph, max_nof_edges);

  return graph;
}

void gt_scaffolder_graph_delete(GtScaffoldGraph *graph)
{
  gt_assert(graph != NULL);

  if (graph->vertices != NULL) {
    gt_free(graph->vertices->header_seq);
    gt_free(graph->vertices->edges[0]);
    gt_free(graph->vertices->edges);
  }

  if (graph->edges != NULL) {
    gt_free(graph->edges->end);
    gt_free(graph->edges->start);
  }

  gt_free(graph->edges);
  gt_free(graph->vertices);
  gt_free(graph);
}

/* Initialize a new vertex in <*graph>. Each vertex represents a contig and
   contains information about the sequence header <*header_seq>, sequence
   length <seq_len>, A-statistics <astat> and estimated copy number <copy_num>*/
void gt_scaffolder_graph_add_vertex(GtScaffoldGraph *graph,
                                    const char *header_seq,
                                    GtUword seq_len,
                                    float astat,
                                    float copy_num)
{
  gt_assert(graph != NULL);
  gt_assert(graph->nof_vertices < graph->max_nof_vertices);

  GtUword nextfree = graph->nof_vertices;

  /* Initialize vertex */
  graph->vertices[nextfree].id = nextfree; /* SD: remove without breaking */
  graph->vertices[nextfree].seq_len = seq_len;
  graph->vertices[nextfree].astat = astat;
  graph->vertices[nextfree].copy_num = copy_num;
  graph->vertices[nextfree].nof_edges = 0;
  if (header_seq != NULL) {
    graph->vertices[nextfree].header_seq =
      gt_malloc( sizeof (char) * strlen(header_seq) );
    strncpy( graph->vertices[nextfree].header_seq, header_seq,
             strlen(header_seq) );
  }
  /* SD: Initialize state? graph->vertices[nextfree].state = GIS_UNVISITED; */

  /* Allocate initial space for pointer to outgoing edges */
  graph->vertices[nextfree].edges = gt_malloc(sizeof(*graph->edges));

  graph->nof_vertices++;
}

/* Initialize a new, directed edge in <*graph>. Each edge between two contig
   vertices <vstartID> and <vendID> contains information about the distance
   <dist>, standard deviation <std_dev>, number of pairs <num_pairs> and the
   direction of <vstartID> <dir> and corresponding <vendID> <same> */
void gt_scaffolder_graph_add_edge(GtScaffoldGraph *graph,
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

  /* Assign edge to start vertice */
  gt_assert(vstartID < graph->nof_vertices && vendID < graph->nof_vertices);
  /* Allocate new space for pointer to this edge */
  if (graph->vertices[vstartID].nof_edges > 0) {
    graph->vertices[vstartID].edges =
      /* SK: realloc zu teuer? Besser: DistEst parsen und gezielt allokieren */
      gt_realloc( graph->vertices[vstartID].edges, sizeof (*graph->edges) *
                  (graph->vertices[vstartID].nof_edges + 1) );
  }
  /* Assign adress of this edge to the pointer */
  graph->vertices[vstartID].edges[graph->vertices[vstartID].nof_edges] =
    &graph->edges[nextfree];

  graph->vertices[vstartID].nof_edges++;

  graph->nof_edges++;
}

static GtScaffoldGraphEdge *gt_scaffolder_graph_find_edge(
                                    GtScaffoldGraph *graph,
                                    GtUword vertexid_1,
                                    GtUword vertexid_2)
{
  GtScaffoldGraphEdge *edge;

  for (edge = (graph->vertices + vertexid_1)->edges[0];
       edge < ((graph->vertices + vertexid_1)->edges[0] +
          (graph->vertices + vertexid_1)->nof_edges); edge++) {
    if (edge->end->id == vertexid_2)
      return edge;
  }
  return NULL;
}

/* determines corresponding vertex id to contig header */
static GtUword gt_scaffolder_graph_get_vertex_id(GtScaffoldGraph *graph,
                                             const char* header_seq,
                                             GtError *err)
{
  GtScaffoldGraphVertex *vertex;
  GtUword had_err;

  had_err = -1;
  for (vertex = graph->vertices;
       vertex < (graph->vertices + graph->nof_vertices); vertex++) {
    if (strcmp(vertex->header_seq, header_seq) == 0) {
      had_err = 0;
      break;
    }
  }

  /* contig header was not found */
  if (had_err == -1)
    gt_error_set(err, " distance and contig file inconsistent ");

  return vertex->id;
}

/* assign edge <*edge> new attributes */
static void gt_scaffolder_graph_alter_edge(GtScaffoldGraphEdge *edge,
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
int gt_scaffolder_graph_print(const GtScaffoldGraph *g,
                              const char *filename,
                              GtError *err)
{
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
                                       GtFile *f)
{
  GtScaffoldGraphVertex *v;
  GtScaffoldGraphEdge *e;
  /* 0: GIS_UNVISITED, 1: GIS_POLYMORPHIC, 2: GIS_INCONSISTENT,
     3: GIS_VISITED, 4: GIS_PROCESSED */
  const char *color_array[5] = {"black", "gray", "gray", "red", "green"};

  /* print first line into f */
  gt_file_xprintf(f, "graph {\n");

  /* iterate over all vertices and print them. add attribute color according
     to the current state */
  for (v = g->vertices; v < (g->vertices + g->nof_vertices); v++) {
    gt_file_xprintf(f, "%lu [color=\"%s\"];\n", v->id, color_array[v->state]);
  }

  /* iterate over all edges and print them. add attribute color according to
     the current state and label the edge with the distance*/
  for (e = g->edges; e < (g->edges + g->nof_edges); e++) {
    gt_file_xprintf(f, "%lu -- %lu [color=\"%s\" label=\"%ld\"];\n",
                    e->start->id, e->end->id,
                    color_array[e->state], e->dist);
  }

  /* print the last line into f */
  gt_file_xprintf(f, "}\n");
}

/* parse distance information of contigs in abyss-dist-format and
   save them as edges of scaffold graph */
/* LG: check for "mate-flag"? */
static int gt_scaffolder_graph_read_distances(const char *filename,
                                              GtScaffoldGraph *graph,
                                              bool ismatepair,
                                              GtError *err)
{
  FILE *infile;
  char *buffer, *c, *field;
  GtUword pos, bufferlen, result, num_pairs, fieldsize, rootctgid, ctgid;
  GtScaffoldGraphEdge *edge;
  GtWord dist;
  float std_dev;
  bool firstfield, nextfirstfield, sense, same;
  int had_err;

  had_err = 0;
  edge = NULL;
  infile = fopen(filename, "rb");
  if (infile == NULL) {
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

  if (result != bufferlen) {
    gt_error_set (err , " incomplete read of file %s ", filename);
    had_err = -1;
  }

  pos = 0;
  firstfield = true;
  nextfirstfield = true;
  /* sense direction as default */
  sense = true;

  /* iterate over file character by character */
  for (c = buffer; c != buffer + bufferlen; c++) {
    field[pos] = *c;
    pos++;

    if (fieldsize == pos) {
      fieldsize *= 2;
      field = gt_realloc(field, sizeof (*field) * fieldsize);
    }

    if (*c == ' ' || *c == '\n') {
      field[pos-1] = '\0';
      pos = 0;

      /* parse root contig ID */
      if (firstfield) {
        firstfield = false;

        rootctgid = gt_scaffolder_graph_get_vertex_id(graph, field, err);
        /* exit if distance and contig file inconsistent */
        gt_error_check(err);

        /* Debbuging:
          printf("rootctgid: %s\n",field);*/
      }
      else {
        /* parse contig-record with attributes distance, number of pairs,
           standard deviation of distance */
        if (sscanf(field,"%ld,%lu,%f",&dist,&num_pairs,&std_dev) == 3) {
          /* check if edge between vertices already exists */
          edge = gt_scaffolder_graph_find_edge(graph, rootctgid, ctgid);
          if (edge != NULL) {
            /*  LG: laut SGA edge->std_dev < std_dev,  korrekt? */
            if (ismatepair == false && edge->std_dev < std_dev) {
              /* LG: Ueberpruefung Kantenrichtung notwendig? */
              gt_scaffolder_graph_alter_edge(edge, dist, std_dev, num_pairs,
              sense, same);
            }
            else {
              /*Conflicting-Flag?*/
            }
          }
          else
            gt_scaffolder_graph_add_edge(graph, rootctgid, ctgid, dist, std_dev,
                                         num_pairs, sense, same);
          /* Debbuging:
             printf("dist: %ld\n num_pairs: %lu\n std_dev:"
                 "%f\n num5: %lu\n sense: %d\n\n",dist, num_pairs, std_dev,
                 num5,sense);*/

        }
        /* switch direction if semicolon occurs */
        else if (*field == ';')
          sense = !sense;
        nextfirstfield = true;
      }
    }
    /* parse contig ID */
    else if (*c == ',' && nextfirstfield) {
      nextfirstfield = false;
      field[pos-1] = '\0';
      pos = 0;
      /* parsing composition,
         '+' indicates same strand and '-' reverse strand */
      same = field[pos-2] == '+' ? true : false;

      ctgid = gt_scaffolder_graph_get_vertex_id(graph, field, err);
      /* exit if distance and contig file inconsistent */
      gt_error_check(err);

      /* Debbuging:
         printf("ctgid: %s\n",field);*/
    }
    if (*c == '\n') {
      firstfield = true;
      /* sense direction as default */
      sense = true;
    }
  }
  return had_err;
}

/* counts contigs with minimum length in callback data
   (fasta reader callback function, gets called after fasta entry
   has been read) */
static int gt_scaffolder_graph_count_ctg(GtUword length,
                                         void *data,
                                         GtError* err)
{
  int had_err;
  GtScaffoldGraphCallbackData *callback_data =
  (GtScaffoldGraphCallbackData*) data;

  had_err = 0;
  if (length >= callback_data->min_ctg_len)
    callback_data->nof_valid_ctg++;
  if (length == 0) {
    gt_error_set (err , " invalid sequence length ");
    had_err = -1;
  }
  return had_err;
}

/* saves header to callback data
   (fasta reader callback function, gets called for each description
    of fasta entry) */
static int gt_scaffolder_graph_save_header(const char *description,
                                           GtUword length,
                                           void *data, GtError *err)
{
  int had_err;
  GtScaffoldGraphCallbackData *callback_data =
  (GtScaffoldGraphCallbackData*) data;

  had_err = 0;
  callback_data->header_seq = description;
  callback_data->header_len = length;
  if (length == 0) {
    gt_error_set (err , " invalid header length ");
    had_err = -1;
  }
  return had_err;
}

/* saves header, sequence length of contig to scaffolder graph
   (fasta reader callback function, gets called after fasta entry
   has been read) */
static int gt_scaffolder_graph_save_ctg(GtUword seq_length,
                                        void *data,
                                        GtError* err)
{
  int had_err;
  GtScaffoldGraphCallbackData *callback_data =
  (GtScaffoldGraphCallbackData*) data;

  had_err = 0;
  if (seq_length > callback_data->min_ctg_len)
    gt_scaffolder_graph_add_vertex(callback_data->graph,
    callback_data->header_seq, seq_length, 0.0, 0.0);

  if (seq_length == 0) {
    gt_error_set (err , " invalid sequence length ");
    had_err = -1;
  }
  return had_err;
}

/* create scaffold graph from file */
GtScaffoldGraph *gt_scaffolder_graph_new_from_file(const char *ctg_filename,
                                                   GtUword min_ctg_len,
                                                   const char *dist_filename,
                                                   GtError *err)
{
  GtScaffoldGraph *graph;
  GtFastaReader* reader;
  GtStr *str_filename;
  int had_err;
  GtScaffoldGraphCallbackData *callback_data;

  had_err = 0;
  str_filename = gt_str_new();
  gt_str_set(str_filename, ctg_filename);
  callback_data = gt_malloc(sizeof(*callback_data));
  callback_data->nof_valid_ctg = 0;
  callback_data->min_ctg_len = min_ctg_len;
  /* parse contigs in FASTA-format and save them as vertices of
     scaffold graph */
  reader = gt_fasta_reader_rec_new(str_filename);
  had_err = gt_fasta_reader_run(reader, NULL, NULL,
            gt_scaffolder_graph_count_ctg, callback_data, err);

  graph = gt_malloc(sizeof(*graph));
  gt_scaffolder_graph_create_vertices(graph, callback_data->nof_valid_ctg);

  callback_data->graph = graph;
  had_err = gt_fasta_reader_run(reader, gt_scaffolder_graph_save_header, NULL,
            gt_scaffolder_graph_save_ctg, callback_data, err);
  gt_fasta_reader_delete(reader);
  gt_str_delete(str_filename);

  /* parse distance information of contigs in abyss-dist-format and
     save them as edges of scaffold graph */
  had_err =
    gt_scaffolder_graph_read_distances(dist_filename, graph, false, err);

  if (had_err != 0)
    gt_error_check(err);

  return graph;
}

/* check if unique order of edges <*edge1>, <*edge2> with probability
   <cutoff> exists */
static bool gt_scaffolder_graph_ambiguousorder(const GtScaffoldGraphEdge *edge1,
                                               const GtScaffoldGraphEdge *edge2,
                                               float cutoff)
{
  float expval, variance, interval, prob12, prob21;

  expval = edge1->dist - edge2->dist;
  /* sichere Multiplikation, Division? */
  variance = 2 * (edge1->std_dev * edge1->std_dev) -
             (edge2->std_dev * edge2->std_dev);
  interval = -expval / sqrt(variance);
  prob12 = 0.5 * (1 + erf(interval) );
  prob21 = 1.0 - prob12;

  return (prob12 <= cutoff && prob21 <= cutoff);
}

static GtUword calculate_overlap (GtScaffoldGraphEdge *edge1,
                                  GtScaffoldGraphEdge *edge2)
{
  GtUword overlap = 0;
  GtWord intersect_start, intersect_end;
  if (edge2->dist > (edge1->end->seq_len - 1) ||
      edge1->dist > (edge2->end->seq_len - 1)) {
    intersect_start = MAX(edge1->dist, edge2->dist);
    intersect_end = MAX(edge1->end->seq_len - 1, edge2->end->seq_len - 1);
    overlap = intersect_end - intersect_start + 1;
  }
  return overlap;
}

static void
gt_scaffolder_graph_check_mark_polymorphic(GtScaffoldGraphEdge *edge1,
                                           GtScaffoldGraphEdge *edge2,
                                           float pcutoff,
                                           float cncutoff)
{
  GtScaffoldGraphVertex *poly_vertex;
  GtScaffoldGraphEdge *edge;

  if (gt_scaffolder_graph_ambiguousorder(edge1, edge2, pcutoff) &&
      (edge1->end->copy_num + edge2->end->copy_num) < cncutoff) {
    /* mark vertex with lower copy number as polymorphic */
    if (edge1->end->copy_num < edge2->end->copy_num)
      poly_vertex = edge1->end;
    else
      poly_vertex = edge2->end;
    /* mark all edges of the polymorphic vertex as polymorphic */
    if (poly_vertex->state != GIS_POLYMORPHIC) {
      for (edge = poly_vertex->edges[0];
           edge < (poly_vertex->edges[0] + poly_vertex->nof_edges); edge++)
        edge->state = GIS_POLYMORPHIC;
      poly_vertex->state = GIS_POLYMORPHIC;
    }
  }
}

/* mark polymorphic edges/vertices and inconsistent edges in scaffold graph */
int gt_scaffolder_graph_filter(GtScaffoldGraph *graph,
                               float pcutoff,
                               float cncutoff,
                               GtUword ocutoff)
{
  GtScaffoldGraphVertex *vertex;
  GtScaffoldGraphEdge *edge1, *edge2;
  GtUword overlap;
  GtUword maxoverlap = 0;
  unsigned int dir; /* int statt bool, weil Iteration bislang nicht möglich */
  int had_err = 0;

  /* iterate over all vertices */
  for (vertex = graph->vertices;
       vertex < (graph->vertices + graph->nof_vertices); vertex++) {
    /* iterate over directions (sense/antisense) */
    for (dir = 0; dir < 2; dir++) {
      /* iterate over all pairs of edges */
      for (edge1 = vertex->edges[0];
           edge1 < (vertex->edges[0] + vertex->nof_edges); edge1++) {
        for (edge2 = edge1 + 1; edge2 < (vertex->edges[0] + vertex->nof_edges);
             edge2++) {
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
      for (edge1 = vertex->edges[0];
           edge1 < (vertex->edges[0] + vertex->nof_edges); edge1++) {
        for (edge2 = edge1 + 1; edge2 < (vertex->edges[0] + vertex->nof_edges);
             edge2++) {
          if (edge1->sense == dir && edge2->sense == dir &&
              edge1->state != GIS_POLYMORPHIC &&
              edge2->state != GIS_POLYMORPHIC) {
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
             edge1 < (vertex->edges[0] + vertex->nof_edges); edge1++)
          edge1->state = GIS_INCONSISTENT;
      }
    }
  }
  return had_err;
}

/* check if vertex holds just sense or antisense edges */
static bool gt_scaffolder_graph_isterminal(const GtScaffoldGraphVertex *vertex)
{
  GtScaffoldGraphEdge *edge;
  bool dir;

  dir = vertex->edges[0]->sense;
  for (edge = (vertex->edges[0] + 1); edge < (vertex->edges[0] +
       vertex->nof_edges); edge++) {
    if (edge->sense != dir)
      return false;
  }

  return true;
}

/* remove cycles
static void gt_scaffolder_removecycles(GtScaffoldGraph *graph) {
}*/

/* create new walk */
static GtScaffoldGraphWalk *gt_scaffolder_walk_new(void)
{
  GtScaffoldGraphWalk *walk;

  walk = gt_malloc(sizeof(*walk));
  walk->total_contig_len = 0;
  walk->size = 0;
  walk->nof_edges = 0;
  return walk;
}

/* remove walk <*walk> */
static void gt_scaffolder_walk_delete(GtScaffoldGraphWalk *walk)
{
  if (walk != NULL)
    gt_free(walk->edges);
  gt_free(walk);
}

/* add edge <*edge> to walk <*walk> */
static void gt_scaffolder_walk_addegde(GtScaffoldGraphWalk *walk,
                                       GtScaffoldGraphEdge *edge)
{
  if (walk->size == walk->nof_edges) {
    walk->size += INCREMENT_SIZE;
    walk->edges = gt_realloc(walk->edges, walk->size*sizeof(*walk->edges));
  }
  walk->edges[walk->nof_edges] = edge;
  walk->total_contig_len += edge->end->seq_len;
  walk->nof_edges++;
}

GtScaffoldGraphWalk *gt_scaffolder_create_walk(GtScaffoldGraph *graph,
                 GtScaffoldGraphVertex *start)
{
  /* BFS-Traversierung innerhalb aktueller Zusammenhangskomponente
     ausgehend von terminalen Knoten zu terminalen Knoten */
  GtQueue *wqueue;
  GtArray *terminal_vertices;
  GtScaffoldGraphEdge *edge, *reverseedge, *nextedge, **edgemap;
  GtScaffoldGraphVertex *endvertex, *currentvertex, *nextendvertex;
  GtUword lengthbestwalk, lengthcwalk;
  GtScaffoldGraphWalk *bestwalk, *currentwalk;
  float distance, *distancemap;
  bool dir;
  lengthbestwalk = 0;
  bestwalk = gt_scaffolder_walk_new();

  wqueue = gt_queue_new();
  terminal_vertices = gt_array_new(sizeof(start));

  /* SK: Mit GtWord_Min / erwarteter Genomlaenge statt 0 initialisieren */
  distancemap = calloc(graph->nof_vertices, sizeof (*distancemap));
  edgemap = gt_malloc(sizeof (*edgemap)*graph->nof_vertices);

  dir = start->edges[0]->sense;
  for (edge = start->edges[0];
       edge < (start->edges[0] + start->nof_edges); edge++) {
    endvertex = edge->end;

    /* SK: genometools hashes verwenden, Dichte evaluieren
       SK: DistEst beim Einlesen prüfen */
    distancemap[endvertex->id] = edge->dist;
    edgemap[endvertex->id] = edge;

    gt_queue_add(wqueue, edge);
  }

  while (gt_queue_size(wqueue) != 0) {
    edge = (GtScaffoldGraphEdge*)gt_queue_get(wqueue);
    endvertex = edge->end;

    /* store all terminal vertices */
    if (gt_scaffolder_graph_isterminal(endvertex))
      gt_array_add(terminal_vertices, endvertex);

    for (nextedge = endvertex->edges[0];
    nextedge < (endvertex->edges[0] + endvertex->nof_edges); nextedge++) {
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

  /* Ruecktraversierung durch EdgeMap für alle terminalen Knoten
     Konstruktion des Walks  */
  while (gt_array_size(terminal_vertices) != 0) {
    currentvertex = gt_array_pop(terminal_vertices);

    currentwalk = gt_scaffolder_walk_new();
    while (currentvertex->id != start->id) {
      reverseedge = edgemap[currentvertex->id];
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
void gt_scaffolder_makescaffold(GtScaffoldGraph *graph)
{
  GtScaffoldGraphVertex *vertex, *currentvertex, *nextvertex, *start;
  GtScaffoldGraphEdge *edge;
  GtScaffoldGraphWalk *walk;
  GtUword ccnumber;
  GtQueue *vqueue;
  GtArray *terminal_vertices, *cc_walks;

  /* Entfernung von Zyklen
     gt_scaffolder_removecycles(graph); */

  /* Iteration ueber alle Knoten, Makierung aller Knoten als unbesucht */
  for (vertex = graph->vertices;
       vertex < (graph->vertices + graph->nof_vertices); vertex++) {
    if (vertex->state != GIS_POLYMORPHIC)
      vertex->state = GIS_UNVISITED;
  }

 /* BFS-Traversierung durch Zusammenhangskomponenten des Graphen,
    siehe GraphSearchTree.h */
  ccnumber = 0;
  vqueue = gt_queue_new();
  terminal_vertices = gt_array_new(sizeof(vertex));
  cc_walks = gt_array_new(sizeof(walk));

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
      currentvertex = (GtScaffoldGraphVertex*)gt_queue_get(vqueue);
      /*currentvertex.cc = ccnumber;*/

      /* store all terminal vertices to calculate all paths between them */
      if (gt_scaffolder_graph_isterminal(currentvertex))
  gt_array_add(terminal_vertices, currentvertex);

      currentvertex->state = GIS_VISITED;
      for (edge = currentvertex->edges[0];
           edge < (currentvertex->edges[0] + currentvertex->nof_edges);
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

    /* calculate all paths between terminal vertices in this cc */
    while (gt_array_size(terminal_vertices) != 0) {
      start = gt_array_pop(terminal_vertices);
      walk = gt_scaffolder_create_walk(graph, start);
      gt_array_add(cc_walks, walk);
    }

    /* TODO: the best walk in this cc has to be chosen */
  }

  gt_array_delete(terminal_vertices);
  gt_array_delete(cc_walks);
  gt_queue_delete(vqueue);
}
