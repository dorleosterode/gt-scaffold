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
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/minmax.h"
#include "core/fasta_reader_rec.h"
#include "core/error.h"

/* increment size for realloc of walk */
#define INCREMENT_SIZE 32

typedef enum { GIS_UNVISITED, GIS_POLYMORPHIC, GIS_INCONSISTENT,
               GIS_REPEAT, GIS_VISITED, GIS_PROCESSED } GraphItemState;

/* vertex of scaffold graph (describes one contig) */
typedef struct GtScaffolderGraphVertex {
  /* unique vertex ID */
  GtUword id;
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
struct GtScaffolderGraph {
  GtScaffolderGraphVertex *vertices;
  GtUword nof_vertices;
  GtUword max_nof_vertices;
  GtScaffolderGraphEdge *edges;
  GtUword nof_edges;
  GtUword max_nof_edges;
};

/* linear scaffold */
typedef struct GtScaffolderGraphWalk {
  GtUword nof_edges;
  GtUword size;
  GtUword total_contig_len;
  GtScaffolderGraphEdge **edges;
}GtScaffolderGraphWalk;

/* for parsing valid contigs,
   e.g. contigs with minimum length <min_ctg_len> */
typedef struct {
  GtUword nof_valid_ctg;
  GtUword min_ctg_len;
  GtStr *header_seq;
  GtScaffolderGraph *graph;
} GtScaffolderGraphFastaReaderData;

/* DistanceMap */
/* EdgeMap */

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

  if (graph->vertices != NULL) {
    /* If existent, free header_seq and pointer to outgoing edges first */
    for ( vertex = graph->vertices;
          vertex < (graph->vertices + graph->nof_vertices);
          vertex++
        )
    {
      if (vertex->header_seq != NULL)
        gt_str_delete(vertex->header_seq);
      if (vertex->edges != NULL)
        gt_free(vertex->edges);
    }
    gt_free(graph->vertices);
  }

  if (graph->edges != NULL)
    gt_free(graph->edges);

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
  graph->vertices[nextfree].id = nextfree; /* SD: remove without breaking */
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

/* Initialize a new, directed edge in <*graph>. Each edge between two contig
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

static GtScaffolderGraphEdge *gt_scaffolder_graph_find_edge(
                                    const GtScaffolderGraph *graph,
                                    GtUword vertexid_1,
                                    GtUword vertexid_2)
{
  GtScaffolderGraphEdge **edge;
  GtScaffolderGraphVertex *v1 = graph->vertices + vertexid_1;

  for (edge = v1->edges; edge < (v1->edges + v1->nof_edges); edge++) {
    if ((*edge)->end->id == vertexid_2)
      return *edge;
  }
  return NULL;
}

/* determines corresponding vertex id to contig header */
/* SD: Binaersuche separat testen */
static int gt_scaffolder_graph_get_vertex_id(const GtScaffolderGraph *graph,
                                             GtUword *vertex_id,
                                             const GtStr *header_seq,
                                             GtError *err)
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
        *vertex_id = mid_vertex->id;
        break;
      }
      else if (cmp < 0)
        min_vertex = mid_vertex + 1;
      else
        max_vertex = mid_vertex - 1;
    }

  /* contig header was not found */
  if (found == false)
  {
    had_err = -1;
    gt_error_set(err, " distance and contig file inconsistent ");
  }
  return had_err;
}

/* assign edge <*edge> new attributes */
static void gt_scaffolder_graph_alter_edge(GtScaffolderGraphEdge *edge,
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
    gt_file_xprintf(f, GT_WU " [color=\"%s\" label=\"%s\"];\n", v->id,
                    color_array[v->state], gt_str_get(v->header_seq));
  }

  /* iterate over all edges and print them. add attribute color according to
     the current state and label the edge with the distance*/
  for (e = g->edges; e < (g->edges + g->nof_edges); e++) {
    gt_file_xprintf(f,
                    GT_WU " -- " GT_WU " [color=\"%s\" label=\"" GT_WD "\"];\n",
                    e->start->id, e->end->id,
                    color_array[e->state], e->dist);
  }

  /* print the last line into f */
  gt_file_xprintf(f, "}\n");
}


/* count records */
static int gt_scaffolder_graph_count_distances(const char *file_name,
                                              GtUword *nof_distances,
                                              GtError *err)
{
  FILE *file;
  char line[1024], *field, ctg_header[1024];
  GtUword num_pairs, record_counter;
  GtWord dist;
  float std_dev;
  int had_err;

  had_err = 0;
  record_counter = 0;
  file = fopen(file_name, "rb");
  if (file == NULL){
    had_err = -1;
    gt_error_set(err, " can not read file %s ",file_name);
  }

  if (had_err != -1)
  {
    /* iterate over each line of file until eof (contig record) */
    while (fgets(line, 1024, file) != NULL)
    {
      /* split line by first space delimiter */
      field = strtok(line," ");

      /* iterate over space delimited records */
      while (field != NULL)
      {
        /* count records */
        if (sscanf(field,"%[^<,],%ld,%lu,%f", ctg_header, &dist, &num_pairs,
            &std_dev) == 4)
          record_counter++;

        /* split line by next space delimiter */
        field = strtok(NULL," ");
      }
    }
  }
  *nof_distances = record_counter;

  fclose(file);
  return had_err;
}

/* sort by lexicographic ascending order */
static int gt_scaffolder_graph_vertices_compare(const void *a, const void *b)
{
  GtScaffolderGraphVertex *vertex_a =
  (GtScaffolderGraphVertex*) a;
  GtScaffolderGraphVertex *vertex_b =
  (GtScaffolderGraphVertex*) b;
  return gt_str_cmp(vertex_a->header_seq, vertex_b->header_seq);
}


/* parse distance information of contigs in abyss-dist-format and
   save them as edges of scaffold graph, PRECONDITION: header contains
   no commas and spaces */
/* LG: check for "mate-flag"? */
/* SK: DistParser in eigenes Modul auslagern? */
/* Bsp.: Ctg1 Ctg2+,15,10,5.1 ; Ctg3-,65,10,5.1 */
static int gt_scaffolder_graph_read_distances(const char *filename,
                                              GtScaffolderGraph *graph,
                                              bool ismatepair,
                                              GtError *err)
{
  FILE *file;
  /* SD: Konstante setzen? */
  char line[1024], *field, ctg_header[1024];
  GtUword num_pairs, root_ctg_id, ctg_id, ctg_header_len;
  GtWord dist;
  float std_dev;
  bool same, sense;
  GtScaffolderGraphEdge *edge;
  int had_err;
  GtStr *gt_str_ctg_header, *gt_str_field;
  GtScaffolderGraphVertex *v;

  had_err = 0;
  root_ctg_id = 0;
  ctg_id = 0;

  file = fopen(filename, "rb");
  if (file == NULL){
    had_err = -1;
    gt_error_set(err, " can not read file %s ",filename);
  }

  /* sort by lexicographic ascending order */
  qsort(graph->vertices, graph->nof_vertices, sizeof(*graph->vertices),
       gt_scaffolder_graph_vertices_compare);
  /* update vertex ID
     LG: Knoten ID eigentlich redundant? SD: Ja. */
  for (v = graph->vertices; v < (graph->vertices + graph->nof_vertices); v++)
    v->id = v - graph->vertices;

  if (had_err != -1)
  {
    /* iterate over each line of file until eof (contig record) */
    while (fgets(line, 1024, file) != NULL)
    {
      /* remove '\n' from end of line */
      line[strlen(line)-1] = '\0';
      /* set sense direction as default */
      sense = true;
      /* split line by first space delimiter */
      field = strtok(line," ");

      /* get vertex id corresponding to root contig header */
      gt_str_field = gt_str_new_cstr(field);
      had_err = gt_scaffolder_graph_get_vertex_id(graph, &root_ctg_id,
                gt_str_field, err);
      gt_str_delete(gt_str_field);

      /* exit if distance and contig file inconsistent */
      gt_error_check(err);
      /* Debbuging: printf("rootctgid: %s\n",field);*/

      /* iterate over space delimited records */
      while (field != NULL)
      {
        /* parse record consisting of contig header, distance,
           number of pairs, std. dev. */
        /* SD: %[^<,] ist eine negierte Zeichenklasse (Workaround weil %s nicht
               funktioniert
        */
        if (sscanf(field,"%[^<,],%ld,%lu,%f", ctg_header, &dist, &num_pairs,
            &std_dev) == 4)
        {
          /* parsing composition,
           '+' indicates same strand and '-' reverse strand */
          ctg_header_len = strlen(ctg_header);
          same = ctg_header[ctg_header_len - 1] == '+' ? true : false;

          /* cut composition sign */
          ctg_header[ctg_header_len - 1] = '\0';
          gt_str_ctg_header = gt_str_new_cstr(ctg_header);
          /* get vertex id corresponding to contig header */
          had_err = gt_scaffolder_graph_get_vertex_id(graph, &ctg_id,
                    gt_str_ctg_header, err);
          gt_str_delete(gt_str_ctg_header);

          /* exit if distance and contig file inconsistent */
          gt_error_check(err);

          /* check if edge between vertices already exists */
          edge = gt_scaffolder_graph_find_edge(graph, root_ctg_id, ctg_id);
          if (edge != NULL)
          {
            /*  LG: laut SGA edge->std_dev < std_dev,  korrekt? */
            if (ismatepair == false && edge->std_dev < std_dev)
            {
              /* LG: Ueberpruefung Kantenrichtung notwendig? */
              gt_scaffolder_graph_alter_edge(edge, dist, std_dev, num_pairs,
              sense, same);
            }
            /*else { Conflicting-Flag? }*/
          }
          else
            gt_scaffolder_graph_add_edge(graph, root_ctg_id, ctg_id, dist, std_dev,
                                       num_pairs, sense, same);
          /* Debbuging:
             printf("ctgid: %s\n",field);
             printf("dist: " GT_WD "\n num_pairs: " GT_WU "\n std_dev:"
                "%f\n num5: " GT_WU "\n sense: %d\n\n",dist, num_pairs, std_dev,
                num5,sense);
          */
        }
        /* switch direction */
        else if (*field == ';')
          sense = !sense;

        /* split line by next space delimiter */
        field = strtok(NULL," ");
      }
    }
  }

  fclose(file);
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
  GtScaffolderGraphFastaReaderData *fasta_reader_data =
  (GtScaffolderGraphFastaReaderData*) data;

  had_err = 0;
  if (length >= fasta_reader_data->min_ctg_len)
    fasta_reader_data->nof_valid_ctg++;
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
  GtScaffolderGraphFastaReaderData *fasta_reader_data =
  (GtScaffolderGraphFastaReaderData*) data;
  GtStr *gt_str_description;
  char *new_description, *space_ptr;

  had_err = 0;

  gt_str_description = gt_str_new_cstr(description);
  new_description = gt_str_get(gt_str_description);
  /* cut header sequence after first space */
  space_ptr = strchr(new_description, ' ');
  if (space_ptr != NULL)
    *space_ptr = '\0';
  gt_str_set(gt_str_description, new_description);

  fasta_reader_data->header_seq = gt_str_description;
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
  GtScaffolderGraphFastaReaderData *fasta_reader_data =
  (GtScaffolderGraphFastaReaderData*) data;

  had_err = 0;
  if (seq_length > fasta_reader_data->min_ctg_len)
  {
    gt_scaffolder_graph_add_vertex(fasta_reader_data->graph,
    fasta_reader_data->header_seq, seq_length, 0.0, 0.0);

    gt_str_delete(fasta_reader_data->header_seq);
  }


  if (seq_length == 0) {
    gt_error_set (err , " invalid sequence length ");
    had_err = -1;
  }
  return had_err;
}


/* load astatics and copy number of every contig and mark repeated contigs */
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
  int had_err = 0;
  GtStr *gt_str_field;
  GtScaffolderGraphVertex *vertex;

  file = fopen(filename, "rb");
  if (file == NULL){
    had_err = -1;
    gt_error_set(err, " can not read file %s ", filename);
  }

  if (had_err != -1)
  {
    /* iterate over each line of file until eof (contig record) */
    while (fgets(line, 1024, file) != NULL)
    {
      /* remove '\n' from end of line */
      line[strlen(line)-1] = '\0';

      /* split line by first tab delimiter */
      field = strtok(line,"\t");

      /* get vertex id corresponding to root contig header */
      gt_str_field = gt_str_new_cstr(field);
      had_err = gt_scaffolder_graph_get_vertex_id(graph, &root_ctg_id,
                gt_str_field, err);
      gt_str_delete(gt_str_field);
      /* exit if distance and contig file inconsistent */
      gt_error_check(err);

      field_counter = 0;
      copy_num = 0.0;
      astat = 0.0;
      copy_num_found = false;
      astat_found = false;
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

/* create scaffold graph from file */
/* TODO: include a-statistics, copy number */
GtScaffolderGraph *gt_scaffolder_graph_new_from_file(const char *ctg_filename,
                                                     GtUword min_ctg_len,
                                                     const char *dist_filename,
                                                     GtError *err)
{
  GtScaffolderGraph *graph;
  GtFastaReader* reader;
  GtStr *str_filename;
  int had_err;
  GtScaffolderGraphFastaReaderData fasta_reader_data;
  GtUword nof_distances;

  had_err = 0;
  graph = NULL;
  str_filename = gt_str_new_cstr(ctg_filename);
  fasta_reader_data.nof_valid_ctg = 0;
  fasta_reader_data.min_ctg_len = min_ctg_len;
  /* count contigs */
  reader = gt_fasta_reader_rec_new(str_filename);
  had_err = gt_fasta_reader_run(reader, NULL, NULL,
            gt_scaffolder_graph_count_ctg, &fasta_reader_data, err);
  gt_fasta_reader_delete(reader);

  if (had_err == 0)
  {
    /* allocate memory for vertices of scaffolder graph */
    graph = gt_malloc(sizeof(*graph));
    graph->vertices = NULL;
    graph->edges = NULL;
    gt_scaffolder_graph_init_vertices(graph, fasta_reader_data.nof_valid_ctg);

    fasta_reader_data.graph = graph;
    /* parse contigs in FASTA-format and save them as vertices of
       scaffold graph */
    reader = gt_fasta_reader_rec_new(str_filename);
    had_err = gt_fasta_reader_run(reader, gt_scaffolder_graph_save_header,
              NULL, gt_scaffolder_graph_save_ctg, &fasta_reader_data, err);
    gt_fasta_reader_delete(reader);

    if (had_err == 0)
    {
      /* count distance information */
      nof_distances = 0;
      had_err = gt_scaffolder_graph_count_distances(dist_filename,
              &nof_distances, err);

      if (had_err == 0)
      {
        /* allocate memory for edges of scaffolder graph */
        gt_scaffolder_graph_init_edges(graph, nof_distances);
        /* parse distance information of contigs in abyss-dist-format and
           save them as edges of scaffold graph */
        had_err =
          gt_scaffolder_graph_read_distances(dist_filename, graph, false, err);
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

  gt_str_delete(str_filename);
  return graph;
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

static GtUword gt_scaffolder_calculate_overlap(GtScaffolderGraphEdge *edge1,
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
static bool
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

/*  remove cycles */
/* static void gt_scaffolder_removecycles(GtScaffolderGraph *graph) { */

/*     /\* found one cc first remove all cycles starting at terminal */
/*        vertices. keep the array with terminal vertices clean. *\/ */
/*     while (gt_array_size(terminal_vertices) != 0) { */
/*       /\* search for cycles starting at every terminal vertex. if a */
/* 	 backedge is found, mark both vertices connected with this */
/* 	 backedge. *\/ */
/*     } */

/* } */

/* DFS to detect Cycles given a starting vertex */
/* TODO: remember all visited vertices to change the state to
   GIS_UNVISITED before search with the next starting vertex. */

/* SD: Commented this function because compiler complaines about it being unused */

/*static GtScaffolderGraphEdge
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
	 point and return a boolean or something like that *//*
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
static GtScaffolderGraphWalk *gt_scaffolder_walk_new(void)
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
static void gt_scaffolder_walk_delete(GtScaffolderGraphWalk *walk)
{
  if (walk != NULL)
    gt_free(walk->edges);
  gt_free(walk);
}

/* add edge <*edge> to walk <*walk> */
static void gt_scaffolder_walk_addegde(GtScaffolderGraphWalk *walk,
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
    distancemap[endvertex->id] = edge->dist;
    edgemap[endvertex->id] = edge;

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
    currentvertex = *(GtScaffolderGraphVertex **) gt_array_pop(terminal_vertices);
    gt_assert(currentvertex < graph->vertices + graph->nof_vertices);
    gt_assert(currentvertex >= graph->vertices);

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
void gt_scaffolder_makescaffold(GtScaffolderGraph *graph)
{
  gt_assert(graph != NULL);

  GtScaffolderGraphVertex *vertex, *currentvertex, *nextvertex, *start;
  GtUword eid;
  GtScaffolderGraphWalk *walk;
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
	gt_assert((*(GtScaffolderGraphVertex **)gt_array_get_last(terminal_vertices)) == currentvertex);
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

    /* TODO: the best walk in this cc has to be chosen */
  }

  gt_array_delete(terminal_vertices);
  gt_array_delete(cc_walks);
  gt_queue_delete(vqueue);
}

/* Function to test basic graph functionality on different scenarios:
- Create graph and allocate space for <max_nof_vertices> vertices and
  <max_nof_edges> edges.
- Initialize vertices if <init_vertices> is true and create <nof_vertices>
  vertices.
- Initialize edges if <init_edges> is true and create <nof_edges> edges.
- Delete graph. */
int gt_scaffolder_test_graph(GtUword max_nof_vertices,
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
    /* SD: Muss dennoch gefreed werden */
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
    char outfile[] = "gt_scaffolder_test.dot";
    err = gt_error_new();
    gt_scaffolder_graph_print(graph, outfile, err);
  }

  gt_scaffolder_graph_delete(graph);
  if (graph != NULL)
    had_err = -1;

  return had_err;
}
