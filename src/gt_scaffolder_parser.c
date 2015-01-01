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

#include "core/fasta_reader_rec.h"
#include "core/ma_api.h"

#include "gt_scaffolder_graph.h"
#include "gt_scaffolder_parser.h"

const GtUword BUFSIZE = 1024;

/* for parsing valid contigs,
   e.g. contigs with minimum length <min_ctg_len> */
typedef struct {
  GtUword nof_valid_ctg;
  GtUword min_ctg_len;
  GtStr *header_seq;
  GtScaffolderGraph *graph;
} GtScaffolderGraphFastaReaderData;

/* sort by lexicographic ascending order */
static int gt_scaffolder_graph_vertices_compare(const void *a, const void *b)
{
  GtScaffolderGraphVertex *vertex_a =
  (GtScaffolderGraphVertex*) a;
  GtScaffolderGraphVertex *vertex_b =
  (GtScaffolderGraphVertex*) b;
  return gt_str_cmp(vertex_a->header_seq, vertex_b->header_seq);
}

/* test parsing distance records */
int gt_scaffolder_parser_read_distances_test(const char *filename,
                                             char *output_filename,
                                             GtError *err)
{
  FILE *file;
  GtFile *f;
  char line[BUFSIZE+1], *field, ctg_header[BUFSIZE+1], sign;
  GtUword ctg_header_len;
  GtWord dist, num_pairs;
  float std_dev;
  bool same, sense, first_antisense;
  int had_err;

  had_err = 0;
  file = fopen(filename, "rb");
  if (file == NULL) {
    had_err = -1;
    gt_error_set(err, "can not read file %s",filename);
  }

  if (had_err != -1) {
    f = gt_file_new(output_filename, "w", err);
    if (f == NULL)
      had_err = -1;
  }

  if (had_err != -1)
  {
    /* iterate over each line of file until eof (contig record) */
    while (fgets(line, BUFSIZE, file) != NULL)
    {
      /* remove '\n' from end of line */
      line[strlen(line)-1] = '\0';
      /* set sense direction as default */
      sense = true;
      /* split line by first space delimiter */
      field = strtok(line," ");

      /* write parsed distance information to file */
      gt_file_xprintf(f, "%s", field);
      first_antisense = true;

      /* iterate over space delimited records */
      while (field != NULL)
      {
        /* parse record consisting of contig header, distance,
           number of pairs, std. dev. */
        if (sscanf(field,"%[^>,]," GT_WD "," GT_WD ",%f", ctg_header, &dist,
            &num_pairs, &std_dev) == 4)
        {
          /* ignore invalid records */
          if (num_pairs < 0) {
            had_err = -1;
            gt_error_set(err, "Invalid value for number of pairs");
            break;
          }

          /* parsing composition,
           '+' indicates same strand and '-' reverse strand */
          ctg_header_len = strlen(ctg_header);
          same = ctg_header[ctg_header_len - 1] == '+' ? true : false;

          /* cut composition sign */
          ctg_header[ctg_header_len - 1] = '\0';

          /* write parsed distance information to file */
          sign = same == true ? '+' : '-';
          gt_file_xprintf(f, " %s%c," GT_WD "," GT_WD ",%.1f", ctg_header,
                  sign, dist, num_pairs, std_dev);
        }
        /* switch direction */
        else if (*field == ';')
          sense = sense ? false : true;

        /* split line by next space delimiter */
        field = strtok(NULL," ");

        if (!sense && first_antisense) {
          gt_file_xprintf(f, " ;");
          first_antisense = false;
        }

      }
      if (had_err == -1)
        break;

      if (sense)
        gt_file_xprintf(f, " ;");
      gt_file_xprintf(f, "\n");
    }
  }

  if (file != NULL)
    gt_file_delete(f);

  fclose(file);
  return had_err;
}

/* count records and check integrity of abyss-dist-format */
int gt_scaffolder_parser_count_distances(const GtScaffolderGraph *graph,
                                               const char *file_name,
                                               GtUword *nof_distances,
                                               GtError *err)
{
  FILE *file;
  char line[BUFSIZE+1], *field, ctg_header[BUFSIZE+1], comp_sign;
  GtUword record_counter, *edge_counter,
          line_record_counter;
  GtWord dist, num_pairs;
  float std_dev;
  int had_err;
  bool valid_contig;
  GtStr *gt_str_field;
  GtScaffolderGraphVertex *v, *ctg, *root_ctg;

  had_err = 0;
  record_counter = 0;
  gt_str_field = gt_str_new();

  /* sort by lexicographic ascending order */
  qsort(graph->vertices, graph->nof_vertices, sizeof (*graph->vertices),
       gt_scaffolder_graph_vertices_compare);
  /* update vertex ID
     LG: Knoten ID eigentlich redundant? SD: Ja. LG: OK?!
  for (v = graph->vertices; v < (graph->vertices + graph->nof_vertices); v++)
    v->index = v - graph->vertices;*/

  file = fopen(file_name, "rb");
  if (file == NULL) {
    had_err = -1;
    gt_error_set(err, "can not read file %s",file_name);
  }

  edge_counter = gt_calloc(graph->nof_vertices, sizeof (*edge_counter));

  if (had_err != -1)
  {
    /* iterate over each line of file until eof (contig record) */
    while (fgets(line, BUFSIZE, file) != NULL)
    {
      /* split line by first space delimiter */
      field = strtok(line," ");

      gt_str_set(gt_str_field, field);
      valid_contig = gt_scaffolder_graph_get_vertex(graph, &root_ctg,
                gt_str_field);

      /* split line by next space delimiter */
      field = strtok(NULL," ");

      /* if no records exist */
      if (field == NULL) {
        had_err = -1;
        gt_error_set(err, "Invalid record in dist file %s",
                           file_name);
        break;
      }

      if (valid_contig) {

        line_record_counter = 0;

        /* iterate over space delimited records */
        while (field != NULL)
        {
          /* count records */
          /* SD: Keep an eye on negated string, might fail, did before */
          if (sscanf(field,"%[^>,]," GT_WD "," GT_WD ",%f", ctg_header,
              &dist, &num_pairs, &std_dev) == 4) {

            /* detect invalid records */
            if (num_pairs < 0) {
              had_err = -1;
              gt_error_set(err, "Invalid value for number of pairs in dist "
                                "file %s", file_name);
              break;
            }

            comp_sign = ctg_header[strlen(ctg_header) - 1];
            if (comp_sign != '+' && comp_sign != '-') {
              had_err = -1;
              gt_error_set(err, "Invalid composition sign in dist file %s",
                                 file_name);
              break;
            }

            /* cut composition sign */
            ctg_header[strlen(ctg_header) - 1] = '\0';

            gt_str_set(gt_str_field, ctg_header);
            /* get vertex id corresponding to contig header */
            valid_contig = gt_scaffolder_graph_get_vertex(graph, &ctg,
                           gt_str_field);

            if (valid_contig) {

              edge_counter[ctg-graph->vertices] += 1;
              line_record_counter++;
            }

          }
          /* detect invalid record */
          else if (*field != ';') {
            had_err = -1;
            gt_error_set(err, "Invalid record in dist file %s",
                               file_name);
            break;
          }

          /* split line by next space delimiter */
          field = strtok(NULL," ");
        }
        if (had_err == -1)
          break;
        edge_counter[root_ctg-graph->vertices] += line_record_counter;
        record_counter += line_record_counter;
      }
    }
    fclose(file);
  }

  if (had_err != -1) {
    /* allocate memory for edges of vertices */
    for (v = graph->vertices; v < (graph->vertices + graph->nof_vertices);
         v++) {
      if (edge_counter[gt_scaffolder_graph_get_vertex_id(graph, v)] != 0)
        v->edges = gt_malloc(sizeof(*v->edges) *
         edge_counter[gt_scaffolder_graph_get_vertex_id(graph, v)]);
    }
    *nof_distances = record_counter;
  }

  gt_free(edge_counter);
  gt_str_delete(gt_str_field);

  return had_err;
}

/* parse distance information of contigs in abyss-dist-format and
   save them as edges of scaffold graph
   PRECONDITION: header contains no commas and spaces */
/* LG: check for "mate-flag"? */
/* Bsp.: Ctg1 Ctg2+,15,10,5.1 ; Ctg3-,65,10,5.1 */
int gt_scaffolder_parser_read_distances(const char *filename,
                                              GtScaffolderGraph *graph,
                                              bool ismatepair,
                                              GtError *err)
{
  FILE *file;
  /* SD: Konstante setzen? */
  char line[BUFSIZE+1], *field, ctg_header[BUFSIZE+1];
  GtUword ctg_header_len;
  GtWord dist, num_pairs;
  float std_dev;
  bool same, sense, valid_contig;
  GtScaffolderGraphEdge *edge;
  GtScaffolderGraphVertex *root_ctg, *ctg;
  int had_err;
  GtStr *gt_str_field;

  had_err = 0;
  gt_str_field = gt_str_new();

  file = fopen(filename, "rb");
  if (file == NULL) {
    had_err = -1;
    gt_error_set(err, " can not read file %s ",filename);
  }

  if (had_err != -1)
  {
    /* iterate over each line of file until eof (contig record) */
    while (fgets(line, BUFSIZE, file) != NULL)
    {
      /* remove '\n' from end of line */
      line[strlen(line)-1] = '\0';
      /* set sense direction as default */
      sense = true;
      /* split line by first space delimiter */
      field = strtok(line," ");

      /* get vertex id corresponding to root contig header */
      /* SK: Moeglichkeit einer set Funktion evaluieren */
      gt_str_set(gt_str_field, field);
      valid_contig = gt_scaffolder_graph_get_vertex(graph, &root_ctg,
                gt_str_field);

      /* Debbuging: printf("rootctgid: %s\n",field);*/

      if (valid_contig) {
        /* iterate over space delimited records */
        while (field != NULL)
        {
          /* parse record consisting of contig header, distance,
             number of pairs, std. dev. */
          /* SD: %[^>,] ist eine negierte Zeichenklasse (Workaround weil %s
                 nicht funktioniert)
          */
          if (sscanf(field,"%[^>,]," GT_WD "," GT_WD ",%f", ctg_header, &dist,
              &num_pairs, &std_dev) == 4)
          {

            /* parsing composition,
             '+' indicates same strand and '-' reverse strand */
            ctg_header_len = strlen(ctg_header);
            same = ctg_header[ctg_header_len - 1] == '+' ? true : false;

            /* cut composition sign */
            ctg_header[ctg_header_len - 1] = '\0';

            gt_str_set(gt_str_field, ctg_header);
            /* get vertex id corresponding to contig header */
            valid_contig = gt_scaffolder_graph_get_vertex(graph, &ctg,
                      gt_str_field);

            if (valid_contig) {
              /* check if edge between vertices already exists */
              edge = gt_scaffolder_graph_find_edge(root_ctg, ctg);
              if (edge != NULL)
              {
                /*  LG: laut SGA edge->std_dev < std_dev,  korrekt? */
                if (!ismatepair && edge->std_dev < std_dev)
                {
                  /* LG: Ueberpruefung Kantenrichtung notwendig? */
                  gt_scaffolder_graph_alter_edge(edge, dist, std_dev,
                                               num_pairs,sense, same);
                }
                /*else { Conflicting-Flag? }*/
              }
              else
                gt_scaffolder_graph_add_edge(graph, root_ctg, ctg, dist,
                                              std_dev,num_pairs, sense, same);
            }
          }
          /* switch direction */
          else if (*field == ';')
            sense = sense ? false : true;

          /* split line by next space delimiter */
          field = strtok(NULL," ");
        }
      }
    }
    fclose(file);
  }

  gt_str_delete(gt_str_field);
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

/* str_dup with gt_malloc */
char *gt_strdup (const char *source)
{
  char *dest;
  dest = gt_malloc (strlen (source) + 1);
  if (dest == NULL) return NULL;
  strcpy (dest, source);
  return dest;
}

/* saves header to callback data
   (fasta reader callback function, gets called for each description
    of fasta entry) */
static int gt_scaffolder_graph_save_header(const char *description,
                                           GtUword length,
                                           void *data, GtError *err)
{
  int had_err;
  char *writeable_description, *space_ptr;
  GtScaffolderGraphFastaReaderData *fasta_reader_data =
  (GtScaffolderGraphFastaReaderData*) data;

  had_err = 0;
  writeable_description = gt_strdup(description);
  /* cut header sequence after first space */
  space_ptr = strchr(writeable_description, ' ');
  if (space_ptr != NULL)
    *space_ptr = '\0';
  /* overwrite current GtString */
  gt_str_set(fasta_reader_data->header_seq, writeable_description);
  gt_free(writeable_description);

  if (length == 0) {
    gt_error_set (err , "Invalid header length");
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
  GtStr *cloned_gt_str;
  GtScaffolderGraphFastaReaderData *fasta_reader_data =
  (GtScaffolderGraphFastaReaderData*) data;

  had_err = 0;
  if (seq_length > fasta_reader_data->min_ctg_len)
  {
    cloned_gt_str = gt_str_clone(fasta_reader_data->header_seq);
    gt_scaffolder_graph_add_vertex(fasta_reader_data->graph,
    cloned_gt_str, seq_length, 0.0, 0.0);
  }

  if (seq_length == 0) {
    gt_error_set (err , "Invalid sequence length");
    had_err = -1;
  }
  return had_err;
}

/* count contigs */
int gt_scaffolder_parser_count_contigs(const char *filename,
                                       GtUword min_ctg_len,
                                       GtUword *nof_contigs,
                                       GtError *err)
{
  GtFastaReader* reader;
  GtStr *str_filename;
  GtScaffolderGraphFastaReaderData fasta_reader_data;
  int had_err;

  str_filename = gt_str_new_cstr(filename);
  fasta_reader_data.nof_valid_ctg = 0;
  fasta_reader_data.min_ctg_len = min_ctg_len;

  reader = gt_fasta_reader_rec_new(str_filename);
  had_err = gt_fasta_reader_run(reader, NULL, NULL,
            gt_scaffolder_graph_count_ctg, &fasta_reader_data, err);
  gt_fasta_reader_delete(reader);
  gt_str_delete(str_filename);

  *nof_contigs = fasta_reader_data.nof_valid_ctg;

  return had_err;
}

/* parse contigs in FASTA-format and save them as vertices of
   scaffold graph */
int gt_scaffolder_parser_read_contigs(GtScaffolderGraph *graph,
                                      const char *filename,
                                      GtUword min_ctg_len,
                                      GtError *err)
{
  GtFastaReader* reader;
  GtStr *str_filename;
  GtScaffolderGraphFastaReaderData fasta_reader_data;
  int had_err;

  str_filename = gt_str_new_cstr(filename);
  fasta_reader_data.header_seq = gt_str_new();
  fasta_reader_data.nof_valid_ctg = 0;
  fasta_reader_data.min_ctg_len = min_ctg_len;
  fasta_reader_data.graph = graph;

  reader = gt_fasta_reader_rec_new(str_filename);
  had_err = gt_fasta_reader_run(reader, gt_scaffolder_graph_save_header,
            NULL, gt_scaffolder_graph_save_ctg, &fasta_reader_data, err);
  gt_fasta_reader_delete(reader);
  gt_str_delete(str_filename);
  gt_str_delete(fasta_reader_data.header_seq);
  return had_err;
}
