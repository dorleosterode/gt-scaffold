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

#include "core/init_api.h"
#include "core/types_api.h"

#include "gt_scaffolder_graph.h"
#include "gt_scaffolder_algorithms.h"
#include "gt_scaffolder_parser.h"

/* adapted from SGA examples */
#define MIN_CONTIG_LEN 200

#define COPY_NUM_CUTOFF 0.3
#define ASTAT_NUM_CUTOFF 20.0

#define PROBABILITY_CUTOFF 0.01
#define COPY_NUM_CUTOFF_2 1.5
#define OVERLAP_CUTOFF 400

int main(int argc, char **argv)
{
  GtError *err;
  GtScaffolderGraph *graph;
  GtHashmap *hashmap;
  char module[32], *contig_filename, *dist_filename, *astat_filename,
       *hash_filename;
  int had_err = 0;

  if (argc == 1 || sscanf(argv[1], "%s", module) != 1) {
    fprintf(stderr,"Usage: %s <module> <arguments>\n" ,argv[0]);
    exit(EXIT_FAILURE);
  }

  /* initialize */
  gt_lib_init();
  /* create error object */
  err = gt_error_new();

  if (strcmp(module, "graph") == 0) {
    GtUword max_nof_vertices, max_nof_edges, nof_vertices, nof_edges;
    int init_vertices_tmp, init_edges_tmp, print_graph_tmp;
    bool init_vertices, init_edges, print_graph;

    if (argc != 9 ||
        sscanf(argv[2], GT_WD, &max_nof_vertices) != 1 ||
        sscanf(argv[3], GT_WD, &max_nof_edges) != 1 ||
        sscanf(argv[4], "%d", &init_vertices_tmp) != 1 ||
        sscanf(argv[5], GT_WD, &nof_vertices) != 1 ||
        sscanf(argv[6], "%i", &init_edges_tmp) != 1 ||
        sscanf(argv[7], GT_WD, &nof_edges) != 1 ||
        sscanf(argv[8], "%i", &print_graph_tmp) != 1)
    {
      fprintf(stderr, "USAGE: <max_nof_vertices> <max_nof_edges> "
              "<init_vertices> <nof_vertices> <init_edges> <nof_edges> "
               "<print_graph> %s\n", argv[0]);
      return EXIT_FAILURE;
    } else {
      init_vertices = init_vertices_tmp;
      init_edges = init_edges_tmp;
      print_graph = print_graph_tmp;

      /* Create graph with wrapper construction function and delete it */
      had_err = gt_scaffolder_graph_test(max_nof_vertices, max_nof_edges,
        init_vertices, nof_vertices, init_edges, nof_edges, print_graph, err);

      if (had_err != 0)
          fprintf(stderr,"ERROR: %s\n",gt_error_get(err));
    }
  }

  else if (strcmp(module, "parser") == 0) {
    if (argc != 3) {
      fprintf(stderr, "Usage: <DistEst file>\n");
      return EXIT_FAILURE;
    } else {
      dist_filename = argv[2];
      had_err = gt_scaffolder_parser_read_distances_test(dist_filename,
                "gt_scaffolder_parser_test_read_distances.de", err);
    }
  }

  else if (strcmp(module, "scaffold") == 0) {
    if (argc != 6) {
      fprintf(stderr, "Usage:<FASTA-file with contigs> <DistEst file> "
                      "<astat file> <hash outfile>\n");
      return EXIT_FAILURE;
    } else {
      graph = NULL;
      hashmap = NULL;
      contig_filename = argv[2];
      dist_filename = argv[3];
      astat_filename = argv[4];
      hash_filename = argv[5];

      had_err = gt_scaffolder_graph_new_from_file(&graph, contig_filename,
                MIN_CONTIG_LEN, dist_filename, &hashmap, err);

      if (had_err == 0) {
        gt_scaffolder_parser_hashmap_test(hashmap, hash_filename, err);
      }

      if (had_err == 0) {
        /* load astatistics and copy number from file */
        had_err = gt_scaffolder_graph_mark_repeats(astat_filename, graph,
                  COPY_NUM_CUTOFF, ASTAT_NUM_CUTOFF, err);
      }

      if (had_err == 0) {

        gt_scaffolder_graph_print(graph,
              "gt_scaffolder_algorithms_test_mark_repeats.dot", err);

        /* mark polymorphic vertices, edges and inconsistent edges */
        gt_scaffolder_graph_filter(graph, PROBABILITY_CUTOFF,
                  COPY_NUM_CUTOFF_2, OVERLAP_CUTOFF);
        gt_scaffolder_graph_print(graph,
              "gt_scaffolder_algorithms_test_filter.dot", err);

        /* SD: makescaffold includes _removecycles(), perform on copy? */
        gt_scaffolder_removecycles(graph);
        gt_scaffolder_graph_print(graph,
              "gt_scaffolder_algorithms_test_removecycles.dot", err);

        gt_scaffolder_makescaffold(graph);
        gt_scaffolder_graph_print(graph,
              "gt_scaffolder_algorithms_test_makescaffold.dot", err);
      }

      if (had_err == 0) {
        GtUword i;
        GtStr *ids, *seq;
        GtScaffolderGraphRecord *rec;
        GtArray *recs = gt_scaffolder_graph_iterate_scaffolds(graph);

        gt_scaffolder_graph_write_scaffold(recs, "gt_scaffolder_new_write.scaf",
          err);

        ids = gt_str_new();
        for (i = 0; i < gt_array_size(recs); i++) {
          rec = *(GtScaffolderGraphRecord **) gt_array_get(recs, i);
          gt_str_reset(ids);
          seq = gt_scaffolder_graph_generate_string(rec, ids);
          gt_str_delete(seq);
          /* deleting all recs after stringgeneration */
          gt_scaffolder_graph_record_delete(rec);
        }
        gt_array_delete(recs);
        gt_str_delete(ids);
      }

      if (had_err != 0)
        fprintf(stderr,"ERROR: %s\n",gt_error_get(err));

      gt_hashmap_delete(hashmap);
      gt_scaffolder_graph_delete(graph);
    }
  }
  else {
    fprintf(stderr,"Usage: %s <module> <arguments>\n" ,argv[0]);
    exit(EXIT_FAILURE);
  }

  gt_error_delete(err);
  gt_lib_clean();
  return had_err;
}
