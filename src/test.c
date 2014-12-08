#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/init_api.h"
#include "core/types_api.h"

#include "gt_scaffolder_graph.h"
#include "gt_scaffolder_algorithms.h"

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
  char *contig_filename, *dist_filename, *astat_filename;
  int had_err = 0;

  /* initialize */
  gt_lib_init();
  /* create error object */
  err = gt_error_new();

  /* LG: gt_scaffolder_graph_test has to be adapted to new
         gt_scaffolder_graph_add_edge_ptr_to_vertex function (allocation of
         space for edges of every vertex in advance, instead of reallocation)
  */
  /* SK: err auch uebergeben */
  /* Create graph with wrapper construction function and delete it */
  gt_scaffolder_graph_test(5, 8, false, 0, false, 0, false);

  /* Create graph, only initialize vertices, delete it */
  gt_scaffolder_graph_test(5, 8, true, 0, false, 0, false);

  /* Create graph, initialize and create vertices, delete it */
  gt_scaffolder_graph_test(5, 8, true, 5, false, 0, false);

  /* Create graph, initialize vertices and edges, create vertices, delete it */
  gt_scaffolder_graph_test(5, 8, true, 5, true, 0, false);

  /* Create graph, initialize and create vertices and edges, delete it */
  gt_scaffolder_graph_test(5, 8, true, 5, true, 8, false);

  /* Create graph, initialize and create vertices and edges, print it,
     delete it */
  gt_scaffolder_graph_test(5, 8, true, 5, true, 8, true);

  if (argc == 4) {
    graph = NULL;
    contig_filename = argv[1];
    dist_filename = argv[2];
    astat_filename = argv[3];

    /* load contigs and distance information from file */
    /* LG: return had_err by function ...new_from_file? */
    graph = gt_scaffolder_graph_new_from_file(contig_filename, MIN_CONTIG_LEN,
            dist_filename, err);
    gt_scaffolder_graph_print(graph, "gt_scaffolder_parser_test_complete.dot",
                              err);
    /* load astatistics and copy number from file */
    had_err = gt_scaffolder_graph_mark_repeats(astat_filename, graph,
              COPY_NUM_CUTOFF, ASTAT_NUM_CUTOFF, err);

    if (had_err == 0) {
      gt_scaffolder_graph_print(graph,
            "gt_scaffolder_algorithms_test_filter_repeats.dot", err);
      /* mark polymorphic vertices, edges and inconsistent edges */
      had_err = gt_scaffolder_graph_filter(graph, PROBABILITY_CUTOFF,
                COPY_NUM_CUTOFF_2, OVERLAP_CUTOFF);
      if (had_err == 0)
        gt_scaffolder_graph_print(graph,
            "gt_scaffolder_algorithms_test_filter_polymorphism.dot", err);
    }
    if (had_err != 0)
      fprintf(stderr,"ERROR: %s\n",gt_error_get(err));

    gt_scaffolder_graph_delete(graph);
  }
  else
    fprintf(stderr, "Usage:<FASTA-file with contigs> <distance-file with"
           " est. distances between contigs> <astat file>\n");

  gt_error_delete(err);
  gt_lib_clean();
  return EXIT_SUCCESS;
}
