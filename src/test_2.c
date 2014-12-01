#include <stdio.h>
#include "core/types_api.h"
#include "core/error.h"
#include "gt_scaffolder_graph.h"

/* adapted from SGA examples */
const MIN_CONTIG_LEN = 200;

int main(int argc, char **argv)
{
  GtError *err;
  GtScaffolderGraph *graph;
  GtUword min_ctg_len;

  graph = NULL;
  min_ctg_len = 10;
  /* initialize */
  gt_lib_init();
  /* create error object */
  err = gt_error_new();

  if (argc == 3)
  {
    /* load contigs and distance information from file */
    graph = gt_scaffolder_graph_new_from_file(argv[1], MIN_CONTIG_LEN,
            argv[2], err);
    gt_scaffolder_graph_print(graph, "output.dot", err);
  }
  else
    printf("Usage:<FASTA-file with contigs> <distance-file with"
           " est. distances between contigs>\n");
  
  gt_scaffolder_graph_delete(graph);
  gt_error_delete(err);
  gt_lib_clean();

  return 0;
}
