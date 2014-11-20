#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gt_scaffolder_graph.h"
#include "core/types_api.h"


int main(void){
  GtScaffoldGraph *graph;
  char outfile[] = "foo.dot";

  graph = new_graph();
  write_graph(graph, outfile);

  return EXIT_SUCCESS;
}
