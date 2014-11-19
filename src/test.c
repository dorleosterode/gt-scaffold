#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gt_scaffolder_graph.h"
#include "core/types_api.h"


int main(void){
  GtScaffoldGraph *graph;
 
  graph = new_graph();

  graph_show(graph);

  return EXIT_SUCCESS;
}
