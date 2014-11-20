#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gt_scaffolder_graph.h"
#include "core/types_api.h"

/* Datenstruktur Scaffold-Graph */

#define ANTISENSE 1
#define SENSE 0
#define REVERSE 1
#define SAME 0


int main(void){
  GtScaffoldGraph *graph;
  char outfile[] = "foo.dot";

  /* Baue das Haus vom Nikolaus, 5 Knoten, 8 Kanten: */
  graph = new_graph(5, 8);

  graph_add_vertex(graph, 100, 20, 40);
  graph_add_vertex(graph, 100, 20, 40);
  graph_add_vertex(graph, 100, 20, 40);
  graph_add_vertex(graph, 100, 20, 40);
  graph_add_vertex(graph, 100, 20, 40);

  graph_add_edge(graph, 0, 1, 2, 1.5, 4, SENSE, SAME);
  graph_add_edge(graph, 0, 2, 2, 1.5, 4, SENSE, SAME);
  graph_add_edge(graph, 0, 3, 2, 1.5, 4, SENSE, SAME);
  graph_add_edge(graph, 1, 2, 2, 1.5, 4, SENSE, SAME);
  graph_add_edge(graph, 1, 3, 2, 1.5, 4, SENSE, SAME);
  graph_add_edge(graph, 2, 3, 2, 1.5, 4, SENSE, SAME);
  graph_add_edge(graph, 2, 4, 2, 1.5, 4, SENSE, SAME);
  graph_add_edge(graph, 3, 4, 2, 1.5, 4, SENSE, SAME);

  write_graph(graph, outfile);

  return EXIT_SUCCESS;
}
