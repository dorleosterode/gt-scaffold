#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gt_scaffolder_graph.h"
#include "core/types_api.h"

int main(void){
  GtError *err;

  /* initialize */
  gt_lib_init();
  /* create error object */
  err = gt_error_new();

  /* Create graph and delete it */
  gt_scaffolder_test_graph(5, 8, false, 0, false, 0, false);

  /* Create graph, only initialize vertices, delete it */
  gt_scaffolder_test_graph(5, 8, true, 0, false, 0, false);

  /* Create graph, initialize and create vertices, delete it */
  gt_scaffolder_test_graph(5, 8, true, 5, false, 0, false);

  /* Create graph, initialize vertices and edges, create vertices, delete it */
  gt_scaffolder_test_graph(5, 8, true, 5, true, 0, false);

  /* Create graph, initialize and create vertices and edges, delete it */
  gt_scaffolder_test_graph(5, 8, true, 5, true, 8, true);

  gt_error_delete(err);
  return EXIT_SUCCESS;
}
