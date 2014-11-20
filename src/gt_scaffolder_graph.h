/* Gt-Scaffolder auf Basis des SGA-Algorithmus
   written by Dorle Osterode, Lukas Goetz, Stefan Dang
   <stefan.dang@studium.uni-hamburg.de>
   Encoding: UTF-8
*/
#include "core/types_api.h"


#ifndef GT_SCAFFOLDER_GRAPH_H
#define GT_SCAFFOLDER_GRAPH_H

typedef struct GtScaffoldGraph GtScaffoldGraph;
typedef struct GtScaffoldGraphVertex GtScaffoldGraphVertex;
typedef struct GtScaffoldGraphEdge GtScaffoldGraphEdge;


/* Testfunktionen fuer Graph-Datenstruktur */
GtScaffoldGraph *new_graph(void);
int write_graph(struct GtScaffoldGraph *g, char *filename);
void print_graph(struct GtScaffoldGraph *g, FILE *f);

/* the scaffolder_graph tool */
GtScaffoldGraph *gt_scaffolder_graph_new_from_file(const char *filename, int *had_err);
int gt_scaffolder_graph_filtering(GtScaffoldGraph *graph, float pcutoff,
    float cncutoff, GtUword ocutoff);
void gt_scaffolder_makescaffold(GtScaffoldGraph *graph);


#endif
