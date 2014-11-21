/* Gt-Scaffolder auf Basis des SGA-Algorithmus
   written by Dorle Osterode, Lukas Goetz, Stefan Dang
   <stefan.dang@studium.uni-hamburg.de>
   Encoding: UTF-8
*/
#include <stdio.h>
#include "core/types_api.h"


#ifndef GT_SCAFFOLDER_GRAPH_H
#define GT_SCAFFOLDER_GRAPH_H

/* SK: GraphItemState oder GtScaffolderGraphName
       Gt-Namenskonvention für Zustände einhalten (docs/ oder manuals/developermanual)
       Automatische Prüfung durch scripts/src_check */
typedef enum { GS_UNVISITED, GS_REPEAT, GS_POLYMORPHIC, GS_INCONSISTENT,
               GS_VISITED, GS_PROCESSED } GraphState;

/* Knoten */
typedef struct GtScaffoldGraphVertex {
  /* eindeutige ID fuer den Knoten */
  GtUword id;
  /* Headersequenz des zugehoerigen Contigs */
  char *header;
  /* Laenge der Sequenz, die der Contig darstellt */
  GtUword seqlen;
  /* Wert der A-Statistik, um Contigs als REPEAT oder UNIQUE
     klassifizieren zu koennen;
     in Genom-Tools vom Typ float */
  float astat;
  /* abgeschaetzte Anzahl an Vorkommen des Contigs im Genom */
  float copynum;
  GtUword nofedges;
  /* Sammlung von Kanten, die von dem Contig abgehen */
  /* Speicherung der ID der Kanten statt Ptr auf Kanten? */
  struct GtScaffoldGraphEdge **edges;
  /* Zustand des Knotens, außer GS_INCONSISTENT */
  GraphState state;
} GtScaffoldGraphVertex;

/* Kante */
typedef struct GtScaffoldGraphEdge {
  /* eindeutige ID fuer die Kante */
  GtUword id;
  /* Pointer zu dem Knoten, zu dem die Kante fuehrt */
  struct GtScaffoldGraphVertex *end;
  /* Pointer zu dem Knoten, von dem die Kante kommt */
  struct GtScaffoldGraphVertex *start;
  /* Abschaetzung der Entfernung der verbundenen Contigs */
  GtWord dist;
  /* Standardabweichung von der abgeschaetzten Entfernung */
  float stddev;
  /* Anzahl der Distanzinformationen, die ein Anzeichen fuer die
  Verbindung der Contigs geben */
  GtUword numpairs;
  /* Zustand der Kante */
  GraphState state;
  /* enthaelt die Richtung (Sense, Antisense) und welche
     Straenge die paired-Information enthalten (die gleiche
     Richtung oder das Reverse) */
  bool dir; /* SK: forward */
  bool comp; /* SK: sense */
} GtScaffoldGraphEdge;

/* Graph */
typedef struct GtScaffoldGraph {
  struct GtScaffoldGraphVertex *vertices;
  GtUword nofvertices;
  GtUword maxnofvertices;
  struct GtScaffoldGraphEdge *edges;
  GtUword nofedges;
  GtUword maxnofedges;
} GtScaffoldGraph;

/* Datenstruktur fuer Queue-Element */
/* SD: Benötigt? Redundant. */
typedef struct Pair {
  float dist;
  struct GtScaffoldGraphEdge *edge;
} Pair;

/* Datenstruktur fuer Walk
SK: Umbenennen! */
typedef struct Walk {
  GtUword nofedges;
  GtUword size;
  GtUword totalcontiglen;
  struct GtScaffoldGraphEdge **edges;
} Walk;

/* Datenstruktur fuer DistanceMap */
/* Datenstruktur fuer EdgeMap */


/* Grundlegende Graphfunktionen
SK: gt_scaffold_graph_new, const auf nicht-Pointer entfernen */
GtScaffoldGraph *new_graph(const GtUword nofvertices, const GtUword nofedges);
void graph_add_vertex(GtScaffoldGraph *graph, const GtUword seqlen,
  const float astat, const float copynum);
void graph_add_edge(GtScaffoldGraph *graph, const GtUword vstartID,
  const GtUword vendID, const GtWord dist, const float stddev,
  const GtUword numpairs, const bool dir, const bool comp);
/* SK: gt_scaffold_graph_delete */

/* Darstellung des Graphen
SK: gt_file benutzen, print_graph_generic / print_graph
SK: sga Format unterstützen asgq (?) */
int write_graph(const struct GtScaffoldGraph *g, const char *filename);
void print_graph(const struct GtScaffoldGraph *g, FILE *f);

/* the scaffolder_graph tool
SK: Fasta-Iterator (core/fastareaderseqiterator) */
GtScaffoldGraph *gt_scaffolder_graph_new_from_file(const char *ctgfilename,
              GtUword minctglen); /* SK: GtError Objekt */
int gt_scaffolder_graph_filtering(GtScaffoldGraph *graph, float pcutoff,
    float cncutoff, GtUword ocutoff); /* SK: GtError Objekt; nicht benötigt */
void gt_scaffolder_makescaffold(GtScaffoldGraph *graph); /* SK: nicht-öffentlich? */


#endif
