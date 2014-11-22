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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h> /* SK: gt_file*/
#include <stdbool.h>
#include <math.h>
#include "gt_scaffolder_graph.h"
#include <core/assert_api.h>
#include "core/queue_api.h"
#include "core/types_api.h"
#include "core/minmax.h"
#include "core/ma_api.h"

/* data structure of the scaffolder graph */

GtScaffoldGraph *new_graph(GtUword nofvertices, GtUword nofedges) {
  GtScaffoldGraph *graph;

  gt_assert(nofvertices > 0);
  gt_assert(nofedges > 0);

  /* SK: Parameter besser benennen */
  graph = gt_malloc(sizeof(*graph));
  graph->vertices = gt_malloc(sizeof(*graph->vertices) * nofvertices);
  graph->edges = gt_malloc(sizeof(*graph->edges) * nofedges);
  graph->nofvertices = 0;
  graph->maxnofvertices = nofvertices;
  graph->nofedges = 0;
  graph->maxnofedges = nofedges;

  return graph;
}

void graph_add_vertex(GtScaffoldGraph *graph, const GtUword seqlen,
  const float astat, const float copynum) {
  gt_assert(graph != NULL);
  gt_assert(graph->nofvertices < graph->maxnofvertices);

  /* Knoten erstellen */
  /* SK: ID redundant */
  graph->vertices[graph->nofvertices].id = graph->nofvertices;
  graph->vertices[graph->nofvertices].seqlen = seqlen;
  graph->vertices[graph->nofvertices].astat = astat;
  graph->vertices[graph->nofvertices].copynum = copynum;
  graph->vertices[graph->nofvertices].nofedges = 0;
  /* SD: Standardstatus einfügen? Spart evtl einen Initialisierungs-Durchlauf.
  graph->vertices[graph->nofvertices].state = GS_UNVISITED; */

  graph->nofvertices++;
}

void graph_add_edge(GtScaffoldGraph *graph, const GtUword vstartID,
  const GtUword vendID, const GtWord dist, const float stddev,
  const GtUword numpairs, const bool dir, const bool comp) {

  gt_assert(graph != NULL);
  gt_assert(graph->nofedges < graph->maxnofedges);

  /* Kante erstellen */
  /* SK: ID redundant */
  graph->edges[graph->nofedges].id = graph->nofedges;
  graph->edges[graph->nofedges].start = graph->vertices + vstartID;
  graph->edges[graph->nofedges].end = graph->vertices + vendID;
  graph->edges[graph->nofedges].dist = dist;
  graph->edges[graph->nofedges].stddev = stddev;
  graph->edges[graph->nofedges].numpairs = numpairs;
  graph->edges[graph->nofedges].dir = dir;
  graph->edges[graph->nofedges].comp = comp;

  /* Kante im Startknoten eintragen */
  gt_assert(vstartID < graph->nofvertices && vendID < graph->nofvertices);
  /* SK: Erstes Malloc auslagern */
  if(graph->vertices[vstartID].nofedges == 0) {
    graph->vertices[vstartID].edges = gt_malloc( sizeof(*graph->edges) );
  }
  /* SK: realloc zu teuer? Besser: DistEst parsen und gezielt allokieren */
  else {
    graph->vertices[vstartID].edges = gt_realloc(graph->vertices[vstartID].edges,
                                              sizeof(*graph->edges) *
                                              (graph->vertices[vstartID].nofedges + 1));
  }
  graph->vertices[vstartID].edges[graph->vertices[vstartID].nofedges] =
    &graph->edges[graph->nofedges];

  graph->vertices[vstartID].nofedges++;

  graph->nofedges++;
}

/* dotfile in datei filename ausgeben */
int write_graph(const struct GtScaffoldGraph *g, const char *filename) {
  int err = 0;

  FILE *f = fopen(filename, "w");
  if (f == NULL)
    return errno;

  /* TODO: errorhandling einfuehren */
  print_graph(g, f);

  fclose(f);

  return err;
}

/* dotfile ausgeben */
void print_graph(const struct GtScaffoldGraph *g, FILE *f) {
  GtUword i;
  struct GtScaffoldGraphVertex *v;
  struct GtScaffoldGraphEdge *e;
  char *color;

  /* die erste Zeile in der Dot-Datei schreiben */
  fprintf(f, "graph {\n");

  /* alle Knoten durchgehen und schreiben */
  for (i = 0; i < g->nofvertices; i++) {
    v = &g->vertices[i];

    /* SK: const char Array für Farben verwenden, enum als Indizes */
    if (v->state == GS_REPEAT || v->state == GS_POLYMORPHIC || v->state == GS_INCONSISTENT) {
      color = "gray";
    } else if (v->state == GS_VISITED) {
      color = "red";
    } else if (v->state == GS_PROCESSED) {
      color = "green";
    } else {
      color = "black";
    }

    fprintf(f, "%lu [color=\"%s\"];\n", v->id, color);
  }

  /* alle Kanten durchgehen und schreiben */
  for (i = 0; i < g->nofedges; i++) {
    e = &g->edges[i];

    if (e->state == GS_REPEAT || e->state == GS_POLYMORPHIC || e->state == GS_INCONSISTENT) {
      color = "gray";
    } else if (e->state == GS_VISITED) {
      color = "red";
    } else if (e->state == GS_PROCESSED) {
      color = "green";
    } else {
      color = "black";
    }

    fprintf(f, "%lu -- %lu [color=\"%s\" label=\"%lu\"];\n", e->start->id, e->end->id, color, e->dist);
  }

  /* die letzte Zeile schreiben */
  fprintf(f, "}\n");
}

/* Einlesen der Contigs im FASTA-Format und Speicherung der Contigs
   mit deren Header und Laenge als Knoten des Scaffold Graphen */
static void read_contigs(const char *filename, GtScaffoldGraph *graph,
                         GtUword minctglen) {
  GtUword num_of_sequences = 0, entryid;
  char currentchar = '\0';
  FILE *infile;
  bool firstseq = true, firstline = true;
  char buffer[BUFSIZ], itembuf[BUFSIZ];

  infile = fopen(filename, "rb");
  if (infile == NULL) {
     fprintf(stderr,"ERROR: Can not open file %s !\n",filename);
     exit(EXIT_FAILURE);
  }
  /* Zaehlen der Sequenzen */
  while (EOF != (currentchar = fgetc(infile))) {
    if (currentchar == '>')
    num_of_sequences++;
  }

  if (num_of_sequences == 0) {
     fprintf(stderr,"ERROR: No sequences !\n");
     exit(EXIT_FAILURE);
  }

  entryid = 0;
  graph->nofvertices = num_of_sequences;
  graph->vertices = gt_malloc(sizeof(*graph->vertices) * num_of_sequences);

  /* Lesen der Datei zeilenweise */
  rewind(infile);
  while (fgets(buffer, BUFSIZ, infile) != NULL) {
    if (sscanf(buffer, ">%s", itembuf) == 1) {
      /* Lesen des Header */
      if (!firstseq) {
        /* Wenn Laenge des aktuellen Contigs unter Schwellenwert minctglen,
           zugehoerigen Header loeschen, ansonsten Contig-Zaehler erhoehen */
        if (graph->vertices[entryid].seqlen < minctglen)
          free(graph->vertices[entryid].header);
        else
          entryid++;

        firstline = true;
      }
      else
        firstseq = false;
      /* Speicherung neuen Header */
      graph->vertices[entryid].header = gt_malloc(sizeof
            (*graph->vertices[entryid].header) * strlen(itembuf) + 1);
      strncpy(graph->vertices[entryid].header, itembuf, strlen(itembuf));
    }
    else if (sscanf(buffer, "%s", itembuf) == 1) {
      /* Lesen eines Kommentars */
      if (itembuf[0] == ';')
        continue;
      /* Lesen einer Sequenz */
      if (firstline) {
        /* erste Zeile der Sequenz */
        graph->vertices[entryid].seqlen = strlen(itembuf);
        firstline = false;
      }
      else
        graph->vertices[entryid].seqlen += strlen(itembuf);
    }
  }

  fclose(infile);
}


/* Erzeugung des Scaffold Graphen */
GtScaffoldGraph *gt_scaffolder_graph_new_from_file(const char *ctgfilename,
                 GtUword minctglen) {
  GtScaffoldGraph *graph;

  graph = gt_malloc(sizeof(*graph));
  if (graph == NULL) {
     fprintf(stderr,"ERROR: Memory allocation failed!\n");
     exit(EXIT_FAILURE);
  }
  /* Einlesen der Contigs im FASTA-Format und Speicherung der Contigs
     mit deren Header und Laenge als Knoten des Scaffold Graphen */
  read_contigs(ctgfilename, graph, minctglen);
  /* Einlesen der Distanzinformationen der Contigs im Abyss-Dist-Format
     und Speicherung als Kanten des Scaffold Graphen */
  //read_distances(distfilename, graph);

  return graph;
}



/* Pruefung auf eindeutige Ordnung der Kanten edge1, edge 2 mit Wahrscheinlichkeit
   cutoff */
static bool gt_scaffolder_graph_ambiguousorder(const GtScaffoldGraphEdge *edge1,
      const GtScaffoldGraphEdge *edge2, float cutoff) {

  float expval, variance, interval, prob12, prob21;

  expval = edge1->dist - edge2->dist;
  /* sichere Multiplikation, Division */
  /* SK: pow nicht verwenden, stattdessen mit sich selbst multiplizieren */
  variance = 2 * pow(edge1->stddev,2) - pow(edge2->stddev,2);
  interval = -expval / sqrt(variance);
  prob12 = 0.5 * (1 + erf(interval) );
  prob21 = 1.0 - prob12;

  return (prob12 <= cutoff && prob21 <= cutoff);
}

/* Makierung polymorpher Kanten/Knoten und inkonsistenter Kanten im
   Scaffold Graphen */
int gt_scaffolder_graph_filtering(GtScaffoldGraph *graph, float pcutoff,
    float cncutoff, GtUword ocutoff) {
  GtScaffoldGraphVertex *vertex, *polymorphic_vertex;
  GtScaffoldGraphEdge *edge1, *edge2;
  GtUword vid, eid_1, eid_2, eid_3, overlap;
  GtUword maxoverlap = 0;
  float sum_copynum;
  unsigned int dir; /* int statt bool, weil Iteration bislang nicht möglich */
  GtWord intersect_start, intersect_end;
  int had_err = 0;

  /* Iteration ueber alle Knoten */
  for (vid = 0; vid < graph->nofvertices; vid++) {
    vertex = &graph->vertices[vid];
    /* Iteration ueber alle Richtungen (Sense/Antisense) */
    for (dir = 0; dir < 2; dir++) {
      /* Iteration ueber alle Kantenpaare */
      for (eid_1 = 0; eid_1 < vertex->nofedges; eid_1++) {
        for (eid_2 = eid_1+1; eid_2 < vertex->nofedges; eid_2++) {
          edge1 = vertex->edges[eid_1];
          edge2 = vertex->edges[eid_2];
          if (edge1->dir == dir && edge2->dir == dir) {
            /* Pruefung des Kantenpaares edge1, edge2 auf polymorphe Merkmale */
            sum_copynum = edge1->end->copynum + edge2->end->copynum;
            if (gt_scaffolder_graph_ambiguousorder(edge1, edge2, pcutoff) &&
                sum_copynum < cncutoff) {
              /* Markierung Endknoten mit kleinerer estCopyNum als polymorph */
              if (edge1->end->copynum < edge2->end->copynum)
                polymorphic_vertex = edge1->end;
              else
                polymorphic_vertex = edge2->end;
              /* Markierung aller Sense- /Antisensekanten des polymorphen Knoten
                 als polymorph
              SK: Prüfen ob Knoten polymorph ist */
              for (eid_3 = 0; eid_3 < polymorphic_vertex->nofedges; eid_3++)
                polymorphic_vertex->edges[eid_3]->state = GS_POLYMORPHIC;
              polymorphic_vertex->state = GS_POLYMORPHIC;
            }
            /* SD: Nur das erste Paar polymoprh markieren? */
          }
        }
      }

      /* keine Pruefung auf inkonsistente Kanten fuer polymorphe Knoten
         notwendig */
      if (vertex->state == GS_POLYMORPHIC)
        break;
      /* Iteration ueber alle nicht-polymorphen Kantenpaare in derselben Richtung */
      for (eid_1 = 0; eid_1 < vertex->nofedges; eid_1++) {
        for (eid_2 = eid_1+1; eid_2 < vertex->nofedges; eid_2++) {
          edge1 = vertex->edges[eid_1];
          edge2 = vertex->edges[eid_2];
          if (edge1->dir == dir && edge2->dir == dir &&
            edge1->state != GS_POLYMORPHIC && edge2->state != GS_POLYMORPHIC) {
            /* TODO: calculate overlapp
            SK: In eigene Funktion auslagern */
            if (edge2->dist > (edge1->end->seqlen - 1) ||
            edge1->dist > (edge2->end->seqlen - 1)) {
              intersect_start = MAX(edge1->dist, edge2->dist);
              intersect_end = MAX(edge1->end->seqlen - 1, edge2->end->seqlen - 1);
              overlap = intersect_end - intersect_start + 1;
              if (overlap > maxoverlap) {
                maxoverlap = overlap;
              }
            }
          }
        }
      }

     /* Pruefung aller Kantenpaare auf maximale Ueberlappung > ocutoff */
      if (maxoverlap > ocutoff) {
        for (eid_1 = 0; eid_1 < vertex->nofedges; eid_1++)
         vertex->edges[eid_1]->state = GS_INCONSISTENT;
      }
    }
  }
  /* SD: Knoten & Kanten werden nicht gelöscht, sondern der Einfachheit halber später geprüft*/
  return had_err;
}

/* Ueberpruefung ob Knoten terminal ist, d.h. nur sense oder antisense-Kanten
   vorliegen */
static bool gt_scaffolder_graph_isterminal(const GtScaffoldGraphVertex *vertex) {
  GtUword sense = 0, antisense = 0, eid;
  GtScaffoldGraphEdge *edge;

  /* Nicht zählen, Schleife abbrechen ueber != prev_sense */
  for (eid = 0; eid < vertex->nofedges; eid++) {
    edge = vertex->edges[eid];
    if (edge->dir == SENSE)
      sense++;
    else
      antisense++;
  }

  return ((sense == 0 && antisense != 0) || (sense != 0 && antisense == 0));
}

/* Entfernung von Zyklen
static void gt_scaffolder_removecycles(GtScaffoldGraph *graph) {
}*/

/* Erstellung eines neuen Walks */
static Walk *gt_scaffolder_walk_new(void) {
  Walk *walk;

  walk = gt_malloc(sizeof(*walk));
  walk->totalcontiglen = 0;
  walk->size = 0;
  walk->nofedges = 0;
  return walk;
}

/* Loeschen eines Walks */
static void gt_scaffolder_walk_delete(Walk *walk) {
  if (walk != NULL)
    free(walk->edges);
  free(walk);
}

/* Ausgabe der Contig-Gesamtlaenge eines Walks */
static GtUword gt_scaffolder_walk_getlength(Walk *walk) {
  return walk->totalcontiglen;
}

/* Hinzufuegen einer Kante zum Walk */
static void gt_scaffolder_walk_addegde(Walk *walk, GtScaffoldGraphEdge *edge) {
  if (walk->size == walk->nofedges) {
    /* SK: 10 als Konstante definieren, const unsigned long increment_size 2er Potenz */
    walk->size += 10;
    walk->edges = gt_realloc(walk->edges, walk->size*sizeof(*walk->edges));
  }
  walk->edges[walk->nofedges] = edge;
  walk->totalcontiglen += edge->end->seqlen;
  walk->nofedges++;
}

/* Konstruktion des Scaffolds mit groesster Contig-Gesamtlaenge */
void gt_scaffolder_makescaffold(GtScaffoldGraph *graph) {
  GtScaffoldGraphVertex *vertex, *currentvertex, *nextvertex, *nextendvertex,
                        *endvertex;
  GtScaffoldGraphEdge *edge, *nextedge, *reverseedge, **edgemap;
  GtUword vid, eid, ccnumber, lengthcwalk, lengthbestwalk;
  GtQueue *vqueue, *wqueue;
  float distance, *distancemap;
  Pair *pair, *updatepair;
  bool dir;
  Walk *bestwalk, *currentwalk;

  /* Entfernung von Zyklen
     gt_scaffolder_removecycles(graph); */

  /* Iteration ueber alle Knoten, Makierung aller Knoten als unbesucht */
  for (vid = 0; vid < graph->nofvertices; vid++) {
    vertex = &graph->vertices[vid];
    /* SD: Existieren Repeat-Knoten, nach Prozessierung der AStatistik? */
    if (vertex->state == GS_REPEAT || vertex->state == GS_POLYMORPHIC)
      continue; /* SK: negieren statt continue */
    vertex->state = GS_UNVISITED;
  }

 /* BFS-Traversierung durch Zusammenhangskomponenten des Graphen,
    siehe GraphSearchTree.h */
  ccnumber = 0;
  vqueue = gt_queue_new();
  pair = gt_malloc(sizeof(*pair));
  updatepair = gt_malloc(sizeof(*updatepair));
  /* SK: Mit GtWord_Min / erwarteter Genomlaenge statt 0 initialisieren */
  distancemap = calloc(graph->nofvertices, sizeof(*distancemap));
  edgemap = gt_malloc(sizeof(*edgemap)*graph->nofvertices);

  for (vid = 0; vid < graph->nofvertices; vid++) {
    vertex = &graph->vertices[vid];
    if (vertex->state == GS_REPEAT || vertex->state == GS_POLYMORPHIC ||
        vertex->state == GS_VISITED)
      continue;
    ccnumber += 1;
    vertex->state = GS_PROCESSED;
    gt_queue_add(vqueue, vertex); // TODO: hier moechte evtl nicht die adresse von dem pointer uebergeben werden

    while (gt_queue_size(vqueue) != 0) {
      currentvertex = (GtScaffoldGraphVertex*)gt_queue_get(vqueue);
      //currentvertex.cc = ccnumber;

      /* BFS-Traversierung innerhalb aktueller Zusammenhangskomponente
         ausgehend von terminalen Knoten zu terminalen Knoten */
      lengthbestwalk = 0;
      bestwalk = gt_scaffolder_walk_new();
      /* SK: In weitere Funktion mit eigener Queue auslagern */
      wqueue = gt_queue_new();

      if (gt_scaffolder_graph_isterminal(vertex)) {
        dir = vertex->edges[0]->dir;
        for (eid = 0; eid < vertex->nofedges; eid++) {
          edge = vertex->edges[eid];
          endvertex = edge->end;
          pair->edge = edge;
          pair->dist = edge->dist;

          /* SK: genometools hashes verwenden, Dichte evaluieren
             SK: DistEst beim Einlesen prüfen */
          distancemap[endvertex->id] = edge->dist;
          edgemap[endvertex->id] = edge;

          gt_queue_add(wqueue, pair);
        }
        while(gt_queue_size(wqueue) != 0) {
          pair = (Pair*)gt_queue_get(wqueue);
          edge = pair->edge;
          endvertex = edge->end;

          /* Ruecktraversierung durch EdgeMap wenn terminaler Knoten erreicht,
             Konstruktion des Walks  */
          if (gt_scaffolder_graph_isterminal(endvertex)) {

            currentwalk = gt_scaffolder_walk_new();
            currentvertex = endvertex;
            while (currentvertex->id != vertex->id) {
              reverseedge = edgemap[currentvertex->id];
             /* Start NICHT end */
              currentvertex = reverseedge->end;
             /* Speicherung des aktuellen Walks */
             /* Kante vorher duplizieren */
              gt_scaffolder_walk_addegde(currentwalk, reverseedge);
            }

            /* Ermittelung Contig-Gesamtlaenge des aktuellen Walks  */
            lengthcwalk = gt_scaffolder_walk_getlength(currentwalk);
            if ( lengthcwalk > lengthbestwalk) {
              gt_scaffolder_walk_delete(bestwalk);
              bestwalk = currentwalk;
              lengthbestwalk = lengthcwalk;
            }
            else
              gt_scaffolder_walk_delete(currentwalk);
          }

          /* SD: Terminal Set implementieren, bestWalk über Rücktraversierung */
          for (eid = 0; eid < endvertex->nofedges; eid++) {
            nextedge = endvertex->edges[eid];
            if (nextedge->dir == dir) {
              nextendvertex = nextedge->end;
              distance = pair->dist + nextedge->dist;

              if (distancemap[nextendvertex->id] == 0 ||
                  distancemap[nextendvertex->id] > distance) {
                distancemap[nextendvertex->id] = distance;
                edgemap[nextendvertex->id] = nextedge;
                updatepair->edge = nextedge;
                updatepair->dist = distance;
                gt_queue_add(wqueue, updatepair);
              }
            }
          }
        }
      }


      currentvertex->state = GS_VISITED;
      for (eid = 0; eid < currentvertex->nofedges; eid++) {
        edge = currentvertex->edges[eid];
        nextvertex = edge->end;
        if (vertex->state == GS_REPEAT || vertex->state == GS_POLYMORPHIC)
          continue;
        if (nextvertex->state == GS_UNVISITED) {
          nextvertex->state = GS_PROCESSED;
          gt_queue_add(vqueue, nextvertex);
        }
      }
    }
  }
}



