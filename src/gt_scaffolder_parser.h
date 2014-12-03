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

#include "gt_scaffolder_graph.h"

#ifndef GT_SCAFFOLDER_PARSER_H
#define GT_SCAFFOLDER_PARSER_H

/* sort by lexicographic ascending order */
int gt_scaffolder_graph_vertices_compare(const void *a, const void *b);

/* count records */
int gt_scaffolder_graph_count_distances(const GtScaffolderGraph *graph,
                                               const char *file_name,
                                               GtUword *nof_distances,
                                               GtError *err);

/* parse distance information of contigs in abyss-dist-format and
   save them as edges of scaffold graph, PRECONDITION: header contains
   no commas and spaces */
/* LG: check for "mate-flag"? */
/* SK: DistParser in eigenes Modul auslagern? */
/* Bsp.: Ctg1 Ctg2+,15,10,5.1 ; Ctg3-,65,10,5.1 */
int gt_scaffolder_graph_read_distances(const char *filename,
                                              GtScaffolderGraph *graph,
                                              bool ismatepair,
                                              GtError *err);

/* counts contigs with minimum length in callback data
   (fasta reader callback function, gets called after fasta entry
   has been read) */
int gt_scaffolder_graph_count_ctg(GtUword length,
                                         void *data,
                                         GtError* err);

/* saves header to callback data
   (fasta reader callback function, gets called for each description
    of fasta entry) */
int gt_scaffolder_graph_save_header(const char *description,
                                           GtUword length,
                                           void *data, GtError *err);

/* saves header, sequence length of contig to scaffolder graph
   (fasta reader callback function, gets called after fasta entry
   has been read) */
int gt_scaffolder_graph_save_ctg(GtUword seq_length,
                                        void *data,
                                        GtError* err);

#endif
