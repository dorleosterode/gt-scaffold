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

/* test parsing distance records */
int gt_scaffolder_parser_read_distances_test(const char *filename,
                                             char *output_filename,
                                             GtError *err);

/* count records and check integrity of abyss-dist-format */
int gt_scaffolder_parser_count_distances(const GtScaffolderGraph *graph,
                                               const char *file_name,
                                               GtUword *nof_distances,
                                               GtError *err);

/* parse distance information of contigs in abyss-dist-format and
   save them as edges of scaffold graph */
int gt_scaffolder_parser_read_distances(const char *filename,
                                              GtScaffolderGraph *graph,
                                              bool ismatepair,
                                              GtError *err);

/* count contigs */
int gt_scaffolder_parser_count_contigs(const char *filename,
                                       GtUword min_ctg_len,
                                       GtUword *nof_contigs,
                                       GtError *err);

/* parse contigs in FASTA-format and save them as vertices of
   scaffold graph */
int gt_scaffolder_parser_read_contigs(GtScaffolderGraph *graph,
                                      const char *filename,
                                      GtUword min_ctg_len,
                                      bool astat_is_annotated,
                                      GtError *err);
#endif
