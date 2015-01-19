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

#ifndef GT_SCAFFOLDER_BAMPARSER_H
#define GT_SCAFFOLDER_BAMPARSER_H

/* for saving distance information */
typedef struct {
  char *id;
  GtWord dist;
  GtUword nof_pairs;
  float std_dev;
  bool sense;
  bool same;
} Ctg;

typedef struct {
  char *root_ctg_id;
  Ctg *ctg;
  GtUword nof_ctg;
  GtUword size;
} Record;

typedef struct {
  Record *record;
  GtUword nof_record;
  GtUword size;
} DistRecords;

void init_dist_records(DistRecords *dist);
void write_dist_records(DistRecords dist);
void delete_dist_records(DistRecords dist);

/* read paired information from bam file and corresponding hist file */
int gt_scaffolder_bamparser_read_paired_information(DistRecords *dist,
                                                    const char *bam_filename,
                                                    const char *hist_filename,
                                                    GtWord min_dist,
                                                    GtWord max_dist,
                                                    GtUword min_qual,
                                                    GtUword min_nof_pairs,
                                                    GtUword min_ref_length,
                                                    GtUword min_align,
                                                    GtError *err);
#endif