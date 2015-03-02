/*
  Copyright (c) 2015 Lukas Götz
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

/* data type for saving paired information
   between contigs according to abyss dist format */
typedef struct {
  GtStr *id;
  GtWord dist;
  GtUword nof_pairs;
  float std_dev;
  bool sense;
  bool same;
} Ctg;

typedef struct {
  GtStr *root_ctg_id;
  Ctg *ctg;
  GtUword nof_ctg;
  GtUword size;
} Record;

typedef struct {
  Record *record;
  GtUword nof_record;
  GtUword size;
} DistRecords;

/* initialize distance records */
DistRecords *gt_scaffolder_bamparser_init_dist_records(void);

/* print distance records in abyss dist format into file filename */
int gt_scaffolder_bamparser_print_dist_records(const DistRecords *dist,
                                               const char *filename,
                                               GtError *err);

/* delete distance records */
void gt_scaffolder_bamparser_delete_dist_records(DistRecords *dist);

/* read paired information from bam file */
int gt_scaffolder_bamparser_read_paired_information(DistRecords *dist,
                                                    const char *bam_filename,
                                                    GtWord min_dist,
                                                    GtWord max_dist,
                                                    GtUword min_qual,
                                                    GtUword min_nof_pairs,
                                                    GtUword min_ref_length,
                                                    GtUword min_align,
                                                    GtError *err);
#endif
