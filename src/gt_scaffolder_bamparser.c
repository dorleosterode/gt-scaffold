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

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "core/str_api.h"
#include "core/unused_api.h"
#include "core/alphabet_api.h"
#include "core/minmax.h"
#include "core/error.h"
#include "core/ma_api.h"
#include "core/types_api.h"
#include "core/hashmap_api.h"
#include "extended/sam_alignment.h"
#include "extended/samfile_iterator.h"

/* for definition of functions missing in genometools */
#include "external/samtools-0.1.18/bam.h"
#include "external/samtools-0.1.18/sam.h"

#include "gt_scaffolder_bamparser.h"
#include "gt_scaffolder_parser.h"

#define INCREMENT_SIZE 1024
#define INCREMENT_SIZE_2 64

struct GtSamAlignment{
  bam1_t       *s_alignment;
  GtAlphabet   *alphabet;
  GtUchar      *seq_buffer,
               *qual_buffer;
  GtUword s_bufsize,
                q_bufsize,
                rightmost;
};

struct GtSamfileIterator {
  GtAlphabet     *alphabet;
  GtSamAlignment *current_alignment;
  char           *filename,
                 *mode;
  samfile_t      *samfile;
  void           *aux;
  GtUword   ref_count;
};

/* data type for saving information about
   reads aligned to some contig  */
typedef struct Read {
  bool is_reverse;
  bool is_mreverse;
  bool is_unmapped;
  bool is_munmapped;
  bool is_valid;
  int tid;
  int mtid;
  GtUword pos;
  GtWord ref_read_start;
  GtWord mref_read_start;
  GtWord isize;
  GtStr *query_name;
  GtStr *ref_name;
  GtStr *mref_name;
  GtUword ref_len;
  GtUword mref_len;
} Read;

typedef struct ReadSet {
  Read *read;
  GtUword nof;
  GtUword size;
} ReadSet;

typedef struct FragmentData {
  GtWord *frag_pos;
  GtUword *frag_size;
  GtUword nof_frag_pos;
  GtUword size_frag_pos;
  GtUword nof_frag_size;
  GtUword size_frag_size;
  GtUword ma;
} FragmentData;

typedef struct HistogramData {
  GtHashmap *hash_map;
  GtWord *keys;
  GtUword *values;
  GtUword value_sum;
  GtUword nof_pairs;
  GtUword size;
  double mean;
  double std_dev;
  GtWord min;
  GtWord max;
  bool lib_rf;
} HistogramData;

typedef struct PmfData {
  double *dist;
  GtUword nof;
  double mean;
  double std_dev;
  double minp;
} PmfData;

typedef struct CompareData {
  GtWord dist;
  double likelihood;
  GtUword nof_pairs;
  PmfData pmf_data;
} CompareData;

/* initialize distance records */
DistRecords *gt_scaffolder_bamparser_init_dist_records(void) {
  DistRecords *dist;

  dist = gt_malloc(sizeof (*dist));
  dist->nof_record = 0;
  dist->size = 0;
  dist->record = NULL;
  return dist;
}

/* print distance records in abyss dist format into file filename */
int gt_scaffolder_bamparser_print_dist_records(const DistRecords *dist,
                                               const char *filename,
                                               GtError *err) {
  GtUword index, index_2;
  bool set_dir;
  GtFile *file;
  int had_err = 0;

  file = gt_file_new(filename, "w", err);
  if (file == NULL) {
    gt_error_set (err , "distance file can not be created");
    had_err = -1;
  }

  if (!had_err) {
    for (index = 0; index < dist->nof_record; index++) {
      gt_file_xprintf(file, "%s",gt_str_get(dist->record[index].root_ctg_id));
      set_dir = false;
      for (index_2 = 0; index_2 < dist->record[index].nof_ctg; index_2++) {
        if (dist->record[index].ctg[index_2].same && !set_dir) {
          gt_file_xprintf(file," ;");
          set_dir = true;
        }

        gt_file_xprintf(file," %s%c," GT_WD "," GT_WU ",%.1f",
                       gt_str_get(dist->record[index].ctg[index_2].id),
                       dist->record[index].ctg[index_2].sense ? '+' : '-',
                       dist->record[index].ctg[index_2].dist,
                       dist->record[index].ctg[index_2].nof_pairs,
                       dist->record[index].ctg[index_2].std_dev);
      }
      if (!set_dir)
        gt_file_xprintf(file," ;");
      gt_file_xprintf(file,"\n");
    }
    gt_file_delete(file);
  }
  return had_err;
}

/* initialize distance records */
void gt_scaffolder_bamparser_delete_dist_records(DistRecords *dist) {
 GtUword index, index_2;

  for (index = 0; index < dist->size; index++) {
    for (index_2 = 0; index_2 < dist->record[index].nof_ctg; index_2++)
      gt_free(dist->record[index].ctg[index_2].id);
      gt_free(dist->record[index].ctg);
      gt_free(dist->record[index].root_ctg_id);
    }
  gt_free(dist->record);
  gt_free(dist);
}

/* create new distance record */
static void create_dist_record(DistRecords *dist_records,
                               GtStr *root_ctg_id) {
  GtUword index;

  /* resize */
  if (dist_records->nof_record == dist_records->size) {
    dist_records->size += INCREMENT_SIZE_2;
    dist_records->record = gt_realloc(dist_records->record,
           sizeof (*(dist_records->record)) * dist_records->size);
    for (index = dist_records->nof_record; index < dist_records->size;
         index++) {
      dist_records->record[index].size = 0;
      dist_records->record[index].nof_ctg = 0;
      dist_records->record[index].ctg = NULL;
      dist_records->record[index].root_ctg_id = NULL;
    }
  }

  dist_records->record[dist_records->nof_record].root_ctg_id = root_ctg_id;
  dist_records->nof_record++;
}

/* add contig to existing distance record */
static void add_contig_dist_record(DistRecords *dist_records,
                                  GtStr *ctg_id,
                                  double std_dev,
                                  GtWord dist,
                                  GtUword nof_pairs,
                                  bool sense,
                                  bool same) {
  GtUword current_record, next_free_ctg;

  current_record = dist_records->nof_record-1;
  next_free_ctg = dist_records->record[current_record].nof_ctg;

  /* resize */
  if (dist_records->record[current_record].nof_ctg ==
      dist_records->record[current_record].size) {
    dist_records->record[current_record].size += INCREMENT_SIZE_2;
    dist_records->record[current_record].ctg = gt_realloc(
           dist_records->record[current_record].ctg,
           sizeof (*(dist_records->record[current_record].ctg)) *
           dist_records->record[current_record].size);
  }

  dist_records->record[current_record].ctg[next_free_ctg].id = ctg_id;
  dist_records->record[current_record].ctg[next_free_ctg].std_dev =
                                              std_dev;
  dist_records->record[current_record].ctg[next_free_ctg].dist = dist;
  dist_records->record[current_record].ctg[next_free_ctg].nof_pairs =
                                              nof_pairs;
  dist_records->record[current_record].ctg[next_free_ctg].sense = sense;
  dist_records->record[current_record].ctg[next_free_ctg].same = same;
  dist_records->record[current_record].nof_ctg++;
}

static void calc_stat_of_histogram(HistogramData *histogram_data) {
  GtUword index, value, value_sum;
  GtWord key, squares, total, min, max;
  double variance;

  min = GT_WORD_MAX;
  max = GT_WORD_MIN;

  squares = 0;
  total = 0;
  value_sum = 0;

  for (index = 0; index < histogram_data->nof_pairs; index++) {
    key = histogram_data->keys[index];
    value = histogram_data->values[index];
    /* ignore removed keys */
    if (key != GT_WORD_MAX) {
      value_sum += value;
      total += key * value;
      squares += key * key * value;
      if (min > key)
        min = key;
      if (max < key)
        max = key;
    }
  }
  histogram_data->value_sum = value_sum;
  histogram_data->mean = (double)total / value_sum;
  variance = (squares - ((double)total * total / value_sum)) / value_sum;
  histogram_data->std_dev = sqrt(variance);
  histogram_data->min = min;
  histogram_data->max = max;
}

/* remove outliers in histogram by assigning GT_WORD_MAX value */
static void remove_outliers_in_histogram(HistogramData *histogram_data) {
  GtUword index, value_sum;
  GtWord q1, q3, l_value, u_value, *key;

  value_sum = 0;
  for (index = 0; index < histogram_data->nof_pairs; index++) {
    if (histogram_data->keys[index] != GT_WORD_MAX) {
      value_sum += histogram_data->values[index];
      if (value_sum >= ceil(0.25 * histogram_data->value_sum)) {
        index++;
        break;
      }
    }
  }
  q1 = histogram_data->keys[index-1];

  value_sum = 0;
  for (index = 0; index < histogram_data->nof_pairs; index++) {
    if (histogram_data->keys[index] != GT_WORD_MAX) {
      value_sum += histogram_data->values[index];
      if (value_sum >= ceil(0.75 * histogram_data->value_sum)) {
        index++;
        break;
      }
    }
  }
  q3 = histogram_data->keys[index-1];

  l_value = q1 - 20 * (q3 - q1);
  u_value = q3 + 20 * (q3 - q1);

  for (key = histogram_data->keys;
       key < histogram_data->keys + histogram_data->nof_pairs; key++) {
    if (*key != GT_WORD_MAX && (*key < l_value || *key > u_value))
      *key = GT_WORD_MAX;
  }
}

/* trim off the bottom fraction/2 and top fraction/2 data points by
   assigning GT_WORD_MAX value */
static void trim_fraction_of_histogram(HistogramData *histogram_data,
                                       double fraction) {
  GtUword index;
  double low_cutoff = fraction / 2;
  double high_cutoff = 1.0 - fraction / 2;
  double cumulative = 0, temp_total;

  for (index = 0; index < histogram_data->nof_pairs; index++) {
    if (histogram_data->keys[index] != GT_WORD_MAX) {
      temp_total = cumulative + (double)histogram_data->values[index]
                              / histogram_data->value_sum;
      if (temp_total <= low_cutoff || cumulative >= high_cutoff)
        histogram_data->keys[index] = GT_WORD_MAX;
      cumulative = temp_total;
    }
  }
}

/* remove noise from the histogram (noise is defined as a
   sample x where h[x-1] == NULL && h[x+1] == NULL),
   two elements have at least to remain */
/* precondition: keys are sorted */
static void remove_noise_in_histogram(HistogramData *histogram_data) {
  GtWord *key;
  GtUword index, nof_removed, nof;

  nof = histogram_data->nof_pairs;
  nof_removed = 0;

  gt_assert(nof > 1);
  for (index = 0; index < nof; index++) {
    key = histogram_data->keys;
    /* left interval bound */
    if (index == 0) {
      if (key[index] < key[index + 1]-1 &&
        nof_removed < nof-2) {
        key[index] = GT_WORD_MAX;
        nof_removed++;
      }
    }
    /* right interval bound */
    else if (index == nof-1) {
      if (key[index] > key[index - 1]+1 &&
          nof_removed < nof-2)
        key[index] = GT_WORD_MAX;
    }
    /* within interval */
    else {
      if (key[index] > key[index - 1]+1 &&
          key[index] < key[index + 1]-1 &&
          nof_removed < nof-2) {
        key[index] = GT_WORD_MAX;
        nof_removed++;
      }
    }
  }
}

/* create probability mass function (pmf) based on histogram */
static void create_pmf(PmfData *pmf_data, HistogramData histogram_data) {
  GtUword key, *nof;

  pmf_data->mean = histogram_data.mean;
  pmf_data->std_dev = histogram_data.std_dev;
  pmf_data->minp = 1.0 / histogram_data.value_sum;
  pmf_data->dist = gt_malloc(sizeof (*pmf_data->dist) * (histogram_data.max + 1));
  pmf_data->nof = histogram_data.max + 1;
  for (key = 0; key <= histogram_data.max; key++) {
    nof = gt_hashmap_get(histogram_data.hash_map, &key);
    if (nof != NULL && *nof > 0)
      pmf_data->dist[key] = (double)*nof / histogram_data.value_sum;
    else
      pmf_data->dist[key] = pmf_data->minp;
  }
}

static double window(GtWord x1,
                     GtWord x2,
                     GtWord x) {
  GtWord return_val, x3;

  x3 = x1 + x2;
  if (x <= 0)
    return_val = 1;
  else if (x < x1)
    return_val = x;
  else if (x < x2)
    return_val = x1;
  else if (x < x3)
    return_val = x3 - x;
  else
    return_val = 1;

  return (double)return_val / x1;
}

/* sort fragment sizes */
static int compare_fragments(const void *a,
                             const void *b) {
  GtUword *fragment_a = (GtUword*) a;
  GtUword *fragment_b = (GtUword*) b;
  int return_val;

  if (*(fragment_a+1) - *fragment_a < *(fragment_b+1) - *fragment_b)
    return_val = -1;
  else if (*(fragment_a+1) - *fragment_a > *(fragment_b+1) - *fragment_b)
    return_val = 1;
  else
    return_val = 0;

  return return_val;
}

/* calculate fragment size distribution */
static void calculate_fragment_dist(FragmentData *fragment_data,
                                    GtWord factor) {
  GtUword size, index, next_free;

  qsort(fragment_data->frag_pos, fragment_data->nof_frag_pos,
        sizeof (GtUword)*2, compare_fragments);

  fragment_data->nof_frag_size = 0;

  for (index = 0; index < fragment_data->nof_frag_pos; index++) {
    size = fragment_data->frag_pos[(index*2)+1] -
           fragment_data->frag_pos[index*2];

    /* resize fragment array */
    if (fragment_data->nof_frag_size == fragment_data->size_frag_size) {
      fragment_data->size_frag_size += INCREMENT_SIZE;
      fragment_data->frag_size = gt_realloc(fragment_data->frag_size,
                                 sizeof (*(fragment_data->frag_size)) *
                                 fragment_data->size_frag_size * 2);
    }

    if (index == 0 ||
       (index > 0 && (size != (fragment_data->frag_pos[((index-1)*2)+1] -
                               fragment_data->frag_pos[(index-1)*2])))) {
      next_free = fragment_data->nof_frag_size * 2;
      fragment_data->frag_size[next_free] = size - factor;
      fragment_data->frag_size[next_free+1] = 1;
      fragment_data->nof_frag_size++;
    }
    else
      fragment_data->frag_size[next_free+1]++;
  }
}

/* compute the log likelihood that these samples came from the
   specified distribution shifted by the parameter theta */
static int compute_likelihood(double *likelihood,
                              GtUword *nof_pairs,
                              GtWord theta,
                              FragmentData *fragment_data,
                              PmfData pmf_data,
                              GtError *err)
{
  int had_err = 0;
  GtUword index;
  GtWord new_frag_size;
  double pmf_prob;

  /* iterate over fragment sizes */
  for (index = 0; index < fragment_data->nof_frag_size; index++) {
    new_frag_size = fragment_data->frag_size[index*2] + theta;

    /* OBSTACLE: behavior for negative index ambiguously described */
    if (new_frag_size >= 0 && new_frag_size < pmf_data.nof)
      pmf_prob = pmf_data.dist[new_frag_size];
    else
      pmf_prob = pmf_data.minp;

    *likelihood += fragment_data->frag_size[(index*2)+1] * log(pmf_prob);

    if (pmf_prob > pmf_data.minp)
      *nof_pairs += fragment_data->frag_size[(index*2)+1];

    if (pmf_prob < 0) {
      gt_error_set (err , "negative probability");
      had_err = -1;
      break;
    }
  }

  return had_err;
}

static int maximum_likelihood_estimate(GtWord *dist,
                                       GtUword *nof_pairs,
                                       GtWord min_dist,
                                       GtWord max_dist,
                                       FragmentData *fragment_data,
                                       PmfData pmf_data,
                                       GtUword len_ref,
                                       GtUword len_mref,
                                       GtError *err) {
  int had_err = 0;
  GtUword best_n = 0, n, min_frag_size, max_frag_size, index;
  GtWord best_theta = min_dist, theta;
  double pmf_prob, best_likelihood = (double)GT_WORD_MIN, c, likelihood;

  min_frag_size = fragment_data->frag_size[0];
  max_frag_size = fragment_data->frag_size[(fragment_data->nof_frag_size-1)*2];

  min_dist = MAX(min_dist, 0 - min_frag_size);
  max_dist = MIN(max_dist, pmf_data.nof - max_frag_size - 1);

  for (theta = min_dist; theta <= max_dist; theta++) {
  /* Calculate the normalizing constant of the PMF, f_theta(x) */
    c = 0;

    for (index = 0; index <  pmf_data.nof; index++) {

      if (index < pmf_data.nof)
        pmf_prob = pmf_data.dist[index];
      else
        pmf_prob = pmf_data.minp;

      c += pmf_prob * window(len_ref, len_mref, index - theta);
    }

    likelihood = 0;
    n = 0;
    had_err = compute_likelihood(&likelihood, &n, theta, fragment_data,
              pmf_data, err);
    if (had_err != 0)
      break;

    likelihood -= fragment_data->nof_frag_pos * log(c);

    if (n > 0 && likelihood > best_likelihood) {
      best_likelihood = likelihood;
      best_theta = theta;
      best_n = n;
    }
  }
  *nof_pairs = best_n;
  *dist = best_theta;

  return had_err;
}

/* estimate the distance between two contigs
   using maximum likelihood estimator */
static int estimate_dist_using_mle(GtWord *dist,
                                   GtUword *nof_pairs,
                                   GtWord min_dist,
                                   GtWord max_dist,
                                   FragmentData *fragment_data,
                                   PmfData pmf_data,
                                   GtUword len_ref,
                                   GtUword len_mref,
                                   bool rf,
                                   GtError *err) {

  int had_err = 0;
  GtUword temp;

  len_ref -= fragment_data->ma - 1;
  len_mref -= fragment_data->ma - 1;

  /* swap reference and mate reference length */
  if (len_ref > len_mref) {
    temp = len_ref;
    len_ref = len_mref;
    len_mref = temp;
  }

  /* library is oriented reverse-forward */
  if (rf) {
    calculate_fragment_dist(fragment_data, 0);
    had_err = maximum_likelihood_estimate(dist, nof_pairs, min_dist, max_dist,
                            fragment_data, pmf_data, len_ref, len_mref, err);

  }
  /* library is oriented forward-reverse */
  /* Subtract 2*(l-1) from each sample */
  else {
    calculate_fragment_dist(fragment_data, 2 * (fragment_data->ma - 1));
    had_err = maximum_likelihood_estimate(dist, nof_pairs, min_dist, max_dist,
                             fragment_data, pmf_data, len_ref, len_mref, err);
    *dist = MAX(min_dist, *dist - 2 * (GtWord)(fragment_data->ma - 1));
  }

  return had_err;
}

/* callback function
int compare_pmf_histogram(void *key,
                          void *value,
                          void *data,
                          GtError *err) {
  GtWord *key_2 = (GtUword*) key;
  GtUword *value_2 = (GtUword*) value;
  CompareData *compare_data = (CompareData*) data;
  double pmf_prob;

  if ((*key_2 + compare_data->dist) < compare_data->pmf_data.nof)
    pmf_prob = compare_data->pmf_data.dist[*key_2 + compare_data->dist];
  else
    pmf_prob = compare_data->pmf_data.minp;

  if (pmf_prob > compare_data->pmf_data.minp)
    compare_data->nof_pairs += *value_2;
  return 0;
}
*/

/* Estimate the distance between two contigs using the difference of
   the population mean and the sample mean
int estimate_dist_using_mean(FragmentData fragment_data,
                             PmfData pmf_data,
                             GtUword *nof_pairs,
                             GtWord *dist,
                             GtError *err) {
  int had_err;
  HistogramData histogram_data;
  CompareData compare_data;

  had_err = 0;

  histogram_data = create_histogram_from_dist(fragment_data);

  *dist = roundl(pmf_data.mean - histogram_data.mean);

 Count the number of fragments that agree with the distribution
  compare_data.nof_pairs = 0;
  compare_data.pmf_data = pmf_data;
  compare_data.dist = *dist;
  had_err = gt_hashmap_foreach(histogram_data.hash_map,
                               compare_pmf_histogram, &compare_data, err);
  if (had_err == 0)
    *nof_pairs = compare_data.nof_pairs;

  return had_err;
}*/

/* calculate provisional fragment size as if the contigs
   were perfectly adjacent with no overlap or gap and
   save only fragment size if combination of start and
   end point are unique */
static void calculate_fragment(GtUword start_read,
                               bool read_reverse,
                               GtUword start_mread,
                               bool mread_reverse,
                               GtUword len_ref,
                               GtUword len_mref,
                               FragmentData *fragment_data,
                               bool rf) {

  GtUword start_frag, end_frag, align;
  GtWord size, *pair;
  bool found;

  /* correct start position of read, mate read
     if necessary */
  if (read_reverse)
    start_read = len_ref - start_read;
  if (!mread_reverse)
    start_mread = len_mref - start_mread;

  /* determine start, end position of fragment */
  if (rf) {
    start_frag = start_mread;
    end_frag = len_mref + start_read;
  }
  else {
    start_frag = start_read;
    end_frag = len_ref + start_mread;
  }

  /* resize fragment array */
  if (fragment_data->nof_frag_pos == fragment_data->size_frag_pos) {
    fragment_data->size_frag_pos += INCREMENT_SIZE;
    fragment_data->frag_pos = gt_realloc(fragment_data->frag_pos,
           sizeof (*(fragment_data->frag_pos)) *
                  fragment_data->size_frag_pos * 2);
  }

  /* check if fragment already exists */
  found = false;
  for (pair = fragment_data->frag_pos;
       pair < fragment_data->frag_pos + (fragment_data->nof_frag_pos * 2);
       pair += 2) {
    if (*pair == start_frag && *(pair+1) == end_frag) {
      found = true;
      break;
    }
  }

  /* save fragment size with unique start and end point */
  if (!found) {
    size = end_frag - start_frag;
    fragment_data->frag_pos[fragment_data->nof_frag_pos * 2] = start_frag;
    fragment_data->frag_pos[(fragment_data->nof_frag_pos * 2) + 1] = end_frag;
    fragment_data->nof_frag_pos++;

    if (!rf && size <= 2 * (GtWord)(fragment_data->ma - 1)) {
      align = size / 2;
      fragment_data->ma = MIN(fragment_data->ma, align);
    }
  }
}

/* !functions missing in genometools! */

/* determine length of reference */
GtUword gt_samfile_iterator_reference_length(const GtSamfileIterator *s_iter,
                                             int32_t reference_num) {
  gt_assert(reference_num >= 0);
  gt_assert(reference_num < s_iter->samfile->header->n_targets);
  return s_iter->samfile->header->target_len[reference_num];
}

/* determine id of mate reference */
int32_t gt_sam_alignment_mate_ref_num(GtSamAlignment *bam_align)
{
  gt_assert(bam_align != NULL);
  return bam_align->s_alignment->core.mtid;
}

/* determine insert size of read pair */
GtUword gt_sam_alignment_insert_size(GtSamAlignment *bam_align) {
  gt_assert(bam_align != NULL);
  return (GtUword) bam_align->s_alignment->core.isize;
}

GtUword gt_sam_alignment_query_length(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  return (GtUword) sam_alignment->s_alignment->core.l_qname;
}

/* calculate reference and mate reference position at read start */
static int calc_read_start(GtSamAlignment *bam_align,
                           GtWord *ref_read_start,
                           GtError *err) {
  int had_err = 0;
  GtUchar op;
  bool first = true;
  GtUword index, len;
  GtUword qstart = 0;
  GtUword qlen = 0;
  GtUword qspan = 0;
  GtUword tspan = 0;

  for (index = 0; index < gt_sam_alignment_cigar_length(bam_align); index++) {
    op = gt_sam_alignment_cigar_i_operation(bam_align, index);
    len = gt_sam_alignment_cigar_i_length(bam_align, index);

    /* soft clipping, hard clipping */
    if (op == 'S' || op == 'H') {
      if (first)
        qstart = len;
      qlen += len;
    }
    /* alignment match, sequence mismatch, sequence match */
    else if (op == 'M' || op == 'X' || op == '=') {
      qlen += len;
      qspan += len;
      tspan += len;
    }
    /* insertion to the reference */
    else if (op == 'I') {
      qlen += len;
      qspan += len;
    }
    /* deletion from the reference, skipped region from the reference,
       padding (silent deletion from padded reference) */
    else if (op == 'D' || op == 'N' || op == 'P')
      tspan += len;
    /* invalid cigar operation occurs */
    else {
      had_err = -1;
      gt_error_set(err, "invalid cigar operation");
      break;
    }
    first = false;
  }

  gt_assert(qstart + qspan <= qlen);
  if (gt_sam_alignment_is_reverse(bam_align))
    *ref_read_start = gt_sam_alignment_pos(bam_align) +
                      tspan + (qlen - qspan - qstart);
  else
    *ref_read_start = gt_sam_alignment_pos(bam_align) - qstart;

  return had_err;
}

/* analyze read set */
static int analyze_read_set(DistRecords *dist_records,
                            ReadSet read_set,
                            GtUword min_nof_pairs,
                            GtWord min_dist,
                            GtWord max_dist,
                            FragmentData *fragment_data,
                            PmfData pmf_data,
                            bool rf,
                            GtError *err) {
  GtUword nof_pairs;
  GtWord dist;
  double std_dev;
  bool sense, same, first_contig, new_contig;
  GtStr *ctg_id;
  Read *read;
  int had_err = 0;

  first_contig = true;
  new_contig = false;

  /* iterate over reads */
  for (read = read_set.read; read < read_set.read+read_set.nof; read++) {

    /* analyze read subset */
    if (read != read_set.read && fragment_data->nof_frag_pos != 0 &&
        (read->tid != (read-1)->tid ||
         read->mtid != (read-1)->mtid ||
         read->is_reverse != (read-1)->is_reverse)) {

      /* start new record if switch of contig (reference) occurred or
         first contig (reference) is present */
      if (first_contig || new_contig) {
        create_dist_record(dist_records, gt_str_clone((read-1)->ref_name));
        first_contig = false;
        new_contig = false;
      }

      nof_pairs = 0;
      /* calculate distance and nof pairs */
      if (fragment_data->nof_frag_pos >= min_nof_pairs)
        had_err = estimate_dist_using_mle(&dist, &nof_pairs, min_dist,
                                         max_dist, fragment_data, pmf_data,
                                      read->ref_len, read->mref_len, rf, err);

      /*  append new entry to current record */
      if (nof_pairs >= min_nof_pairs) {
        ctg_id = gt_str_clone((read-1)->mref_name);
        std_dev = pmf_data.std_dev / sqrt(nof_pairs);
        sense = (read-1)->is_reverse != (read-1)->is_mreverse;
        same = (read-1)->is_reverse;
        add_contig_dist_record(dist_records, ctg_id, std_dev, dist, nof_pairs,
                               sense, same);
      }

      /* reset fragment array */
      fragment_data->nof_frag_pos = 0;
    }

    /* check if switch of contig (reference) occurred */
    if (read != read_set.read && read->tid != (read-1)->tid)
      new_contig = true;

    /* check if read is valid */
    if (!read->is_munmapped && read->is_valid) {

      /* calculate fragment size and save it if start and
         end point are unique */
      calculate_fragment(read->ref_read_start, read->is_reverse,
                         read->ref_read_start+read->isize, read->is_mreverse,
                         read->ref_len, read->mref_len, fragment_data, rf);

    }
  }

  nof_pairs = 0;
  /* calculate distance and nof pairs */
  if (fragment_data->nof_frag_pos >= min_nof_pairs)
    had_err = estimate_dist_using_mle(&dist, &nof_pairs, min_dist,
                                      max_dist, fragment_data, pmf_data,
                                      read->ref_len, read->mref_len, rf, err);

  /*  append new entry to current record */
  if (nof_pairs >= min_nof_pairs) {
    ctg_id = gt_str_clone((read-1)->mref_name);
    std_dev = pmf_data.std_dev / sqrt(nof_pairs);
    sense = (read-1)->is_reverse != (read-1)->is_mreverse;
    same = (read-1)->is_reverse;
    add_contig_dist_record(dist_records, ctg_id, std_dev, dist, nof_pairs,
                           sense, same);
  }

  return had_err;
}

/* load reads from bam-file */
static int load_read_set(ReadSet *readset,
                         const char *bam_filename,
                         GtUword min_ref_length,
                         GtUword min_qual,
                         GtError *err) {

  int had_err;
  GtAlphabet *alpha;
  GtSamfileIterator *bam_iterator;
  GtSamAlignment *bam_align;
  GtWord ref_read_start;
  GtUword next_free;
  GtStr *name;

  had_err = 0;
  alpha = gt_alphabet_new_dna();
  name = gt_str_new();
  bam_iterator = gt_samfile_iterator_new_bam(bam_filename, alpha, err);

  if (bam_iterator == NULL)
    had_err = -1;

  if (had_err == 0) {
    /* iterate over reads and save them */
    while (gt_samfile_iterator_next(bam_iterator, &bam_align) > 0) {

      had_err = calc_read_start(bam_align, &ref_read_start, err);
      if (had_err == 0) {

        /* resize read set */
        if (readset->nof == readset->size) {
          readset->size += INCREMENT_SIZE;
          readset->read = gt_realloc(readset->read,
                                sizeof (*readset->read) * readset->size);
        }

        next_free = readset->nof;

        /* save query name of read */
        gt_str_set(name, gt_sam_alignment_identifier(bam_align));
        readset->read[next_free].query_name = gt_str_clone(name);
        readset->read[next_free].is_unmapped =
                              gt_sam_alignment_is_unmapped(bam_align);
        readset->read[next_free].is_reverse =
                              gt_sam_alignment_is_reverse(bam_align);
        readset->read[next_free].ref_read_start = ref_read_start;
        readset->read[next_free].pos = gt_sam_alignment_pos(bam_align);

        /* save id of reference and mate reference */
        readset->read[next_free].tid = gt_sam_alignment_ref_num(bam_align);
        readset->read[next_free].mtid =
                                 gt_sam_alignment_mate_ref_num(bam_align);

        /* save name und length of reference */
        if (readset->read[next_free].tid >= 0) {
          gt_str_set(name, gt_samfile_iterator_reference_name(
                         bam_iterator, readset->read[next_free].tid));
          readset->read[next_free].ref_name = gt_str_clone(name);
          readset->read[next_free].ref_len =
                         gt_samfile_iterator_reference_length(
                         bam_iterator, readset->read[next_free].tid);

        }
        else {
          readset->read[next_free].ref_name = NULL;
          readset->read[next_free].ref_len = 0;
        }

        /* save name und length of mate reference */
        if (readset->read[next_free].mtid >= 0) {
          gt_str_set(name, gt_samfile_iterator_reference_name(
                         bam_iterator, readset->read[next_free].mtid));
          readset->read[next_free].mref_name = gt_str_clone(name);
          readset->read[next_free].mref_len =
                 gt_samfile_iterator_reference_length(
                  bam_iterator, readset->read[next_free].mtid);

        }
        else {
          readset->read[next_free].mref_name = NULL;
          readset->read[next_free].mref_len = 0;
        }

        /* check if read is valid */
        if (readset->read[next_free].ref_len < min_ref_length ||
            gt_sam_alignment_mapping_quality(bam_align) < min_qual ||
            !gt_sam_alignment_is_paired(bam_align) ||
            readset->read[next_free].tid == readset->read[next_free].mtid ||
            readset->read[next_free].is_unmapped)
          readset->read[next_free].is_valid = false;
        else
          readset->read[next_free].is_valid = true;

        readset->nof++;
      }
    }
    gt_samfile_iterator_delete(bam_iterator);
  }

  gt_alphabet_delete(alpha);
  gt_str_delete(name);

  return had_err;
}

static int compare_fragment(const void *a,
                            const void *b)
{
  GtUword *fragement_a = (GtUword*) a;
  GtUword *fragement_b = (GtUword*) b;

  if (*fragement_a > *fragement_b)
    return 1;
  else if (*fragement_a < *fragement_b)
    return -1;
  else
    return 0;
}

/* correct strand orientation and mapping status of mate read,
   set insert size of read pair */
static void fixmate(Read *read_0,
                    Read *read_1) {

  read_0->is_munmapped = false;
  read_0->is_mreverse = false;

  if (read_1->is_unmapped)
    read_0->is_munmapped = true;
  if (read_1->is_reverse)
    read_0->is_mreverse = true;

  if (read_0->is_munmapped)
    read_0->isize = 0;
  else
    read_0->isize = read_1->ref_read_start - read_0->ref_read_start;
}

/* create histogram based on fragment sizes of reads */
static void create_histogram(HistogramData *histogram_data,
                             ReadSet readset) {
  GtUword *fragment_set, fragment_set_nof, fragment_set_size,
           next_free, *fragment, insert_size, *value,
           nof_rf, nof_fr, value_sum;
  GtWord *key;
  Read *read;

  fragment_set = NULL;
  fragment_set_nof = 0;
  fragment_set_size = 0;

  /* iterate over reads pairs and save their fragment sizes */
  for (read = readset.read; read < readset.read+readset.nof; read++) {

    /* check if sequential reads are a read pair */
    if (read+1 < readset.read+readset.nof &&
        gt_str_cmp(read->query_name, (read+1)->query_name) == 0) {

      fixmate(read, read+1);
      fixmate(read+1, read);

      /* check if read pair is valid */
      if (!read->is_unmapped && !(read+1)->is_unmapped &&
           read->tid == (read+1)->tid &&
           read->is_reverse != (read+1)->is_reverse) {

        if (read->is_reverse)
          insert_size = (read+1)->isize;
        else
          insert_size = read->isize;

        /* resize fragment array */
        if (fragment_set_nof == fragment_set_size) {
          fragment_set_size += INCREMENT_SIZE;
          fragment_set = gt_realloc(fragment_set,
                                 sizeof (*fragment_set) * fragment_set_size);
        }

        /* save insert size of read pair */
        fragment_set[fragment_set_nof] = insert_size;
        fragment_set_nof++;
      }
    }
  }

  /* sort fragments by their sizes */
  qsort(fragment_set, fragment_set_nof, sizeof (*fragment_set),
        compare_fragment);

  next_free = 0;

  /* iterate over sorted fragments and count their frequencies */
  for (fragment = fragment_set; fragment < fragment_set+fragment_set_nof;
       fragment++) {

    /* resize histogram hashmap */
    if (histogram_data->nof_pairs == histogram_data->size) {
      histogram_data->size += INCREMENT_SIZE;
      histogram_data->values = gt_realloc(histogram_data->values,
                    sizeof (*histogram_data->values) * histogram_data->size);
      histogram_data->keys = gt_realloc(histogram_data->keys,
                    sizeof (*histogram_data->keys) * histogram_data->size);
    }

    /* check if first fragment is present or fragment changed */
    if (fragment == fragment_set || *fragment != *(fragment-1)) {

      /* ignore fragments with frequency lower than 2 */
      if (fragment != fragment_set && *fragment != *(fragment-1) &&
          histogram_data->values[next_free] < 3)
        histogram_data->nof_pairs--;

      next_free = histogram_data->nof_pairs;
      histogram_data->keys[next_free] = *fragment;
      histogram_data->values[next_free] = 1;
      histogram_data->nof_pairs++;
    }
    else
      histogram_data->values[next_free]++;
  }

  /* ignore fragments with frequency lower than 2 */
  if (histogram_data->values[next_free] < 3)
    histogram_data->nof_pairs--;

  nof_rf = 0;
  nof_fr = 0;
  value_sum = 0;
  /* save fragment frequencies in hashmap */
  for (key = histogram_data->keys, value = histogram_data->values;
       key < histogram_data->keys+histogram_data->nof_pairs; key++, value++) {

    /* count keys in respective range */
    if (*key >= INT_MIN && *key <= 0)
      nof_rf++;
    if (*key >= 1 && *key <= INT_MAX)
      nof_fr++;
    value_sum += *value;

    gt_hashmap_add(histogram_data->hash_map, key, value);
    /* printf("" GT_WD "\t" GT_WU " \n", *key, *value); */
  }

  histogram_data->lib_rf = nof_fr < nof_rf;
  histogram_data->value_sum = value_sum;

  /* remove noise */
  remove_noise_in_histogram(histogram_data);
  /* calculate statistics */
  calc_stat_of_histogram(histogram_data);

  /* remove outlier */
  remove_outliers_in_histogram(histogram_data);
  /* calculate statistics */
  calc_stat_of_histogram(histogram_data);

  /* trim off fraction */
  trim_fraction_of_histogram(histogram_data, 0.0001);
  /* calculate statistics */
  calc_stat_of_histogram(histogram_data);

  gt_free(fragment_set);
}

/* callback function to sort reads by query name */
static int compare_read_2(const void *a,
                          const void *b)
{
  Read *read_a = (Read*) a;
  Read *read_b = (Read*) b;

  return gt_str_cmp(read_a->query_name, read_b->query_name);
}

/* callback function to sort reads by reference id and orientation
   of read relative to reference, reference id of mate read,
   orientation of read relative to his mate read, alignment pos */
static int compare_read_3(const void *a,
                          const void *b)
{
  Read *read_a = (Read*) a;
  Read *read_b = (Read*) b;
  int return_val;

  if (read_a->tid < read_b->tid)
    return_val = -1;
  else if (read_a->tid > read_b->tid)
    return_val = 1;
  else {
    if (read_a->is_reverse && !read_b->is_reverse)
      return_val = 1;
    else if (!read_a->is_reverse && read_b->is_reverse)
      return_val = -1;
    else if (read_a->mtid > read_b->mtid)
      return_val = 1;
    else if (read_a->mtid < read_b->mtid)
      return_val = -1;
    else if (read_a->is_reverse != read_a->is_mreverse)
      return_val = -1;
    else if (read_a->is_reverse == read_a->is_mreverse)
      return_val = 1;
    else if (read_a->pos < read_b->pos)
      return_val = -1;
    else
      return_val = 1;
  }

  return return_val;
}

/* read paired information from bam file and corresponding hist file */
int gt_scaffolder_bamparser_read_paired_information(DistRecords *dist,
                                                    const char *bam_filename,
                                                    GtWord min_dist,
                                                    GtWord max_dist,
                                                    GtUword min_qual,
                                                    GtUword min_nof_pairs,
                                                    GtUword min_ref_length,
                                                    GtUword min_align,
                                                    GtError *err) {
  int had_err = 0;
  bool rf;
  FragmentData fragment_data;
  HistogramData histogram_data;
  PmfData pmf_data;
  ReadSet read_set;
  Read *read;

  /* initialize read set */
  read_set.size = 0;
  read_set.nof = 0;
  read_set.read = NULL;

  /* initialize fragment array */
  fragment_data.nof_frag_pos = 0;
  fragment_data.size_frag_pos = 0;
  fragment_data.frag_pos = NULL;
  fragment_data.nof_frag_size = 0;
  fragment_data.size_frag_size = 0;
  fragment_data.frag_size = NULL;

  /* initialize hashmap */
  histogram_data.hash_map = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  histogram_data.nof_pairs = 0;
  histogram_data.size = 0;
  histogram_data.values = NULL;
  histogram_data.keys = NULL;

  /* initialize pmf */
  pmf_data.dist = NULL;
  pmf_data.nof = 0;

  /* load reads from bam-file */
  had_err = load_read_set(&read_set, bam_filename, min_ref_length,
                          min_qual, err);

  if (had_err == 0) {
    /* sort reads by query name */
    qsort(read_set.read, read_set.nof, sizeof (*read_set.read),
    compare_read_2);

    /* create histogram based on fragment sizes of reads */
    create_histogram(&histogram_data, read_set);

    /* create probability mass function (pmf) based on histogram */
    create_pmf(&pmf_data, histogram_data);

    /* default cutoff values */
    max_dist = pmf_data.nof-1;
    fragment_data.ma = min_align;

    /* determine the orientation of the library */
    rf = histogram_data.lib_rf;
    /* TODO
       if (rf)
       negate each element of histogram
    */

    /* sort reads by reference id and orientation of read relative to
       reference, reference id of mate read, orientation of read relative
       to his mate read, alignment pos */
    qsort(read_set.read, read_set.nof, sizeof (*read_set.read),
                                       compare_read_3);

    analyze_read_set(dist, read_set, min_nof_pairs, min_dist, max_dist,
                   &fragment_data, pmf_data, rf, err);

  }
  /* clean up */
  for (read = read_set.read; read < read_set.read+read_set.nof; read++) {
    gt_free(read->query_name);
    gt_free(read->ref_name);
    gt_free(read->mref_name);
  }
  gt_free(read_set.read);

  gt_free(histogram_data.keys);
  gt_free(histogram_data.values);
  gt_hashmap_delete(histogram_data.hash_map);
  gt_free(fragment_data.frag_pos);
  gt_free(fragment_data.frag_size);
  gt_free(pmf_data.dist);

  return had_err;
}
