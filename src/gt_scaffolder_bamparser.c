#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

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

/* BAM reads */
typedef struct Read {
  bool reverse;
  bool mreverse;
  GtUword tid;
  GtUword mtid;
  GtWord ref_read_start;
  GtWord mref_read_start;
} Read;

typedef struct ReadSet {
  Read *read;
  GtUword nof;
  GtUword size;
} ReadSet;

typedef struct FragmentData {
  GtWord *frag_pos;
  GtWord *frag_size;
  GtUword nof;
  GtUword size;
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

void init_dist_records(DistRecords *dist) {
  GtUword index;

  dist->nof_record = 0;
  dist->size = 50;
  dist->record = gt_malloc(sizeof (*(dist->record)) * dist->size);
  for (index = 0; index < dist->size; index++) {
    dist->record[index].size = 20;
    dist->record[index].nof_ctg = 0;
    dist->record[index].root_ctg_id = NULL;
    dist->record[index].ctg = gt_malloc(sizeof (*(dist->record[index].ctg))*
                            dist->record[index].size);
  }
}

void write_dist_records(DistRecords dist) {
  GtUword index, index_2;
  bool set_dir;

  for (index = 0; index < dist.nof_record; index++) {
    printf("%s",dist.record[index].root_ctg_id);
    set_dir = false;
    for (index_2 = 0; index_2 < dist.record[index].nof_ctg; index_2++) {
      if (dist.record[index].ctg[index_2].same && !set_dir) {
        printf(" ;");
        set_dir = true;
      }

      printf(" %s%c," GT_WD "," GT_WU ",%.1f",
                       dist.record[index].ctg[index_2].id,
                       dist.record[index].ctg[index_2].sense ? '+' : '-',
                       dist.record[index].ctg[index_2].dist,
                       dist.record[index].ctg[index_2].nof_pairs,
                       dist.record[index].ctg[index_2].std_dev);
     }
     if (!set_dir)
       printf(" ;");
     printf("\n");
  }
}

void delete_dist_records(DistRecords dist) {
 GtUword index, index_2;

  for (index = 0; index < dist.size; index++) {
    for (index_2 = 0; index_2 < dist.record[index].nof_ctg; index_2++)
      gt_free(dist.record[index].ctg[index_2].id);
      gt_free(dist.record[index].ctg);
      gt_free(dist.record[index].root_ctg_id);
    }
  gt_free(dist.record);
}

/* sort read set by orientation of read relative to reference,
   reference id of mate read and orientation of read relative
   to his mate read */
static int compare_read_order(const void *a,
                              const void *b)
{
  int return_val;
  Read *read_a = (Read*) a;
  Read *read_b = (Read*) b;

  if (read_a->reverse && !read_b->reverse)
    return_val = 1;
  else if (!read_a->reverse && read_b->reverse)
    return_val = -1;
  else if (read_a->mtid > read_b->mtid)
    return_val = 1;
  else if (read_a->mtid < read_b->mtid)
    return_val = -1;
  else if (read_a->reverse != read_a->mreverse)
    return_val = -1;
  else if (read_a->reverse == read_a->mreverse)
    return_val = 1;

  return return_val;
}

/* initialize histogram */
static void init_histogram(HistogramData *histogram_data) {
  histogram_data->hash_map = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  histogram_data->nof_pairs = 0;
  histogram_data->size = 1000;
  histogram_data->value_sum = 0;
  histogram_data->keys = gt_malloc(sizeof (*histogram_data->keys)
                        * histogram_data->size);
  histogram_data->values = gt_malloc(sizeof (*histogram_data->values)
                        * histogram_data->size);
}

/* create histogram based on fragment distribution */
static void create_histogram_from_dist(FragmentData fragment_data,
                                       GtWord factor,
                                       HistogramData *histogram_data) {
  GtWord sum, *size, new_size,min, max;
  GtUword *value;

  sum = 0;
  min = GT_WORD_MAX;
  max = GT_WORD_MIN;

  /* iterate over fragment sizes */
  for (size = fragment_data.frag_size;
       size < fragment_data.frag_size + fragment_data.nof; size++) {

    /* resize histogram hashmap  */
    if (histogram_data->nof_pairs == histogram_data->size) {
      histogram_data->size += 100;
      histogram_data->values = gt_realloc(histogram_data->values,
                     sizeof (*histogram_data->values) * histogram_data->size);
      histogram_data->keys = gt_realloc(histogram_data->keys,
                     sizeof (*histogram_data->keys) * histogram_data->size);
    }

    /* fill histogram hash map */
    new_size = *size - factor;
    value = gt_hashmap_get(histogram_data->hash_map, &new_size);
    if (value == NULL) {
      histogram_data->values[histogram_data->nof_pairs] = 1;
      histogram_data->keys[histogram_data->nof_pairs] = new_size;
      gt_hashmap_add(histogram_data->hash_map,
                     histogram_data->keys + histogram_data->nof_pairs,
                     histogram_data->values + histogram_data->nof_pairs);
      histogram_data->nof_pairs++;
    }
    else
      (*value)++;

    sum += *size;

    if (min > new_size)
      min = new_size;
    if (max < new_size)
      max = new_size;

  }

  /* save mean */
  histogram_data->mean = (double)sum / fragment_data.nof;
  histogram_data->value_sum = fragment_data.nof;
  histogram_data->min = min;
  histogram_data->max = max;
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

/* create histogram based on histogram file */
static int create_histogram_from_file(const char *file_name,
                                      HistogramData *histogram_data,
                                      GtError *err) {
  FILE *file;
  int had_err;
  char line[1024];
  GtWord key, value;
  GtUword value_sum, nof_rf, nof_fr;

  had_err = 0;
  file = fopen(file_name, "rb");
  if (file == NULL) {
    had_err = -1;
    gt_error_set(err, "can not read file %s",file_name);
  }

  nof_rf = 0;
  nof_fr = 0;
  value_sum = 0;

  if (had_err != -1) {
    /* iterate over lines of hist-file */
    while (fgets(line, 1024, file) != NULL) {
      /* resize histogram hashmap */
      if (histogram_data->nof_pairs == histogram_data->size) {
        histogram_data->size += 100;
        histogram_data->values = gt_realloc(histogram_data->values,
                      sizeof (*histogram_data->values) * histogram_data->size);
        histogram_data->keys = gt_realloc(histogram_data->keys,
                      sizeof (*histogram_data->keys) * histogram_data->size);
      }
      /* parse line of hist-file and save key/value pair in hashmap */
      if (sscanf(line,"" GT_WD "\t" GT_WD "", &key, &value) == 2 &&
          value >= 0) {
        if (key >= 0) {
          histogram_data->values[histogram_data->nof_pairs] = (GtUword) value;
          histogram_data->keys[histogram_data->nof_pairs] = key;
          gt_hashmap_add(histogram_data->hash_map,
                           histogram_data->keys + histogram_data->nof_pairs,
                           histogram_data->values + histogram_data->nof_pairs);
          histogram_data->nof_pairs++;

          value_sum += value;
        }

        /* count keys in range*/
        if (key >= INT_MIN && key <= 0)
          nof_rf++;
        if (key >= 1 && key <= INT_MAX)
          nof_fr++;

      }
      else {
        had_err = -1;
        gt_error_set(err, "negative count value in histogram file %s",
                     file_name);
        break;
      }
    }
  }
  fclose(file);

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

  return had_err;
}

/* create probability mass function (pmf) based on histogram */
static PmfData create_pmf(HistogramData histogram_data) {
  PmfData pmf_data;
  GtUword key, *nof;

  pmf_data.mean = histogram_data.mean;
  pmf_data.std_dev = histogram_data.std_dev;
  pmf_data.minp = 1.0 / histogram_data.value_sum;
  pmf_data.dist = gt_malloc(sizeof (*pmf_data.dist) * (histogram_data.max + 1));
  pmf_data.nof = histogram_data.max + 1;
  for (key = 0; key <= histogram_data.max; key++) {
    nof = gt_hashmap_get(histogram_data.hash_map, &key);
    if (nof != NULL && *nof > 0)
      pmf_data.dist[key] = (double)*nof / histogram_data.value_sum;
    else
      pmf_data.dist[key] = pmf_data.minp;
  }

  return pmf_data;
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

/* callback function */
static int compare_pmf_histogram(void *key,
                                 void *value,
                                 void *data,
                                 GtError *err) {
  GtWord *key_2 = (GtWord*) key;
  GtUword *value_2 = (GtUword*) value;
  CompareData *compare_data = (CompareData*) data;
  double pmf_prob;
  int had_err = 0;
  GtWord index = *key_2 + compare_data->dist;

  /* OBSTACLE: behavior for negative index ambiguously described */
  if (index >= 0 && index < compare_data->pmf_data.nof)
    pmf_prob = compare_data->pmf_data.dist[index];
  else
    pmf_prob = compare_data->pmf_data.minp;

  compare_data->likelihood += *value_2 * log(pmf_prob);

  if (pmf_prob > compare_data->pmf_data.minp)
    compare_data->nof_pairs += *value_2;

  if (pmf_prob < 0) {
    gt_error_set (err , " negative probability ");
    had_err = -1;
  }

  return had_err;
}

/* callback function */
static int compare_order(const void *a,
                         const void *b) {
  GtWord *key_1 = (GtWord*) a;
  GtWord *key_2 = (GtWord*) b;
  int return_val;

  if (*key_1 > *key_2)
    return_val = 1;
  else if (*key_1 < *key_2)
    return_val = -1;
  else
    return_val = 0;

  return return_val;
}

/* compute the log likelihood that these samples came from the
   specified distribution shifted by the parameter theta */
static int compute_likelihood(double *likelihood,
                              GtUword *nof_pairs,
                              GtWord theta,
                              HistogramData histogram_data,
                              PmfData pmf_data,
                              GtError *err)
{
  CompareData compare_data;
  int had_err = 0;

  compare_data.likelihood = 0;
  compare_data.nof_pairs = 0;
  compare_data.pmf_data = pmf_data;
  compare_data.dist = theta;
  had_err = gt_hashmap_foreach_ordered(histogram_data.hash_map,
            compare_pmf_histogram, &compare_data, compare_order, err);
  *likelihood = compare_data.likelihood;
  *nof_pairs = compare_data.nof_pairs;

  return had_err;
}

static int maximum_likelihood_estimate(GtWord *dist,
                                       GtUword *nof_pairs,
                                       GtWord min_dist,
                                       GtWord max_dist,
                                       HistogramData histogram_data,
                                       PmfData pmf_data,
                                       GtUword len_ref,
                                       GtUword len_mref,
                                       GtError *err) {
  int had_err = 0;
  GtUword best_n = 0, n;
  GtWord best_theta = min_dist, index, theta;
  double pmf_prob, best_likelihood = (double)GT_WORD_MIN, c, likelihood;

  min_dist = MAX(min_dist, 0 - histogram_data.max);
  max_dist = MIN(max_dist, pmf_data.nof - histogram_data.min - 1);

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

    had_err = compute_likelihood(&likelihood, &n, theta, histogram_data,
              pmf_data, err);
    if (had_err != 0)
      break;

    likelihood -= histogram_data.value_sum * log(c);
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
                                   FragmentData fragment_data,
                                   PmfData pmf_data,
                                   GtUword len_ref,
                                   GtUword len_mref,
                                   bool rf,
                                   HistogramData histogram_data,
                                   GtError *err) {

  int had_err = 0;
  GtUword temp;

  len_ref -= fragment_data.ma - 1;
  len_mref -= fragment_data.ma - 1;

  /* swap reference and mate reference length */
  if (len_ref > len_mref) {
    temp = len_ref;
    len_ref = len_mref;
    len_mref = temp;
  }

  /* library is oriented reverse-forward */
  if (rf) {
    create_histogram_from_dist(fragment_data, 0, &histogram_data);
    had_err = maximum_likelihood_estimate(dist, nof_pairs, min_dist, max_dist,
                            histogram_data, pmf_data, len_ref, len_mref, err);

  }
  /* library is oriented forward-reverse */
  /* Subtract 2*(l-1) from each sample */
  else {
    create_histogram_from_dist(fragment_data, 2 * (fragment_data.ma - 1),
                               &histogram_data);
    had_err = maximum_likelihood_estimate(dist, nof_pairs, min_dist, max_dist,
                             histogram_data, pmf_data, len_ref, len_mref, err);
    *dist = MAX(min_dist, *dist - 2 * (GtWord)(fragment_data.ma - 1));
  }

  /* reset hashmap */
  gt_hashmap_reset(histogram_data.hash_map);
  histogram_data.nof_pairs = 0;
  histogram_data.value_sum = 0;

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

/* calculate fragments and save unique fragments */
static void calculate_fragments(GtUword start_read,
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
  if (fragment_data->nof == fragment_data->size) {
    fragment_data->size += 100;
    fragment_data->frag_pos = gt_realloc(fragment_data->frag_pos,
           sizeof (*(fragment_data->frag_pos)) * fragment_data->size * 2);
    fragment_data->frag_size = gt_realloc(fragment_data->frag_size,
           sizeof (*(fragment_data->frag_size)) * fragment_data->size);
  }

  /* check if fragment already exists */
  found = false;
  for (pair = fragment_data->frag_pos;
       pair < fragment_data->frag_pos + (fragment_data->nof * 2); pair += 2) {
    if (*pair == start_frag && *(pair+1) == end_frag) {
      found = true;
      break;
    }
  }

  /* save unique fragments and their sizes */
  if (!found) {
    size = end_frag - start_frag;
    fragment_data->frag_pos[fragment_data->nof * 2] = start_frag;
    fragment_data->frag_pos[(fragment_data->nof * 2) + 1] = end_frag;
    fragment_data->frag_size[fragment_data->nof] = size;

    fragment_data->nof++;

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

/* determine input size of read */
GtUword gt_sam_alignment_input_size(GtSamAlignment *bam_align) {
  gt_assert(bam_align != NULL);
  return (GtUword) bam_align->s_alignment->core.isize;
}

/* calculate reference and mate reference position at read start */
static int calc_read_start(GtSamAlignment *bam_align,
                           GtWord *ref_read_start,
                           GtWord *mref_read_start,
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

  *mref_read_start = *ref_read_start + gt_sam_alignment_input_size(bam_align);

  return had_err;
}

/* check if read pair is valid */
static bool is_read_pair_valid(const GtSamfileIterator *bam_iterator,
                               GtSamAlignment *bam_align,
                               GtUword min_ref_length,
                               GtUword min_qual) {
  bool result;
  if (gt_sam_alignment_is_unmapped(bam_align) ||
      gt_sam_alignment_mate_is_unmapped(bam_align) ||
      (gt_sam_alignment_ref_num(bam_align) ==
       gt_sam_alignment_mate_ref_num(bam_align)) ||
      !gt_sam_alignment_is_paired(bam_align) ||
      (gt_samfile_iterator_reference_length(bam_iterator,
       gt_sam_alignment_ref_num(bam_align)) < min_ref_length) ||
      (gt_sam_alignment_mapping_quality(bam_align) < min_qual))
    result = false;
  else
    result = true;

  return result;
}

static int analyze_read_set(DistRecords *dist_records,
                           ReadSet readset,
                           GtUword min_nof_pairs,
                           GtWord min_dist,
                           GtWord max_dist,
                           FragmentData fragment_data,
                           PmfData pmf_data,
                           bool rf,
                           GtSamfileIterator *bam_iterator,
                           HistogramData histogram_data,
                           GtError *err) {
  GtUword last_mref_id, nof_pairs, len_ref, len_mref = 0;
  GtWord dist;
  bool last_reverse;
  Read *read;
  int had_err = 0;
  Record *record;

  /* sort read set by orientation of read relative to reference,
     reference id of mate read and orientation of read relative
     to his mate read */
  qsort(readset.read, readset.nof, sizeof (*readset.read), compare_read_order);
  last_mref_id = readset.read[0].mtid;
  last_reverse = readset.read[0].reverse;
  len_ref = gt_samfile_iterator_reference_length(bam_iterator,
                                                 readset.read[0].tid);

  /* save reference name */
  /* resize */
  record = dist_records->record+dist_records->nof_record;
  record->root_ctg_id = gt_strdup(gt_samfile_iterator_reference_name(
                        bam_iterator,readset.read[0].tid));
  dist_records->nof_record++;

  /* iterate over (sorted) read set (characterized by same reference id) */
  for (read = readset.read; read < readset.read + readset.nof; read++) {

    len_mref = gt_samfile_iterator_reference_length(bam_iterator, read->mtid);
    /* analyze sub read set (characterized by same reference id of
       mate read and same orientation of read relative to his mate read) */
    if (last_mref_id != read->mtid || last_reverse != read->reverse) {

      nof_pairs = 0;
      /* calculate distance and nof pairs */
      if (fragment_data.nof >= min_nof_pairs) {
        had_err = estimate_dist_using_mle(&dist, &nof_pairs, min_dist,
                                         max_dist, fragment_data, pmf_data,
                                         len_ref, len_mref, rf,
                                         histogram_data, err);
      }

      /* save record */
      if (nof_pairs >= min_nof_pairs) {
        /* resize */
        record->ctg[record->nof_ctg].id =
                   gt_strdup(gt_samfile_iterator_reference_name(
                             bam_iterator, last_mref_id));
        record->ctg[record->nof_ctg].std_dev = pmf_data.std_dev /
                                               sqrt(nof_pairs);
        record->ctg[record->nof_ctg].dist = dist;
        record->ctg[record->nof_ctg].nof_pairs = nof_pairs;
        record->ctg[record->nof_ctg].sense = (last_reverse != read->reverse);
        record->ctg[record->nof_ctg].same = last_reverse;
        record->nof_ctg++;
      }

      /* reset fragment array */
      fragment_data.nof = 0;
    }

    /* calculate fragments and save unique fragments */
    calculate_fragments(read->ref_read_start, read->reverse,
                        read->mref_read_start, read->mreverse,
                        len_ref, len_mref, &fragment_data, rf);

    last_mref_id = read->mtid;
    last_reverse = read->reverse;
  }

  nof_pairs = 0;
  /* calculate distance and nof pairs */
  if (fragment_data.nof >= min_nof_pairs ) {
    had_err = estimate_dist_using_mle(&dist, &nof_pairs, min_dist,
                                    max_dist, fragment_data, pmf_data,
                                    len_ref, len_mref, rf, histogram_data, err);
  }

  /* save record */
  if (nof_pairs >= min_nof_pairs) {
    /* resize */
    record->ctg[record->nof_ctg].id =
                  gt_strdup(gt_samfile_iterator_reference_name(
                            bam_iterator, last_mref_id));
    record->ctg[record->nof_ctg].std_dev = pmf_data.std_dev / sqrt(nof_pairs);
    record->ctg[record->nof_ctg].dist = dist;
    record->ctg[record->nof_ctg].nof_pairs = nof_pairs;
    record->ctg[record->nof_ctg].sense = (!last_reverse);
    record->ctg[record->nof_ctg].same = last_reverse;
    record->nof_ctg++;
  }

  /* reset fragment array */
  fragment_data.nof = 0;

  return had_err;
}

int gt_scaffolder_bamparser_read_paired_information(DistRecords *dist,
                                                    const char *bam_filename,
                                                    const char *hist_filename,
                                                    GtWord min_dist,
                                                    GtWord max_dist,
                                                    GtUword min_qual,
                                                    GtUword min_nof_pairs,
                                                    GtUword min_ref_length,
                                                    GtUword min_align,
                                                    GtError *err) {
  int read_bytes, had_err = 0;
  GtUword last_ref_id;
  GtWord ref_read_start, mref_read_start;
  bool rf;
  FragmentData fragment_data;
  HistogramData histogram_data, histogram_data_2;
  PmfData pmf_data;
  ReadSet readset;
  GtSamfileIterator *bam_iterator;
  GtSamAlignment *bam_align;
  GtAlphabet *alpha;

  alpha = gt_alphabet_new_dna();

  /* create read set */
  readset.size = 1000;
  readset.nof = 0;
  readset.read = gt_malloc(sizeof (*readset.read) * readset.size);

  /* create fragment array */
  fragment_data.nof = 0;
  fragment_data.size = 1000;
  fragment_data.frag_pos = gt_malloc(sizeof (*(fragment_data.frag_pos))
                        * fragment_data.size * 2);
  fragment_data.frag_size = gt_malloc(sizeof (*(fragment_data.frag_pos))
                        * fragment_data.size);

  /* initialize histogram and create histogram based on hist-file */
  init_histogram(&histogram_data);
  init_histogram(&histogram_data_2);
  had_err = create_histogram_from_file(hist_filename, &histogram_data, err);
  /* create probability mass function (pmf) based on histogram */
  pmf_data = create_pmf(histogram_data);

  /* default cutoff values */
  max_dist = pmf_data.nof-1;
  fragment_data.ma = min_align;

  /* determine the orientation of the library */
  rf = histogram_data.lib_rf;
  /* if (rf)
       negate each element of histogram
  */

  /* open bam-file, read first alignment and determine its reference id */
  bam_iterator = gt_samfile_iterator_new_bam(bam_filename, alpha, NULL);
  read_bytes = gt_samfile_iterator_next(bam_iterator, &bam_align);
  last_ref_id =  gt_sam_alignment_ref_num(bam_align);

  /* iterate over reads sorted by reference id */
  while (read_bytes > 0) {

    /* analyze read set (characterized by same reference id) */
    if (last_ref_id != gt_sam_alignment_ref_num(bam_align) &&
        readset.nof != 0) {
      had_err = analyze_read_set(dist, readset, min_nof_pairs,
                min_dist, max_dist, fragment_data, pmf_data, rf,
                bam_iterator, histogram_data_2, err);
      readset.nof = 0;
    }

    /* save valid reads to read set (characterized by same reference id) */
    if (is_read_pair_valid(bam_iterator, bam_align,
                           min_ref_length, min_qual)) {

      /* resize read set */
      if (readset.nof == readset.size) {
        readset.size += 100;
        readset.read = gt_realloc(readset.read,
                           sizeof (*readset.read) * readset.size);
      }

      /* calculate reference and mate reference position at read start */
      had_err = calc_read_start(bam_align, &ref_read_start,
                                &mref_read_start, err);
      if (had_err != 0)
        break;

      readset.read[readset.nof].tid = gt_sam_alignment_ref_num(bam_align);
      readset.read[readset.nof].reverse =
                                       gt_sam_alignment_is_reverse(bam_align);
      readset.read[readset.nof].mtid = gt_sam_alignment_mate_ref_num(bam_align);
      readset.read[readset.nof].mreverse =
                                    gt_sam_alignment_mate_is_reverse(bam_align);
      readset.read[readset.nof].ref_read_start = ref_read_start;
      readset.read[readset.nof].mref_read_start = mref_read_start;
      readset.nof++;
    }

    last_ref_id =  gt_sam_alignment_ref_num(bam_align);
    read_bytes = gt_samfile_iterator_next(bam_iterator, &bam_align);
  }

  gt_sam_alignment_delete(bam_align);
  gt_samfile_iterator_delete(bam_iterator);
  gt_hashmap_delete(histogram_data.hash_map);
  gt_hashmap_delete(histogram_data_2.hash_map);
  gt_alphabet_delete(alpha);

  gt_free(readset.read);
  gt_free(histogram_data.keys);
  gt_free(histogram_data.values);
  gt_free(histogram_data_2.keys);
  gt_free(histogram_data_2.values);
  gt_free(fragment_data.frag_pos);
  gt_free(fragment_data.frag_size);
  gt_free(pmf_data.dist);

  return had_err;
}
