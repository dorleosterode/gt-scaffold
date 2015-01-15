#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "core/minmax.h"
#include "core/error.h"
#include "core/ma_api.h"
#include "core/types_api.h"
#include "external/samtools-0.1.18/bam.h"
#include "core/hashmap_api.h"

/* BAM reads */
typedef struct Read {
  bool reverse;
  bool mreverse;
  GtUword tid;
  GtUword mtid;
  GtWord target_query_start;
  GtWord mtarget_query_start;
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

/* sort read set by orientation of read relative to reference,
   reference id of mate read and orientation of read relative
   to his mate read */
static int read_compare(const void *a, const void *b)
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

/* calculate statistics of cigar string of alignment */
static int calc_cigar_stats(bam1_t *align,
                            GtWord *target_query_start,
                            GtWord *mtarget_query_start,
                            GtError *err) {
  uint32_t *cigar;
  int op, had_err = 0;
  GtUword len, index;
  bool first = true;
  GtUword qstart = 0;
  GtUword qlen = 0;
  GtUword qspan = 0;
  GtUword tspan = 0;
  cigar = bam1_cigar(align);

  for (index = 0; index < align->core.n_cigar; index++) {
    op = cigar[index] & 0xf;
    len = cigar[index] >> 4;
    if (op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP) {
      if (first)
        qstart = len;
      qlen += len;
    }
    else if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) {
      qlen += len;
      qspan += len;
      tspan += len;
    }
    else if (op == BAM_CINS) {
      qlen += len;
      qspan += len;
    }
    else if (op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CPAD)
      tspan += len;
    else {
      had_err = -1;
      gt_error_set(err, "invalid cigar string");
      break;
    }
    first = false;
  }

  gt_assert(qstart + qspan <= qlen);
  if (bam1_strand(align))
    *target_query_start = align->core.pos + tspan + (qlen - qspan - qstart);
  else
    *target_query_start = align->core.pos - qstart;

  *mtarget_query_start = *target_query_start + align->core.isize;
  

  return had_err;
}


/* test callback function */
int show(void *key,
                     void *value,
                     void *data,
                     GtError *err) {
  GtWord diff, *key_2 = (GtUword*) key;
  GtUword *value_2 = (GtUword*) value;
  printf("%ld %ld\n",*key_2,*value_2);
  return 0;
}

/* initialize histogram */
void init_histogram(HistogramData *histogram_data){
  histogram_data->hash_map = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  histogram_data->nof_pairs = 0;
  histogram_data->size = 500;
  histogram_data->value_sum = 0;
  histogram_data->keys = gt_malloc(sizeof (*histogram_data->keys)
                        * histogram_data->size);
  histogram_data->values = gt_malloc(sizeof (*histogram_data->values)
                        * histogram_data->size);
}

/* create histogram based on fragment distribution */
HistogramData create_histogram_from_dist(FragmentData fragment_data,
                                         GtWord factor){
  HistogramData histogram_data;
  GtWord sum, *size, new_size,min, max;
  GtUword index, *value;

  /* initialize histogram  */
  init_histogram(&histogram_data);
  sum = 0;
  min = GT_WORD_MAX;
  max = GT_WORD_MIN;

  /* iterate over fragment sizes */
  for (size = fragment_data.frag_size;
       size < fragment_data.frag_size + fragment_data.nof; size++){

    /* resize histogram hashmap  */
    if (histogram_data.nof_pairs == histogram_data.size) {
      histogram_data.size += 100;
      histogram_data.values = gt_realloc(histogram_data.values,
                          sizeof(*histogram_data.values) * histogram_data.size);
      histogram_data.keys = gt_realloc(histogram_data.keys,
                          sizeof(*histogram_data.keys) * histogram_data.size);
    }

//    printf("SIZE:%ld\n",*size);

    /* fill histogram hash map */
    new_size = *size - factor;
    value = gt_hashmap_get(histogram_data.hash_map, &new_size);
    if (value == NULL) {
      histogram_data.values[histogram_data.nof_pairs] = 1;
      histogram_data.keys[histogram_data.nof_pairs] = new_size;
      gt_hashmap_add(histogram_data.hash_map,
                     histogram_data.keys + histogram_data.nof_pairs,
                     histogram_data.values + histogram_data.nof_pairs);
      histogram_data.nof_pairs++;
    }
    else
      (*value)++;

//   printf("\nKEY:%ld Value:%lu\n",*size,value==NULL?0:*value);

    sum += *size;

    if (min > new_size)
      min = new_size;
    if (max < new_size)
      max = new_size;

  }
//  gt_hashmap_foreach(histogram_data.hash_map, show, NULL, NULL);

  /* save mean */
  histogram_data.mean = (double)sum / fragment_data.nof;
  histogram_data.value_sum = fragment_data.nof;
  histogram_data.min = min;
  histogram_data.max = max;

  return histogram_data;
}

/* create histogram based on hist-file */
int create_histogram_from_file(const char *file_name,
                               HistogramData *histogram_data,
                               GtError *err){
  FILE *file;
  int had_err;
  char line[1024];
  GtWord key, value, squares, total, min, max;
  GtUword value_sum, nof_rf, nof_fr;
  double variance;

  had_err = 0;
  file = fopen(file_name, "rb");
  if (file == NULL) {
    had_err = -1;
    gt_error_set(err, "can not read file %s",file_name);
  }

  nof_rf = 0;
  nof_fr = 0;
  squares = 0;
  total = 0;
  value_sum = 0;
  min = GT_WORD_MAX;
  max = GT_WORD_MIN;

  if (had_err != -1) {
    /* iterate over lines of hist-file */
    while (fgets(line, 1024, file) != NULL) {
      /* resize histogram hashmap */
      if (histogram_data->nof_pairs == histogram_data->size) {
        histogram_data->size += 100;
        histogram_data->values = gt_realloc(histogram_data->values,
                       sizeof(*histogram_data->values) * histogram_data->size);
        histogram_data->keys = gt_realloc(histogram_data->keys,
                       sizeof(*histogram_data->keys) * histogram_data->size);
      }
      /* remove noise, outliers, apply trim_fraction */

      /* parse line of hist-file and save key/value pair in hashmap */
      if (sscanf(line,"" GT_WD "\t" GT_WD "", &key, &value) == 2 && value >= 0) {
        if (key >= 0) {
          histogram_data->values[histogram_data->nof_pairs] = (GtUword) value;
          histogram_data->keys[histogram_data->nof_pairs] = key;
          gt_hashmap_add(histogram_data->hash_map,
                           histogram_data->keys + histogram_data->nof_pairs,
                           histogram_data->values + histogram_data->nof_pairs);
          histogram_data->nof_pairs++;

          value_sum += value;
          total += key * value;
          squares += key * key * value;
          if (min > key)
            min = key;
          if (max < key)
            max = key;          
        }

        /* count keys in range*/
        if (key >= INT_MIN && key <= 0)
          nof_rf++;
        if (key >= 1 && key <= INT_MAX)
          nof_fr++;

      }
      else {
        had_err = -1;
        gt_error_set(err, "negative count value in histogram file %s",file_name);
        break;
      }
    }
  }

  histogram_data->lib_rf = nof_fr < nof_rf;
  histogram_data->mean = (double)total / value_sum;
  variance = (squares - ((double)total * total / value_sum)) / value_sum;
  histogram_data->std_dev = sqrt(variance);
  histogram_data->min = min;
  histogram_data->max = max;
  histogram_data->value_sum = value_sum;

//   gt_hashmap_foreach(histogram_data->hash_map, show, NULL, NULL);

  return had_err;
}

/* create probability mass function (pmf) based on histogram */
PmfData create_pmf(HistogramData histogram_data){
  PmfData pmf_data;
  GtUword key, *nof;
  double *dist;

  pmf_data.mean = histogram_data.mean;
  pmf_data.std_dev = histogram_data.std_dev;
  pmf_data.minp = 1.0 / histogram_data.value_sum;
  pmf_data.dist = gt_malloc(sizeof (*pmf_data.dist) * (histogram_data.max + 1));
  pmf_data.nof = histogram_data.max + 1;
//  printf("PMF: (counts:%ld)\n",histogram_data.value_sum);
  for (key = 0; key <= histogram_data.max; key++){
    nof = gt_hashmap_get(histogram_data.hash_map, &key);
    if (nof != NULL && *nof > 0)
      pmf_data.dist[key] = (double)*nof / histogram_data.value_sum;
    else
      pmf_data.dist[key] = pmf_data.minp;

//    printf("%f\n",pmf_data.dist[key]);
  }

  return pmf_data;
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
   Fehler setzen ?!
  return 0;
}
*/

double window(GtWord x1, GtWord x2, GtWord x){
  GtWord return_val, x3 = x1 + x2;

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
int compare_pmf_histogram(void *key,
                          void *value,
                          void *data,
                          GtError *err) {
  GtWord *key_2 = (GtWord*) key;
  GtUword *value_2 = (GtUword*) value;
  CompareData *compare_data = (CompareData*) data;
  double pmf_prob;
  GtWord index = *key_2 + compare_data->dist;

  if (index < compare_data->pmf_data.nof)
    pmf_prob = compare_data->pmf_data.dist[index];
  else
    pmf_prob = compare_data->pmf_data.minp;

  compare_data->likelihood += *value_2 * log(pmf_prob);

//  printf("PROB:%f LIKELIHOOD:%f KEY:%ld+%ld VAL:%lu\n", pmf_prob, compare_data->likelihood, *key_2,compare_data->dist, *value_2);
  if (pmf_prob > compare_data->pmf_data.minp)
    compare_data->nof_pairs += *value_2;

  /* Fehler setzen ?! */
  return 0;
}

int compare(const void *a, const void *b){
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
double compute_likelihood(GtWord theta,
                          HistogramData histogram_data,
                          PmfData pmf_data,
                          GtUword *nof_pairs)
{
  CompareData compare_data;

  compare_data.likelihood = 0;
  compare_data.nof_pairs = 0;
  compare_data.pmf_data = pmf_data;
  compare_data.dist = theta;
  gt_hashmap_foreach_ordered(histogram_data.hash_map,
                     compare_pmf_histogram, &compare_data, compare, NULL);
  *nof_pairs = compare_data.nof_pairs;

  return  compare_data.likelihood;
}

GtWord maximum_likelihood_estimate(GtUword *nof_pairs,
                                   GtWord min_dist,
                                   GtWord max_dist,
                                   HistogramData histogram_data,
                                   PmfData pmf_data,
                                   GtUword len_ref,
                                   GtUword len_mref){

  /* When randomly selecting fragments that span a given point,
   * longer fragments are more likely to be selected than
   * shorter fragments.
  */

  GtUword best_n = 0, n;
  GtWord best_theta = min_dist, index, theta;
  double pmf_prob, best_likelihood = FLT_MIN, c = 0, likelihood;

  min_dist = MAX(min_dist, 0 - histogram_data.max);
  max_dist = MIN(max_dist, pmf_data.nof - histogram_data.min - 1);


  printf("==========================================================\n");
  for (theta = min_dist; theta <= max_dist; theta++) {
  /* Calculate the normalizing constant of the PMF, f_theta(x) */

    printf("lenref:%lu lenmref:%lu pmf_data.nof:%lu\n",len_ref, len_mref,pmf_data.nof);
    for (index = 0; index <  pmf_data.nof; index++) {

      if (index < pmf_data.nof)
        pmf_prob = pmf_data.dist[index];
      else
        pmf_prob = pmf_data.minp;

      c += pmf_prob * window(len_ref, len_mref, index - theta);
    }
    printf("pmf[i]: %f i: %lu theta:%ld window():%f \n",pmf_prob, index-1, theta, window(len_ref, len_mref, index - theta-1));

    likelihood = compute_likelihood(theta, histogram_data, pmf_data, &n);
    printf("LIKELIHOOD before:%f\n",likelihood);
    likelihood -= histogram_data.value_sum * log(c);
    printf("(hist-val: %lu logc: %f)\n",histogram_data.value_sum, log(c));
    printf("LIKELIHOOD after:%f\n",likelihood);
    if (n > 0 && likelihood > best_likelihood) {
      best_likelihood = likelihood;
      best_theta = theta;
      best_n = n;
    }
  }

  *nof_pairs = best_n;
  return best_theta;
}

/* Estimate the distance between two contigs
   using maximum likelihood estimator */
GtWord estimate_dist_using_mle(GtUword *nof_pairs,
                               GtWord min_dist,
                               GtWord max_dist,
                               FragmentData fragment_data,
                               PmfData pmf_data,
                               GtUword len_ref,
                               GtUword len_mref,
                               bool rf){

  int had_err;
  HistogramData histogram_data;
  GtUword temp;
  GtWord dist;

  had_err = 0;

  len_ref -= fragment_data.ma - 1;
  len_mref -= fragment_data.ma - 1;

//  printf("MA:%ld\n",fragment_data.ma);

  /* swap reference and mate reference length */
  if (len_ref > len_mref) {
    temp = len_ref;
    len_ref = len_mref;
    len_mref = temp;
  }

  /* library is oriented reverse-forward */
  if (rf) {
    histogram_data = create_histogram_from_dist(fragment_data, 0);
    dist = maximum_likelihood_estimate(nof_pairs, min_dist, max_dist, histogram_data,
                                pmf_data, len_ref, len_mref);

  }
  /* library is oriented forward-reverse */
  /* Subtract 2*(l-1) from each sample */
  else {
    histogram_data = create_histogram_from_dist(fragment_data, 2 *
                     (fragment_data.ma - 1));
    dist = maximum_likelihood_estimate(nof_pairs, min_dist, max_dist, histogram_data,
                     pmf_data, len_ref, len_mref);
    dist = MAX(min_dist, dist - 2 * (GtWord)(fragment_data.ma - 1));
  }

  return dist;
}

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
void calculate_fragments(GtUword start_read,
                         bool read_reverse,
                         GtUword start_mread,
                         bool mread_reverse,
                         GtUword len_ref,
                         GtUword len_mref,
                         FragmentData *fragment_data,
                         bool rf){

  GtUword start_frag, end_frag, align;
  GtWord size;
  bool found;
  GtWord *pair;

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
       pair < fragment_data->frag_pos + (fragment_data->nof * 2); pair += 2){
//    printf("%ld==%ld AND %ld==%ld\n",*pair, start_frag,*(pair+1),end_frag);
    if (*pair == start_frag && *(pair+1) == end_frag) {
      found = true;
//      printf("FOUND\n");
      break;
    }
  }

  /* save unique fragments and their sizes */
  if (!found) {
    printf("start:%ld end:%ld size:%ld ",start_frag,end_frag, end_frag - start_frag);
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
//  printf("start:%ld end:%ld",start_frag,end_frag);
}

void parse_paired_info(char *bam_filename, char *hist_filename) {
  int read_bytes, had_err;
  bam_header_t *header;
  bamFile file;
  bam1_t *align;
  GtUword last_tid, last_mtid, len_ref, len_mref,
          nof_pairs, *value, min_nof_pairs, min_ref_length, min_align,counter;
  GtWord dist, target_query_start, mtarget_query_start, min_dist, max_dist;
  uint32_t min_qual;
  bool last_same, found, rf;
  double std_dev;
  FragmentData fragment_data;
  ReadSet readset;
  GtError *err;
  HistogramData histogram_data;
  PmfData pmf_data;
  Read *read;

  /* initialize */
  gt_lib_init();

  /* create error object */
  err = gt_error_new();

  /* create read set */
  align = gt_malloc(sizeof (*align));
  readset.size = 100;
  readset.nof = 0;
  readset.read = gt_malloc(sizeof (*readset.read)
                      * readset.size);

  /* create fragment array */
  fragment_data.nof = 0;
  fragment_data.size = 100;
  fragment_data.frag_pos = gt_malloc(sizeof (*(fragment_data.frag_pos))
                        * fragment_data.size * 2);
  fragment_data.frag_size = gt_malloc(sizeof (*(fragment_data.frag_pos))
                        * fragment_data.size);

  /* initialize histogram and create histogram based on hist-file */
  init_histogram(&histogram_data);
  had_err = create_histogram_from_file(hist_filename, &histogram_data, err);
  /* create probability mass function (pmf) based on histogram */
  pmf_data = create_pmf(histogram_data);

  /* default cutoff values */
  min_qual = 10;
  min_nof_pairs = 5;
  min_ref_length = 200;
  min_dist = -99;
  max_dist = pmf_data.nof-1;
  min_align = 100;
  fragment_data.ma = min_align;
  counter = 0;

  /* determine the orientation of the library */
  rf = histogram_data.lib_rf;
  /* if (rf)
       negate each element of histogram
  */

  /* open bam-file, read header and first alignment */
  file = bam_open(bam_filename,"r");
  header = bam_header_read(file);
  read_bytes = bam_read1(file, align);

  last_tid = align->core.tid;
  /* iterate over reads sorted by reference id */
  while (read_bytes > 0) {
    counter++;

    /* analyze read set (characterized by same reference id) */
    if (last_tid != align->core.tid && readset.nof != 0) {
      /* sort read set by orientation of read relative to reference,
         reference id of mate read and orientation of read relative
         to his mate read */
      qsort(readset.read, readset.nof, sizeof (*readset.read),
            read_compare);
      last_mtid = readset.read[0].mtid;
      last_same = (readset.read[0].reverse == readset.read[0].mreverse);
      printf("%s\n",header->target_name[last_tid]);

      /* iterate over (sorted) read set (characterized by same reference id) */
      for (read = readset.read; read < readset.read + readset.nof; read++) {
        len_ref = header->target_len[read->tid];
        len_mref = header->target_len[read->mtid];

        /* analyze sub read set (characterized by same reference id of
           mate read and same orientation of read relative to his mate read) */
        if (last_mtid != read->mtid ||
            last_same != (read->reverse == read->mreverse)) {
        

          nof_pairs = 0;
          /* calculate distance and nof pairs */
          if (fragment_data.nof >= min_nof_pairs && counter < 100) {

            /* another method to esimate distance
               estimate_dist_using_mean(fragment_data, pmf_data, &nof_pairs,
                                    &dist, err); */
            dist = estimate_dist_using_mle(&nof_pairs, min_dist, max_dist,
                           fragment_data, pmf_data, len_ref, len_mref, rf);
          }

          /* write reads */
          if (nof_pairs >= min_nof_pairs) {
              std_dev = pmf_data.std_dev / sqrt(nof_pairs);
              printf("  %s,%ld,%lu,%.1f",header->target_name[last_mtid],
                                         dist, nof_pairs, std_dev);
          }

          /* reset fragment array */
          fragment_data.nof = 0;
        }

        /* calculate fragments and save unique fragments */
        calculate_fragments(read->target_query_start, read->reverse,
                            read->mtarget_query_start, read->mreverse,
                            len_ref, len_mref, &fragment_data, rf);

        last_mtid = read->mtid;
        last_same = (read->reverse == read->mreverse);
      }

      if (fragment_data.nof >= min_nof_pairs && counter < 100) {

        /* another method to esimate distance
           estimate_dist_using_mean(fragment_data, pmf_data, &nof_pairs,
                                    &dist, err); */
        dist = estimate_dist_using_mle(&nof_pairs, min_dist, max_dist,
                            fragment_data, pmf_data, len_ref, len_mref, rf);
      }

      /* write reads */
      if (nof_pairs >= min_nof_pairs) {
        std_dev = pmf_data.std_dev / sqrt(nof_pairs);
        printf("  %s,%ld,%lu,%.1f\n",header->target_name[last_mtid],
                                         dist, nof_pairs, std_dev);
      }

      /* reset fragment array */
      fragment_data.nof = 0;

      readset.nof = 0;
      printf("\n");
    }

    /* ignore if read or mate of read is unmapped,
                 read and mate of read map to same ref,
                 read is not paired,
                 reference length of read is smaller than cutoff */
    if (!(align->core.flag & BAM_FUNMAP) &&
        !(align->core.flag & BAM_FMUNMAP) &&
          align->core.tid != align->core.mtid &&
         (align->core.flag & BAM_FPAIRED) &&
          header->target_len[align->core.tid] >= min_ref_length
         /*&& align->core.qual < min_qual*/) {

      /* resize read set */
      if (readset.nof == readset.size) {
        readset.size += 100;
        readset.read = gt_realloc(readset.read,
                           sizeof (*readset.read) * readset.size);
      }

      /* calculate statistics of cigar string of alignment */
      had_err = calc_cigar_stats(align, &target_query_start,
                                 &mtarget_query_start, err);
      if (had_err != 0)
        break;

      /* save current reads to read set (characterized by same reference id) */
      readset.read[readset.nof].tid = align->core.tid;
      readset.read[readset.nof].reverse = bam1_strand(align);
      readset.read[readset.nof].mtid = align->core.mtid;
      readset.read[readset.nof].mreverse = bam1_mstrand(align);
      readset.read[readset.nof].target_query_start = target_query_start;
      readset.read[readset.nof].mtarget_query_start = mtarget_query_start;
      readset.nof++;
    }

    last_tid = align->core.tid;
    read_bytes = bam_read1(file, align);
  }

  return EXIT_SUCCESS;
}
