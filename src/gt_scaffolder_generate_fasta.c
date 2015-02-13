/*
  Copyright (c) 2015 Dorle Osterode, Stefan Dang, Lukas Götz
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

#include "core/array_api.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/queue_api.h"
#include "core/arraydef.h"
#include "core/cstr_api.h"

#include "match/rdj-ovlfind-dp.h"
#include "match/rdj-strgraph.h"
#include "match/rdj-spmlist.h"

#include "gt_scaffolder_graph.h"
#include "gt_scaffolder_algorithms.h"

struct GtScaffolderGraphResolveStats {
  GtUword graph_resolved;
  GtUword num_solutions;
  GtUword no_path;
  GtUword num_gaps;
  GtUword overlap_resolved;
  GtUword overlap_try;
  GtUword unresolved;
  GtUword singletons;
};

void gt_scaffolder_graph_reverse_gt_str(GtStr *str) {
  GtUword len;
  gt_assert(str != NULL);
  len = gt_str_length(str);

  if (len > 0) {
    GtUword i;
    char *rev = gt_malloc(sizeof (char) * (len + 1));
    char *cstr = gt_str_get(str);

    for (i = 0; i < len; i++) {
      rev[i] = cstr[len - i - 1];
    }
    rev[len] = '\0';

    gt_str_set(str, rev);
    gt_free(rev);
  }
}

static char complement_base(char c) {
  switch (c) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'a': return 't';
    case 't': return 'a';
    case 'c': return 'g';
    case 'g': return 'c';
    default: return 'N';
  }
}

void gt_scaffolder_graph_reverse_complement_gt_str(GtStr *str) {
  GtUword len;
  gt_assert(str != NULL);
  len = gt_str_length(str);

  if (len > 0) {
    GtUword i;
    char *rev = gt_malloc(sizeof (char) * (len + 1));
    char *cstr = gt_str_get(str);

    for (i = 0; i < len; i++) {
      rev[i] = complement_base(cstr[len - i - 1]);
    }
    rev[len] = '\0';

    gt_str_set(str, rev);
    gt_free(rev);
  }
}

static bool gt_scaffolder_graph_graph_resolve(GtScaffolderGraphEdge *edge,
                                              GtStr *resv_seq,
                                              GtStrgraph *strgraph,
                                              GtEncseq *encseq,
                                              GtHashmap *contigs) {
  GtUword i, j;
  bool found;

  gt_assert(edge != NULL);
  gt_assert(resv_seq != NULL);

  i = (GtUword) gt_hashmap_get(contigs, gt_str_get(edge->start->header_seq));
  j = (GtUword) gt_hashmap_get(contigs, gt_str_get(edge->end->header_seq));

  found = gt_strgraph_traverse_from_to(strgraph, encseq, i, j,
          edge->dist, edge->sense, resv_seq);

  return found;
}

typedef struct GtScaffolderGraphAlignmentData {
  GtUword best_dist;
  GtUword u_len;
} GtScaffolderGraphAlignmentData;

/* callback function. is called for every alignment found by gt_ovlfind_dp() */
static void gt_scaffolder_graph_best_alignment_by_distance(GtUword u_len,
                                                        GT_UNUSED GtUword v_len,
                                                           GtUword dist,
                                                           bool suff,
                                                           void *data) {
  if (suff) {
    GtScaffolderGraphAlignmentData *d = (GtScaffolderGraphAlignmentData *) data;
    /* TODO:? SGA uses similarity scores instead of distances */
    if (dist < d->best_dist) {
      /* found new minimal distance */
      d->best_dist = dist;
      d->u_len = u_len;
    }
  }
}

static bool gt_scaffolder_graph_overlap_resolve(GtScaffolderGraphEdge *edge,
                                                GtStr *seq,
                                                GtStr *next_seq,
                                                GtStr *resv_seq,
                                                GT_UNUSED GtUword max_edist,
                                                GtUword min_length) {
  GtUword upper_bound;
  /* 500 is used in sga (Algorithms/OverlapTools)*/
  GtUword max_alignment_length = 500;
  GtUword align_len;
  double max_error;
  GtScaffolderGraphAlignmentData data;

  gt_assert(edge != NULL);
  gt_assert(seq != NULL);
  gt_assert(next_seq != NULL);
  gt_assert(resv_seq != NULL);

  if (gt_str_length(seq) != 0 && gt_str_length(next_seq) != 0) {

    /* set upper bound for overlap */
    upper_bound = (GtUword) -1 * edge->dist + 3.0f * edge->std_dev;
    if (upper_bound > max_alignment_length)
      return false;

    /* calculate the length of the strings to align */
    align_len = MIN3(upper_bound, gt_str_length(seq), gt_str_length(next_seq));
    if (align_len > max_alignment_length)
      return false;

    /* max_error = max_edist / */
    /*               (float) MAX( gt_str_length(seq),
                     gt_str_length(next_seq) ); */

    /* TODO: test the max error rate of sga */
    max_error = 0.05;

    /* initialize data for callback */
    data.best_dist = GT_UWORD_MAX;
    data.u_len = 0;

    /* calculate all alignments of the to strings */
    gt_ovlfind_dp(gt_str_get(seq) + (gt_str_length(seq) - align_len),
                  align_len,
                  gt_str_get(next_seq),
                  align_len,
                  max_error,
                  GT_OVLFIND_SPM,
                  min_length,
                  false,
                  gt_scaffolder_graph_best_alignment_by_distance,
                  &data);

    if (data.best_dist < GT_UWORD_MAX) {
      /* found good overlap */
      char *str;
      GtUword resv_len = gt_str_length(next_seq) - data.u_len;

      str = gt_malloc((resv_len + 1) * sizeof (*str));
      memcpy(str, gt_str_get(next_seq) + data.u_len, resv_len);
      str[resv_len] = '\0';
      gt_str_set(resv_seq, str);

      gt_free(str);

      return true;
    }

  }
  return false;
}

static void gt_scaffolder_graph_introduce_gap(GtScaffolderGraphEdge *edge,
                                              GtUword min_gap_length,
                                              GtStr *next_seq,
                                              GtStr *resv_seq) {
  gt_assert(edge != NULL);
  gt_assert(next_seq != NULL);
  gt_assert(resv_seq != NULL);

  char *seq = NULL;
  char *next_cseq = gt_str_get(next_seq);
  GtUword next_seq_len = gt_str_length(next_seq);

  /* overlap couldn't be resolved */
  if (edge->dist < 0) {
    GtWord overlap = edge->dist * -1;
    GtUword rest_len = next_seq_len - overlap;
    GtUword len = min_gap_length + rest_len + 1;
    seq = gt_malloc(len * sizeof (*seq));

    memset(seq, 'N', min_gap_length);

    memcpy(seq + min_gap_length, next_cseq + overlap, rest_len + 1);

    gt_assert(seq[len-1] == '\0');
  }
  else {
    GtUword gap_len = MAX(edge->dist, min_gap_length);
    GtUword len = gap_len + next_seq_len + 1;
    seq = gt_malloc(len * sizeof (*seq));

    memset(seq, 'N', gap_len);

    memcpy(seq + gap_len, next_cseq, next_seq_len + 1);

    gt_assert(seq[len-1] == '\0');
  }

  gt_str_set(resv_seq, seq);
  gt_free(seq);
}

static void gt_scaffolder_graph_get_sequence(GtEncseq *encseq,
                                             GtUword seqnum,
                                             GtStr *seq) {

  GtUword pos, nof_chars, l;
  GtUword inc = 16384;
  GtArraychar contig_seq;
  char *seq_cstr;

  pos = gt_encseq_seqstartpos(encseq, seqnum);
  nof_chars = gt_encseq_seqlength(encseq, seqnum);
  GT_INITARRAY(&contig_seq, char);
  for (l = 0; l < nof_chars; l++, pos++) {
    char *c;
    GT_GETNEXTFREEINARRAY(c, &contig_seq, char, inc);
    *c = gt_encseq_get_decoded_char(encseq, pos,
                                    GT_READMODE_FORWARD);
  }

  seq_cstr = gt_malloc(sizeof (*seq_cstr) * (contig_seq.nextfreechar + 1));
  memcpy(seq_cstr, contig_seq.spacechar, contig_seq.nextfreechar);
  seq_cstr[contig_seq.nextfreechar] = '\0';
  gt_str_set(seq, seq_cstr);
  gt_free(seq_cstr);
  GT_FREEARRAY(&contig_seq, char);
}

GtStr *gt_scaffolder_graph_generate_string(GtScaffolderGraphRecord *rec,
                                           GtStr *ids,
                                           GtStrgraph *strgraph,
                                           GtEncseq *encseq,
                                           GtHashmap *contigs,
                                           struct GtScaffolderGraphResolveStats
                                           *stats) {
  GtStr *seq, *root_id;
  GtArray *id_array;
  GtUword seqnum;

  /* initialize seq with gt_str of the root-node of rec. we need
     the sequence for that */
  seqnum = (GtUword) gt_hashmap_get(contigs, gt_str_get(rec->root->header_seq));

  seq = gt_str_new();

  gt_scaffolder_graph_get_sequence(encseq, seqnum, seq);

  id_array = gt_array_new(sizeof (GtStr *));
  root_id = gt_str_clone(rec->root->header_seq);
  gt_str_append_char(root_id, '+');
  gt_array_add(id_array, root_id);

  if (gt_array_size(rec->edges) > 0) {
    GtUword i;
    GtScaffolderGraphEdge *edge;
    GtStr *resv_seq = gt_str_new();
    GtStr *out_id;
    bool resolved;
    bool root_dir;
    bool rel_comp = true;
    bool prev_comp = true;

    /* set root direction */
    edge = *(GtScaffolderGraphEdge **)gt_array_get(rec->edges, 0);
    root_dir = edge->sense;

    if (!root_dir)
      gt_scaffolder_graph_reverse_gt_str(seq);

    /* iterate over all edges in the scaffold */
    for (i = 0; i < gt_array_size(rec->edges); i++) {
      stats->num_gaps++;
      gt_str_reset(resv_seq);
      edge = *(GtScaffolderGraphEdge **)gt_array_get(rec->edges, i);

      /* store relative composition to root-contig */
      if (!edge->same)
        rel_comp = rel_comp ? false : true;

      /* try to find unique walk through graph to resolve the gap */
      resolved = gt_scaffolder_graph_graph_resolve(edge, resv_seq,
                                                   strgraph, encseq, contigs);

      if (resolved) {
        stats->graph_resolved++;
        /* check if we have to reverse complement the sequence */
        if (!prev_comp) {
          /* reverse complement the sequence */
          gt_scaffolder_graph_reverse_complement_gt_str(resv_seq);
        }

        if (!root_dir)
          gt_scaffolder_graph_reverse_gt_str(resv_seq);
      }

      /* check if the contigs overlap and resolve the overlap */
      if (!resolved) {
        GtStr *next_seq = gt_str_new();

        /* get the sequence of edge->end.  maybe this initialization
           is not needed! */
        seqnum = (GtUword) gt_hashmap_get(contigs, edge->end->header_seq);

        gt_scaffolder_graph_get_sequence(encseq, seqnum, next_seq);

        /* reverse complement the sequence */
        if (!rel_comp)
          gt_scaffolder_graph_reverse_complement_gt_str(next_seq);

        if (!root_dir)
          gt_scaffolder_graph_reverse_gt_str(next_seq);

        if (edge->dist < 0) {
          /* TODO: determine what values should be used for max_error
             and min_overlap_length.
             max_edist = max_error * MAX(|seq|,|next_seq|) */
          stats->overlap_try++;
          resolved = gt_scaffolder_graph_overlap_resolve(edge, seq,
                                                         next_seq, resv_seq,
                                                         0, 1);
        }

        /* introduce a gap between the contigs */
        if (!resolved) {
          /* TODO: calculate min_gap_length */
          gt_scaffolder_graph_introduce_gap(edge, 25, next_seq, resv_seq);
          stats->unresolved++;
        }
        else
          stats->overlap_resolved++;

        gt_str_delete(next_seq);
      }

      /* get the header of the current end-vertex and add the sense
         information */
      gt_str_append_str(seq, resv_seq);
      out_id = gt_str_clone(edge->end->header_seq);
      gt_str_append_char(out_id, rel_comp ? '+' : '-');
      gt_array_add(id_array, out_id);

      prev_comp = rel_comp;
    }

    gt_str_delete(resv_seq);

    if (!root_dir) {
      GtUword i, num_ids;
      GtStr *out_id;

      gt_scaffolder_graph_reverse_gt_str(seq);
      /* reverse the order of the ids */
      num_ids = gt_array_size(id_array);
      for (i = 0; i < num_ids; i++) {
        out_id = *(GtStr **) gt_array_get(id_array, num_ids - i - 1);
        gt_str_append_str(ids, out_id);
        gt_str_delete(out_id);
      }
    }
    else {
      GtUword i;
      GtStr *out_id;

      for (i = 0; i < gt_array_size(id_array); i++) {
        out_id = *(GtStr **) gt_array_get(id_array, i);
        gt_str_append_str(ids, out_id);
        gt_str_delete(out_id);
      }
    }
  }

  /* singleton scaffold */
  if (gt_array_size(id_array) == 1) {
    stats->singletons++;
    GtStr *out_id = *(GtStr **) gt_array_get(id_array, 0);
    gt_str_append_str(ids, out_id);
    gt_str_delete(out_id);
  }

  gt_array_delete(id_array);

  return seq;
}

/*
  needed parameters:
  - contig_filename
  - err
  - spm_filename
  - spm_suffix // maybe just assume .spm
  - scaffold_records
  - output_filename

  needed outupt:
  - had_err
 */

int gt_scaffolder_graph_generate_fasta(char *contig_file,
                                       char *spm_file,
                                       char *fasta_file,
                                       GtArray *recs,
                                       GtError *err) {

  int had_err;

  /* prepare the whole strgraph and stuff */
  GtStrArray *files = gt_str_array_new();
  GtEncseqEncoder *enc = gt_encseq_encoder_new();

  gt_str_array_add_cstr(files, contig_file);

  gt_encseq_encoder_set_input_dna(enc);
  gt_encseq_encoder_clip_desc(enc);
  had_err = gt_encseq_encoder_encode(enc, files, contig_file, err);

  gt_encseq_encoder_delete(enc);
  gt_str_array_delete(files);

  if (had_err == 0) {
    GtEncseq *encseq;
    GtStrgraph *strgraph;
    char spm_complete[1025];
    GtEncseqLoader *load = gt_encseq_loader_new();
    encseq = gt_encseq_loader_load(load, contig_file, err);

    gt_encseq_loader_delete(load);
    if (encseq == NULL) {
      return -1;
    }

    strgraph = gt_strgraph_new(gt_encseq_num_of_sequences(encseq));

    memcpy(spm_complete, spm_file, strlen(spm_file));
    memcpy(spm_complete + strlen(spm_file), ".0.spm", 6);

    gt_assert(spm_complete[strlen(spm_complete)] == '\0');

    had_err = gt_spmlist_parse(spm_complete, 0, gt_spmproc_strgraph_count,
                               (void *)strgraph, err);

    if (had_err == 0) {

      gt_strgraph_allocate_graph(strgraph, 0, encseq);

      had_err = gt_strgraph_load_spm_from_file(strgraph, 0,
                                               false, NULL,
                                               spm_file, 1, ".spm", err);
    }

    if (had_err == 0) {
      const char *desc;
      GtFile *out;
      GtStr *ids, *seq;
      GtUword i;
      GtScaffolderGraphRecord *rec;
      GtUword desc_len;
      char contig[1025];
      GtWord seq_num;
      struct GtScaffolderGraphResolveStats stats;
      GtHashmap *contigs = gt_hashmap_new(GT_HASH_STRING, gt_free_func, NULL);

      for (seq_num = 0; seq_num < gt_encseq_num_of_sequences(encseq);
           seq_num ++) {
        desc = gt_encseq_description(encseq, &desc_len, seq_num);

        memcpy(contig, desc, desc_len);
        contig[desc_len] = '\0';

        gt_hashmap_add(contigs, gt_cstr_dup(contig), (void *) seq_num);
      }

      /* prepare stats */
      stats.graph_resolved = 0;
      stats.overlap_resolved = 0;
      stats.overlap_try = 0;
      stats.unresolved = 0;
      stats.singletons = 0;
      stats.num_gaps = 0;

      out = gt_file_new(fasta_file, "w", err);

      if (out == NULL)
        return -1;

      ids = gt_str_new();
      /* need the mirrored sequence to generate the sequence of the
         traversed strgraph */
      had_err = gt_encseq_mirror(encseq, err);

      /* TODO: clean-up before return! */
      if (had_err != 0)
        return had_err;

      for (i = 0; i < gt_array_size(recs); i++) {
        rec = *(GtScaffolderGraphRecord **) gt_array_get(recs, i);
        gt_str_set(ids, "> ");
        seq = gt_scaffolder_graph_generate_string(rec, ids, strgraph,
                                                  encseq, contigs, &stats);
        /* write seq to fasta file */
        gt_file_xfputs(gt_str_get(ids), out);
        gt_file_xfputs("\n", out);
        /* TODO: print 80 chars per line */
        gt_file_xfputs(gt_str_get(seq), out);
        gt_file_xfputs("\n", out);

        gt_str_delete(seq);
        /* deleting all recs after stringgeneration */
      }
      gt_str_delete(ids);

      gt_file_delete(out);

      gt_hashmap_delete(contigs);

      /* print stats */
      printf("number of gaps attempted: %lu\n", stats.num_gaps);
      printf("resolved with graph: %lu\n", stats.graph_resolved);
      printf("tried to solve with overlap: %lu\n", stats.overlap_try);
      printf("resolved with overlap: %lu\n", stats.overlap_resolved);
      printf("unresolved: %lu\n", stats.unresolved);
      printf("singletons: %lu\n", stats.singletons);
    }
    gt_encseq_delete(encseq);
  }

  return had_err;
}
