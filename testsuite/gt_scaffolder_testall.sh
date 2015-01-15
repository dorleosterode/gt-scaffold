#!/usr/bin/env bash

# Copyright (c) 2014 Dorle Osterode, Stefan Dang, Lukas Götz
# Copyright (c) 2014 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

set -e -x

TESTDATA=../testdata/gt-scaffold-assets
TESTDIRS="e_coli e_lambda m_tuberculosis s_aureus"
TESTSUITE=../testsuite

# DOWNLOAD / UPDATE TESTDATA
if [[ ! -d $TESTDATA ]]; then
  git clone https://github.com/stepf/gt-scaffold-assets $TESTDATA
else
  tmp=$(pwd); cd $TESTDATA && git pull && cd "$tmp"
fi

# TEST SCAFFOLD MODULE ON ALL TESTDATA
# usage: test.x graph <FASTA-file with contigs> <DistEst file> <astat file>
for dir in $TESTDIRS; do
  full_dir="$TESTDATA/$dir"
  ./test.x scaffold "$full_dir/primary-contigs.fa" "$full_dir/libPE.de" \
    "$full_dir/libPE.astat" hashmap_out.fa
  $TESTSUITE/diff_graph_files.rb "$full_dir/07_sga_makeScaffolds.dot" \
    gt_scaffolder_algorithms_test_makescaffold.dot
done

# test graph with 1 vertex and 0 edges
full_dir="$TESTDATA/phix"
! ./test.x scaffold "$full_dir/primary-contigs.fa" "$full_dir/libPE.de" \
   "$full_dir/libPE.astat" hashmap_out.fa
