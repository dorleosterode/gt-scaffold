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

TESTDATA=../testdata
FALSE=0
TRUE=1

# TEST BASIC MODULE
# usage: test.x graph <max_nof_vertices> <max_nof_edges> <init_vertices>
#                     <nof_vertices> <init_edges> <nof_edges> <print_graph>
./test.x graph 5 8 $FALSE 0 $FALSE 0 $FALSE
./test.x graph 5 8 $TRUE 0 $FALSE 0 $FALSE
./test.x graph 5 8 $TRUE 5 $FALSE 0 $FALSE
! ./test.x graph 5 8 $TRUE 6 $FALSE 0 $FALSE
./test.x graph 5 8 $TRUE 5 $TRUE 0 $FALSE
./test.x graph 5 8 $TRUE 5 $TRUE 8 $FALSE
! ./test.x graph 5 8 $TRUE 0 $TRUE 8 $FALSE
! ./test.x graph 5 8 $TRUE 5 $TRUE 9 $FALSE
./test.x graph 5 8 $TRUE 5 $TRUE 8 $TRUE
diff gt_scaffolder_graph_test.dot $TESTDATA/gt_scaffolder_graph_test_expected.dot

# TEST PARSER MODULE
# usage: test.x graph <DistEst file>
./test.x parser $TESTDATA/wrong_libPE_1.de
./test.x parser $TESTDATA/wrong_libPE_2.de
./test.x parser $TESTDATA/libPE.de
diff gt_scaffolder_parser_test_read_distances.de $TESTDATA/libPE.de

# TEST FILTER MODULE
# usage: test.x graph <FASTA-file with contigs> <DistEst file> <astat file>
./test.x scaffold $TESTDATA/primary-contigs.fa $TESTDATA/libPE.de $TESTDATA/libPE.astat
diff gt_scaffolder_algorithms_test_mark_repeats.dot \
  $TESTDATA/gt_scaffolder_algorithms_test_mark_repeats_expected.dot
diff gt_scaffolder_algorithms_test_filter.dot \
  $TESTDATA/gt_scaffolder_algorithms_test_filter_expected.dot
diff gt_scaffolder_algorithms_test_removecycles.dot \
  $TESTDATA/gt_scaffolder_algorithms_test_removecycles_expected.dot
diff gt_scaffolder_algorithms_test_makescaffold.dot \
  $TESTDATA/gt_scaffolder_algorithms_test_makescaffold_expected.dot
