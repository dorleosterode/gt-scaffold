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

# Create graph with wrapper construction function and delete it
./test.x graph 5 8 $FALSE 0 $FALSE 0 $FALSE
# Create graph, only initialize vertices, delete it
./test.x graph 5 8 $TRUE 0 $FALSE 0 $FALSE
# Create graph, initialize and create vertices, delete it
./test.x graph 5 8 $TRUE 5 $FALSE 0 $FALSE
# Create graph, initialize vertices and edges, create vertices, delete it
./test.x graph 5 8 $TRUE 5 $TRUE 0 $FALSE
# Create graph, initialize and create vertices and edges, delete it
./test.x graph 5 8 $TRUE 5 $TRUE 8 $FALSE
# Create graph, initialize and create vertices and edges, print it, delete it
./test.x graph 5 8 $TRUE 5 $TRUE 8 $TRUE
diff gt_scaffolder_graph_test.dot \
  $TESTDATA/gt_scaffolder_graph_test_expected.dot

./test.x parser $TESTDATA/libPE.de
diff gt_scaffolder_parser_test_read_distances.de \
  $TESTDATA/libPE.de

./test.x filter $TESTDATA/primary-contigs.fa $TESTDATA/libPE.de $TESTDATA/libPE.astat
diff gt_scaffolder_algorithms_test_filter_repeats.dot \
  $TESTDATA/gt_scaffolder_algorithms_test_filter_repeats_expected.dot
