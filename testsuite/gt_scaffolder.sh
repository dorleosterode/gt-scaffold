#!/usr/bin/env bash

set -e -x

TESTDATA=../testdata

./test.x $TESTDATA/primary-contigs.fa $TESTDATA/libPE.de $TESTDATA/libPE.astat

diff gt_scaffolder_graph_test.dot \
  $TESTDATA/gt_scaffolder_graph_test_expected.dot

diff gt_scaffolder_parser_test_complete.dot \
  $TESTDATA/gt_scaffolder_parser_test_complete_expected.dot

diff gt_scaffolder_algorithms_test_filter_repeats.dot \
  $TESTDATA/gt_scaffolder_algorithms_test_filter_repeats_expected.dot
