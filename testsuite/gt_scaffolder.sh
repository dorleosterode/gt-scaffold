#!/usr/bin/env bash

set -e -x

./test.x
diff gt_scaffolder_test.dot ../testdata/gt_scaffolder_test_expected_result.dot
