Name "gt scaffolder graph: init graph"
Keywords "graph init"
Test do
  run("#{$bin}test.x graph 5 8 0 0 0 0 0", :retval => 0)
end

Name "gt scaffolder graph: init vertices"
Keywords "graph vertex"
Test do
  run("#{$bin}test.x graph 5 8 1 0 0 0 0", :retval => 0)
end

Name "gt scaffolder graph: add vertices"
Keywords "graph vertex"
Test do
  run("#{$bin}test.x graph 5 8 1 5 0 0 0", :retval => 0)
end

Name "gt scaffolder graph: add too many vertices"
Keywords "graph vertex"
Test do
  run("#{$bin}test.x graph 5 8 1 6 0 0 0", :retval => 2)
end

Name "gt scaffolder graph: init edges"
Keywords "graph vertex edge"
Test do
  run("#{$bin}test.x graph 5 8 1 5 1 0 0", :retval => 0)
end

Name "gt scaffolder graph: add vertices and edges"
Keywords "graph vertex edge"
Test do
  run("#{$bin}test.x graph 5 8 1 5 1 8 0", :retval => 0)
end

# Name "gt scaffolder graph: add vertex to unitialized struct"
# Keywords "graph vertex"
# Test do
#   run("#{$bin}test.x graph 5 8 1 0 1 8 0", :retval => 2)
# end

Name "gt scaffolder graph: add too many edges"
Keywords "graph vertex"
Test do
  run("#{$bin}test.x graph 5 8 1 5 1 9 0", :retval => 2)
end

Name "gt scaffolder graph: draw toy graph"
Keywords "graph dot"
Test do
  run("#{$bin}test.x graph 5 8 1 5 1 8 1", :retval => 0)
  run("diff $(pwd)/gt_scaffolder_graph_test.dot #{$testdata}/gt_scaffolder_graph_test_expected.dot", :retval => 0)
end

Name "gt scaffolder parser: read erroneous .de file (1/2)"
Keywords "parser DistEst"
Test do
  run("#{$bin}test.x parser #{$testdata}/wrong_libPE_1.de", :retval => 0)
end

Name "gt scaffolder parser: read erroneous .de file (2/2)"
Keywords "parser DistEst"
Test do
  run("#{$bin}test.x parser #{$testdata}/wrong_libPE_2.de", :retval => 0)
end

Name "gt scaffolder parser: parse .de file"
Keywords "graph dot"
Test do
  run("#{$bin}test.x parser #{$testdata}/libPE.de", :retval => 0)
  run("diff $(pwd)/gt_scaffolder_parser_test_read_distances.de #{$testdata}/libPE.de", :retval => 0)
end

# TODO: Decide whether to combine all tests of the scaffold module into a single
#       one. At the moment the module is executed for every single test,
#       effectively generating the same results several times (it's fast!).
#       Advantage: No redundant computation of scaffolds.
#       Disadvantage: Harder to see, where an potential error occured.
Name "gt scaffolder scaffold: make scaffold"
Keywords "scaffold"
Test do
  run("#{$bin}test.x scaffold #{$testdata}/primary-contigs.fa #{$testdata}/libPE.de #{$testdata}/libPE.astat false", :retval => 0)
end

Name "gt scaffolder scaffold: mark repeats"
Keywords "scaffold repeat"
Test do
  run("#{$bin}test.x scaffold #{$testdata}/primary-contigs.fa #{$testdata}/libPE.de #{$testdata}/libPE.astat false", :retval => 0)
  run("diff $(pwd)/gt_scaffolder_algorithms_test_mark_repeats.dot #{$testdata}gt_scaffolder_algorithms_test_mark_repeats_expected.dot", :retval => 0)
end

Name "gt scaffolder scaffold: mark inconsistent and polymorphic"
Keywords "scaffold inconsistent polymorphic"
Test do
  run("#{$bin}test.x scaffold #{$testdata}/primary-contigs.fa #{$testdata}/libPE.de #{$testdata}/libPE.astat false", :retval => 0)
  run("diff $(pwd)/gt_scaffolder_algorithms_test_filter.dot #{$testdata}gt_scaffolder_algorithms_test_filter_expected.dot", :retval => 0)
end

Name "gt scaffolder scaffold: detect directed cycles"
Keywords "scaffold cycle"
Test do
  run("#{$bin}test.x scaffold #{$testdata}/primary-contigs.fa #{$testdata}/libPE.de #{$testdata}/libPE.astat false", :retval => 0)
  run("diff $(pwd)/gt_scaffolder_algorithms_test_removecycles.dot #{$testdata}gt_scaffolder_algorithms_test_removecycles_expected.dot", :retval => 0)
end

Name "gt scaffolder scaffold: final scaffold"
Keywords "scaffold cycle"
Test do
  run("#{$bin}test.x scaffold #{$testdata}/primary-contigs.fa #{$testdata}/libPE.de #{$testdata}/libPE.astat false", :retval => 0)
  run("diff $(pwd)/gt_scaffolder_algorithms_test_makescaffold.dot #{$testdata}gt_scaffolder_algorithms_test_makescaffold_expected.dot", :retval => 0)
  run("#{$testsuite}diff_graph_files.rb #{$testdata}/sga_makeScaffolds.dot \
  gt_scaffolder_algorithms_test_makescaffold.dot", :retval => 0)
end
