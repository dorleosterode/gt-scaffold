Name "scaffold: make scaffold"
Keywords "scaffold"
Test do
  run("#{$bin}test.x scaffold #{$testdata}/primary-contigs.fa #{$testdata}/libPE.de #{$testdata}/libPE.astat", :retval => 0)
end

Name "scaffold: mark repeats"
Keywords "scaffold repeat"
Test do
  run("diff #{$bin}gt_scaffolder_algorithms_test_mark_repeats.dot #{$testdata}gt_scaffolder_algorithms_test_mark_repeats_expected.dot", :retval => 0)
end

Name "scaffold: mark inconsistent and polymorphic"
Keywords "scaffold inconsistent polymorphic"
Test do
  run("diff #{$bin}gt_scaffolder_algorithms_test_filter.dot #{$testdata}gt_scaffolder_algorithms_test_filter_expected.dot", :retval => 0)
end

Name "scaffold: detect directed cycles"
Keywords "scaffold cycle"
Test do
  run("diff #{$bin}gt_scaffolder_algorithms_test_removecycles.dot #{$testdata}gt_scaffolder_algorithms_test_removecycles_expected.dot", :retval => 0)
end

Name "scaffold: final scaffold"
Keywords "scaffold cycle"
Test do
  run("diff #{$bin}gt_scaffolder_algorithms_test_makescaffold.dot #{$testdata}gt_scaffolder_algorithms_test_makescaffold_expected.dot", :retval => 0)
  run("#{$testsuite}diff_graph_files.rb #{$testdata}/sga_makeScaffolds.dot \
  gt_scaffolder_algorithms_test_makescaffold.dot", :retval => 0)
end
