Name "parser: read erroneous .de file (1/2)"
Keywords "parser DistEst"
Test do
  run("#{$bin}test.x parser #{$testdata}/wrong_libPE_1.de", :retval => 0)
end

Name "parser: read erroneous .de file (2/2)"
Keywords "parser DistEst"
Test do
  run("#{$bin}test.x parser #{$testdata}/wrong_libPE_2.de", :retval => 0)
end

Name "parser: parse .de file"
Keywords "graph dot"
Test do
  run("#{$bin}test.x parser #{$testdata}/libPE.de", :retval => 0)
  run("diff #{$bin}gt_scaffolder_parser_test_read_distances.de #{$testdata}/libPE.de", :retval => 0)
end
