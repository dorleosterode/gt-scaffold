Name "graph: init graph"
Keywords "graph init"
Test do
  run("#{$bin}test.x graph 5 8 0 0 0 0 0", :retval => 0)
end

Name "graph: init vertices"
Keywords "graph vertex"
Test do
  run("#{$bin}test.x graph 5 8 1 0 0 0 0", :retval => 0)
end

Name "graph: add vertices"
Keywords "graph vertex"
Test do
  run("#{$bin}test.x graph 5 8 1 5 0 0 0", :retval => 0)
end

Name "graph: add too many vertices"
Keywords "graph vertex"
Test do
  run("#{$bin}test.x graph 5 8 1 6 0 0 0", :retval => 2)
end

Name "graph: init edges"
Keywords "graph vertex edge"
Test do
  run("#{$bin}test.x graph 5 8 1 5 1 0 0", :retval => 0)
end

Name "graph: add vertices and edges"
Keywords "graph vertex edge"
Test do
  run("#{$bin}test.x graph 5 8 1 5 1 8 0", :retval => 0)
end

# Name "graph: add vertex to unitialized struct"
# Keywords "graph vertex"
# Test do
#   run("#{$bin}test.x graph 5 8 1 0 1 8 0", :retval => 2)
# end

Name "graph: add too many edges"
Keywords "graph vertex"
Test do
  run("#{$bin}test.x graph 5 8 1 5 1 9 0", :retval => 2)
end

Name "graph: draw toy graph"
Keywords "graph dot"
Test do
  run("#{$bin}test.x graph 5 8 1 5 1 8 1", :retval => 0)
  run("diff #{$bin}gt_scaffolder_graph_test.dot #{$testdata}/gt_scaffolder_graph_test_expected.dot", :retval => 0)
end
