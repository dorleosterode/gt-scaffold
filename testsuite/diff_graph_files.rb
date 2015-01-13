#!/usr/bin/env ruby

#compare hash tables with vertex and edge information of two graphs
def compare_graphs(hash_table_1, hash_table_2)

  has_diff = 0

  hash_table_1.keys.each do |vertex|
    if hash_table_2.key?(vertex)
      diff_1 = hash_table_1[vertex] - hash_table_2[vertex]
      diff_2 = hash_table_2[vertex] - hash_table_1[vertex]
      #show edge existing only in sga graph file
      if !diff_1.empty?
        print("< ",vertex," to ",diff_1.join(",")," (edge)\n")
        has_diff = 1
      end
      #show edge existing only in other graph file
      if !diff_2.empty?
        print("> ",vertex," to ",diff_2.join(",")," (edge)\n")
        has_diff = 1
      end
    #show vertex existing only in sga graph file
    else
      print("< ",vertex," (vertex)\n")
      has_diff = 1
    end
  end

  hash_table_2.keys.each do |vertex|
    #show vertex existing only in graph file
    if !hash_table_1.key?(vertex)
      print("> ",vertex," (vertex)\n")
      has_diff = 1
    end
  end

  exit has_diff

end

#parse edges and vertices of graph in SGA format
def read_graph_sga(graph)
  hash_table = Hash.new

  graph.each do |line|
    #parse edge
    if line.match(/\"(.+)\" -> \"(.+)\" \[.+\];/) != nil
      if !hash_table.key?($1)
        STDERR.puts "ERROR: edge contains unknown vertex"
        exit 1
      end
      edge_array = hash_table[$1]
      edge_array.push($2)
      hash_table[$1] = edge_array
    #parse vertex
    elsif line.match(/\"(.+)\" \[.+\];/) != nil
      hash_table[$1] = Array.new
    end
  end

  return hash_table
end

#parse edges and vertices of graph in gt-scaffolder format
def read_graph_gt_scaffolder(graph)
  hash_table = Hash.new
  id_to_header_seq = Hash.new
  visible_state = ["black", "magenta", "red", "green"]

  graph.each do |line|
    #parse edge
    if (line.match(/(\d+) -> (\d+) \[color=\"(.+?)\".+\];/) != nil &&
        visible_state.include?($3))
      if !id_to_header_seq.key?($1) || !id_to_header_seq.key?($2)
        STDERR.puts "ERROR: edge contains unknown vertex"
        exit 1
      end
      edge_array = hash_table[id_to_header_seq[$1]]
      edge_array.push(id_to_header_seq[$2])
      hash_table[id_to_header_seq[$1]] = edge_array
    #parse vertex
    elsif (line.match(/(\d+) \[color=\"(.+)\" label=\"(.+)\"\];/) != nil &&
           visible_state.include?($2))
      hash_table[$3] = Array.new
      id_to_header_seq[$1] = $3
    end
  end

  return hash_table
end

if ARGV.length != 2
  STDERR.puts "USAGE: #{$0} <sga graph file> <graph file>"
  exit 1
end

file_name_1 = ARGV[0]
file_name_2 = ARGV[1]
begin
  graph_file_1 = File.new(file_name_1,"r")
  graph_file_2 = File.new(file_name_2,"r")
rescue => err
  STDERR.puts "ERROR: Can not open file #{file_name_1}/#{file_name_2}: #{err}"
  exit 1
end

graph_1 = graph_file_1.readlines
graph_2 = graph_file_2.readlines

graph_file_1.close
graph_file_2.close

hash_table_1 = read_graph_sga(graph_1)
hash_table_2 = read_graph_gt_scaffolder(graph_2)
compare_graphs(hash_table_1, hash_table_2)
