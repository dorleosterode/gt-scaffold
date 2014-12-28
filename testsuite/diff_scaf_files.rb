#!/usr/bin/env ruby

#compare scaffold arrays
def compare_scaffolds(scaffold_arrays_1, scaffold_arrays_2)

  diff = []
  scaffold_arrays_1.each_with_index do |scaffold_array_1, num|
    scaffold_arrays_2.each do |scaffold_array_2|
      diff = scaffold_array_1 - scaffold_array_2
      #scaffold_array_1 shares with scaffold_array_2
      #at least one common element
      if diff.size < scaffold_array_1.size
        break
      end
    end
    #show contigs existing only in first scaffold file
    if !diff.empty?
      print("< ",diff.join(", ")," in line ",num+1,"\n")
    end
  end

  scaffold_arrays_2.each_with_index do |scaffold_array_2, num|
    scaffold_arrays_1.each do |scaffold_array_1|
      diff = scaffold_array_2 - scaffold_array_1
      if diff.size < scaffold_array_2.size
        break
      end
    end
    #show contigs existing only in second scaffold file
    if !diff.empty?
      print("> ",diff.join(", ")," in line ",num+1,"\n")
    end
  end

end

#read scaffolds of scaffold file and save them as arrays
def read_scaffolds(scaffolds)
  scaffold_array = []
  scaffold_arrays = []
  scaffolds.each do |scaffold|
    scaffold_array = scaffold.split(/\t/)
    scaffold_array.each do |element|
      element.chomp!
      element.sub!(/,-?\d+,\d+\.\d+,\d+,\d+,\D*$/,"")
    end
    scaffold_arrays.push(scaffold_array)
  end
  return scaffold_arrays
end

#display total number of scaffolds, number of single scaffolds
#saved in scaffold file
def show_stat(scaffold_arrays)
  nof_single_scaffolds = 0
  nof_scaffolds = 0
  scaffold_arrays.each do |scaffold_array|
    if scaffold_array.size == 1
      nof_single_scaffolds += 1
    end
    nof_scaffolds += 1
  end
  print("Total number of scaffolds: ",nof_scaffolds,"\n")
  print("Number of single scaffolds: ",nof_single_scaffolds,"\n\n")
end

if ARGV.length != 2
  STDERR.puts "USAGE: #{$0} <scaf file> <scaf file>"
  exit 1
end

file_name_1 = ARGV[0]
file_name_2 = ARGV[1]
begin
  scaffold_file_1 = File.new(file_name_1,"r")
  scaffold_file_2 = File.new(file_name_2,"r")
rescue => err
  STDERR.puts "ERROR: Can not open file #{file_name_1}/#{file_name_2}: #{err}"
  exit 1
end

scaffolds_1 = scaffold_file_1.readlines
scaffolds_2 = scaffold_file_2.readlines

scaffold_arrays_1 = read_scaffolds(scaffolds_1)
scaffold_arrays_2 = read_scaffolds(scaffolds_2)
compare_scaffolds(scaffold_arrays_1, scaffold_arrays_2)
show_stat(scaffold_arrays_1)
show_stat(scaffold_arrays_2)

scaffold_file_1.close
scaffold_file_2.close
