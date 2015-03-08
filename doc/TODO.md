## TODO

## OBSTACLES
- [insert line with standard deviation]
- mark all edges of end vertex of inconsistent edge in twin direction ?
  (ScaffoldVisitors.cpp: ScaffoldLinkValidator::visit: line 188)
  (...algorithms.c: ...filter: line 314)
- mark vertex as repeat if copy number below cutoff ? -> untrusted
  (ScaffoldVisitors.cpp: ScaffoldAStatisticVisitor::visit: line 86)
  (...algorithms.c: ...mark_repeats: line 83)

## DONE
- output: convert scaffolds to fasta [Dorle]
  - evaluate and test the implementation
  - quast evaluation [Stefan]
- copy number is missing in contigs of readjoiner, astat of readjoiner different
  to astat of SGA [Lukas]
  - copy number:
    - gt readjoiner prefilter -copynum (> *.rcn)
    - gt readjoiner assembly -copynum -astat
- creation of bam-file from aligned reads of readjoiner [Lukas]
  - alternative: custom format containing fragment sizes, paired-end information
- finish testsuite [Stefan]
- tidy code [Dorle, Lukas, Stefan]
- documentation [Dorle, Lukas, Stefan], 15-20 pages
  - UML diagram
  - bamparser
  - up-to-date results
  - potential discussion: C vs script languages
    - possibilities: traverse string graph etc

--
- build tool architecture [Dorle]
- allow min contig length as parameter [dorle]
- algorithms: mark twin edge as GIS_SCAFFOLD -> (implement
  removeMultiEdgeScaffold function not needed) [Dorle]
- create gt scaffold tool [dorle]
- parser: parse astat from readjoiner contig file [Lukas]
  - readjoiner: -astat flag
- parser / new module: generate DistEst file from .bam [Lukas]
  - look into abyss estimation (paper) (does not exist)

--
- ported test to gt testsuite [Stefan]
- finalize presentation [Dorle / Lukas / Stefan]
- N50 [Stefan]
  - verify SGA: gt seqstat -contigs <fasta>
  - own: extended/assembly_stats_calculator.h
- output: convert scaffolds to fasta [Dorle]
  - traverse string graph
- parser: tranform asqg format to readjoiner-compatible format [Dorle]

--
- output: convert scaffolds to fasta [Dorle]
  - resolve overlaps
  - fill gaps
  - global alignment (not needed)
- test: automate all tests [Stefan]
- test: test on e coli k12 [Stefan]
- graph: exit, when graph only contains 1 vertex and 0 edges [Stefan]

--
- test: diff-script: sga-dot and gt-scaffolder dot, sga-scaf and gt-scaffolder scaf[Lukas]
- test: extensively test on several test sets [Lukas]
- scan-build: fix dead init of variables [Stefan]

--
- SGA: print .dot for step-by-step comparison [Stefan]
- test: test makeScaffold and removeCycles [Stefan]
- test pipeline: prettify and shellcheck for review [Stefan]
- test: ensure expected failures – found workaround [Stefan]
- test: refactor basic module [Stefan]
- test: modularized everything [Stefan]
- removed attribute index of vertices, adapted graph functions [Lukas]
- added test function for distance parser [Lukas]
- test: parser scenarios (e.g. negative contig length, wrong input) [Lukas]
- implement new script to filter DistEst-file [Dorle]
- extend graph_print-function to print the dir-information as arrowhead and tail [Dorle]
- implement set_vertex_status(v, status) and set_edge_status(e, status)
   to correctly mark all vertices and edges [Dorle]
- algorithms: mark single contigs as scaffolds [Dorle]
- graph: print scaffold only [Dorle]
