## TODO
- output: convert scaffolds to fasta [Dorle]
  - traverse string graph
- parser: tranform asqc format to readjoiner-compatible format [Dorle]
- parser / new module: generate DistEst file from .bam [Lukas]
- parser: parse astat from readjoiner contig file [Lukas / Stefan]
- parser: build hashmap contig_id (headerseq after first whitespace) -> seq [Stefan]
  - filepointer?

## OBSTACLES
- algorithms: implement removeMultiEdgeScaffold function [Dorle]
- [insert line with standard deviation]
- mark all edges of end vertex of inconsistent edge in twin direction ?
  (ScaffoldVisitors.cpp: ScaffoldLinkValidator::visit: line 188)
  (...algorithms.c: ...filter: line 314)
- mark vertex as repeat if copy number below cutoff ?
  (ScaffoldVisitors.cpp: ScaffoldAStatisticVisitor::visit: line 86)
  (...algorithms.c: ...mark_repeats: line 83)

## DONE
- output: convert scaffolds to fasta [Dorle]
  - resolve overlaps
  - fill gaps
  - global alignment [Kurtz] (not needed)
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
