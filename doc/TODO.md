## TODO
- output: convert scaffolds to fasta [Dorle]
  - traverse string graph
- parser: tranform asqg format to readjoiner-compatible format [Dorle]
- parser / new module: generate DistEst file from .bam [Lukas]
- parser: parse astat from readjoiner contig file [Lukas / Stefan]
- parser: build hashmap contig_id (headerseq after first whitespace) -> seq [Stefan]
  - evaluate md5 implementation using GtEncSeq
- phix, no error [Stefan]

--
## Presentation
  Motivation
    - only brief reminder about assembly problem
  Methods
    - in-depth review of scaffolding step (NP something)
      - intuition about single subproblems, pseudocode not necessary
      - illustrations, e.g. subgraph from dot files
  Results
    - compare gt scaffold and sga
      - time / memory (cluster, 8 GB)
        - script/rdj-spacepeak.sh <binary>
        - bookkeeping in sga -v?
      - complete genomes
        - S. cerevisiae (reason: small?)
        - C. elegans (reason: assembly competitions? sga paper)
  Discussion / Conclusion
    - evaluate reasons for "another" implementation of scaffolding algorithm
      - fewer depencencies
      - potential future ideas / implementations

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
