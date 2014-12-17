## TODO
- algorithms: implement removeMultiEdgeScaffold function [Dorle]
- algorithms: mark single contigs as scaffolds [Dorle]
- graph: print scaffold only [Dorle]
- output: evaluate possibility to convert scaffolds to fasta [Lukas]
- test: parser scenarios (e.g. negative contig length, wrong input) [Lukas]
- SGA: print .dot for step-by-step comparison [Stefan]
- test: test makeScaffold and removeCycles [Stefan]

## OBSTACLES
- scan-build: fix dead init of variables or not? [Stefan]
- ScaffoldLinkValidator::visit: line 188, mark all edges in twin dir?

## DONE
- test pipeline: prettify and shellcheck for review [Stefan]
- test: ensure expected exits – found workaround [Stefan]
- test: refactor basic module [Stefan]
- test: modularized everything [Stefan]
- removed attribute index of vertices, adapted graph functions [Lukas]
- added test function for distance parser [Lukas]
- implement new script to filter DistEst-file [Dorle]
- extend graph_print-function to print the dir-information as arrowhead and tail
  [Dorle]
- implement set_vertex_status(v, status) and set_edge_status(e, status)
   to correclty mark all vertices and edges [Dorle]
