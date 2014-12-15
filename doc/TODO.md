## TODO
- SGA: print .dot for step-by-step comparison [Stefan]
- test pipeline: prettify and shellcheck for review [Stefan]
- implement set_vertex_status(v, status) and set_edge_status(e, status)
   to correclty mark all vertices and edges [Dorle]
- extend graph_print-function to print the dir-information as arrowhead and tail
  [Dorle]

## OBSTACLES
- scan-build: fix dead init of variables or not? [Stefan]

## DONE
- test: ensure expected exits – found workaround [Stefan]
- test: refactor basic module [Stefan]
- test: modularized everything [Stefan]
- removed attribute index of vertices, adapted graph functions [Lukas]
- added test function for distance parser [Lukas]
- implement new script to filter DistEst-file [Dorle]

