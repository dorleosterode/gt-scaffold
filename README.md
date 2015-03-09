# gt Scaffolder

gt Scaffolder is a free scaffolding software and was developed as a student
project (in the research group for genome informatics) at the Center for
Bioinformatics, Hamburg. It is based on the [GenomeTools](https://github.com/genometools/genometools)-Library and is
intended to be used in a *de novo* assembly pipeline with [Readjoiner](http://www.zbh.uni-hamburg.de/en/research/genome-informatics/software/readjoiner.html)
(an assembly software available in the GenomeTools).

gt Scaffolder tries to solve the *Contig Scaffolding Problem* in the style
of [SGA Scaffold](https://github.com/jts/sga) by using a scaffold graph and thus achieves competitive
results at a low time and memory footprint. It is also fully compatible with
the SGA assembly pipeline.

Please note that gt Scaffolder is not ready for production yet. The finished
software will be included in GenomeTools. As of now it can successfully
compute scaffolds using SGA-compatible input files. Both the interface to
Readjoiner as well as the possibility to reconstruct the target sequence are
still immature and need further development.

### Building

gt Scaffolder should run on every POSIX compliant UNIX system, for example,
Linux, Mac OS X, and OpenBSD. This will build gt Scaffolder both as a stand-alone
binary and as a gt Tool.

#### Build as a part of GenomeTools
```bash
# Clone gt Scaffolder
git clone https://github.com/dorleosterode/gt-scaffold.git

# Clone and switch to gt_scaffolder branch of our GenomeTools fork
git clone https://github.com/dorleosterode/genometools.git
cd genometools
git checkout gt_scaffolder

# Set GTDIR
GTDIR=$(pwd)

# Simlink gt Scaffolder files into GenomeTools
cd src/match/
ln -s ../../../gt-scaffold/src/gt*.c .
ln -s ../../../gt-scaffold/src/gt*.h .
cd -

# Build GenomeTools with gt scaffolder as a tool
make -j4 cairo=no pango=no
cd ..

# Build gt Scaffolder stand-alone binary
cd gt-scaffold/src
make
```
