# CARNAVAL: Cellular Automata for RNA Virtual A-Life

DNA and RNA are fundamental to biogenesis, nanotechnology, and artificial life.
The [RNA World](https://en.wikipedia.org/wiki/RNA_world) hypothesis posits that self-replicating RNA emerged early from the primordial soup, concentrated within amphiphilic micelles or hydrothermal vents.
Many naturally-occurring ribozymes, riboswitches, and regulatory sequences provide circumstantial evidence for the RNA world.
Further, experiments _in vitro_ demonstrate the plausibility of RNA-like [self-replicators](https://www.ncbi.nlm.nih.gov/pubmed/24388759),
while simulations _in silico_ recapitulate aspects of micelles and [ribo-cells](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3303737/).
Nucleic acids are also powerful tools for nanotechnological AI,
having been used to build molecular [logic circuits](https://en.wikipedia.org/wiki/Toehold_mediated_strand_displacement)
including a basic perceptron,
data storage [devices](https://en.wikipedia.org/wiki/DNA_digital_data_storage),
and programmable self-assembling [materials](https://en.wikipedia.org/wiki/DNA_origami).

In bioinformatics, simulations of RNA folding kinetics typically model
secondary structure as an [abstract graph](https://en.wikipedia.org/wiki/Nucleic_acid_structure_prediction) of base-pair contacts.
Many biophysical parameters of these simulations can be measured experimentally
or calculated by statistical mechanics arguments.
In general however these approaches generally rely on approximate treatments of spatial phenomena such as the entropic cost of loop closure.
Since coordinates are not included, any local concentration must explicitly be built into the model.

An alternative, more amenable to spatial structure while avoiding the cost of an all-atom simulation,
is a coarse-grained lattice-based molecular dynamics that tracks the positions of individual bases.
Several such approaches model an RNA molecule as a ``two-tolerant random walk'': bases occupy individual sites on the lattice, with two bases allowed to occupy the same site if (and only if) they form a basepair.
Previous work has investigated [square](https://doi.org/10.1103/PhysRevE.68.051904), [triangular](https://www.ncbi.nlm.nih.gov/pubmed/18020662) and [face-centered cubic](https://www.ncbi.nlm.nih.gov/pubmed/20210413) lattices.

This software implements a model of template-directed RNA self-replication based on a simplified two-tolerant model of RNA folding kinetics on a cubic lattice, allowing diagonal steps.
Reaction rules extend this diffusion model allowing template-driven polymerization of nucleic acid sequences.
Simulation tests of the model reproduce experimentally determined RNA structure and demonstrate self-replication.

## Installation

Requirements:

- A C++ compiler (C++11 or later), e.g. [clang](https://clang.llvm.org/)
- The [BOOST](https://www.boost.org/) library
- GNU Make, or equivalent e.g. [BioMake](https://github.com/evoldoers/biomake)

Type `make` and then `bin/carnaval -h`, and off you go.
