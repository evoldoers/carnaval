# CARNAVAL: Cellular Automata for RNA Virtual A-Life

DNA and RNA are fundamental to biogenesis, nanotechnology, and artificial life.
The RNA World hypothesis posits that self-replicating RNA emerged early from the primordial soup, concentrated within amphiphilic micelles or hydrothermal vents.
Many naturally-occurring ribozymes, riboswitches, and regulatory sequences provide circumstantial evidence for the RNA world.
Further, experiments in vitro demonstrate the plausibility of RNA-like self-replicators,
while simulations in silico recapitulate aspects of micelles and ribo-cells.
Nucleic acids are also powerful tools for nanotechnological AI,
having been used to build molecular logic circuits
including a basic perceptron,
data storage devices,
and programmable self-assembling materials.

In bioinformatics, simulations of RNA folding kinetics typically model
secondary structure as an abstract graph of base-pair contacts.
Many biophysical parameters of these simulations can be measured experimentally
or calculated by statistical mechanics arguments.
In general however these approaches generally rely on approximate treatments of spatial phenomena such as the entropic cost of loop closure.
Since coordinates are not included, any local concentration must explicitly be built into the model.

An alternative, more amenable to spatial structure while avoiding the cost of an all-atom simulation,
is a coarse-grained lattice-based molecular dynamics that tracks the positions of individual bases.
Several such approaches model an RNA molecule as a ``two-tolerant random walk'': bases occupy individual sites on the lattice, with two bases allowed to occupy the same site if (and only if) they form a basepair.
Previous work has investigated square, triangular and face-centered cubic lattices.

This software implements a model of template-directed RNA self-replication based on a simplified two-tolerant model of RNA folding kinetics on a cubic lattice, allowing diagonal steps.
Reaction rules extend this diffusion model allowing template-driven polymerization of nucleic acid sequences.
Simulation tests of the model reproduce experimentally determined RNA structure and demonstrate self-replication.

## Installation

Requirements:

- A C++ compiler (C++11 or later)
- [BOOST](https://www.boost.org/)

Type `make` and then `bin/carnaval -h`, and off you go.
