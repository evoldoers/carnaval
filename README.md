# CARNAVAL: Cellular Automata RNA for Virtual A-Life

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

## Folding kinetics

The most basic way to run CARNAVAL is as a simulation of RNA folding kinetics.
Use `--init` to specify an RNA sequence, `--unit-moves` to specify the number of random steps per monomeric unit,
and `--folds` to log progress.

For example:

~~~~
bin/carnaval --init AAAAAAAAAGGGGGGGGGUUUUUUUUUCCCCC --unit-moves 100000000 --folds
~~~~

You can use `--csv`, `--json`, or `--bitmap` to save the posterior base-pairing probabilities in various formats.

## Template-directed polymerization

You can seed the space with monomers using `--density` and watch for the formation of sequences using `--seqs`:

~~~~
bin/carnaval --init AACCUUGG --density 0.1 --unit-moves 1000000 --seqs
~~~~

To experiment with the energy model and its effect on replication fidelity, use the options (e.g. `--gu` to change the wobble basepair energy)
or edit the JSON file representing the state of the world, which you can read and write using `--load` and `--save`.

## Folding kinetics on a surface

CARNAVAL can simulate strands tethered to an inert surface, which is the
natural geometry for a microarray or a DNA origami tile.

An `--surface` axis fills the plane at coordinate 0 normal to that axis with
sites that never diffuse and never base-pair. They are pure excluded volume.
`--tether` then places copies of the `--init` sequence on a grid of the given
spacing, each covalently linked at its 5' end to a site in that plane, extending
away from the surface:

~~~~
bin/carnaval -x 32 -y 32 -z 16 --dna --surface z --tether 4 --init GGCCAATTGC --unit-moves 100000
~~~~

That builds a 2D array of 64 tethered 10-mers on a 32x32 surface. Strands extend
along z by default on a three-dimensional board, or along y on a flat one;
`--tether-axis` overrides this.

`--pin-first` is the lighter-weight alternative when you want a single strand
held at one end without building a wall: it pins the first base of `--init` to
its lattice site, leaving the rest of the chain free.

`--dna` reads and writes T in place of U and forbids wobble pairs. `--no-wobble`
forbids wobble pairs without changing the alphabet. Note that the energy model
remains the simplified one described in the paper, with a single stacking term;
it is not a nearest-neighbour parameterization, so it should be recalibrated
before quantitative DNA work.

## Polymerization is off by default

`--bond` sets the probability that two adjacent template-bound monomers form a
covalent bond. It now defaults to 0, so a plain folding simulation contains no
ligase activity. Earlier versions parsed this parameter but never read it, and
ligated unconditionally whenever the geometry allowed. To reproduce the
template-directed replication described in the paper, set it explicitly:

~~~~
bin/carnaval --init AACCUUGG --density 0.1 --bond 1 --unit-moves 1000000 --seqs
~~~~

## Simulated time

The simulator's native time unit is the Monte Carlo sweep: one attempted move
per mobile unit. Sweeps are reported in the log and in saved board files, and
are the only quantity that can be mapped onto an experiment.

`--seconds-per-sweep` records such a mapping and makes the log and the JSON
output carry an elapsed time in seconds alongside the sweep count. The value is
not predicted by the model; it has to be fitted, by matching a simulated
relaxation to a measured one. Until it is fitted, treat sweep counts as an
arbitrary but internally consistent clock.

## Performance notes

Folding the 76-base tRNA for 76 million attempted moves, on one core:

| build | wall time | accepted moves/s |
|---|---|---|
| previous release | 5.0s | 0.87M |
| current | 1.5s | 2.83M |

Three changes account for this. Board dimensions that are powers of two let the
periodic wrap use a bit mask rather than an integer division, which matters
because that wrap sits under every cell access. The chain-adjacency test, which
rejected about 78% of all attempted moves and cost roughly a quarter of all
retired instructions, is now a single bit test against a cached per-unit mask of
legal directions; the mask depends only on where a unit's chain neighbours sit,
so it is two table lookups to rebuild and is refreshed only for the few units a
move can affect. The random number generator is xoshiro256++ rather than
Mersenne Twister, which was itself over a quarter of retired instructions.

Non-power-of-two boards still work and print a note suggesting otherwise.

`--sort-units` renumbers units in lattice order for cache locality. It is a pure
relabelling and leaves the Markov chain unchanged, but measured at only a few
percent, so it is off by default.

`--legal-moves-only` samples directions from the legal-move mask instead of
rejecting illegal ones, with the Metropolis-Hastings correction for the altered
proposal density. It preserves the equilibrium ensemble but changes the relative
mobility of units, so it distorts kinetics. It also buys only about 20% more
simulated time per second of compute, because each attempt becomes much more
expensive. It is off by default and is not recommended unless you want
equilibrium structure and do not care about rates.

The samplers were checked against each other by comparing base-pairing
probability matrices for a 12-base hairpin over 24 seeds at a temperature where
the chain mixes freely. All three agree on every pair to within Monte Carlo
error. At low temperature the chain mixes poorly enough that rare, long-lived
conformations are not reliably estimated from a single trajectory, whichever
sampler is used; run many seeds and compare matched constructs.

## Tests and measurements

~~~~
make test            # correctness regression suite
make bench           # tRNA throughput benchmark, checked against a recorded baseline
make foldtime        # folding-time measurement for the tRNA and a hammerhead
make foldtime-long   # the same at the run length used in the paper
~~~~

`make test` runs `test/regress.sh`, which checks that the board stays
internally consistent across every feature combination, that a seed reproduces
its trajectory, that a tethered board survives a JSON round trip with its
surface sites and pinning intact, that `--bond` really gates ligation, that
inert sites never leave their plane, and that `--dna` changes both the alphabet
and the wobble rule. Board validity is checked by saving: `assertValid` runs
inside `toJson`, so a successful save means cell occupancy and every chain link
agree with the unit list.

`make bench` folds the tRNA and reports attempted moves, accepted moves and
sweeps per second. It compares against `test/baseline.json` and fails if
accepted moves per second has regressed by more than a factor of two. Record a
new baseline for a different machine with `python3 test/benchmark.py --update`.

`make foldtime` measures two quantities that are easy to conflate. The first is
the first-passage time to a reference structure: how many sweeps until each
helix, and then all of them together, first appears. That is a distribution over
seeds rather than a single number, and some seeds never arrive, so it is
reported as a median and quartiles over the seeds that did alongside the
fraction that did. No mean is given, because with censored runs a mean needs an
assumed distribution. The second is the equilibration of the fold energy: how
many sweeps until the seed-averaged energy stops drifting. The energy settles
long before the reference structure is found, so it is the weaker claim.

Reference constructs live in `test/constructs.py`. The tRNA is PDB entry 1ehz,
with the sequence and the secondary structure taken from the paper's data
(`carnaval-paper/data/1ehz.seq`, whose second line is the PDB annotation);
parsing it gives the four canonical helices plus one tertiary pair. The
hammerhead is a hammerhead-topology construct rather than a natural isolate: a
three-way junction carrying the two published core consensus motifs as its
junction linkers, a GNRA tetraloop on stem II as in the viral hammerheads, and a
GUC cleavage triplet in the stem III loop. Its stems are designed, so the target
structure is known by construction instead of relying on an annotation.
Substituting a natural sequence means replacing the sequence and the helix pair
lists, and wants an annotated structure for that isolate.

### Measured results

On one core of an Intel Xeon at 2.10GHz, folding the tRNA for 76 million
attempted moves takes 1.33s, at 3.25M accepted moves per second and 753,000
sweeps per second. That is the recorded baseline in `test/baseline.json`.

Folding times, 16 seeds each, at a budget of 100 million sweeps: 7.6 x 10^9
moves for the 76-base tRNA and 5.1 x 10^9 for the 51-base hammerhead. That is
ten times the run length used in the paper, which is what it takes to make the
medians estimable rather than censored. Quantiles are over the full sample with
non-arrivals ordered last.

| target | seeds reaching | median (sweeps) | median (moves) | q25 | q75 |
|---|---|---|---|---|---|
| tRNA T-arm | 16/16 | 462,500 | 3.5 x 10^7 | 118,750 | 3,281,250 |
| tRNA D-arm | 16/16 | 525,000 | 4.0 x 10^7 | 150,000 | 2,568,750 |
| tRNA anticodon arm | 16/16 | 1,687,500 | 1.3 x 10^8 | 262,500 | 5,800,000 |
| tRNA acceptor stem | 16/16 | 2,162,500 | 1.6 x 10^8 | 787,500 | 7,518,750 |
| tRNA all four helices | 15/16 | 13,512,500 | 1.0 x 10^9 | 2,337,500 | 43,143,750 |
| tRNA tertiary pair | 11/16 | 61,387,500 | 4.7 x 10^9 | 28,200,000 | >budget |
| hammerhead stem II | 16/16 | 100,000 | 5.1 x 10^6 | 43,750 | 300,000 |
| hammerhead stem III | 16/16 | 512,500 | 2.6 x 10^7 | 25,000 | 17,618,750 |
| hammerhead stem I | 16/16 | 1,250,000 | 6.4 x 10^7 | 318,750 | 2,125,000 |
| hammerhead all three stems | 16/16 | 2,687,500 | 1.4 x 10^8 | 400,000 | 17,618,750 |
| hammerhead stems + core unpaired | 16/16 | 6,312,500 | 3.2 x 10^8 | 1,543,750 | 26,750,000 |

Four things are worth taking from this.

The native tRNA fold takes about 10^9 moves to find, with a median of 13.5
million sweeps. That is roughly 1.4 times the entire run length used in the
paper, which lands near the 25th percentile of the distribution. The paper's
figure of a tRNA structure recovered in 7.6 x 10^8 steps was therefore a
better-than-median run, not a typical one.

The distribution is extremely broad, and that matters more than the median.
Arrivals at the full tRNA fold spanned 1.15 to 93.4 million sweeps, a factor of
eighty, with an interquartile ratio of eighteen. Any single run tells you almost
nothing. Quote a median over many seeds with its censoring fraction.

The energy plateau is itself budget-dependent, so "equilibrated" is not a
property this measurement establishes. Raising the budget tenfold moved the
tRNA plateau from 69.8 to 73.0 and the time to reach 95% of it from 700,000 to
1,400,000 sweeps. The longer runs kept finding better structures, so the
plateau is a lower bound that grows with the budget rather than a converged
value. Treat the energy equilibration time as a floor, not an estimate.

The long-range closing helix is the bottleneck, and the slow cases are
bimodal rather than merely variable. The tRNA acceptor stem, which spans the
whole chain, has a median first passage four times that of the T-arm. The
hammerhead's stem III shows eight seeds arriving at the 25,000-sweep sampling
floor and eight more spread from 975,000 to 74.8 million, which is kinetic
trapping rather than a spread around a central value. That is the behaviour to
expect from a model with no bond rigidity term, where a long-range contact has
to be found by undirected diffusion.

Even once folded, the target structure is not the whole ensemble. Over the
whole run the tRNA spends 19.4% +- 3.7% of sampled time with all four helices
present, and the hammerhead 76.6% +- 4.8% with all three stems. The tRNA figure
is low partly because much of each run precedes first arrival.
