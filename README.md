# RNA Structure Aligner

Information about the project, results and included algorithms can be found here: [Thesis](https://drive.google.com/open?id=17tl4rKGbjy284WpHZ5hWgyhpB6T08-vA)


## Usage
RNA 3D structure aligner usage: <method> <reference.pdb> <target.pdb> <aligning-mode> <rmsd-threshold> <output-directory> (optional)</br>
Provided methods: geometric, genetic</br>
Provided aligning modes: sequence-dependent (seq-dep), sequence-independent (seq-indep). By default sequence-dependent is used.</br>
RMSD threshold: maximal RMSD score of resultant alignment. By default 3.5A is used.

Example input parameters:</br>
genetic example/1/1_solution_0_rpr.pdb example/1/testA_rpr.pdb seq-indep 3.5 example/1/output

**Warning - currently included .pdb parser was written just for the sake of showing basic functionality and should not be used  outside of this project!.**

### Methods
Both methods require coarse-grained structures as an input for alignment!</br>
Required RMSD and number of threads can be specified in AlignerConfig structure which also contains method-specific configurations.

- genetic - GeneticAligner, uses genetic metaheuristic, multithreaded. 
- geometric - GeometricAligner. Strong mathematical roots with greedy expansion of the kernel.

### Output 
Both methods output AlignerOutput structure.

- **aligned** - How many nucleotides were aligned within RMSD limit.
- **referenceIndexes** - Indexes of the reference structure that were used for the alignment.
- **targetMapping** - Indexes of the target structure to which were reference structure was mapped. Index *i-th* nucleotide in **referenceIndexes** is linked to the *i-th* nucleotide in **targetMapping**.
- **superimposer** - Superimposer structure. It contains final shift and rotation matrices that can be used to align the whole structure. Main structure and required functions were rewritten from BioJava SVDSuperimposer implementation to enable more lightweight 'Coordinates' structure instead of 'Atom'.
