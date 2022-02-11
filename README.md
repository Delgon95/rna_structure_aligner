# RNA Structure Aligner

Information about the project, results and included algorithms can be found here: [Thesis](https://drive.google.com/open?id=17tl4rKGbjy284WpHZ5hWgyhpB6T08-vA)

# Build

Please run:

    ./build.sh
    
or

    mvn package

`Java` and `Maven` are required to build the project.

Tested with:
Java 11 + Maven 3.6.3

# Usage

RNA Hugs is a command line application, which requires to be provided
with:

-   `--reference` and the path to reference PDB file.
-   `--target` and the path to target PDB file.


```
java -jar target/rna-hugs-1.1-jar-with-dependencies.jar --reference reference.pdb --target target.pdb --method genetic --rmsd 2

usage: RNA Hugs
 -m,--method <method>             (optional) Method used to align input structures.
                                  Available: geometric, genetic
                                  Default: geometric
    --mode <aligning-mode>        (optional) Aligning mode used with each method. Can be either
                                  sequence independent or sequence dependent
                                  Available: seq-indep, seq-dep
                                  Default: seq-indep
 -o,--output <path>               (optional) Output directory for all results and alignements.
                                  Default: use directory of target structure.
 -r,--reference <reference.pdb>   Reference structure in .pdb format.
    --rmsd <rmsd>                 (optional) Maximum RMSD (in Ångström) that the aligned fragment
                                  will not exceed.
                                  Default: 3.5
 -t,--target <target.pdb>         Target structure in .pdb format.
    --threads <threads>           (optional) Number of threads used by algoritm. Easy way to speedup
                                  the processing.
                                  Default: all system threads
    --time-limit <seconds>        (optional) Maximum execution time in seconds before program
                                  returns. Execution time can be well below set value if no
                                  improvements are found for long time.
                                  Default: 300
```

**Warning - currently included .pdb parser was written just for the sake of showing basic functionality and should not be used  outside of this project!.**

# Methods

Both methods require coarse-grained structures as an input for alignment!</br>
Required RMSD and number of threads can be specified in AlignerConfig structure which also contains method-specific configurations.

- genetic - GeneticAligner, uses genetic metaheuristic, multithreaded. 
- geometric - GeometricAligner. Strong mathematical roots with greedy expansion of the kernel.

# Configuration

Both algorithms be further configured and fine adjusted. Configs are located in `AlignerConfig.java` file / `AlignerConfig` class.

Explanation of more important config parameters.

### General
    threads              - Number of threads used by algoritm. Easy way to speedup the processing.
    rmsdLimit            - Maximum RMSD of the aligned structure fragment. 
                           Same as 'rmsd-threshold' when using the program.
    returnTime           - Maximum execution time in seconds before program returns.
                           Execution time can be well below set value if no improvements 
                           are found for up to 'waitBuffer' time.
    waitBufferPercentage - Maximum time allowed without improvement to the result.
    waitBufferFlat       - Minimum allowed time without improvement to the result.

### Geometric
    pairRmsdLimit   - How similar (RMSD) should the dual nucleotide alignments be.
    tripleRmsdLimit - How similar should the triples (full cores) of nucleotide alignments be. 
    dualCoreBatches - How many batches of 2 nucleotide cores should be made.
                      High value means cores will have to be recalculated more times
                      but the overall RAM requirement will decrease.

### Genetic
    geometricPopulation - Use Geometric method to generate initial population for Genetic.
                          This can (it can also slow it down) speed up the final result generation.
                          This removes a large portion of randomness which can be detrimental to 
                          generating different results.
    resetThreadTime     - Time without improvement to the population for thread (one population size) to restart
                          and regenerate the population and start again.
    bestPercentage      - How many best specimens from each generation should be moved to the next one.
                          Increase can improve speed of improving the result but can limit the total search 
                          scope which can lead to lower final alignment.
    populationSize      - Number of speciments in the population. The higher the value the longer it takes for
                          each epoch but should provide better search scope and mixing of specimens.

# Examples

To compare results and performance of different methods and aligning modes please run:

    ./run_examples.sh
    
Binaries obtained from running `./build.sh` are required to run examples script.


# Output 
Both methods output AlignerOutput structure.

- **aligned** - How many nucleotides were aligned within RMSD limit.
- **referenceIndexes** - Indexes of the reference structure that were used for the alignment.
- **targetMapping** - Indexes of the target structure to which were reference structure was mapped. Index *i-th* nucleotide in **referenceIndexes** is linked to the *i-th* nucleotide in **targetMapping**.
- **superimposer** - Superimposer structure. It contains final shift and rotation matrices that can be used to align the whole structure. Main structure and required functions were rewritten from BioJava SVDSuperimposer implementation to enable more lightweight 'Coordinates' structure instead of 'Atom'.
