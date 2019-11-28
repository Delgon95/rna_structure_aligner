package edu.put.ma.rna_aligner;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Semaphore;

public class GeometricAligner {
  private final AlignerConfig config;
  private final ArrayList<Nucleotide> referenceStructure;
  private final ArrayList<Nucleotide> targetStructure;
  private long startTime;

  private Semaphore semaphore = new Semaphore(1);
  private ArrayList<Integer> bestChainReference = new ArrayList<Integer>();
  private ArrayList<Integer> bestChainTarget = new ArrayList<Integer>();
  private double bestRmsd = Double.MAX_VALUE;

  GeometricAligner(final AlignerConfig _config,
      final ArrayList<Nucleotide> _referenceStructure,
      final ArrayList<Nucleotide> _targetStructure) {
    super();
    config = _config;
    referenceStructure = _referenceStructure;
    targetStructure = _targetStructure;
    startTime = System.currentTimeMillis();
  }

  public AlignerOutput calculate() {
    startTime = System.currentTimeMillis();
    // Precompute phase.
    // Calculate Distances: Computational: n^2 + m^2 | Memory n^2 + m^2
    Dist[][] referenceDistances = CalulateDistances(referenceStructure);
    Dist[][] targetDistances = CalulateDistances(targetStructure);

    // Create reference RNA Pairs: Computational: n^2 | Memory n^2
    ArrayList<NData> referencePairs = new ArrayList<NData>();
    referencePairs.ensureCapacity(referenceStructure.size() * referenceStructure.size());
    for (int i = 0; i < referenceStructure.size(); i++) {
      for (int j = i + 1; j < referenceStructure.size(); j++) {
        referencePairs.add(new NData(i, j));
      }
    }
    Collections.shuffle(referencePairs);

    ArrayList<ArrayList<NData>> validPairs = new ArrayList<ArrayList<NData>>();
    for (int i = 0; i < referencePairs.size(); i++) {
      validPairs.add(new ArrayList<NData>());
      referencePairs.get(i).rmsd = i;
    }

    for (int batch = 1; batch <= config.dualCoreBatches; batch++) {
      if (ShouldTerminate()) {
        return CreateAlignerOutput();
      }
      // Changes validPairs structure to contain n-th batch.
      CalculatePairCoresBatch(validPairs, referencePairs, referenceDistances,  targetDistances, batch);

      // Get max valid pair candidates for validPairs. 
      int maxValidCandidates = 0;
      for (final ArrayList<NData> valid : validPairs) {
        maxValidCandidates = Math.max(maxValidCandidates, valid.size());
      }

      // Result stored in the global variable.
      for (int tripleBatch = 1; tripleBatch <= 2; tripleBatch++) {
        FindTripleCoresAndCalculate(
            referenceDistances,
            targetDistances,
            referencePairs,
            validPairs,
            maxValidCandidates,
            tripleBatch);
      }
    }
    return CreateAlignerOutput();
  }

  private boolean ShouldTerminate() {
    return (System.currentTimeMillis() - startTime > config.returnTime * 1000);
  }

  private AlignerOutput CreateAlignerOutput() {
    ArrayList<Integer> sortedReference = new ArrayList<Integer>();
    ArrayList<Integer> sortedTarget = new ArrayList<Integer>();
    ArrayList<Nucleotide> nucleotidesReference = new ArrayList<Nucleotide>();
    ArrayList<Nucleotide> nucleotidesTarget = new ArrayList<Nucleotide>();

    // Sort nucleotide indexes.
    int index = 0;
    for (int i = 0; i < referenceStructure.size(); i++) {
      if (index == bestChainReference.size()) {
        break;
      }
      for (int j = 0; j < bestChainReference.size(); j++) {
        if (bestChainReference.get(j) == i) {
          sortedReference.add(bestChainReference.get(j));
          nucleotidesReference.add(referenceStructure.get(bestChainReference.get(j)));
          sortedTarget.add(bestChainTarget.get(j));
          nucleotidesTarget.add(targetStructure.get(bestChainTarget.get(j)));
          index++;
          break;
        }
      }
    }
    final Superimposer superimposer = Calculations.FitForRMSD(Nucleotide.NucleotidesToList(nucleotidesReference),
        Nucleotide.NucleotidesToList(nucleotidesTarget));

    System.out.println(System.currentTimeMillis() - startTime);

    return new AlignerOutput(sortedReference.size(), sortedReference, sortedTarget, superimposer);
  }

  private Dist[][] CalulateDistances(final ArrayList<Nucleotide> nucleotides) {
    Dist[][] distances = new Dist[nucleotides.size()][nucleotides.size()];
    for (int i = 0; i < nucleotides.size(); i++) {
      for (int j = 0; j < nucleotides.size(); j++) {
        distances[i][j] = new Dist(nucleotides.get(i), nucleotides.get(j));
      }
    }
    return distances;
  }

  // Search whole result space for cores with 2 nucleotides with RMSD within batch.
  private void CalculatePairCoresBatch(ArrayList<ArrayList<NData>> validPairs,
      ArrayList<NData> referencePairs,
      final Dist[][] referenceDistances,
      final Dist[][] targetDistances,
      int batch) {
    final double minimumRmsd = (Double.valueOf((batch - 1)) / Double.valueOf(config.dualCoreBatches)) * config.rmsdLimit;
    final double maximumRmsd = (Double.valueOf(batch) / Double.valueOf(config.dualCoreBatches)) * config.rmsdLimit;

    ForkJoinPool threadPool = null;

    try {
      threadPool = new ForkJoinPool(config.threads);

      threadPool.submit(() ->
      referencePairs.stream().parallel().forEach(ndata -> {
        if (ShouldTerminate()) {
          return;
        }
        ArrayList<Coordinates> referenceAtoms = new ArrayList<Coordinates>(Arrays.asList(new Coordinates[targetStructure.get(0).representatives.size() * 2]));
        NucleotidesToAtoms(referenceAtoms, ndata.index1, ndata.index2, referenceStructure);
        // Calculate pair candidates.
        ArrayList<NData> pairCandidates = CalculatePairCandidates(
            referenceDistances[ndata.index1][ndata.index2], 
            referenceAtoms,
            targetDistances,
            minimumRmsd, maximumRmsd);

        Collections.sort(pairCandidates, new NDataComparator());
        // Here rmsd is pair index!
        validPairs.set((int) ndata.rmsd, pairCandidates);
      })
          ).get();
    } catch (InterruptedException | ExecutionException e) {
      throw new RuntimeException(e);
    } finally {
      if (threadPool != null) {
        threadPool.shutdown();
      }
    }
  }

  private ArrayList<NData> CalculatePairCandidates(final Dist referenceDistance,
      final ArrayList<Coordinates> referenceAtoms,
      final Dist[][] targetDistances,
      final double minimumRmsd,
      final double maximumRmsd) {
    ArrayList<NData> pairCandidates = new ArrayList<NData>();

    final double similarityMax = Math.pow(maximumRmsd * 2, 2) * referenceStructure.get(0).representatives.size();

    ArrayList<Coordinates> targetAtoms = new ArrayList<Coordinates>(Arrays.asList(new Coordinates[targetStructure.get(0).representatives.size() * 2]));
    for (int i = 0; i < targetDistances.length; i++) {
      for (int j = 0; j < targetDistances.length; j++) {
        if (i != j) {
          double similarity = Dist.Similarity(referenceDistance, targetDistances[i][j]);
          // Here rmsdLimit != similarity metric as some calculations were omitted from calculations for speed
          // and moved to similarityMax variable.
          if (similarity <= similarityMax) {
            NucleotidesToAtoms(targetAtoms, i, j, targetStructure);
            // Calculate real RMSD with best rotation and shift.
            // Does not change targetAtoms.
            Superimposer superimposer = Calculations.FitForRMSD(referenceAtoms, targetAtoms);
            double rmsd = Calculations.CalculateRMSD(referenceAtoms, targetAtoms, superimposer);

            if (rmsd >= minimumRmsd && rmsd <= maximumRmsd) {
              pairCandidates.add(new NData(i, j, rmsd, superimposer));
            }
          }
        }
      }				
    }
    return pairCandidates;
  }

  // Does not create new Atoms structures, just references to the original structure.
  private void NucleotidesToAtoms(ArrayList<Coordinates> atoms, int first, int second, final ArrayList<Nucleotide> nucleotides) {
    int index = 0;
    final Nucleotide nucl1 = nucleotides.get(first);
    for (int i = 0; i < nucl1.representatives.size(); i++, index++) {
      atoms.set(index, nucl1.representatives.get(i));
    }
    final Nucleotide nucl2 = nucleotides.get(second);
    for (int i = 0; i < nucl2.representatives.size(); i++, index++) {
      atoms.set(index, nucl2.representatives.get(i));		
    }
  }

  // Does not create new Atoms structures, just references to the original structure.
  private void NucleotidesToAtoms(ArrayList<Coordinates> atoms, int third, final ArrayList<Nucleotide> nucleotides) {
    final Nucleotide nucl3 = nucleotides.get(third);
    int index = nucl3.representatives.size() * 2;
    for (int i = 0; i < nucl3.representatives.size(); i++, index++) {
      atoms.set(index, nucl3.representatives.get(i));
    }
  }

  private ArrayList<NData> FindTriplesCores(
      final Dist[][] referenceDistances,
      final Dist[][] targetDistances,
      ArrayList<Coordinates> referenceAtoms,
      ArrayList<Coordinates> targetAtoms,
      NData candidate,
      final NData ndata,
      final int index,
      final int batch) {
    ArrayList<NData> triplesCandidates = new ArrayList<NData>();
    Superimposer superimposer = null;

    final double similarityMax = Math.pow(config.rmsdLimit * 3, 2) * referenceStructure.get(0).representatives.size();

    for (int i = ndata.index2 + 1; i < referenceStructure.size(); i++) {
      // Fill last (third) nucleotide worth of atoms.
      NucleotidesToAtoms(referenceAtoms, i, referenceStructure);
      for (int j = 0; j < targetStructure.size(); j++) {
        if (j != candidate.index1 && j != candidate.index2) {
          // Fill last (third) nucleotide worth of atoms.
          double similarity = Dist.Similarity(referenceDistances[ndata.index1][ndata.index2],
              targetDistances[candidate.index1][candidate.index2],
              referenceDistances[ndata.index1][i],
              targetDistances[candidate.index1][j],
              referenceDistances[ndata.index2][i],
              targetDistances[candidate.index2][j]);
          if (similarity < similarityMax) {
            NucleotidesToAtoms(targetAtoms, j, targetStructure);
            superimposer = Calculations.FitForRMSD(referenceAtoms, targetAtoms);
            final ArrayList<Coordinates> fitTargetAtoms = Calculations.MoveAtomsForRMSD(
                targetAtoms,
                superimposer);
            double rmsd = Calculations.CalculateRMSD(referenceAtoms, fitTargetAtoms);

            if (rmsd < config.rmsdLimit) {
              triplesCandidates.add(new NData(i, j, rmsd, superimposer));
            }
          }
        }
      }
    }

    Collections.sort(triplesCandidates, new NDataComparator());

    // Cut cores up to batch preferences.
    // Add some kind of triples check. Can calculate all data again to find out that the second batch is empty.
    ArrayList<NData> triplesCandidatesBatch = new ArrayList<NData>();
    int best = (int) Math.max(config.tripleCoreBatchMinimum, triplesCandidates.size() * config.tripleCoreBestPercentage);
    int limit = Math.max(0, triplesCandidates.size() - best);
    if (batch == 1) {
      for (int c = triplesCandidates.size() - 1; c >= limit; c--) {
        triplesCandidatesBatch.add((NData) triplesCandidates.get(c).clone());
      }
    } else if (batch == 2) {
      for (int c = limit - 1; c >= 0; c--) {
        triplesCandidatesBatch.add((NData) triplesCandidates.get(c).clone());
      }
    } else {
      return null;
    }

    // Used all candidates in the first batch.
    if (batch == 1 && limit == 0) {
      // Use rmsd variable as substitute for a flag that this candidate should not be calculated again.
      candidate.rmsd = -1.0;
    }

    return triplesCandidatesBatch;
  }

  private void FindTripleCoresAndCalculate(
      final Dist[][] referenceDistances,
      final Dist[][] targetDistances,
      final ArrayList<NData> referencePairs,
      final ArrayList<ArrayList<NData>> validPairs,
      final int maxValidCandidates,
      final int batch) {
    // Check candidates from the first ones (best dual RMSD).
    // Parallelism can be moved one step up (here)
    for (int index = 0; index < maxValidCandidates; index++) {
      ForkJoinPool threadPool = null;
      final int finalIndex = index;
      try {
        threadPool = new ForkJoinPool(config.threads);

        threadPool.submit(() ->
        referencePairs.stream().parallel().forEach(ndata -> {
          if (ShouldTerminate()) {
            return;
          }

          // validPairs rmsd value represent nucleotide index.
          final ArrayList<NData> candidates = validPairs.get((int) ndata.rmsd);
          int candidateIndex = candidates.size() - 1 - finalIndex;

          if (candidateIndex < 0) {
            return;
          }

          NData candidate = candidates.get(candidateIndex);

          // All possible cores were calculated during first batch.
          if (candidate.rmsd < 0.0) {
            return;
          }

          ArrayList<Coordinates> referenceAtoms = new ArrayList<Coordinates>(Arrays.asList(new Coordinates[targetStructure.get(0).representatives.size() * 3]));
          ArrayList<Coordinates> targetAtoms = new ArrayList<Coordinates>(Arrays.asList(new Coordinates[targetStructure.get(0).representatives.size() * 3]));

          // Filled 2 nucleotides worth of Atoms
          NucleotidesToAtoms(referenceAtoms, ndata.index1, ndata.index2, referenceStructure);
          NucleotidesToAtoms(targetAtoms, candidate.index1, candidate.index2, targetStructure);

          // Find tripleCores already batch trimmed.
          final ArrayList<NData> triplesCoresCandidates = FindTriplesCores(
              referenceDistances,
              targetDistances,
              referenceAtoms,
              targetAtoms,
              candidate,
              ndata, 
              finalIndex, 
              batch);

          // returns null if ShouldTerminate was true while finding cores.
          if (triplesCoresCandidates == null) {
            return;
          }

          // Find final result from found triple cores.
          for (final NData tripleCore : triplesCoresCandidates) {
            NucleotidesToAtoms(referenceAtoms, tripleCore.index1, referenceStructure);
            NucleotidesToAtoms(targetAtoms, tripleCore.index2, targetStructure);
            FindFromTriple(referenceAtoms, targetAtoms, ndata, candidate,
                tripleCore.index1, tripleCore.index2,
                tripleCore.superimposer);
          }
        })
            ).get();
      } catch (InterruptedException | ExecutionException e) {
        throw new RuntimeException(e);
      } finally {
        if (threadPool != null) {
          threadPool.shutdown();
        }
      }
    }
  }

  private void FindFromTriple(final ArrayList<Coordinates> referenceAtoms, 
      final ArrayList<Coordinates> targetAtoms, final NData ndata, final NData candidate, int referenceIndex3,
      int targetIndex3, final Superimposer superimposer) {

    // Here are all triple nucleotide candidates.
    // Prepare for calculations.
    ArrayList<Integer> chainReference = new ArrayList<Integer>(referenceStructure.size());
    ArrayList<Integer> chainTarget = new ArrayList<Integer>(targetStructure.size());
    chainReference.add(ndata.index1);
    chainReference.add(ndata.index2);
    chainReference.add(referenceIndex3);
    chainTarget.add(candidate.index1);
    chainTarget.add(candidate.index2);
    chainTarget.add(targetIndex3);

    // Prepare target structure rotated and shifted to the triple core.
    ArrayList<Nucleotide> targetStructureMoved = new ArrayList<Nucleotide>(targetStructure.size());
    for (final Nucleotide res : targetStructure) {
      targetStructureMoved.add((Nucleotide) res.clone());
    }

    for (Nucleotide nucl : targetStructureMoved) {
      for (Coordinates representative : nucl.representatives) {
        Calculations.rotate(representative, superimposer.getRotation());
        Calculations.shift(representative, superimposer.getTranslation());
      }
    }

    CalculateTriple(chainReference, chainTarget, targetStructureMoved);
  }

  private void CalculateTriple(ArrayList<Integer> chainReference, ArrayList<Integer> chainTarget,
      ArrayList<Nucleotide> targetStructureMoved) {

    ArrayList<Nucleotide> nucleotidesReference = new ArrayList<Nucleotide>(referenceStructure.size());
    ArrayList<Nucleotide> nucleotidesTarget = new ArrayList<Nucleotide>(targetStructure.size());
    HashSet<Integer> usedReference = new HashSet<Integer>();
    HashSet<Integer> usedTarget = new HashSet<Integer>();
    for (int i = 0; i < chainReference.size(); i++) {
      usedReference.add(chainReference.get(i));
      usedTarget.add(chainTarget.get(i));
      nucleotidesReference.add(referenceStructure.get(chainReference.get(i)));
      nucleotidesTarget.add(targetStructure.get(chainTarget.get(i)));
    }

    // Precompute table with RMSD improvements.
    ArrayList<NData> precomputedDistances = PrecomputeRmsdChanges(targetStructureMoved, usedReference, usedTarget);

    NData bestCandidates = FindBestCandidate(usedReference, usedTarget, precomputedDistances, 
        CalculateRmsdSum(chainReference, chainTarget, targetStructureMoved));

    boolean finalShifted = false;

    while (bestCandidates != null) {
      chainReference.add(bestCandidates.index1);
      chainTarget.add(bestCandidates.index2);
      usedReference.add(bestCandidates.index1);
      usedTarget.add(bestCandidates.index2);
      nucleotidesReference.add(referenceStructure.get(bestCandidates.index1));
      nucleotidesTarget.add(targetStructure.get(bestCandidates.index2));

      double rmsdSum = CalculateRmsdSum(chainReference, chainTarget, targetStructureMoved);

      bestCandidates = FindBestCandidate(usedReference, usedTarget, precomputedDistances, rmsdSum);

      // Could not find new candidate. Try shifting structure. Only once. Only if close to current max.
      // Try to shift only if current structure is close to the current best. Do not shift clearly bad alignments.
      if (bestCandidates == null && !finalShifted && chainReference.size() >= bestChainReference.size() - 3) {
        finalShifted = true;

        // Shift and move.
        rmsdSum = CalculateRmsdSum(chainReference, chainTarget, targetStructureMoved);

        Superimposer superimposer = Calculations.FitForRMSD(Nucleotide.NucleotidesToList(nucleotidesReference), Nucleotide.NucleotidesToList(nucleotidesTarget));
        targetStructureMoved = new ArrayList<Nucleotide>(targetStructure.size());
        for (Nucleotide res : targetStructure) {
          targetStructureMoved.add((Nucleotide) res.clone());
        }

        for (Nucleotide nucl : targetStructureMoved) {
          for (Coordinates representative : nucl.representatives) {
            Calculations.rotate(representative, superimposer.getRotation());
            Calculations.shift(representative, superimposer.getTranslation());
          }
        }

        precomputedDistances = PrecomputeRmsdChanges(targetStructureMoved, usedReference, usedTarget);
        bestCandidates = FindBestCandidate(usedReference, usedTarget, precomputedDistances, rmsdSum);
      } 
    }

    double currentRmsd = Calculations.FitAndCalculateRMSD(Nucleotide.NucleotidesToList(nucleotidesReference), Nucleotide.NucleotidesToList(nucleotidesTarget));

    // Finished adding nucleotides to the structure.
    try {
      semaphore.acquire();
      // New alignment clearly better.
      if (bestChainReference.size() < chainReference.size()) {
        bestRmsd = currentRmsd;
        bestChainReference = chainReference;
        bestChainTarget = chainTarget;
      } else if (bestChainReference.size() == chainReference.size()) { // Possible RMSD improvement.
        if (bestRmsd > currentRmsd) {
          bestRmsd = currentRmsd;
          bestChainReference = chainReference;
          bestChainTarget = chainTarget;
        }
      }
      semaphore.release();
    } catch (InterruptedException e) {
      e.printStackTrace();
    }
  }

  private ArrayList<NData> PrecomputeRmsdChanges(final ArrayList<Nucleotide> targetStructureMoved,
      final HashSet<Integer> usedReference, final HashSet<Integer> usedTarget) {
    ArrayList<NData> rmsdSums = new ArrayList<NData>((referenceStructure.size() - usedReference.size()) * (targetStructureMoved.size() - usedTarget.size()) );
    for (int i = 0; i < referenceStructure.size(); i++) {
      if (!usedReference.contains(i)) {
        for (int j = 0; j < targetStructureMoved.size(); j++) {
          if (!usedTarget.contains(j)) {
            rmsdSums.add(new NData(i, j, CheckAddition(referenceStructure.get(i), targetStructureMoved.get(j))));
          }
        }
      }
    }
    Collections.sort(rmsdSums, new NDataComparator());
    return rmsdSums;
  }

  private double CheckAddition(final Nucleotide lhs, final Nucleotide rhs) {
    double sum = 0;
    for (int i = 0; i < lhs.representatives.size(); i++) {
      sum += Math.pow(Coordinates.getDistance(lhs.representatives.get(i),
          rhs.representatives.get(i)), 2);
    }
    return sum;
  }

  private NData FindBestCandidate(final HashSet<Integer> usedReference, final HashSet<Integer> usedTarget,
      ArrayList<NData> precomputedDistances, double calculatedRmsdSum) {
    int chainSize = usedReference.size();
    int index = -1;
    NData bestCandidates = null;

    for (int i = precomputedDistances.size() - 1; i >= 0; --i) {
      NData best = precomputedDistances.get(i);
      if (usedReference.contains(best.index1) || usedTarget.contains(best.index2)) {
        precomputedDistances.remove(i);
      } else {
        bestCandidates = best;
        index = i; // precomputedDistances.remove(i);
        break;
      }
    }

    if (bestCandidates != null) {
      double rmsd = CalculateFastRmsd(chainSize + 1, calculatedRmsdSum + bestCandidates.rmsd);
      if (rmsd < config.rmsdLimit) {
        precomputedDistances.remove(index);
        return bestCandidates;
      }
    }
    return null;
  }

  private double CalculateFastRmsd(int chainSize, double calculatedRmsdSum) {
    return Math.sqrt(calculatedRmsdSum / (chainSize * referenceStructure.get(0).representatives.size()));
  }

  private double CalculateRmsdSum(ArrayList<Integer> chainReference, ArrayList<Integer> chainTarget,
      ArrayList<Nucleotide> targetStructureMoved) {
    double rmsdSum = 0;
    for (int i = 0; i < chainReference.size(); i++) {
      for (int j = 0; j < referenceStructure.get(0).representatives.size(); j++) {
        rmsdSum += Math.pow(Coordinates.getDistance(referenceStructure.get(chainReference.get(i)).representatives.get(j),
            targetStructureMoved.get(chainTarget.get(i)).representatives.get(j)),
            2); 
      }
    }

    return rmsdSum;
  }
}