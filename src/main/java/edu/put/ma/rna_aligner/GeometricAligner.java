package edu.put.ma.rna_aligner;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Semaphore;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.util.Precision;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GeometricAligner {
  private static final Logger LOGGER = LoggerFactory.getLogger(GeometricAligner.class);

  private final AlignerConfig config;
  private final ArrayList<Nucleotide> referenceStructure;
  private final ArrayList<Nucleotide> targetStructure;
  private long startTime;

  private Semaphore semaphore = new Semaphore(1);
  private ArrayList<Integer> bestChainReference = new ArrayList<Integer>();
  private ArrayList<Integer> bestChainTarget = new ArrayList<Integer>();
  private double bestRmsd = Double.MAX_VALUE;
  private boolean isSequenceDependent;
  private double rmsdLimit;
  private long globalStart = System.currentTimeMillis();
  private long stopTime;

  private boolean createPopulation = false;
  private int populationSize = 0;
  ArrayList<Specimen> populationPool = new ArrayList<Specimen>();

  GeometricAligner(final AlignerConfig _config, final ArrayList<Nucleotide> _referenceStructure,
      final ArrayList<Nucleotide> _targetStructure, final boolean _isSequenceDependent,
      final double _rmsdLimit) {
    super();
    config = _config;
    referenceStructure = _referenceStructure;
    targetStructure = _targetStructure;
    isSequenceDependent = _isSequenceDependent;
    //if ((isSequenceDependent) && (referenceStructure.size() != targetStructure.size()))
    //  isSequenceDependent = !isSequenceDependent;
    rmsdLimit = _rmsdLimit;
    stopTime = (long) ((globalStart + (1000 * config.returnTime * 0.25)));
  }

  public AlignerOutput calculate() {
    startTime = System.currentTimeMillis();
    // Precompute phase.
    // Calculate Distances: Computational: n^2 + m^2 | Memory n^2 + m^2
    final Dist[][] referenceDistances = CalulateDistances(referenceStructure);
    final Dist[][] targetDistances = CalulateDistances(targetStructure);

    // Create reference RNA Pairs: Computational: n^2 | Memory n^2
    ArrayList<ArrayList<NData>> validPairs = new ArrayList<ArrayList<NData>>();
    ArrayList<NData> referencePairsTmp = new ArrayList<NData>();
    referencePairsTmp.ensureCapacity(
        (referenceStructure.size() - 1) * (referenceStructure.size() - 2) / 2);
    validPairs.ensureCapacity(
        (referenceStructure.size() - 1) * (referenceStructure.size() - 2) / 2);
    int id = 0;
    for (int i = 0; i < referenceStructure.size(); i++) {
      for (int j = i + 1; j < referenceStructure.size(); j++) {
        referencePairsTmp.add(new NData(i, j, id++));
        validPairs.add(new ArrayList<NData>());
      }
    }
    Collections.shuffle(referencePairsTmp);
    final ArrayList<NData> referencePairs = referencePairsTmp;

    for (int tripleBatch = 1; tripleBatch <= 2; tripleBatch++) {
      for (int batch = 1; batch <= config.dualCoreBatches; batch++) {
        if (ShouldTerminate()) {
          return CreateAlignerOutput();
        }
        // Changes validPairs structure to contain n-th batch.
        CalculatePairCoresBatch(
            validPairs, referencePairs, referenceDistances, targetDistances, batch);

        // Get max valid pair candidates for validPairs.
        int maxValidCandidates = 0;
        for (final ArrayList<NData> valid : validPairs) {
          maxValidCandidates = Math.max(maxValidCandidates, valid.size());
        }

        // Result stored in the global variable.
        FindTripleCoresAndCalculate(referenceDistances, targetDistances, referencePairs, validPairs,
            maxValidCandidates, tripleBatch);
      }
    }
    return CreateAlignerOutput();
  }

  public ArrayList<Specimen> createPopulation(int populationSize) {
    this.populationSize = populationSize;
    createPopulation = true;
    calculate();
    return populationPool;
  }

  private void updateStopTime(Boolean resultImprovement) {
    // Made an improvement.
    long now = System.currentTimeMillis();
    long timeLeft = ((long) config.returnTime * 1000) - (now - globalStart);

    long timeImprovement = 0;
    if (resultImprovement) {
      timeImprovement =
          (long) Math.max(timeLeft * config.imprResultPercentage, config.imprResultFlat * 1000);
    } else {
      timeImprovement =
          (long) Math.max(timeLeft * config.imprRmsdPercentage, config.imprRmsdFlat * 1000);
    }
    long timeBuffer = stopTime + timeImprovement - now;
    if (timeBuffer < config.waitBufferFlat) {
      timeBuffer = config.waitBufferFlat;
    } else {
      timeBuffer = (long) Math.min(timeBuffer, timeLeft * config.waitBufferPercentage);
    }
    stopTime = (long) Math.min(now + timeBuffer, globalStart + (config.returnTime * 1000));
  }

  private boolean ShouldTerminate() {
    return (System.currentTimeMillis() > stopTime
        || (bestChainReference.size() == referenceStructure.size() || bestChainTarget.size() == targetStructure.size())
        || (populationPool.size() > populationSize
            || (createPopulation && (bestChainReference.size() == referenceStructure.size() || bestChainTarget.size() == targetStructure.size())
                && populationPool.size() > config.threads * 3)));
  }

  private AlignerOutput CreateAlignerOutput() {
    final ArrayList<Integer> sortedReference =
        new ArrayList<Integer>(Arrays.asList(new Integer[referenceStructure.size()]));

    final Integer[] sortedTargetTab = new Integer[referenceStructure.size()];
    Arrays.fill(sortedTargetTab, 0, referenceStructure.size(), -1);
    final ArrayList<Integer> sortedTarget = new ArrayList<Integer>(Arrays.asList(sortedTargetTab));

    final ArrayList<Nucleotide> nucleotidesReference = new ArrayList<Nucleotide>();
    final ArrayList<Nucleotide> nucleotidesTarget = new ArrayList<Nucleotide>();

    for (int i = 0; i < referenceStructure.size(); i++) {
      sortedReference.set(i, i);
    }
    int count = 0;
    for (int i = 0; i < bestChainReference.size(); i++) {
      final int referenceIndex = bestChainReference.get(i).intValue();
      sortedTarget.set(referenceIndex, bestChainTarget.get(i));
      count++;
    }

    for (int i = 0; i < referenceStructure.size(); i++) {
      if (sortedTarget.get(i).intValue() > -1) {
        nucleotidesReference.add(referenceStructure.get(i));
        nucleotidesTarget.add(targetStructure.get(sortedTarget.get(i).intValue()));
      }
    }
    final Superimposer superimposer =
        Calculations.FitForRMSD(Nucleotide.NucleotidesToList(nucleotidesReference),
            Nucleotide.NucleotidesToList(nucleotidesTarget));

    return new AlignerOutput(count, sortedReference, sortedTarget, superimposer,
        System.currentTimeMillis() - startTime, bestRmsd);
  }

  private Dist[][] CalulateDistances(final ArrayList<Nucleotide> nucleotides) {
    Dist[][] distances = new Dist[nucleotides.size()][nucleotides.size()];
    for (int i = 0; i < nucleotides.size(); i++) {
      for (int j = i + 1; j < nucleotides.size(); j++) {
        distances[i][j] = new Dist(nucleotides.get(i), nucleotides.get(j));
        distances[j][i] = new Dist(distances[i][j].distances);
      }
    }
    return distances;
  }

  // Search whole result space for cores with 2 nucleotides with RMSD within batch.
  private void CalculatePairCoresBatch(ArrayList<ArrayList<NData>> validPairs,
      final ArrayList<NData> referencePairs, final Dist[][] referenceDistances,
      final Dist[][] targetDistances, int batch) {
    final double minimumRmsd =
        (Double.valueOf((batch - 1)) / Double.valueOf(config.dualCoreBatches))
        * config.pairRmsdLimit;
    final double maximumRmsd =
        (Double.valueOf(batch) / Double.valueOf(config.dualCoreBatches)) * config.pairRmsdLimit;

    ForkJoinPool threadPool = null;

    try {
      threadPool = new ForkJoinPool(config.threads);

      threadPool
          .submit(() -> referencePairs.stream().parallel().forEach(ndata -> {
            if (ShouldTerminate()) {
              return;
            }
            ArrayList<Coordinates> referenceAtoms = new ArrayList<Coordinates>(
                Arrays.asList(new Coordinates[targetStructure.get(0).representatives.size() * 2]));
            NucleotidesToAtoms(referenceAtoms, ndata.index1, ndata.index2, referenceStructure);
            // Calculate pair candidates.
            final ArrayList<NData> pairCandidates =
                CalculatePairCandidates(ndata, referenceDistances[ndata.index1][ndata.index2],
                    referenceAtoms, targetDistances, minimumRmsd, maximumRmsd);

            Collections.sort(pairCandidates, new NDataComparator());
            // Here rmsd is pair index!
            validPairs.set((int) ndata.rmsd, pairCandidates);
          }))
          .get();
    } catch (InterruptedException | ExecutionException e) {
      throw new RuntimeException(e);
    } finally {
      if (threadPool != null) {
        threadPool.shutdown();
      }
    }
  }

  private ArrayList<NData> CalculatePairCandidates(final NData ndata, final Dist referenceDistance,
      final ArrayList<Coordinates> referenceAtoms,
      final Dist[][] targetDistances, final double minimumRmsd, final double maximumRmsd) {
    ArrayList<NData> pairCandidates = new ArrayList<NData>();

    final double similarityMax =
        Math.pow(maximumRmsd * 2, 2) * referenceStructure.get(0).representatives.size();
    ArrayList<Coordinates> targetAtoms1 = new ArrayList<Coordinates>(
        Arrays.asList(new Coordinates[targetStructure.get(0).representatives.size() * 2]));
    ArrayList<Coordinates> targetAtoms2 = new ArrayList<Coordinates>(
        Arrays.asList(new Coordinates[targetStructure.get(0).representatives.size() * 2]));

      for (int i = 0; i < targetDistances.length; i++) {
        for (int j = i + 1; j < targetDistances.length; j++) {
         final boolean ij = ((StringUtils.equalsIgnoreCase(targetStructure.get(i).getCode(),
             referenceStructure.get(ndata.index1).getCode())) &&
            				(StringUtils.equalsIgnoreCase(targetStructure.get(j).getCode(),
             referenceStructure.get(ndata.index2).getCode())));
         final boolean ji = ((StringUtils.equalsIgnoreCase(targetStructure.get(j).getCode(),
             referenceStructure.get(ndata.index1).getCode())) &&
            				(StringUtils.equalsIgnoreCase(targetStructure.get(i).getCode(),
             referenceStructure.get(ndata.index2).getCode())));
         if (!isSequenceDependent || (ij || ji)) {

            final double similarity = Dist.Similarity(referenceDistance, targetDistances[i][j]);
            // Here rmsdLimit != similarity metric as some calculations were omitted from calculations
            // for speed and moved to similarityMax variable.
            if (similarity <= similarityMax) {
              NucleotidesToAtoms(targetAtoms1, i, j, targetStructure);
              // Calculate real RMSD with best rotation and shift.
              // Does not change targetAtoms.
              Superimposer superimposer = Calculations.FitForRMSD(referenceAtoms, targetAtoms1);
              final double rmsd =
                  Calculations.CalculateRMSD(referenceAtoms, targetAtoms1, superimposer);

              if (rmsd >= minimumRmsd && rmsd <= maximumRmsd && (!isSequenceDependent || ij)) {
                pairCandidates.add(new NData(i, j, rmsd, superimposer));
              }
              NucleotidesToAtoms(targetAtoms2, j, i, targetStructure);
              // Calculate real RMSD with best rotation and shift.
              // Does not change targetAtoms.
              superimposer = Calculations.FitForRMSD(referenceAtoms, targetAtoms2);
              final double rmsd2 =
                  Calculations.CalculateRMSD(referenceAtoms, targetAtoms2, superimposer);

              if (rmsd2 >= minimumRmsd && rmsd2 <= maximumRmsd && (!isSequenceDependent || ji)) {
                pairCandidates.add(new NData(j, i, rmsd2, superimposer));
              }
            }
          }
        }
      }
    return pairCandidates;
  }

  // Does not create new Atoms structures, just references to the original structure.
  private void NucleotidesToAtoms(ArrayList<Coordinates> atoms, int first, int second,
      final ArrayList<Nucleotide> nucleotides) {
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
  private void NucleotidesToAtoms(
      ArrayList<Coordinates> atoms, int third, final ArrayList<Nucleotide> nucleotides) {
    final Nucleotide nucl3 = nucleotides.get(third);
    int index = nucl3.representatives.size() * 2;
    for (int i = 0; i < nucl3.representatives.size(); i++, index++) {
      atoms.set(index, nucl3.representatives.get(i));
    }
  }

  private ArrayList<NData> FindTriplesCores(final Dist[][] referenceDistances,
      final Dist[][] targetDistances, ArrayList<Coordinates> referenceAtoms,
      ArrayList<Coordinates> targetAtoms, NData candidate, final NData ndata, final int batch) {
    ArrayList<NData> triplesCandidates = new ArrayList<NData>();
    Superimposer superimposer = null;

    final double similarityMax =
        Math.pow(config.tripleRmsdLimit * 3, 2) * referenceStructure.get(0).representatives.size();

      for (int i = 0; i < referenceStructure.size(); i++) {
        if (i != ndata.index1 && i != ndata.index2) {
          // Fill last (third) nucleotide worth of atoms.
          for (int j = 0; j < targetStructure.size(); j++) {
            if (j != candidate.index1 && j != candidate.index2 &&
            (!isSequenceDependent || ((StringUtils.equalsIgnoreCase(targetStructure.get(j).getCode(),
             referenceStructure.get(i).getCode()))))) {
              // Fill last (third) nucleotide worth of atoms.
              final double similarity =
                  Dist.Similarity(referenceDistances[ndata.index1][ndata.index2],
                      targetDistances[candidate.index1][candidate.index2],
                      referenceDistances[ndata.index1][i], targetDistances[candidate.index1][j],
                      referenceDistances[ndata.index2][i], targetDistances[candidate.index2][j]);
              if (similarity < similarityMax) {
                NucleotidesToAtoms(referenceAtoms, i, referenceStructure);
                NucleotidesToAtoms(targetAtoms, j, targetStructure);
                superimposer = Calculations.FitforTripleRMSD(referenceAtoms, targetAtoms);
                double miniRmsd =
                    Calculations.CalculateRMSD(referenceAtoms, targetAtoms, superimposer);

                if (miniRmsd < config.tripleRmsdLimit) {
                  triplesCandidates.add(new NData(i, j, miniRmsd, superimposer));
                } else if (miniRmsd < config.tripleRmsdLimit * 1.15) {
                  superimposer = Calculations.FitForRMSD(referenceAtoms, targetAtoms);
                  final ArrayList<Coordinates> fitTargetAtoms =
                      Calculations.MoveAtomsForRMSD(targetAtoms, superimposer);
                  double rmsd = Calculations.CalculateRMSD(referenceAtoms, fitTargetAtoms);
                  if (rmsd < config.tripleRmsdLimit) {
                    triplesCandidates.add(new NData(i, j, rmsd, superimposer));
                  }
                }
              }
            }
          }
        }
      }

    Collections.sort(triplesCandidates, new NDataComparator());

    // Cut cores up to batch preferences.
    // Add some kind of triples check. Can calculate all data again to find out that the second
    // batch is empty.
    ArrayList<NData> triplesCandidatesBatch = new ArrayList<NData>();
    int best = (int) Math.max(
        config.tripleCoreBatchMinimum, triplesCandidates.size() * config.tripleCoreBestPercentage);
    int limit = Math.max(0, triplesCandidates.size() - best);
    if (limit > 0) {
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
    } else {
      triplesCandidatesBatch.addAll(triplesCandidates);
    }

    // Used all candidates in the first batch.
    if (batch == 1 && limit == 0) {
      // Use rmsd variable as substitute for a flag that this candidate should not be calculated
      // again.
      candidate.rmsd = -1.0;
    }

    return triplesCandidatesBatch;
  }

  private void FindTripleCoresAndCalculate(final Dist[][] referenceDistances,
      final Dist[][] targetDistances, final ArrayList<NData> referencePairs,
      final ArrayList<ArrayList<NData>> validPairs, final int maxValidCandidates, final int batch) {
    // Check candidates from the first ones (best dual RMSD).
    // Parallelism can be moved one step up (here)
    for (int index = 0; index < maxValidCandidates; index++) {
      ForkJoinPool threadPool = null;
      final int finalIndex = index;
      try {
        threadPool = new ForkJoinPool(config.threads);

        threadPool
            .submit(() -> referencePairs.stream().parallel().forEach(ndata -> {
              if (ShouldTerminate()) {
                return;
              }

              // validPairs rmsd value represent nucleotide index.
              final ArrayList<NData> candidates = validPairs.get((int) ndata.rmsd);

              if (candidates.size() <= finalIndex) {
                return;
              }

              NData candidate = candidates.get(finalIndex);

              // All possible cores were calculated during first batch.
              if (candidate.rmsd < 0.0) {
                return;
              }

              ArrayList<Coordinates> referenceAtoms = new ArrayList<Coordinates>(Arrays.asList(
                  new Coordinates[targetStructure.get(0).representatives.size() * 3]));
              ArrayList<Coordinates> targetAtoms = new ArrayList<Coordinates>(Arrays.asList(
                  new Coordinates[targetStructure.get(0).representatives.size() * 3]));

              // Filled 2 nucleotides worth of Atoms
              NucleotidesToAtoms(referenceAtoms, ndata.index1, ndata.index2, referenceStructure);
              NucleotidesToAtoms(targetAtoms, candidate.index1, candidate.index2, targetStructure);

              // Find tripleCores already batch trimmed.
              final ArrayList<NData> triplesCoresCandidates = FindTriplesCores(referenceDistances,
                  targetDistances, referenceAtoms, targetAtoms, candidate, ndata, batch);

              // returns null if ShouldTerminate was true while finding cores.
              if (triplesCoresCandidates == null) {
                return;
              }

              // Find final result from found triple cores.
              for (final NData tripleCore : triplesCoresCandidates) {
                NucleotidesToAtoms(referenceAtoms, tripleCore.index1, referenceStructure);
                NucleotidesToAtoms(targetAtoms, tripleCore.index2, targetStructure);
                FindFromTriple(referenceAtoms, targetAtoms, ndata, candidate, tripleCore.index1,
                    tripleCore.index2, tripleCore.superimposer);
              }
            }))
            .get();
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
      final ArrayList<Coordinates> targetAtoms, final NData ndata, final NData candidate,
      int referenceIndex3, int targetIndex3, final Superimposer superimposer) {
    // Here are all triple nucleotide candidates.
    // Prepare for calculations.
    ArrayList<Integer> chainReference = new ArrayList<Integer>(referenceStructure.size());
    ArrayList<Integer> chainTarget = new ArrayList<Integer>(referenceStructure.size());
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
    ArrayList<Nucleotide> nucleotidesReference =
        new ArrayList<Nucleotide>(referenceStructure.size());
    ArrayList<Nucleotide> nucleotidesTarget = new ArrayList<Nucleotide>(referenceStructure.size());
    HashSet<Integer> usedReference = new HashSet<Integer>();
    HashSet<Integer> usedTarget = new HashSet<Integer>();
    for (int i = 0; i < chainReference.size(); i++) {
      usedReference.add(chainReference.get(i));
      usedTarget.add(chainTarget.get(i));
      nucleotidesReference.add(referenceStructure.get(chainReference.get(i)));
      nucleotidesTarget.add(targetStructure.get(chainTarget.get(i)));
    }


    // Precompute table with RMSD improvements.
    ArrayList<NData> precomputedDistances =
        PrecomputeRmsdChanges(targetStructureMoved, usedReference, usedTarget);

    NData bestCandidates = FindBestCandidate(usedReference, usedTarget, precomputedDistances,
        CalculateRmsdSum(chainReference, chainTarget, targetStructureMoved));

    int finalShifted = 2;

    while (bestCandidates != null) {
      chainReference.add(bestCandidates.index1);
      chainTarget.add(bestCandidates.index2);
      usedReference.add(bestCandidates.index1);
      usedTarget.add(bestCandidates.index2);
      nucleotidesReference.add(referenceStructure.get(bestCandidates.index1));
      nucleotidesTarget.add(targetStructure.get(bestCandidates.index2));

      if (chainReference.size() == referenceStructure.size())
        break;

      double rmsdSum = CalculateRmsdSum(chainReference, chainTarget, targetStructureMoved);

      bestCandidates = FindBestCandidate(usedReference, usedTarget, precomputedDistances, rmsdSum);

      // Could not find new candidate. Try shifting structure. Only once. Only if close to current
      // max. Try to shift only if current structure is close to the current best. Do not shift
      // clearly bad alignments.
      if (bestCandidates == null && finalShifted != 0
          && chainReference.size() >= bestChainReference.size() - 3) {
        finalShifted--;

        // Shift and move.
        rmsdSum = CalculateRmsdSum(chainReference, chainTarget, targetStructureMoved);

        final Superimposer superimposer =
            Calculations.FitForRMSD(Nucleotide.NucleotidesToList(nucleotidesReference),
                Nucleotide.NucleotidesToList(nucleotidesTarget));
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

        precomputedDistances =
            PrecomputeRmsdChanges(targetStructureMoved, usedReference, usedTarget);
        bestCandidates =
            FindBestCandidate(usedReference, usedTarget, precomputedDistances, rmsdSum);
      }
    }

    final double currentRmsd = Precision.round(
        Calculations.FitAndCalculateRMSD(Nucleotide.NucleotidesToList(nucleotidesReference),
            Nucleotide.NucleotidesToList(nucleotidesTarget)),
        3);

    // Finished adding nucleotides to the structure.
    try {
      final Specimen spec;
      if (createPopulation) {
        spec = new Specimen(config, referenceStructure, targetStructure, isSequenceDependent);
        spec.initialize(chainReference, chainTarget);
        if (populationSize > populationPool.size()) {
          populationPool.add(spec);
        }
      } else {
        spec = null;
      }
      semaphore.acquire();
      // New alignment clearly better.
      if ((bestChainReference.size() <= chainReference.size())
          || ((bestChainReference.size() == chainReference.size()) && (bestRmsd > currentRmsd))) {
        updateStopTime(bestChainReference.size() < chainReference.size());
        bestRmsd = currentRmsd;
        bestChainReference = chainReference;
        bestChainTarget = chainTarget;
      }
      semaphore.release();
    } catch (InterruptedException e) {
      LOGGER.error(e.getMessage(), e);
    }
  }

  private ArrayList<NData> PrecomputeRmsdChanges(final ArrayList<Nucleotide> targetStructureMoved,
      final HashSet<Integer> usedReference, final HashSet<Integer> usedTarget) {
    ArrayList<NData> rmsdSums =
        new ArrayList<NData>((referenceStructure.size() - usedReference.size())
            * (targetStructureMoved.size() - usedTarget.size()));
    for (int i = 0; i < referenceStructure.size(); i++) {
      if (!usedReference.contains(i)) {
        for (int j = 0; j < targetStructureMoved.size(); j++) {
          if (!usedTarget.contains(j)) {
            rmsdSums.add(new NData(
                i, j, CheckAddition(referenceStructure.get(i), targetStructureMoved.get(j))));
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
      sum += Math.pow(
          Coordinates.getDistance(lhs.representatives.get(i), rhs.representatives.get(i)), 2);
    }
    return sum;
  }

  private NData FindBestCandidate(final HashSet<Integer> usedReference,
      final HashSet<Integer> usedTarget, ArrayList<NData> precomputedDistances,
      double calculatedRmsdSum) {
    int chainSize = usedReference.size();
    int index = -1;
    NData bestCandidates = null;

    for (int i = precomputedDistances.size() - 1; i >= 0; --i) {
      NData best = precomputedDistances.get(i);
      if ((usedReference.contains(best.index1) || usedTarget.contains(best.index2))
          || ((isSequenceDependent && (!StringUtils.equalsIgnoreCase(referenceStructure.get(best.index1).getCode(),
                                                                      targetStructure.get(best.index2).getCode()))))) {
        precomputedDistances.remove(i);
      } else {
        bestCandidates = best;
        index = i; // precomputedDistances.remove(i);
        break;
      }
    }

    if (bestCandidates != null) {
      final int size = chainSize + 1;
      final double rmsd = CalculateFastRmsd(size, calculatedRmsdSum + bestCandidates.rmsd);
      if (rmsd < rmsdLimit) {
        precomputedDistances.remove(index);
        return bestCandidates;
      }
    }
    return null;
  }

  private double CalculateFastRmsd(int chainSize, double calculatedRmsdSum) {
    return Math.sqrt(
        calculatedRmsdSum / (chainSize * referenceStructure.get(0).representatives.size()));
  }

  private double CalculateRmsdSum(ArrayList<Integer> chainReference, ArrayList<Integer> chainTarget,
      ArrayList<Nucleotide> targetStructureMoved) {
    double rmsdSum = 0;
    for (int i = 0; i < chainReference.size(); i++) {
      for (int j = 0; j < referenceStructure.get(0).representatives.size(); j++) {
        rmsdSum +=
            Math.pow(Coordinates.getDistance(
                         referenceStructure.get(chainReference.get(i)).representatives.get(j),
                         targetStructureMoved.get(chainTarget.get(i)).representatives.get(j)),
                2);
      }
    }

    return rmsdSum;
  }
}
