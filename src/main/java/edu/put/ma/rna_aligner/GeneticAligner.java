package edu.put.ma.rna_aligner;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.concurrent.Semaphore;
import java.util.stream.IntStream;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GeneticAligner {
  private static final Logger LOGGER = LoggerFactory.getLogger(GeneticAligner.class);

  private final AlignerConfig config;
  private final ArrayList<Nucleotide> referenceStructure;
  private final ArrayList<Nucleotide> targetStructure;

  private final int crossChance;
  private final int mutationChance;
  private final int newSpecimenChance;

  private Random rand = new Random();
  private long globalStart = System.currentTimeMillis();
  private long stopTime;

  private Semaphore semaphore = new Semaphore(1);
  private int bestAlignmentSize = -1;
  private double bestAlignmentRMSD = 10000;
  private double bestAlignmentIncorrectlyAlignedResiduesRatio = 10000;
  private Specimen bestSpecimen;
  private boolean terminate = false;
  private int numberOfPopulationsWithoutImprovement = 0;
  private boolean isSequenceDependent;
  private double rmsdLimit;

  GeneticAligner(final AlignerConfig _config, final ArrayList<Nucleotide> _referenceStructure,
      final ArrayList<Nucleotide> _targetStructure, final boolean _isSequenceDependent,
      final double _rmsdLimit) {
    super();
    config = _config;
    referenceStructure = _referenceStructure;
    targetStructure = _targetStructure;
    isSequenceDependent = _isSequenceDependent;
    if ((isSequenceDependent) && (referenceStructure.size() != targetStructure.size())) {
      isSequenceDependent = !isSequenceDependent;
    }
    crossChance = config.crossChance;
    mutationChance = config.mutationChance + crossChance;
    newSpecimenChance = config.newSpecimenChance + mutationChance;
    rmsdLimit = _rmsdLimit;
    stopTime = (long) ((globalStart + (1000 * config.returnTime * 0.25)));
  }

  private ArrayList<Specimen> nextGeneration(ArrayList<Specimen> generation) {
    ArrayList<Specimen> new_generation = new ArrayList<Specimen>();
    for (int i = 0; i < (generation.size() * config.bestPercentage); i++) {
      new_generation.add(generation.get(i));
    }

    // Add specimen.
    int size = new_generation.size();
    for (int i = 0; i < config.populationSize - size; i++) {
      int index = rand.nextInt(generation.size());
      Specimen newSpecimen;

      int chance = rand.nextInt(101);
      if (chance <= crossChance) {
        // Cross
        newSpecimen = (Specimen) generation.get(index).clone();
        int sindex = rand.nextInt(generation.size());
        while (sindex == index) {
          sindex = rand.nextInt(generation.size());
        }
        newSpecimen.cross(generation.get(sindex));
      } else if (chance <= mutationChance) {
        // Mutate
        newSpecimen = (Specimen) generation.get(index).clone();
        newSpecimen.mutate();
      } else if (chance <= newSpecimenChance) {
        // Create new completely random specimen.
        do {
          newSpecimen =
              new Specimen(config, referenceStructure, targetStructure, isSequenceDependent);
          newSpecimen.initialize(rand.nextInt(100) + 1);
        } while (!(newSpecimen.getUsedNucleotidesNumber() > 1
            && Double.compare(newSpecimen.calculateRMSD(), rmsdLimit) <= 0));
      } else {
        // Incorrect chances?
        continue;
      }

      newSpecimen.refinement();

      if ((newSpecimen.getUsedNucleotidesNumber() > 0) && (!new_generation.contains(newSpecimen))) {
        new_generation.add(newSpecimen);
      } else {
        // Repeated specimen. Add one more.
        i--;
      }
    }

    Collections.sort(new_generation, new SpecimenComparator(rmsdLimit));
    return new_generation;
  }

  private void updateStopTime(Boolean resultImprovement) {
    // Made an improvement.
    long now = System.currentTimeMillis();
    long timeLeft = ((long) config.returnTime * 1000) - (now - globalStart);

    long timeImprovement = 0;
    if (resultImprovement) {
      timeImprovement =
          (long) Math.max(timeLeft * config.imprResultPercentage, config.imprResultFlat * 1000.);
    } else {
      timeImprovement =
          (long) Math.max(timeLeft * config.imprRmsdPercentage, config.imprRmsdFlat * 1000.);
    }
    long timeBuffer = Math.max(stopTime, now) + timeImprovement - now;
    if (timeBuffer < config.waitBufferFlat) {
      timeBuffer = config.waitBufferFlat;
    } else {
      timeBuffer = (long) Math.min(timeBuffer, timeLeft * config.waitBufferPercentage);
    }
    stopTime = (long) Math.min(now + timeBuffer, globalStart + (config.returnTime * 1000));
  }

  private ArrayList<Specimen> generatePopulationPool(int populationSize) {
    if (config.geometricPopulation) {
      long st = System.currentTimeMillis();
      AlignerConfig conf_tmp = new AlignerConfig();
      conf_tmp.rmsdLimit = config.rmsdLimit;
      conf_tmp.returnTime = 30;
      conf_tmp.waitBufferPercentage = 0.75;
      conf_tmp.waitBufferFlat = 20;
      conf_tmp.imprResultFlat = 15;
      conf_tmp.imprRmsdFlat = 6;

      conf_tmp.threads = config.threads;
      GeometricAligner aligner = new GeometricAligner(
          conf_tmp, referenceStructure, targetStructure, isSequenceDependent, rmsdLimit);
      ArrayList<Specimen> res = aligner.createPopulation(populationSize);
      // stopTime = (long) ((System.currentTimeMillis() + (1000 * config.returnTime * 0.2)));

      updateStopTime(true);
      updateStopTime(true);
      updateStopTime(true);
      return res;
    }
    return null;
  }

  public AlignerOutput calculate() {
    IntStream threads = IntStream.range(0, config.threads);
    numberOfPopulationsWithoutImprovement = 0;
    final ArrayList<Specimen> populationPool =
        generatePopulationPool((int) (Math.ceil(config.threads * 0.50)) * config.populationSize);

    threads.parallel().forEach(index -> {
      // No more than X seconds.

      while (!terminate && (System.currentTimeMillis() < stopTime)) {
        ArrayList<Specimen> population = new ArrayList<Specimen>();
        // Create initial population
        if (config.geometricPopulation) {
          for (int i = index * config.populationSize; i < (index + 1) * config.populationSize;
               ++i) {
            Specimen spec;
            if (populationPool.size() > i) {
              spec = (Specimen) populationPool.get(i).clone();
            } else {
              // Geometric did not produce enough specimens. Select random one.
              spec = new Specimen(config, referenceStructure, targetStructure, isSequenceDependent);
              spec.initialize(rand.nextInt(85) + 5);
              spec.refinement();
              int count = 0;
              while (spec.getUsedNucleotidesNumber() <= 1 && count < 20) {
                spec =
                    new Specimen(config, referenceStructure, targetStructure, isSequenceDependent);
                spec.initialize(rand.nextInt(85) + 5);
                spec.refinement();
              }
            }
            if (!population.contains(spec)) {
              spec.calculateRMSD();
              population.add(spec);
            }
          }
        } else {
          int count = 0;
          while (population.size() < config.populationSize) {
            if (count > 10) {
              break;
            }
            final Specimen spec =
                new Specimen(config, referenceStructure, targetStructure, isSequenceDependent);
            spec.initialize(rand.nextInt(85) + 5);
            spec.refinement();
            if (spec.getUsedNucleotidesNumber() > 1) {
              if (!population.contains(spec)) {
                spec.calculateRMSD();
                population.add(spec);
                count = 0;
              } else {
                count++;
              }
            }
          }
        }

        Collections.sort(population, new SpecimenComparator(rmsdLimit));

        int best_size = population.get(0).getUsedNucleotidesNumber();
        double best_rmsd = population.get(0).calculateRMSD();
        double best_incorrectlyAlignedResiduesRatio =
            population.get(0).getIncorrectlyAlignedResiduesRatio();
        Specimen best_specimen = (Specimen) population.get(0).clone();
        if (config.geometricPopulation) {
          try {
            semaphore.acquire();
            if ((best_rmsd <= rmsdLimit && best_size > bestAlignmentSize)
                || ((Double.compare(best_size, bestAlignmentSize) == 0)
                    && (best_rmsd < bestAlignmentRMSD))) {
              bestAlignmentSize = best_size;
              bestAlignmentRMSD = best_rmsd;
              bestAlignmentIncorrectlyAlignedResiduesRatio = best_incorrectlyAlignedResiduesRatio;
              bestSpecimen = (Specimen) best_specimen.clone();
            }
            if (((bestAlignmentSize == Math.min(referenceStructure.size(), targetStructure.size()))
                    /*&& (Double.compare(bestAlignmentRMSD,rmsdLimit) <= 0)*/)) {
              terminate = true;
              break;
            }
          } catch (InterruptedException e) {
            LOGGER.error(e.getMessage(), e);
          } finally {
            semaphore.release();
          }
        }

        if (!((best_size == Math.min(referenceStructure.size(), targetStructure.size()))
                /*&& (Double.compare(best_rmsd,rmsdLimit) <= 0)*/)) {
          long last_improvement = System.currentTimeMillis();

          // Run while no improvement for config.resetThreadTime
          // seconds. Then restart
          while ((System.currentTimeMillis() - last_improvement < 1000 * config.resetThreadTime
                     || config.geometricPopulation)
              && System.currentTimeMillis() < stopTime) {
            population = nextGeneration(population);

            final int currentSize = population.get(0).getUsedNucleotidesNumber();
            final double currentRmsd = population.get(0).calculateRMSD();
            final double incorrectlyAlignedResiduesRatio =
                population.get(0).getIncorrectlyAlignedResiduesRatio();
            if (((best_incorrectlyAlignedResiduesRatio > incorrectlyAlignedResiduesRatio)
                    && (Double.compare(currentRmsd, rmsdLimit) <= 0))
                || ((Double.compare(
                         best_incorrectlyAlignedResiduesRatio, incorrectlyAlignedResiduesRatio)
                        == 0)
                    && (best_size < currentSize) && (Double.compare(currentRmsd, rmsdLimit) <= 0))
                || ((Double.compare(
                         best_incorrectlyAlignedResiduesRatio, incorrectlyAlignedResiduesRatio)
                        == 0)
                    && (best_size == currentSize) && (best_rmsd > currentRmsd))) {
              last_improvement = System.currentTimeMillis();

              best_size = currentSize;
              best_rmsd = currentRmsd;
              best_incorrectlyAlignedResiduesRatio = incorrectlyAlignedResiduesRatio;
              best_specimen = (Specimen) population.get(0).clone();
              updatePopulationsNumber(true);
              // Is it an improvement to the global result?
              if (bestAlignmentSize < best_size
                  || (best_rmsd < bestAlignmentRMSD && bestAlignmentSize == best_size)) {
                // Update global best result.
                try {
                  semaphore.acquire();
                  if ((best_rmsd <= rmsdLimit && best_size > bestAlignmentSize)
                      || ((Double.compare(best_size, bestAlignmentSize) == 0)
                          && (best_rmsd < bestAlignmentRMSD))) {
                    updateStopTime(bestAlignmentSize < best_size);
                    bestAlignmentSize = best_size;
                    bestAlignmentRMSD = best_rmsd;
                    bestAlignmentIncorrectlyAlignedResiduesRatio =
                        best_incorrectlyAlignedResiduesRatio;
                    bestSpecimen = (Specimen) best_specimen.clone();
                  }
                  if (((bestAlignmentSize
                          == Math.min(referenceStructure.size(), targetStructure.size()))
                          /*&& (Double.compare(bestAlignmentRMSD,rmsdLimit) <= 0)*/)) {
                    terminate = true;
                    break;
                  }
                } catch (InterruptedException e) {
                  LOGGER.error(e.getMessage(), e);
                } finally {
                  semaphore.release();
                }
              }

            } else {
              updatePopulationsNumber(false);
            }

            if (((best_size == Math.min(referenceStructure.size(), targetStructure.size()))
                    /*&& (Double.compare(best_rmsd,rmsdLimit) <= 0)*/)) {
              terminate = true;
              break;
            }
          }
        }

        try {
          semaphore.acquire();
          if ((best_rmsd <= rmsdLimit && best_size > bestAlignmentSize)
              || ((Double.compare(best_size, bestAlignmentSize) == 0)
                  && (best_rmsd < bestAlignmentRMSD))) {
            bestAlignmentSize = best_size;
            bestAlignmentRMSD = best_rmsd;
            bestAlignmentIncorrectlyAlignedResiduesRatio = best_incorrectlyAlignedResiduesRatio;
            bestSpecimen = (Specimen) best_specimen.clone();
          }
          if (((bestAlignmentSize == Math.min(referenceStructure.size(), targetStructure.size()))
                  && (Double.compare(bestAlignmentRMSD, rmsdLimit) <= 0))
              || (numberOfPopulationsWithoutImprovement > 300)) {
            terminate = true;
          }
        } catch (InterruptedException e) {
          LOGGER.error(e.getMessage(), e);
        } finally {
          semaphore.release();
        }
      }
    });

    final ArrayList<Integer> referenceIndexes = new ArrayList<Integer>();
    final ArrayList<Integer> targetMapping = new ArrayList<Integer>();

    for (int i = 0; i < bestSpecimen.primaryNucleotidesUsed.length; ++i) {
      referenceIndexes.add(i);
      if (bestSpecimen.primaryNucleotidesUsed[i] == 1) {
        targetMapping.add(bestSpecimen.secondaryNucleotidesMap[i]);
      } else {
        targetMapping.add(-1);
      }
    }

    return new AlignerOutput(bestSpecimen.getUsedNucleotidesNumber(), referenceIndexes,
        targetMapping,
        Calculations.FitForRMSD(Nucleotide.NucleotidesToListMapped(bestSpecimen.primaryNucleotides,
                                    bestSpecimen.primaryNucleotidesUsed, true),
            Nucleotide.NucleotidesToListMapped(
                bestSpecimen.secondaryNucleotides, bestSpecimen.secondaryNucleotidesMap, false)),
        System.currentTimeMillis() - globalStart, bestAlignmentRMSD);
  }

  private final void updatePopulationsNumber(final boolean init) {
    try {
      semaphore.acquire();
      if (init) {
        numberOfPopulationsWithoutImprovement = 0;
      } else {
        numberOfPopulationsWithoutImprovement++;
      }
    } catch (InterruptedException e) {
      LOGGER.error(e.getMessage(), e);
    } finally {
      semaphore.release();
    }
  }
}
