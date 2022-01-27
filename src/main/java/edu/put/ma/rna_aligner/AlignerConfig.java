package edu.put.ma.rna_aligner;

public class AlignerConfig {
  // General config.
  public double pairRmsdLimit = 0.65;
  public double tripleRmsdLimit = 1.0;
  public double rmsdLimit = 3.5;
  // Number of threads used for calculations. 
  public int threads = 12;
  // Break calculations after (300) seconds and return AlignerOutput.
  public double returnTime = 300;

  // Time buffer that result improvements can accumulate.
  public double waitBufferPercentage = 0.1;
  public long waitBufferFlat = 15;

  // How much buffer time each result improvement provides.
  public double imprResultPercentage = 0.04;
  public long imprResultFlat = 10;
  public double imprRmsdPercentage = 0.01;
  public long imprRmsdFlat = 1;

  public boolean geometricPopulation = true;

  // Geometric
  // How many batches of 2 nucleotide cores should be made;
  // if rmsd limit is 4 and split is 4 batches (25%) - batch up to 1 RMSD, 2, 3 and 4;
  // Small number of batches will require more RAM but will recalculate cores fewer number of times.
  // Batch = 1 will store all possible dual core in the memory. It can be very RAM demanding for larger structures.
  public int dualCoreBatches = 2;

  // tripleCoreBestPercentage how many (percentage wise) best tripleCores should be calculated in the first batch (out of 2).
  // BatchMaximum hard maximum limit of triple cores for candidate to consider in the first batch.
  // Introduced to help with batches that can have huge small number of candidates and should be calculated at once.
  public int tripleCoreBatchMinimum = 25;
  public double tripleCoreBestPercentage = 0.25;

  // Genetic
  // How long can the calculations go without result improvement.
  // After time period, calculations are reset and thread starts again.
  public double resetThreadTime = 60;
  // How many of the best specimen are kept in the next generation.
  public double bestPercentage = 0.15;
  // Generation size / population;
  public int populationSize = 125;
  // Chances in whole percents.
  public int crossChance = 25;
  public int mutationChance = 74;
  public int newSpecimenChance = 1;
  // Multiple mutation chances.
  public int singleMutation = 65;
  public int doubleMutation = 15;
  public int tripleMutation = 10;
  public int quadrupleMutation = 10;

  AlignerConfig () {
    super();
  }
}
