package edu.put.ma.rna_aligner;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.Semaphore;
import java.util.stream.IntStream;

public class GeneticAligner {
	private final AlignerConfig config;
	private final ArrayList<Nucleotide> referenceStructure;
	private final ArrayList<Nucleotide> targetStructure;
	
	private final int crossChance;
	private final int mutationChance;
	private final int newSpecimenChance;
	
	Random rand = new Random();
	long globalStart = System.currentTimeMillis();
	
	GeneticAligner(final AlignerConfig _config,
			final ArrayList<Nucleotide> _referenceStructure,
			final ArrayList<Nucleotide> _targetStructure) {
		super();
		config = _config;
		referenceStructure = _referenceStructure;
		targetStructure = _targetStructure;
		crossChance = config.crossChance;
		mutationChance = config.mutationChance + crossChance;
		newSpecimenChance = config.newSpecimenChance + mutationChance;
	}
	
	
	private ArrayList<Specimen> nextGeneration(ArrayList<Specimen> generation) {
		HashMap<Integer, Boolean> specimen_hashes = new HashMap<Integer, Boolean>();

		ArrayList<Specimen> new_generation = new ArrayList<Specimen>();
		for (int i = 0; i < (generation.size() * config.bestPercentage); i++) {
			new_generation.add(generation.get(i));
			specimen_hashes.put(generation.get(i).getHash(), true);
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
				index = rand.nextInt(generation.size());
				newSpecimen.cross(generation.get(index));
			} else if (chance <= mutationChance) {
				// Mutate
				newSpecimen = (Specimen) generation.get(index).clone();
				newSpecimen.mutate();
			} else if (chance <= newSpecimenChance){
				// Create new completely random specimen.
				newSpecimen = new Specimen(config, referenceStructure, targetStructure);
				newSpecimen.initialize(rand.nextInt(100) + 1);
			} else {
				// Incorrect chances?
				continue;
			}

			if (!specimen_hashes.containsKey(newSpecimen.getHash())) {
				new_generation.add(newSpecimen);
			} else {
				// Repeated specimen. Add one more.
				i--;
			}
		}

		try {
			Collections.sort(new_generation, new SpecimenComparator(config.rmsdLimit));
		} catch (Exception e) {
		}
		return new_generation;
	}
	
	
	Semaphore semaphore = new Semaphore(1);
	int bestAlignmentSize = -1;
	double bestAlignmentRMSD = 10000;
	Specimen bestSpecimen;
	
	public AlignerOutput calculate() {
		IntStream threads = IntStream.range(0, config.threads);
		threads.parallel().forEach(index -> {
			// No more than X seconds.
			while (System.currentTimeMillis() - globalStart < 1000 * config.returnTime) {
				ArrayList<Specimen> population = new ArrayList<Specimen>();
				
				// Create initial population
				for (int i = 0; i < config.populationSize; i++) {
					Specimen spec = new Specimen(config, referenceStructure, targetStructure);
					spec.initialize(rand.nextInt(85) + 5);
					spec.calculateRMSD();
					population.add(spec);
				}
				
				int best_size = 0;
				double best_rmsd = 1000;
				Specimen best_specimen = new Specimen(config, referenceStructure, targetStructure);;
				
				long last_improvement = System.currentTimeMillis();
				
				// Run while no improvement for config.resetThreadTime seconds. Then restart
				while (System.currentTimeMillis() - last_improvement < 1000 * config.resetThreadTime &&
					System.currentTimeMillis() - globalStart < 1000 * config.returnTime) {
				    population = nextGeneration(population);
					long endTime = System.currentTimeMillis();
					
					if (best_size <= population.get(0).getUsedNucleotidesNumber()) {
						last_improvement = endTime;
						best_size = population.get(0).getUsedNucleotidesNumber();
						best_rmsd = population.get(0).calculateRMSD();
						best_specimen = (Specimen) population.get(0).clone();
					} else if (best_size == population.get(0).getUsedNucleotidesNumber()) {
						if (best_rmsd < population.get(0).calculateRMSD()) {
							last_improvement = endTime;
							best_size = population.get(0).getUsedNucleotidesNumber();
							best_rmsd = population.get(0).calculateRMSD();
							best_specimen = (Specimen) population.get(0).clone();
						}
					}
				}
				
				try {
					semaphore.acquire();
					if (best_size > bestAlignmentSize) {
						bestAlignmentSize = best_size;
						bestAlignmentRMSD = best_rmsd;
						bestSpecimen = (Specimen) best_specimen.clone();
					} else if (best_size == bestAlignmentSize) {
						if (best_rmsd < bestAlignmentRMSD) {
							bestAlignmentSize = best_size;
							bestAlignmentRMSD = best_rmsd;
							bestSpecimen = (Specimen) best_specimen.clone();
						}
					}
					semaphore.release();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

			}
		});
		
		ArrayList<Integer> referenceIndexes = new ArrayList<Integer>();
		ArrayList<Integer> targetMapping = new ArrayList<Integer>();
		for (int i = 0; i < bestSpecimen.primaryNucleotidesUsed.length; ++i) {
			if (bestSpecimen.primaryNucleotidesUsed[i] == 1) {
				referenceIndexes.add(i + 1);
				targetMapping.add(bestSpecimen.secondaryNucleotidesMap[i] + 1);
			}
		}
		
		return new AlignerOutput(bestSpecimen.getUsedNucleotidesNumber(), referenceIndexes, targetMapping,
				Calculations.FitForRMSD(Nucleotide.NucleotidesToMappedTable(bestSpecimen.primaryNucleotides, bestSpecimen.primaryNucleotidesUsed, true),
						Nucleotide.NucleotidesToMappedTable(bestSpecimen.secondaryNucleotides, bestSpecimen.secondaryNucleotidesMap, false)));
	}

}
