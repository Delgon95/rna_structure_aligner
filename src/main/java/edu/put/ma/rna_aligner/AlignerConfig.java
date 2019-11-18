package edu.put.ma.rna_aligner;

public class AlignerConfig {

	// General config.
	public double rmsdLimit = 3.5;
	// Number of threads used for calculations. 
	public int threads = 8;
	
	// Genetic
	// How long can the calculations go without result improvement.
	// After time period, calculations are reset and thread starts again.
	public double resetThreadTime = 60;
	// Break calculations after (300) seconds and return AlignerOutput.
	public double returnTime = 300;
	// How many of the best specimen are kept in the next generation.
	public double bestPercentage = 0.15;
	// Generation size / population;
	public int populationSize = 200;
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
