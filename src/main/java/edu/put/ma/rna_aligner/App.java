package edu.put.ma.rna_aligner;

import java.util.ArrayList;


public class App  {
    public static void main(String[] args) {
        if (args.length < 3) {
        	System.out.println("Usage: <method> <reference.pdb> <target.pdb>");
        	System.out.println("Methods: geometric, genetic");
        	System.out.println(".pdb parser is not fully tested!");
        	System.exit(1);
        }
        
        // Parse PDB and create default coarse-grained structures.
    	PdbParser parser = new PdbParser();
    	final ArrayList<Nucleotide> referenceStructure = parser.ParsePdbToCoarseGrained(args[1]);
    	final ArrayList<Nucleotide> targetStructure = parser.ParsePdbToCoarseGrained(args[2]);
		
		AlignerConfig config = new AlignerConfig();

        // Aligner function input: config, coarse-grained RNA structures.
		AlignerOutput output = null;
        if (args[0].equals("geometric")) {

        } else if (args[0].equals("genetic")) {
        	GeneticAligner aligner = new GeneticAligner(config, referenceStructure, targetStructure);
            output = aligner.calculate();
        } else {
        	System.out.println(String.format("Invalid method: %s", args[0]));
        	System.out.println("Methods: geometric, genetic");
        	System.exit(2);
        }
        
        System.out.println(String.format("Mapped %d nucleotides", output.aligned));
        System.out.println("Mapped nucleotides.");
        for (int i = 0; i < output.referenceIndexes.size(); ++i) {
        	System.out.println(String.format("%d\t-> %d", output.referenceIndexes.get(i), output.targetMapping.get(i)));
        }
    }
}
