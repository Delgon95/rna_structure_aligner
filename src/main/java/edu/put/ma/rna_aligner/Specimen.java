package edu.put.ma.rna_aligner;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class Specimen {
  public final ArrayList<Nucleotide> primaryNucleotides;
  public final ArrayList<Nucleotide> secondaryNucleotides;

  public int[] primaryNucleotidesUsed;   // 0 if not used, 1 if used.
  public int[] secondaryNucleotidesMap;  // number of secondary nucleotide that is associated with the same index primary nucleotide.
  public List<Integer> availableNucleotides = new ArrayList<Integer>(); 
  public HashMap<Integer, Integer> secondToFirst = new HashMap<Integer, Integer>();

  private Random rand = new Random();
  private boolean changed = true;
  private double rmsd = Double.MAX_VALUE;
  private boolean hashed = false;
  private int hash = 0;

  private final AlignerConfig config;
  private int singleMutation = 65;
  private int doubleMutation = 15;
  private int tripleMutation = 10;
  private int quadrupleMutation = 10;

  Specimen(final AlignerConfig _config,
      final ArrayList<Nucleotide> _primaryNucleotides,
      final ArrayList<Nucleotide> _secondaryNucleotides) {
    super();
    primaryNucleotides = _primaryNucleotides;
    secondaryNucleotides = _secondaryNucleotides;
    primaryNucleotidesUsed = new int[primaryNucleotides.size()];
    secondaryNucleotidesMap = new int[primaryNucleotides.size()];
    config = _config;
    singleMutation = config.singleMutation;
    doubleMutation = config.singleMutation + config.doubleMutation;
    tripleMutation = config.singleMutation + config.doubleMutation + config.tripleMutation;
    quadrupleMutation = config.singleMutation + config.doubleMutation + config.tripleMutation + config.quadrupleMutation;
  }

  @Override
  public Object clone() {
    Specimen copy = new Specimen(config, primaryNucleotides, secondaryNucleotides);
    for (int i = 0; i < primaryNucleotidesUsed.length; i++) {
      copy.primaryNucleotidesUsed[i] = this.primaryNucleotidesUsed[i];
    }
    for (int i = 0; i < secondaryNucleotidesMap.length; i++) {
      copy.secondaryNucleotidesMap[i] = this.secondaryNucleotidesMap[i];
    }
    for (int val : this.availableNucleotides) {
      copy.availableNucleotides.add(val);
    }
    for (Map.Entry<Integer, Integer> keyVal : this.secondToFirst.entrySet()) {
      copy.secondToFirst.put(keyVal.getKey(), keyVal.getValue());
    }

    copy.changed = this.changed;
    copy.rmsd = this.rmsd;

    return copy;
  }

  public void initialize(int percentage) {
    // Prepare available Residues
    for (int i = 0; i < secondaryNucleotides.size(); i++) {
      availableNucleotides.add(i);
    }
    for (int i = 0; i < primaryNucleotidesUsed.length; i++) {
      // 10-90 chance to use primary residue.
      if (rand.nextInt(100) < percentage) {
        primaryNucleotidesUsed[i] = 1;
        int selected = getRandomAvailable();
        secondaryNucleotidesMap[i] = selected;
        secondToFirst.put(selected, i);
      } else {
        primaryNucleotidesUsed[i] = -1;
        secondaryNucleotidesMap[i] = -1;
      }
    }
  }

  public int getHash() {
    if (!hashed) {
      hashed = true;
      hash = 7;
      for (int i = 0; i < primaryNucleotidesUsed.length; i++) {
        hash += (i * primaryNucleotidesUsed[i]) * i;
      }
      hash %= 95482123;
      for (int i = 0; i < secondaryNucleotidesMap.length; i++) {
        hash += (i * primaryNucleotidesUsed[i]);
      }
      hash %= 95482123;
    }

    return hash;
  }

  private int getRandomAvailable() {
    int index = rand.nextInt(availableNucleotides.size());
    int value = availableNucleotides.get(index);
    availableNucleotides.remove(index);
    return value;
  }

  public void printTables() {
    for (int i = 0; i < primaryNucleotidesUsed.length; i++) {
      System.out.print(String.format("(%d %d)", primaryNucleotidesUsed[i], secondaryNucleotidesMap[i]));
    }
    System.out.println("");
  }

  public double calculateRMSD() {
    if (changed) {
      changed = false;
      if (availableNucleotides.size() == primaryNucleotides.size()) {
        return Double.POSITIVE_INFINITY - 1;
      }
      rmsd = Calculations.FitAndCalculateRMSD(Nucleotide.NucleotidesToListMapped(primaryNucleotides, primaryNucleotidesUsed, true),
          Nucleotide.NucleotidesToListMapped(secondaryNucleotides, secondaryNucleotidesMap, false));
    }
    return rmsd;
  }

  public int getRemainingNucleotidesNumber() {
    return availableNucleotides.size();
  }

  public int getUsedNucleotidesNumber() {
    return primaryNucleotides.size() - availableNucleotides.size();
  }

  private void mutate(int variant) {
    switch (variant) {
    // Flip random
    case 0:
      int index = rand.nextInt(primaryNucleotidesUsed.length);
      if (primaryNucleotidesUsed[index] == 1) {
        primaryNucleotidesUsed[index] = -1;
        availableNucleotides.add(secondaryNucleotidesMap[index]);
        secondToFirst.remove(secondaryNucleotidesMap[index]);
        secondaryNucleotidesMap[index] = -1;
      } else {
        primaryNucleotidesUsed[index] = 1;
        int selected = getRandomAvailable();
        secondaryNucleotidesMap[index] = selected;
        secondToFirst.put(selected, index);
      }
      break;
      // Add Swap 2 random.
    case 1:
      int index1 = rand.nextInt(primaryNucleotidesUsed.length);
      int index2 = rand.nextInt(primaryNucleotidesUsed.length);
      while (index2 == index1) {
        index2 = rand.nextInt(primaryNucleotidesUsed.length);
      }
      int temp_primary = primaryNucleotidesUsed[index1];
      primaryNucleotidesUsed[index1] = primaryNucleotidesUsed[index2];
      primaryNucleotidesUsed[index2] = temp_primary;

      secondToFirst.replace(secondaryNucleotidesMap[index1], index2);
      secondToFirst.replace(secondaryNucleotidesMap[index2], index1);

      temp_primary = secondaryNucleotidesMap[index1];
      secondaryNucleotidesMap[index1] = secondaryNucleotidesMap[index2];
      secondaryNucleotidesMap[index2] = temp_primary;
      break;
    case 2:
      int index3 = rand.nextInt(primaryNucleotidesUsed.length);
      int tries = 0;
      while (primaryNucleotidesUsed[index3] != 1) {
        ++tries;
        if (tries == 5) {
          return;
        }
        index3 = rand.nextInt(primaryNucleotidesUsed.length);
      }
      availableNucleotides.add(secondaryNucleotidesMap[index3]);
      secondToFirst.remove(secondaryNucleotidesMap[index3]);

      int selected = getRandomAvailable();
      secondaryNucleotidesMap[index3] = selected;
      secondToFirst.put(selected, index3);
      break;
    }
  }

  public void mutate() {
    changed = true;
    hashed = false;

    int value = rand.nextInt(100);
    if (value < singleMutation) {
      mutate(rand.nextInt(3));
    } else if (value < doubleMutation) {
      mutate(rand.nextInt(3));
      mutate(rand.nextInt(3));
    } else if (value < tripleMutation) {
      mutate(rand.nextInt(3));
      mutate(rand.nextInt(3));
      mutate(rand.nextInt(3));
    } else if (value < quadrupleMutation){
      mutate(rand.nextInt(3));
      mutate(rand.nextInt(3));
      mutate(rand.nextInt(3));
      mutate(rand.nextInt(3));
    }
  }

  private void cross(final Specimen second, int from, int to) {
    // Take parameters in range from second.
    for (int i = from; i < to; i++) {
      // Take Primary structure.
      this.primaryNucleotidesUsed[i] = second.primaryNucleotidesUsed[i];
      // Remove now not used secondary. Update hashmap of values.
      if (this.secondaryNucleotidesMap[i] >= 0) {
        this.availableNucleotides.add(this.secondaryNucleotidesMap[i]);
        this.secondToFirst.remove(this.secondaryNucleotidesMap[i]);
        this.secondaryNucleotidesMap[i] = -1;
      }
    }

    // Fix empty part of in the range with values from second Specimen
    for (int i = from; i < to; i++) {
      if (this.primaryNucleotidesUsed[i] == 1) {
        int valueToMap = second.secondaryNucleotidesMap[i];
        if (!this.secondToFirst.containsKey(valueToMap)) {
          this.secondaryNucleotidesMap[i] = valueToMap;
          this.secondToFirst.put(valueToMap, i);
          // TODO FIX SLOW
          // find and remove from available.
          for (int j = 0; j < this.availableNucleotides.size(); j++) {
            if (this.availableNucleotides.get(j) == valueToMap) {
              this.availableNucleotides.remove(j);
              break;
            }
          }
        } else {
          // Value inside not added part
          int chance = rand.nextInt(100);
          if (chance < 45) {
            // Do not change outside part.
            this.secondaryNucleotidesMap[i] = -1;
            this.primaryNucleotidesUsed[i] = -1;
          } else if (chance < 101) {
            // Force change in original structure to take this one.
            this.secondaryNucleotidesMap[i] = valueToMap;
            int index = this.secondToFirst.get(valueToMap);
            // Reset outside part. Do not add to available as it is used above.
            this.primaryNucleotidesUsed[index] = -1;
            this.secondaryNucleotidesMap[index] = -1;
          } else {
            // Add random available 
            // TODO
          }
        }
      }
    }
  }

  public void cross(final Specimen second) {
    changed = true;
    hashed = false;

    int value = rand.nextInt(100);
    if (value < 85) {
      //Cross;
      int from = rand.nextInt(primaryNucleotidesUsed.length - 1);
      int to = from + rand.nextInt(primaryNucleotidesUsed.length - from);
      while (from == to) {
        to = from + rand.nextInt(primaryNucleotidesUsed.length - from);
      }
      cross(second, from, to);
    } else {
      // Take at least one, but do not go out of range.
      int range = rand.nextInt(primaryNucleotidesUsed.length - 1) + 1;
      cross(second, range, primaryNucleotidesUsed.length);
    }
  }
}