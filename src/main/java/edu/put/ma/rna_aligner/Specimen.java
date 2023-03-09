package edu.put.ma.rna_aligner;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.util.Precision;

public class Specimen {
  public final ArrayList<Nucleotide> primaryNucleotides;
  public final ArrayList<Nucleotide> secondaryNucleotides;

  public int min_size;
  public int[] primaryNucleotidesUsed; // 0 if not used, 1 if used.
  public int[] secondaryNucleotidesMap; // number of secondary nucleotide that is associated with
                                        // the same index primary nucleotide.
  public List<Integer> availableNucleotides = new ArrayList<Integer>();
  public List<Integer> promisingNucleotides = new ArrayList<Integer>();
  public HashMap<Integer, Integer> secondToFirst = new HashMap<Integer, Integer>();

  private Random rand = new Random();
  private boolean changed = true;
  private double rmsd = Double.MAX_VALUE;
  private double incorrectlyAlignedResiduesRatio = 0;

  private final AlignerConfig config;
  private int singleMutation = 65;
  private int doubleMutation = 15;
  private int tripleMutation = 10;
  private int quadrupleMutation = 10;
  private boolean isSequenceDependent;
  private boolean respectOrder;

  Specimen(final AlignerConfig _config, final ArrayList<Nucleotide> _primaryNucleotides,
      final ArrayList<Nucleotide> _secondaryNucleotides, final boolean _isSequenceDependent,
      final boolean _respectOrder) {
    super();
    if (_primaryNucleotides.size() > _secondaryNucleotides.size()) {
      secondaryNucleotides = _primaryNucleotides;
      primaryNucleotides = _secondaryNucleotides;
    } else {
      primaryNucleotides = _primaryNucleotides;
      secondaryNucleotides = _secondaryNucleotides;
    }
    isSequenceDependent = _isSequenceDependent;
    respectOrder = _respectOrder;
    primaryNucleotidesUsed = new int[primaryNucleotides.size()];
    secondaryNucleotidesMap = new int[primaryNucleotides.size()];
    
    config = _config;
    singleMutation = config.singleMutation;
    doubleMutation = config.singleMutation + config.doubleMutation;
    tripleMutation = config.singleMutation + config.doubleMutation + config.tripleMutation;
    quadrupleMutation = config.singleMutation + config.doubleMutation + config.tripleMutation
        + config.quadrupleMutation;
  }

  @Override
  public Object clone() {
    Specimen copy =
        new Specimen(config, primaryNucleotides, secondaryNucleotides, isSequenceDependent, respectOrder);
    for (int i = 0; i < primaryNucleotidesUsed.length; i++) {
      copy.primaryNucleotidesUsed[i] = this.primaryNucleotidesUsed[i];
    }
    for (int i = 0; i < secondaryNucleotidesMap.length; i++) {
      copy.secondaryNucleotidesMap[i] = this.secondaryNucleotidesMap[i];
    }
    for (int val : this.availableNucleotides) {
      copy.availableNucleotides.add(val);
    }
    for (int val : this.promisingNucleotides) {
      copy.promisingNucleotides.add(val);
    }
    for (Map.Entry<Integer, Integer> keyVal : this.secondToFirst.entrySet()) {
      copy.secondToFirst.put(keyVal.getKey(), keyVal.getValue());
    }

    copy.changed = this.changed;
    copy.rmsd = this.rmsd;

    return copy;
  }

  public void initialize(ArrayList<Integer> chainReference, ArrayList<Integer> chainTarget) {
    // Prepare available Residues
    for (int i = 0; i < secondaryNucleotides.size(); i++) {
      availableNucleotides.add(i);
    }
    for (int i = 0; i < primaryNucleotidesUsed.length; i++) {
      primaryNucleotidesUsed[i] = -1;
      secondaryNucleotidesMap[i] = -1;
    }

    for (int i = 0; i < chainReference.size(); i++) {
      int idxRef = chainReference.get(i);
      int idxTar = chainTarget.get(i);
      primaryNucleotidesUsed[idxRef] = 1;
      int idx = availableNucleotides.indexOf(idxRef);
      // Should always be >= 0.
      if (idx >= 0) {
        availableNucleotides.remove(idx);
      } else {
        System.err.println("idx not in available nucleotidex");
        System.err.println(idxRef);
        System.err.println(idxTar);
        System.err.println(idx);
      }
      secondaryNucleotidesMap[idxRef] = idxTar;
      secondToFirst.put(idxTar, idxRef);
    }
    calculateRMSD();
  }

  public void initialize(int percentage) {
    // Prepare available Residues
    for (int i = 0; i < secondaryNucleotides.size(); i++) {
      availableNucleotides.add(i);
    }
    int selected = -1;
    for (int i = 0; i < primaryNucleotidesUsed.length; i++) {
      // 10-90 chance to use primary residue.
      if (rand.nextInt(100) < percentage) {
        primaryNucleotidesUsed[i] = 1;
        if (selected == -1) {
          selected = getRandomAvailable();
        } else {
          final Integer index = Integer.valueOf(++selected);
          if (availableNucleotides.contains(index)) {
            availableNucleotides.remove(index);
          } else {
            selected = getRandomAvailable();
          }
        }
        secondaryNucleotidesMap[i] = selected;
        if (selected >= 0) {
          secondToFirst.put(selected, i);
        } else {
          primaryNucleotidesUsed[i] = -1;
        }
      } else {
        selected = -1;
        primaryNucleotidesUsed[i] = -1;
        secondaryNucleotidesMap[i] = -1;
      }
    }
    calculateRMSD();
  }

  private int getRandomAvailable() {
    if (promisingNucleotides.size() > 0 && rand.nextInt(100) > 75) {
      final int index = rand.nextInt(promisingNucleotides.size());
      final int value = promisingNucleotides.get(index);
      promisingNucleotides.remove(index);
      availableNucleotides.remove(Integer.valueOf(value));
      return value;
    } else if (availableNucleotides.size() > 0) {
      final int index = rand.nextInt(availableNucleotides.size());
      final int value = availableNucleotides.get(index);
      availableNucleotides.remove(index);
      final Integer leftNeighbour = Integer.valueOf(value - 1);
      if ((!promisingNucleotides.contains(leftNeighbour))
          && (availableNucleotides.contains(leftNeighbour)) && rand.nextInt(100) > 50) {
        promisingNucleotides.add(leftNeighbour);
      }
      final Integer rightNeighbour = Integer.valueOf(value + 1);
      if ((!promisingNucleotides.contains(rightNeighbour))
          && (availableNucleotides.contains(rightNeighbour)) && rand.nextInt(100) > 50) {
        promisingNucleotides.add(rightNeighbour);
      }
      return value;
    }
    return -1;
    /*
    if (availableNucleotides.size() > 0) {
      int index = rand.nextInt(availableNucleotides.size());
      int value = availableNucleotides.get(index);
      availableNucleotides.remove(index);
      return value;
    } else {
      return -1;
    }
    */
    /*
    if (promisingNucleotides.size() > 0) {
      final int index = rand.nextInt(promisingNucleotides.size());
      final int value = promisingNucleotides.get(index);
      promisingNucleotides.remove(index);
      availableNucleotides.remove(Integer.valueOf(value));
      return value;
    } else if (availableNucleotides.size() > 0) {
      final int index = rand.nextInt(availableNucleotides.size());
      final int value = availableNucleotides.get(index);
      availableNucleotides.remove(index);
      final Integer leftNeighbour = Integer.valueOf(value - 1);
      if ((!promisingNucleotides.contains(leftNeighbour))
          && (availableNucleotides.contains(leftNeighbour))) {
        promisingNucleotides.add(leftNeighbour);
      }
      final Integer rightNeighbour = Integer.valueOf(value + 1);
      if ((!promisingNucleotides.contains(rightNeighbour))
          && (availableNucleotides.contains(rightNeighbour))) {
        promisingNucleotides.add(rightNeighbour);
      }
      return value;
    }
    return -1;
    */
  }

  public double calculateRMSD() {
    if (changed) {
      changed = false;
      if (availableNucleotides.size() == secondaryNucleotides.size()) {
        return Double.POSITIVE_INFINITY - 1;
      }
      rmsd = Precision.round(
          Calculations.FitAndCalculateRMSD(
              Nucleotide.NucleotidesToListMapped(primaryNucleotides, primaryNucleotidesUsed, true),
              Nucleotide.NucleotidesToListMapped(
                  secondaryNucleotides, secondaryNucleotidesMap, false)),
          3);
      computeIncorrectlyAlignedResiduesRatio();
    }
    return rmsd;
  }

  public void computeIncorrectlyAlignedResiduesRatio() {
    int incorrectlyAlignedResidues = 0;
    if (isSequenceDependent) {
      for (int i = 0; i < primaryNucleotidesUsed.length; i++) {
        if (primaryNucleotidesUsed[i] == 1) {
          final int targetIndex = secondaryNucleotidesMap[i];
          if ((targetIndex > -1) &&
              (!StringUtils.equalsIgnoreCase(primaryNucleotides.get(i).getCode(),
             secondaryNucleotides.get(targetIndex).getCode()))) {
            incorrectlyAlignedResidues++;
          }
        }
      }
    }
    // Respect Chain Order
    // We want to preserver sequential order of aligned nucleotides.
    // Assuming the first nucleotide is correct and the rest should follow and increase in idx.
    if (respectOrder) {
      int curr_idx = -1;
      for (int i = 0; i < primaryNucleotidesUsed.length; i++) {
        if (primaryNucleotidesUsed[i] == 1) {
          final int targetIndex = secondaryNucleotidesMap[i];
          if (curr_idx == -1 && targetIndex > -1) {
            curr_idx = targetIndex;
            // First one is 100% correct.
            continue;
          }
          if (targetIndex > -1) {
            if (curr_idx > targetIndex) {
              // We moved backwards!
              incorrectlyAlignedResidues++;
            } else {
              // We need to update the current max in sequence.
              curr_idx = targetIndex;
            }
          }
        }
      }
    }

    incorrectlyAlignedResiduesRatio =
        Precision.round(incorrectlyAlignedResidues * 100.0 / getUsedNucleotidesNumber(), 3);
  }

  public void refinement() {
    for (int i = 0; i < secondaryNucleotidesMap.length; i++) {
      if (secondaryNucleotidesMap[i] >= 0) {
        if ((i - 1 >= 0) && (i + 1 < secondaryNucleotidesMap.length)
            && (secondaryNucleotidesMap[i + 1] >= 0) && (secondaryNucleotidesMap[i - 1] >= 0)
            && (Math.abs(secondaryNucleotidesMap[i] - secondaryNucleotidesMap[i + 1]) > 1)
            && (Math.abs(secondaryNucleotidesMap[i - 1] - secondaryNucleotidesMap[i]) > 1)) {
          primaryNucleotidesUsed[i] = -1;
          availableNucleotides.add(secondaryNucleotidesMap[i]);
          secondToFirst.remove(secondaryNucleotidesMap[i]);
          secondaryNucleotidesMap[i] = -1;
          if (!changed) {
            changed = true;
          }
        } else if ((i < secondaryNucleotidesMap.length - 2)
            && (secondaryNucleotidesMap[i + 1] == -1)
            && (secondaryNucleotidesMap[i] + 1 == secondaryNucleotidesMap[i + 2])) {
          final int index1 = i + 1;
          final int index2 = i + 2;
          int temp_primary = primaryNucleotidesUsed[index1];
          primaryNucleotidesUsed[index1] = primaryNucleotidesUsed[index2];
          primaryNucleotidesUsed[index2] = temp_primary;

          secondToFirst.replace(secondaryNucleotidesMap[index2], index1);

          temp_primary = secondaryNucleotidesMap[index1];
          secondaryNucleotidesMap[index1] = secondaryNucleotidesMap[index2];
          secondaryNucleotidesMap[index2] = temp_primary;
          if (!changed) {
            changed = true;
          }
        } else if ((i < secondaryNucleotidesMap.length - 1) && (secondaryNucleotidesMap[i + 1] >= 0)
            && (secondaryNucleotidesMap[i] == secondaryNucleotidesMap[i + 1] + 1)) {
          final int index1 = i;
          final int index2 = i + 1;
          int temp_primary = primaryNucleotidesUsed[index1];
          primaryNucleotidesUsed[index1] = primaryNucleotidesUsed[index2];
          primaryNucleotidesUsed[index2] = temp_primary;

          secondToFirst.replace(secondaryNucleotidesMap[index1], index2);
          secondToFirst.replace(secondaryNucleotidesMap[index2], index1);

          temp_primary = secondaryNucleotidesMap[index1];
          secondaryNucleotidesMap[index1] = secondaryNucleotidesMap[index2];
          secondaryNucleotidesMap[index2] = temp_primary;
          if (!changed) {
            changed = true;
          }
        }
      }
    }
  }

  public int getRemainingNucleotidesNumber() {
    return availableNucleotides.size();
  }

  public int getUsedNucleotidesNumber() {
    return secondaryNucleotides.size() - availableNucleotides.size();
  }

/*
  public int getRemainingNucleotidesNumber() {
    return availableNucleotides.size();
  }

  public int getUsedNucleotidesNumber() {
    return min_size - availableNucleotides.size();
  }
*/

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
          int selected = -1;
          if ((index - 1 >= 0) && (secondaryNucleotidesMap[index - 1] != -1)) {
            selected = secondaryNucleotidesMap[index - 1] + 1;
            final Integer val = Integer.valueOf(selected);
            if (availableNucleotides.contains(val)) {
              availableNucleotides.remove(val);
            } else {
              selected = -1;
              if ((index + 1 < secondaryNucleotidesMap.length)
                  && (secondaryNucleotidesMap[index + 1] != -1)) {
                selected = secondaryNucleotidesMap[index + 1] - 1;
                final Integer val2 = Integer.valueOf(selected);
                if (availableNucleotides.contains(val2)) {
                  availableNucleotides.remove(val2);
                } else {
                  selected = -1;
                }
              }
            }
          }
          if (selected == -1) {
            selected = getRandomAvailable();
          }
          secondaryNucleotidesMap[index] = selected;
          if (selected >= 0) {
            secondToFirst.put(selected, index);
          } else {
            primaryNucleotidesUsed[index] = -1;
          }
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
        if (selected >= 0) {
          secondToFirst.put(selected, index3);
        } else {
          primaryNucleotidesUsed[index3] = -1;
        }
        break;
    }
  }

  public void mutate() {
    changed = true;

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
    } else if (value < quadrupleMutation) {
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
          this.availableNucleotides.remove(Integer.valueOf(valueToMap));
        } else {
          // Value inside not added part
          int chance = rand.nextInt(101);
          if (chance < 45) {
            // Do not change outside part.
            this.secondaryNucleotidesMap[i] = -1;
            this.primaryNucleotidesUsed[i] = -1;
          } else if (chance < 100) {
            // Force change in original structure to take this one.
            this.secondaryNucleotidesMap[i] = valueToMap;
            int index = this.secondToFirst.get(valueToMap);
            this.secondToFirst.remove(valueToMap);

            // Reset outside part. Do not add to available as it is used above.
            this.primaryNucleotidesUsed[index] = -1;
            this.secondaryNucleotidesMap[index] = -1;
          } else {
            // Add random available
            int selected = getRandomAvailable();
            secondaryNucleotidesMap[i] = selected;
            if (selected >= 0) {
              secondToFirst.put(selected, i);
            } else {
              this.primaryNucleotidesUsed[i] = -1;
            }
          }
        }
      }
    }
  }

  public void cross(final Specimen second) {
    changed = true;

    int value = rand.nextInt(100);
    if (value < 85) {
      // Cross;
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

  public final double getRmsd() {
    return rmsd;
  }

  public final double getIncorrectlyAlignedResiduesRatio() {
    return incorrectlyAlignedResiduesRatio;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    Specimen specimen = (Specimen) o;
    return Objects.deepEquals(secondaryNucleotidesMap, specimen.secondaryNucleotidesMap);
  }

  @Override
  public int hashCode() {
    return Objects.hash(secondaryNucleotidesMap);
  }
}
