#!/bin/bash

JAR=`pwd`/target/rna-aligner-1.0-jar-with-dependencies.jar

if [[ ! -f ${JAR} ]]; then
  echo "${JAR} does not exist."
  echo "Build binaries first. Try to build with 'mvn package'"
  exit 1
fi

pushd example/1
REF="1_solution_0A_rpr.pdb"
TARGET="testA_rpr.pdb"
OUTPUT="testA_rpr-output.txt"

echo "Running example 1 - Geometric Sequence independent 3.5A"
rm -r geometric_seq-indep 2>/dev/null ; mkdir geometric_seq-indep 2>/dev/null
java -jar ${JAR} geometric ${REF} ${TARGET} seq-indep 3.5 geometric_seq-indep/
cat geometric_seq-indep/${OUTPUT} ; echo

echo "Running example 1 - Genetic Sequence independent 3.5A"
rm -r genetic_seq-indep 2>/dev/null ; mkdir genetic_seq-indep 2>/dev/null
java -jar ${JAR} genetic ${REF} ${TARGET} seq-indep 3.5 genetic_seq-indep/
cat genetic_seq-indep/${OUTPUT} ; echo

echo "Running example 1 - Geometric Sequence dependent 3.5A"
rm -r geometric_seq-dep 2>/dev/null ; mkdir geometric_seq-dep 2>/dev/null
java -jar ${JAR} geometric ${REF} ${TARGET} seq-dep 3.5 geometric_seq-dep/ 
cat geometric_seq-dep/${OUTPUT} ; echo

echo "Running example 1 - Genetic Sequence dependent 3.5A"
rm -r genetic_seq-dep 2>/dev/null ; mkdir genetic_seq-dep 2>/dev/null
java -jar ${JAR} genetic ${REF} ${TARGET} seq-dep 3.5 genetic_seq-dep/
cat genetic_seq-dep/${OUTPUT} ; echo

popd

pushd example/2
REF="1_solution_0_rpr.pdb"
TARGET="1_das_3_rpr.pdb"
OUTPUT="1_das_3_rpr-output.txt"

echo "Running example 2 - Geometric Sequence independent 3.5A"
rm -r geometric_seq-indep 2>/dev/null ; mkdir geometric_seq-indep 2>/dev/null
java -jar ${JAR} geometric ${REF} ${TARGET} seq-indep 3.5 geometric_seq-indep/
cat geometric_seq-indep/${OUTPUT} ; echo

echo "Running example 2 - Genetic Sequence independent 3.5A"
rm -r genetic_seq-indep 2>/dev/null ; mkdir genetic_seq-indep 2>/dev/null
java -jar ${JAR} genetic ${REF} ${TARGET} seq-indep 3.5 genetic_seq-indep/
cat genetic_seq-indep/${OUTPUT} ; echo

echo "Running example 2 - Geometric Sequence dependent 3.5A"
rm -r geometric_seq-dep 2>/dev/null ; mkdir geometric_seq-dep 2>/dev/null
java -jar ${JAR} geometric ${REF} ${TARGET} seq-dep 3.5 geometric_seq-dep/ 
cat geometric_seq-dep/${OUTPUT} ; echo

echo "Running example 2 - Genetic Sequence dependent 3.5A"
rm -r genetic_seq-dep 2>/dev/null ; mkdir genetic_seq-dep 2>/dev/null
java -jar ${JAR} genetic ${REF} ${TARGET} seq-dep 3.5 genetic_seq-dep/
cat genetic_seq-dep/${OUTPUT} ; echo

popd

pushd example/3
REF="4_0_solution_3V7E_rpr.pdb"
TARGET="4_chen_1_rpr.pdb"
OUTPUT="4_chen_1_rpr-output.txt"

echo "Running example 3 - Geometric Sequence independent 3.5A"
rm -r geometric_seq-indep 2>/dev/null ; mkdir geometric_seq-indep 2>/dev/null
java -jar ${JAR} geometric ${REF} ${TARGET} seq-indep 3.5 geometric_seq-indep/
cat geometric_seq-indep/${OUTPUT} ; echo

echo "Running example 3 - Genetic Sequence independent 3.5A"
rm -r genetic_seq-indep 2>/dev/null ; mkdir genetic_seq-indep 2>/dev/null
java -jar ${JAR} genetic ${REF} ${TARGET} seq-indep 3.5 genetic_seq-indep/
cat genetic_seq-indep/${OUTPUT} ; echo

echo "Running example 3 - Geometric Sequence dependent 3.5A"
rm -r geometric_seq-dep 2>/dev/null ; mkdir geometric_seq-dep 2>/dev/null
java -jar ${JAR} geometric ${REF} ${TARGET} seq-dep 3.5 geometric_seq-dep/ 
cat geometric_seq-dep/${OUTPUT} ; echo

echo "Running example 3 - Genetic Sequence dependent 3.5A"
rm -r genetic_seq-dep 2>/dev/null ; mkdir genetic_seq-dep 2>/dev/null
java -jar ${JAR} genetic ${REF} ${TARGET} seq-dep 3.5 genetic_seq-dep/
cat genetic_seq-dep/${OUTPUT} ; echo

popd

pushd example/4
REF="3_solution_0_rpr.pdb"
TARGET="3_chen_1_rpr.pdb"
OUTPUT="3_chen_1_rpr-output.txt"

echo "Running example 4 - Geometric Sequence independent 3.5A"
rm -r geometric_seq-indep 2>/dev/null ; mkdir geometric_seq-indep 2>/dev/null
java -jar ${JAR} geometric ${REF} ${TARGET} seq-indep 3.5 geometric_seq-indep/
cat geometric_seq-indep/${OUTPUT} ; echo

echo "Running example 4 - Genetic Sequence independent 3.5A"
rm -r genetic_seq-indep 2>/dev/null ; mkdir genetic_seq-indep 2>/dev/null
java -jar ${JAR} genetic ${REF} ${TARGET} seq-indep 3.5 genetic_seq-indep/
cat genetic_seq-indep/${OUTPUT} ; echo

echo "Running example 4 - Geometric Sequence dependent 3.5A"
rm -r geometric_seq-dep 2>/dev/null ; mkdir geometric_seq-dep 2>/dev/null
java -jar ${JAR} geometric ${REF} ${TARGET} seq-dep 3.5 geometric_seq-dep/ 
cat geometric_seq-dep/${OUTPUT} ; echo

echo "Running example 4 - Genetic Sequence dependent 3.5A"
rm -r genetic_seq-dep 2>/dev/null ; mkdir genetic_seq-dep 2>/dev/null
java -jar ${JAR} genetic ${REF} ${TARGET} seq-dep 3.5 genetic_seq-dep/
cat genetic_seq-dep/${OUTPUT} ; echo
