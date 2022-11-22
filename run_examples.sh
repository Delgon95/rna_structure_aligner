#!/bin/bash

JAR=`find $(pwd) -iname "*-with-dependencies*"`

if [[ ! -f ${JAR} ]]; then
  echo "${JAR} does not exist."
  echo "Build binaries first. Try to build with 'mvn package'"
  exit 1
fi

function run_rnahugs() {

  FOLDER="${3}_${4}"
  rm -r "${FOLDER}" 2>/dev/null ; mkdir "${FOLDER}" 2>/dev/null
  java -jar ${JAR} -r ${REF} -t ${TARGET} --method ${3} --mode ${4} --rmsd ${5} -o "${FOLDER}"
  cat ${FOLDER}/*-output.txt ; echo
  cat ${FOLDER}/*sequence-alignment.txt ; echo

}



# Run examples in CIF format
pushd example/cif/2
REF="1_solution_0_rpr.cif"
TARGET="1_das_3_rpr.cif"

echo "Running example 2 cif - Geometric Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-indep" "3.5"

echo "Running example 2 cif - Genetic Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-indep" "3.5"

echo "Running example 2 cif - Geometric Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-dep" "3.5"

echo "Running example 2 cif - Genetic Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-dep" "2.5"

popd


pushd example/cif/3
REF="4_0_solution_3V7E_rpr.cif"
TARGET="4_chen_1_rpr.cif"

echo "Running example 3 cif - Geometric Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-indep" "3.5"

echo "Running example 3 cif - Genetic Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-indep" "3.5"

echo "Running example 3 cif - Geometric Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-dep" "3.5"

echo "Running example 3 cif - Genetic Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-dep" "3.5"

popd

pushd example/cif/4
REF="3_solution_0_rpr.cif"
TARGET="3_chen_1_rpr.cif"

echo "Running example 4 cif - Geometric Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-indep" "3.5"

echo "Running example 4 cif - Genetic Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-indep" "3.5"

echo "Running example 4 cif - Geometric Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-dep" "3.5"

echo "Running example 4 cif - Genetic Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-dep" "3.5"

popd

# Run examples in PDB format
pushd example/pdb/1
REF="1_solution_0A_rpr.pdb"
TARGET="testA_rpr.pdb"

echo "Running example 1 pdb - Geometric Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-indep" "3.5"

echo "Running example 1 pdb - Genetic Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-indep" "3.5"

echo "Running example 1 pdb - Geometric Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-dep" "3.5"

echo "Running example 1 pdb - Genetic Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-dep" "3.5"

popd

pushd example/pdb/2
REF="1_solution_0_rpr.pdb"
TARGET="1_das_3_rpr.pdb"

echo "Running example 2 pdb - Geometric Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-indep" "3.5"

echo "Running example 2 pdb - Genetic Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-indep" "3.5"

echo "Running example 2 pdb - Geometric Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-dep" "3.5"

echo "Running example 2 pdb - Genetic Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-dep" "3.5"

popd


pushd example/pdb/3
REF="4_0_solution_3V7E_rpr.pdb"
TARGET="4_chen_1_rpr.pdb"

echo "Running example 3 pdb - Geometric Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-indep" "3.5"

echo "Running example 3 pdb - Genetic Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-indep" "3.5"

echo "Running example 3 pdb - Geometric Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-dep" "3.5"

echo "Running example 3 pdb - Genetic Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-dep" "3.5"

popd

pushd example/pdb/4
REF="3_solution_0_rpr.pdb"
TARGET="3_chen_1_rpr.pdb"

echo "Running example 4 pdb - Geometric Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-indep" "3.5"

echo "Running example 4 pdb - Genetic Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-indep" "3.5"

echo "Running example 4 pdb - Geometric Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-dep" "3.5"

echo "Running example 4 pdb - Genetic Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-dep" "3.5"

popd

