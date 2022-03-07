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
  java -jar ${JAR} -r ${REF} -t ${TARGET} -method ${3} --mode ${4} --rmsd ${5} -o "${FOLDER}"
  cat ${FOLDER}/*-output.txt ; echo

}

pushd example/1
REF="1_solution_0A_rpr.pdb"
TARGET="testA_rpr.pdb"

echo "Running example 1 - Geometric Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-indep" "3.5"

echo "Running example 1 - Genetic Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-indep" "3.5"

echo "Running example 1 - Geometric Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-dep" "3.5"

echo "Running example 1 - Genetic Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-dep" "3.5"

popd

pushd example/2
REF="1_solution_0_rpr.pdb"
TARGET="1_das_3_rpr.pdb"

echo "Running example 2 - Geometric Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-indep" "3.5"

echo "Running example 2 - Genetic Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-indep" "3.5"

echo "Running example 2 - Geometric Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-dep" "3.5"

echo "Running example 2 - Genetic Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-dep" "3.5"

popd

pushd example/3
REF="4_0_solution_3V7E_rpr.pdb"
TARGET="4_chen_1_rpr.pdb"

echo "Running example 3 - Geometric Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-indep" "3.5"

echo "Running example 3 - Genetic Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-indep" "3.5"

echo "Running example 3 - Geometric Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-dep" "3.5"

echo "Running example 3 - Genetic Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-dep" "3.5"

popd

pushd example/4
REF="3_solution_0_rpr.pdb"
TARGET="3_chen_1_rpr.pdb"

echo "Running example 4 - Geometric Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-indep" "3.5"

echo "Running example 4 - Genetic Sequence independent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-indep" "3.5"

echo "Running example 4 - Geometric Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "geometric" "seq-dep" "3.5"

echo "Running example 4 - Genetic Sequence dependent 3.5A"
run_rnahugs ${REF} ${TARGET} "genetic" "seq-dep" "3.5"
