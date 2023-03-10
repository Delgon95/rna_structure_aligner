#!/bin/bash

JAR=`find $(pwd) -iname "*-with-dependencies*"`

if [[ ! -f ${JAR} ]]; then
  echo "${JAR} does not exist."
  echo "Build binaries first. Try to build with 'mvn package'"
  exit 1
fi

function run_rnahugs() {
  FOLDER="${3}_${4}"
  if [[ ${6} ]]; then
    FOLDER="${FOLDER}_${6}"
  fi
  rm -r "${FOLDER}" 2>/dev/null ; mkdir "${FOLDER}" 2>/dev/null
  COMMAND="java -Xms256m -Xmx8g -jar ${JAR} -r ${REF} -t ${TARGET} --method ${3} --mode ${4} --rmsd ${5} -o ${FOLDER}"
  if [[ ${6} ]]; then
    COMMAND="${COMMAND} --respect-order"
  fi

  echo "Running \"${COMMAND}\""
  ${COMMAND}

  cat ${FOLDER}/*-output.txt ; echo
  cat ${FOLDER}/*sequence-alignment.txt ; echo

}


for TYPE in pdb cif; do 
  for EXAMPLE in {1..4}; do
    # Check if the directory exists. Example 1 in cif is missing.
    EX_FOLDER="example/${TYPE}/${EXAMPLE}"
    if [[ -d ${EX_FOLDER} ]]; then
      FILES=$(find ${EX_FOLDER} -maxdepth 1 -iname "*.${TYPE}" | xargs -L1 basename --)
      REF=$(echo ${FILES} | tr " " "\n" | grep "solution")
      TARGET=$(echo ${FILES} | tr " " "\n" | grep -v "solution")
      
      pushd ${EX_FOLDER}
      for ORDER in "" respect-order; do
        for MODE in seq-indep seq-dep; do
          for METHOD in geometric genetic; do
            echo "Running example ${EXAMPLE} ${TYPE} - ${METHOD} ${MODE} ${ORDER} 3.5A"
            run_rnahugs ${REF} ${TARGET} ${METHOD} ${MODE} "3.5" ${ORDER}
          done
        done
      done
      popd
    fi
  done
done
