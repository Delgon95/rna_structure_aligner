#!/bin/bash

if ! which mvn 1>/dev/null 2>/dev/null ; then
  echo "Did not found required 'mvn' (maven) command."
  exit 1
fi

if ! which java 1>/dev/null 2>/dev/null ; then
  echo "Did not found required 'java' (java) command."
  exit 2
fi

if mvn package; then
  echo "Compile finished succesfully."
  ls -l target
else
  echo "Could not compile RNAHugs package."
  exit 3
fi
