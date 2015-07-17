#!/usr/bin/env bash

# Example Usage: ./docker_helper.sh -j "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m" MarkDuplicates INPUT=fake.bam ...
while getopts "j:" OPTION; do
  case $OPTION in
    j) JVM_ARGS=$OPTARG;;
  esac
done
shift $(expr $OPTIND - 1)
TOOL_WITH_ARGS=$@

# Run Picard with JVM args and program args if present. Assumes picard jar lives in the same directory.
java ${JVM_ARGS} -jar picard.jar ${TOOL_WITH_ARGS}