#!/bin/bash

if [ -f partitions/partition_$1.txt ]; then
    # Do not reconstruct the file if it exists
    exit 0
fi

time (fullgen $1 code 1 stdout logerr | pentagon_partition | filter_valid_clusters > partitions/partition_$1.txt)
