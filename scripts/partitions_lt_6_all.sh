#!/bin/bash

mkdir -p partitions

./partitions_lt_6.sh 20
for i in {24..100..2}
  do
     ./partitions_lt_6.sh $i
  done

echo `cat partitions/partition_*.txt | sort | uniq | wc -l` of 47 partitions appear.
