#! /bin/bash

## fastq-dump
F=$(ls -l ./ |awk '/^d/ {print $NF}')

for i in $F

do
        Fn=${i}

        fastq-dump --split-3 ${Fn}

done

