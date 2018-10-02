#!/bin/bash
while [[ $# > 0 ]]; do
    qsub $1
    shift
done
echo "Finish"