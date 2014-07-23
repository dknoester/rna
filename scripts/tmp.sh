#!/bin/bash

rsync -avz --exclude-from=${HOME}/research/etc/rsync-excludes \
    -e ssh dk@hpc.msu.edu:src/rna/expr/$1/$2/ tmp/
