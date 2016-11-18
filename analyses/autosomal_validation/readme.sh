#!/bin/bash

#bsub -q long -W 24:00 -eo codis.err -oo codis.out \
./run.sh codisparams.sh

#bsub -q long -W 72:00 -eo marsh.err -oo marsh.out -n 5 \
./run.sh marshparams.sh

./calibrate_errors.sh codisparams.sh
./calibrate_errors.sh marshparams.sh
