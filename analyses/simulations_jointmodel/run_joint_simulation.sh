#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

for logmu in $MUTMU
do
    for lencoeff in $LENCOEFF
    do
	# Run simulations
	bsub -q short -W 1:00 -eo ${LOGDIR}/simulate.${logmu}.${lencoeff}.err \
	    -oo ${LOGDIR}/simulate.${logmu}.${lencoeff}.out \
	    -J simulate.${logmu}.${lencoeff} \
	    ./simulate_locus.sh ${logmu} ${lencoeff} ${PARAMFILE}	
	# Get files for estimation
	bsub -q short -W 1:00 -eo ${LOGDIR}/combine.${logmu}.${lencoeff}.err \
	    -oo ${LOGDIR}/combine.${logmu}.${lencoeff}.out \
	    -w "ended(simulate.${logmu}.${lencoeff})" \
	    -J combine.${logmu}.${lencoeff} \
	    ./get_locus_files.sh ${logmu} ${lencoeff} ${PARAMFILE}
	# Run joint estimation
	bsub -q short -W 4:00 -eo ${LOGDIR}/estimate.${logmu}.${lencoeff}.err \
	    -oo ${LOGDIR}/estimate.${logmu}.${lencoeff}.out -n ${NUMPROC} \
	    ./joint_estimation.sh ${logmu} ${lencoeff} ${PARAMFILE}
	    -w "ended(combine.${logmu}.${lencoeff})" \
    done
done