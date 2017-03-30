# Get Y-STR pairs using snpmu=2.1775*10**-8 
# http://www.nature.com/ng/journal/v47/n5/full/ng.3171.html 8.71*10**-10*25

source params_base.sh

# Run Y-STR analysis - SGDP
#bsub -q long -W 200:00 -eo log/sgdpy_ml.err -oo log/sgdpy_ml.out -J sgdpy.ml -n 5 \
sbatch --job-name=sgdpy.ml --time 2000 --partition=compute \
    --account=${ACCOUNT} --get-user-env \
    --error=${BASEDIR}/ystr_validation/log/sgdpy_ml.err \
    --output=${BASEDIR}/ystr_validation/log/sgdpy_ml.out \
    ./run_ystr_estimation_ml.sh sgdp_params.sh

# Run Y-STR analysis - 1000 Genomes
#bsub -q long -W 200:00 -eo log/1kgy_ml.err -oo log/1kgy_ml.out -J 1kgy.ml -n 5 \
sbatch --job-name=1kgy.ml --time 2000 --partition=compute \
    --account=${ACCOUNT} --get-user-env \
    --error=${BASEDIR}/ystr_validation/log/1kgy_ml.err \
    --output=${BASEDIR}/ystr_validation/log/1kgy_ml.out \
    ./run_ystr_estimation_ml.sh 1kg_params.sh

# Calibrate errors
./calibrate_errors.sh sgdp_params.sh
./calibrate_errors.sh 1kg_params.sh


######## Try Y-STRs using SMM
./run_ystr_estimation_smm.sh 1kg_params.sh
./run_ystr_estimation_smm.sh sgdp_params.sh
