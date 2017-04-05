# Get batches
./get_batches_all.sh params.sh

# Run autosomal estimates - maxlik
./run_autosomal_ml.sh params.sh

# Run autosomal estimates - smm
./run_autosomal_smm.sh params.sh

# Combine results while still running so we can peek
./gather_completed.sh params.sh

# Filter and scale
#bsub -q short -W 4:00 -eo filter.err -oo filter.log \
./filter.sh params.sh
./filter_smm.sh params.sh
