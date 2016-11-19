#!/bin/bash
./run_joint_simulation.sh params_joint.sh
./combine_results.sh params_joint.sh | awk '(NF==6)' > simulation_joint_results.txt
