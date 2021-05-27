#!/bin/bash
#
# Script to run an instance of the fitness provider interactively, i.e., this will lock your terminal until the task is finished.
# Note that the RUN-TEST_tmpl folder contains the output of the DFT calculation.
# This allow to manually intervene to shorten the runtime. Do this:
# 1) kill the job submitted to the queue,
# 2) add the 'GAUSSIAN JOB ENDED' string to the *_DFT.log file (see the *_FProvider.log file to verify what the script is waiting for)
# 3) wait for the sleeping process to wake-up and pick-up the results to produce the fitness file *_FIT.sdf
#
cd ..
rm -rf RUN-TEST
cp -r test/RUN-TEST_tmpl RUN-TEST
dir="$(pwd)"
# This is a "small" molecule: it takes ca. 10 minutes to run on 40 cpus at saga.sigma2.no
bash Ru_14-el_fitness_BndLng.sh "$dir/RUN-TEST/Gen000/M00000005_I.sdf" "$dir/RUN-TEST/Gen000/M00000005_FIT.sdf" "$dir/RUN-TEST/Gen000" 5 "$dir/RUN-TEST/MOLUID.txt"
#
# For a larger molecule (takes a bit more than an hour on 40 cpus at saga-sigma2.no)
#bash Ru_14-el_fitness_BndLng.sh "$dir/RUN-TEST/Gen000/M00000002_I.sdf" "$dir/RUN-TEST/Gen000/M00000002_FIT.sdf" "$dir/RUN-TEST/Gen000" 2 "$dir/RUN-TEST/MOLUID.txt"
