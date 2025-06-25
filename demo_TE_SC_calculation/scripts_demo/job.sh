#!/bin/bash

#SBATCH --job-name=demoTE-SC-XCov
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/dev/null
#SBATCH --array=[1-2]%2
#SBATCH --mem=10G
#SBATCH -t 10:00:00

cd $SLURM_SUBMIT_DIR

# if you don't want to use SLURM, just erase lines from 3 to 11 and change manually the variables assignation

data_type='MaxOne' # change to MaxOne to compute EC on the other dataset

if [[ $data_type = 'CRCNS' ]]; then 

    main_path="../Data_CRCNS/"
    array_config="./config_CRCNS.txt"
    path_out="../Output_CRCNS/"

elif [[ $data_type = 'MaxOne' ]]; then 

    main_path="../Data_MaxOne/"
    array_config="./config_MaxOne.txt"
    path_out="../Output_MaxOne/"


else
   echo "Error! Type of data should be CRCNS or MaxOne"

fi


# bin size for data binarization
bs=$(awk -v ArrayTaskID="$SLURM_ARRAY_TASK_ID" '$1==ArrayTaskID {print $2}' $array_config)
Culture=$(awk -v ArrayTaskID="$SLURM_ARRAY_TASK_ID" '$1==ArrayTaskID {print $4}' $array_config)

# max delay for TE, SC, XCov estimation
delay_ms=10.0   # unit: ms (max delay for TE, SC, XCov computation)

CI_tau_window=3 # time window length for CI calculation. 
                # unit: number of binsizes bs
echo $bs
echo $delay_ms

# parameters for significance test
time_jitter=1.   # ms
n_permutations=5 # number of null models

sleep $((RANDOM % 4 + 1))


matlab -nodisplay -nosplash -nodesktop -r " TECN_Calculation('$Culture','$bs','$delay_ms','$CI_tau_window','$main_path','$path_out'); exit;"

wait

matlab -nodisplay -nosplash -nodesktop -r " TECN_Sign('$Culture','$bs','$delay_ms','$CI_tau_window','$time_jitter','$n_permutations','$main_path','$path_out'); exit;"

