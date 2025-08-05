#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=4                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=joycew@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=demulti_part1_R1      #optional: job name
#SBATCH --output=demulti_part1_R1_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=demulti_part1_R1_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


# change sbatch notes to change job and output/error names as I go through each FASTQ file

# /usr/bin/time -v ./part1_script.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -l 101

# /usr/bin/time -v ./part1_script.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -l 101

# /usr/bin/time -v ./part1_script.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -l 8

# /usr/bin/time -v ./part1_script.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -l 8


