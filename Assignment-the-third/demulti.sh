#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=4                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=joycew@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=demulti                  #optional: job name
#SBATCH --output=demulti_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=demulti_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID

# didn't create an environment or activate it here, but I should!

#/usr/bin/time -v ./demulti.py -rp1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
#-rp2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
#-i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
#-i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz

/usr/bin/time -v ./demulti.py -rp1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
-rp2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
-i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
-i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz