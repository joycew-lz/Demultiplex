#!/usr/bin/env python

import bioinfo
import gzip
import matplotlib.pyplot as plt

#--------------
# argparse
#--------------

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="to get file name and length of the list for a script that plots mean q_scores")
    parser.add_argument("-f", help="to specify the file name", type=str, required=True)
    parser.add_argument("-l", help="to specify the length of the list to be created", type=int, required=True)
    return parser.parse_args()

args = get_args()

# Set global variables
file_name = args.f
list_length = args.l

#--------------
# init_list function
#--------------

def init_list(lst: list, value: float=0.0) -> list: # float=0.0 indicates default for value
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    for i in range(list_length):
        lst.append(value)
    return lst # return outside the for loop

# initialize the list init_list()
my_list: list = []
init_list(my_list)

#--------------
# populate list
#--------------
        
def populate_list(file: str): # -> tuple[list, int]
    """Opens a FASTQ file and loops through every record to
    convert the quality score of each line to a numerical Phred quality score with the convert_phred() function.
    Stores each Phred quality score for each base pair to an ongoing sum.
    Returns a sum of the quality scores and the number of lines in the file for each position."""

    # initialize the list init_list()
    my_list: list = []
    init_list(my_list)
    
    num_lines = 0
    
    # open fastq file with gzip and "rt"
    with gzip.open(file_name, "rt") as f:
        for line in f:
            num_lines += 1
            line = line.strip()

            if num_lines % 4 == 0: 
                # NEED AN ONGOING SUM OF THE QUALITY SCORE FOR EVERY POSITION (0-100), use enumerate
                for score_index, score in enumerate(line):
                    my_list[score_index] += bioinfo.convert_phred(score)

    return my_list, num_lines
    
# populate the list
my_list, num_lines = populate_list(file_name)


#--------------
# calculate mean q_score at each base,
# print the list
#--------------

# Iterate over a list when you need the index; use enumerate
print(f'# Base Pair\t Mean Quality Score')
for index, score in enumerate(my_list):
    mean = score/(num_lines/4) # calculate mean for each index of my_list: score/(num_lines/4)
    my_list[index] = mean # store the mean quality value at each base back into my_list at the appropriate position
    print(f'{index}\t{mean}')

#--------------
# plotting
#--------------


x = list(range(len(my_list)))
y = my_list
plt.figure(figsize=(12, 6))
plt.plot(x, y)
plt.grid()
plt.title(f"Mean Phred Quality Score per Base Position for {file_name}")
plt.xlabel("Base Position")
plt.ylabel("Mean Phred Quality Score")

plot_name = "mean_qscores.png"
plt.savefig(plot_name)
print(plot_name)
plt.show()