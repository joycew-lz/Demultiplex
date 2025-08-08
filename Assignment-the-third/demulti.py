#!/usr/bin/env python

#--------------
# import modules
#--------------

import bioinfo
import argparse
import math
import gzip 
import matplotlib.pyplot as plt

#--------------
# argparse
#--------------

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="setting R1, R2, R3, and R4 for a demultiplexing script")
    parser.add_argument("-rp1", help="to specify file name for bioread 1, aka R1", type=str, required=True)
    parser.add_argument("-rp2", help="to specify file name for bioread 2, aka R4", type=str, required=True)
    parser.add_argument("-i1", help="to specify file name for index 1, aka R2", type=str, required=True)
    parser.add_argument("-i2", help="to specify file name for index 2, aka R3", type=str, required=True)
    return parser.parse_args()

args = get_args()

# Set global variables
bioread1 = args.rp1
bioread2 = args.rp2
index1 = args.i1
index2 = args.i2

#--------------
# initialize file names for writing into
#--------------

unknown_R1 = "unknown_R1.fq"
unknown_R2 = "unknown_R2.fq"
hopped_R1 = "hopped_R1.fq"
hopped_R2 = "hopped_R2.fq"


#--------------
# initialize counters/dictionaries for counting, to answer summary questions
#--------------

hopped_counter = 0
unknown_counter = 0
possible_combinations = {} # Initialize a dictionary for counting each possible pair of indexes (both swapped and dual-matched)
dual_matched_only = {} # Initialize a dictionary for counting each pair of properly matched indexes.
hopped_only = {}  # Initialize a dictionary for counting each pair of index-hopped indexes.

#--------------
# 24 known indexes/barcodes
#--------------

known_indexesFile = open("/projects/bgmp/joycew/bioinfo/Bi622/Demultiplex/known_indexes.txt")
known_indexes = set()
for known_index in known_indexesFile:
    known_index = known_index.strip()
    known_indexes.add(known_index)

#--------------
# opening dual-matched files to write into (in the big loop later)
    # loop through each index in the known barcodes set,
    # for each index, open two output files: one for R1 reads and one for R2 reads.
    # store the file handles in a dictionary: the key is the index,
        # and the value is a tuple of (R1_handle, R2_handle)
#--------------

matched_filehandles_dict = {}

for i in known_indexes:
    matched_R1 = open(f"{i}_R1_.fq", "w")
    matched_R2 = open(f"{i}_R2_.fq", "w")
    matched_filehandles_dict[i] = (matched_R1, matched_R2)

# in big loop:

    # for i in range(10):
    # fh.write(f"{i}: Leslie is awesome\n")

#--------------
# reverse complement function
#--------------

def reverse_complement(dna_seq: str):
    '''
    Return the reverse complement of a DNA sequence string.
    Example: 'GTAGCGTA' → 'TACGCTAC'
    '''
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    # Reverse the sequence
    reversed_sequence = dna_seq[::-1]
    # Complement each base
    complement_sequence = "".join([complement_dict[base] for base in reversed_sequence])
    return complement_sequence

#--------------
# opening unknown and hopped files to write into (in the big loop)
# parsing it through four files, one record (four lines) at a time...
# starting the big loop!
    # REMEMBER: writing into unknown and hopped files, writing into fh_R1 and fh_R2 (how to do so is written above)
    # REMEMBER: my counters/dicts for summarizing the results
    # REMEMBER: close all files im writing into when the big loop ends
#--------------


with open(unknown_R1, 'w') as unknown_fh_R1, open(unknown_R2, 'w') as unknown_fh_R2, \
    open(hopped_R1, 'w') as hopped_fh_R1, open(hopped_R2, 'w') as hopped_fh_R2:
    with gzip.open(bioread1, 'rt') as rp1, gzip.open(bioread2, 'rt') as rp2, gzip.open(index1, 'rt') as i1, gzip.open(index2, 'rt') as i2:
        
        # initialize empty lists
        bioread1_lines = []
        bioread2_lines = []
        index1_lines = []
        index2_lines = []

        while True: # remember to strip new lines
            bioread1_lines = [rp1.readline().strip() for i in range(4)] # list of all four lines of a record stored for bioread1
            bioread2_lines = [rp2.readline().strip() for i in range(4)] # index 0 is the header, index 1 is the seq, index 2 is +, index 3 is qs
            index1_lines = [i1.readline().strip() for i in range(4)] # this is using list comprehension.
            index2_lines = [i2.readline().strip() for i in range(4)]

            # exit condition: stop when any of the files are done
            if not bioread1_lines[0] or not bioread2_lines[0] or not index1_lines[0] or not index2_lines[0]:
                break
        
            # grab the sequences from index1 and index2 to get the index-pair:
            index1_seq = index1_lines[1]
            index2_seq = index2_lines[1]
            rev_index2_seq = reverse_complement(index2_seq)
            index_pair = f"{index1_seq}-{rev_index2_seq}"

            # modify header: append index_pair to the end of bioread1 and bioread2's headers
            bioread1_lines[0] = bioread1_lines[0] + f" {index_pair}"
            bioread2_lines[0] = bioread2_lines[0] + f" {index_pair}" 

            # INVALID, unknown category:
            if ("N" in index1_seq or "N" in index2_seq) or (index1_seq not in known_indexes or rev_index2_seq not in known_indexes):
                # write to unknown_R1 and unknown_R2
                for line in bioread1_lines:
                    unknown_fh_R1.write(line + '\n')
                for line in bioread2_lines:
                    unknown_fh_R2.write(line + '\n')
                # increment unknown_counter
                unknown_counter += 1
                continue

            # VALID, meaning index1 and (reverse complement) index2 seq lines found in the 24 known indexes:
            elif index1_seq in known_indexes and rev_index2_seq in known_indexes:
                # Add to possible_combinations dict
                if index_pair in possible_combinations:
                    possible_combinations[index_pair] += 1
                else:
                    possible_combinations[index_pair] = 1

                # VALID, dual-matched category:
                if index1_seq == rev_index2_seq:
                    # write to matched_R1 and matched_R2
                    matched_fh_R1, matched_fh_R2 = matched_filehandles_dict[index1_seq]
                    for line in bioread1_lines:
                        matched_fh_R1.write(line + '\n')
                    for line in bioread2_lines:
                        matched_fh_R2.write(line + '\n')
                    # add to dual-matched dict
                    if index_pair in dual_matched_only:
                        dual_matched_only[index_pair] += 1
                    else:
                        dual_matched_only[index_pair] = 1
                
                # VALID, index-hopped category:
                else:
                    # write to hopped_R1 and hopped_R2
                    for line in bioread1_lines:
                        hopped_fh_R1.write(line + '\n')
                    for line in bioread2_lines:
                        hopped_fh_R2.write(line + '\n')
                    # increment hopped counter
                    hopped_counter += 1
                    # add to hopped_only dict
                    if index_pair in hopped_only:
                        hopped_only[index_pair] += 1
                    else:
                        hopped_only[index_pair] = 1

# remember to close all the files i'm writing into!
# with open automatically closes, so just need to close the matched files.
for matched_R1, matched_R2 in matched_filehandles_dict.values():
        matched_R1.close()
        matched_R2.close()

#--------------
# answering summary questions: summary stats and reporting
#--------------

with open("demulti_summary.md", "w") as out:
    # Getting total number of pairs in each category
    total = sum(possible_combinations.values()) + unknown_counter # TOTAL: matched, hopped, and unknown
    matched = sum(dual_matched_only.values()) # total dual-matched index pairs
    # total index-hopped index pairs is hopped_counter
    # total unknown is unknown_counter

    out.write("\n### Summary Statistics\n")
    out.write(f"Total Number of Reads: {total}  \n")
    out.write(f"Number of Read-Pairs with properly matched indexes: {matched:}  \n")
    out.write(f"Number of Read-Pairs that have index-hopped: {hopped_counter:}  \n")
    out.write(f"Number of Read-Pairs with unknown indexes: {unknown_counter:}  \n")
    out.write(f"Percent Dual-Matched: {matched / total * 100:.3f}%  \n")
    out.write(f"Percent Index-Hopped: {hopped_counter / total * 100:.3f}%  \n")
    out.write(f"Percent Unknown: {unknown_counter / total * 100:.3f}%  \n")

    # report the number and percent of each dual-matched pair:
    out.write("\n### Dual-Matched Index Pairs:\n")
    out.write("| Index Pair | Count | Percentage (%) |\n")
    out.write("|---|---|---|\n")
    for pair, count in dual_matched_only.items():
        percent = (count / matched) * 100
        out.write(f"| {pair} | {count} | {percent:.3f}% |\n")
    
    # report the number and percent of each index-hopped pair:
    out.write("\n### Index-Hopped Index Pairs:\n")
    out.write("| Index Pair | Count | Percentage (%) |\n")
    out.write("|---|---|---|\n")
    for pair, count in hopped_only.items():
        percent = (count / hopped_counter) * 100
        out.write(f"| {pair} | {count} | {percent:.3f}% |\n")

#--------------
# answering summary questions: visualizations
    # Bar chart for visualizing the number of dual-matched index pairs
#--------------

plt.figure(figsize=(12, 6))
plt.bar(list(dual_matched_only.keys()), list(dual_matched_only.values()), color="steelblue")
plt.xlabel("Index Pair")
plt.ylabel("Number of Reads")
plt.title("Number of Dual-Matched Index Pairs")
plt.xticks(rotation=45, fontsize=6)
plt.savefig("matched_barchart.svg")

#--------------
# create a matrix (unfinished idea)
#--------------

# Matrix: want unique values only in row and column of matrix displaying contents of possible_combinations
      # for all possible pairs of indexes:
        # use a matrix? row = index_pair_1 and column = index_pair_2. content is the frequency of each
        # possible pair of indexes.
        
# index_pair_1_set = set()
# index_pair_2_set = set()
# for k in possible_combinations.keys():
#     index_pair_1, index_pair_2 = k.split('-')
#     index_pair_1_set.add(index_pair_1)
#     index_pair_2_set.add(index_pair_2)


