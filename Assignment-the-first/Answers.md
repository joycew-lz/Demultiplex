# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101 bp | phred+33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8 bp | phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8 bp | phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101 bp | phred+33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem:

    We are given four "corresponding" FASTQ files generated from the 2017 BGMP cohort's library preps. The files represent read pairs and their corresponding barcodes:
    * R1 and R4 contain the actual sequencing reads (inserts).
    * R2 and R3 contain the index sequences (barcodes).
    R1 and R4 are a read pair while R2 and R3 are their respective indexes, so their index pair. Each group of four lines from all four files represents one read-pair record with their forward and reverse indexes.
    We must keep in mind that R3 is the reverse complement to R2. 
    We are also given 24 known/dual-matched indexes in a text file. These are the expected, correctly paired barcodes.
    ----
    We have to demultiplex the data by using the barcode information in order to know which sequences came from which samples after they were all sequenced together.
    To demultiplex the data, I need to go through every read (record) of the four files and demultiplex by index-pair (R2 and R3).
    * INVALID (barcodes not found in the list of 24 known indexes): If R2 and R3 index sequences are not in the list of 24 known indexes, or contain "N", they are categorized as "unknown".
    * VALID (barcodes found in the list of 24 known indexes): If R2 and R3 index sequences match each other, they are categorized as "dual-matched".
    * VALID (barcodes found in the list of 24 known indexes): if R2 and R3 index sequences do not match each other, they are categorized as "index-hopping". 


2. Describe output:

    FOR EACH RECORD, look at the two index reads (index-pairs) R2 and R3 and their sequences (the barcodes). SORT into one of the three categories below.

    ---- INVALID
    * If one or both index sequences (the barcodes) are not in the list of 24 known indexes, or contain "N" (unknown/low quality category):
        * Write out the record (as in the record of R1 and R4, the read pair), with modified headers, into:
            * unknown_R1.fq (for read1)
            * unknown_R2.fq (for read2)
        * ALSO, if either of the index reads' mean quality score is lower than the threshold determined in part 1, we'll put it in the same(unknown/low quality category).

    ---- VALID
    * If the index-pair's sequences match each other (consider R3 is reverse complement of R2) (dual-matched):
        * Write out the record (as in the record of R1 and R4, the read pair), with modified headers, into one read1 FASTQ file and one read2 FASTQ file **per matching index-pair**.
            * 48 FASTQ files
            * named like (matchingindex_R1.fq and matchingindex_R2.fq)
                * AAAA_R1.fq
                * AAAA_R2.fq
    * If the index-pair's sequences DO NOT match each other (consider R3 is reverse complement or R2) (index-hopping):
        * Write out the record (as in the record of R1 and R4, the read pair), with modified headers, into:
            * hopped_R1.fq
            * hopped_R2.fq

    At the end, get a summary of:
        
    * the number of read pairs with index-hopping observed,
    * the number of read-pairs with unknown index(es),
    * the numbers for each possible pair of indexes (both swapped and dual matched),
    * the number of read-pairs with properly matched indexes (per index-pair).

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
Done!

4. Pseudocode

    import argparse
    
    import bioinfo (and be sure to update bioinfo.py (and git add/commit/push) before making a copy to my Demultiplexing folder)
    
    import gzip

    argparse:
    - readpair1 = args.rp1 (FILE: R1)
    - readpair2 = args.rp2 (FILE: R4)
    - index1 = args.i1 (FILE: R2)
    - index2 = args.i2 (FILE: R3)
    - qc_cutoff = args.q (DETERMINED IN PART 1)

    Initialize file names for writing into.

    ```
    unknown_r1 = "unknown_R1.fq"
    unknown_r2 = "unknown_R2.fq"
    hopped_r1 = "hopped_R1.fq"
    hopped_r2 "hopped_R2.fq"
    ```

    Answering summary questions: Initialize counters. Initialize a dictionary for counting each possible pair of indexes (both swapped and dual matched). Initialize a dictionary for counting each pair of properly matched indexes.

    ```
    hopped_counter = 0
    unknown_counter = 0
    possible_combinations = {}
    dual_matched_only = {}
    ```

```
remember to strip :)

with gzip.open(readpair1, 'r') as rp1, open('readpair2', 'r') as rp2, open('index1', 'r') as i1, open ('index2', 'r') as i2:
    
    initialize line1, line2, line3, line4 to be empty lists 
    
    parse through each record, looking at all four files at the same time; get the four lines of a record stored on these lists:
        [0th index is from readpair1, 1st is from readpair2, 2nd is from index1, 3rd is from index2]
        line1 list = header
        line2 list = sequence
        line3 list = +
        line4 list = quality

        looking at index1 and index2 (indices 2 and 3 of our lines):
            
            - (*) modify the header for EVERYTHING:
            look at index1 and index2's sequence line: line2[2] and line2[3]
            grab the sequences from index1 and index2 to get the index-pair:
                index1_seq = line2[2].strip()
                index2_seq = line2[3].strip()
                rev_com_index2_seq = reverse_complement(index2_seq)
                
                index_pair = f"{index1_seq}-{rev_com_index2_seq}"
            append index_pair to the end of readpair1 and readpair 2's headers in line1: line1[0] and line1[1]
                line1[0] = line1[0].strip() + f" {index_pair}\n" 
                line1[1] = line1[1].strip() + f" {index_pair}\n"

            - (1) [if statement] do we sort into the unknown/low quality category cuz it doesn't match the 24 known barcodes, or contains an N (aka is this INVALID)?
            look at index1 and index2's sequence line: line2[2] and line2[3]
            if index1 seq line or index2 seq line NOT found in the 24 known matches OR contains "N":
                write this record (all four lines) from readpair1 to unknown_r1
                write this record (all four lines) from readpair2 to unknown_r2
                unknown_counter += 1

            - (2) [if statement] or do we sort into the unknown/low quality category given the threshold qc_cutoff (still INVALID)?
            look at index1 and index2's quality score line: line4[2] and line4[3]
            convert quality score lines to mean phred scores using bioinfo.convert_mean_phred()
            if index1 OR index2's mean phred score is below the qc_cutoff:
                write this record (all four lines) from readpair1 to unknown_r1
                write this record (all four lines) from readpair2 to unknown_r2
                unknown_counter += 1
            
            ------[else statement] after dumping the INVALID reads into 2 files, now I'm looking at VALID reads, meaning index1 and (reverse complement) index2 seq lines found in the 24 known indexes!------

            - (*) add to the possible_combinations dict if the read is VALID:
                key = index_pair
                value = 1, and increment by 1 each time I see it again
            
            - (1) [if statement] do we sort into the dual-matched category (VALID)?
            look at index1 and index2's sequence line: line2[2] and line2[3]
            if index1 seq line == rev_com_index2_seq:
                write this record (all four lines) from readpair1 to a read1 FASTQ file, for each index-pair
                write this record (all four lines) from readpair2 to a read2 FASTQ file, for each index-pair
                (48 of these total for each index-pair)
                add to the dual_matched_only dict:
                    key = index_pair
                    value = 1, and increment by 1 each time I see it again

            - (2) [if statement] or do we sort the unmatched/index hopping category (also VALID)?
            look at index1 and index2's sequence line: line2[2] and line2[3]
            if index1 seq line != rev_com_index2_seq:
                write this record (all four lines) from readpair1 to hopped_r1
                write this record (all four lines) from readpair2 to hopped_r2
                hopped_counter += 1
        
        empty the line lists after looping through one record
```

At the end, get a summary of:
* the number of read pairs with index-hopping observed,
* the number of read-pairs with unknown index(es),
* the numbers for each possible pair of indexes (both swapped and dual matched),
* the number of read-pairs with properly matched indexes (per index-pair).

Questions/more to come:
- How do I initialize and use a dictionary to help me write 52 files? Or should I use append? Any advice?
- For later: need to clarify code for writing into initialized file names

- For later: running this script efficiently: i'll write this as a python script. to run a python script in the terminal, it's `./pythonscript.py [argparse variables]` (bash language :))
    - since i don't want to run this in a login node, i'm going to make a bash script
        - in the bash script, i'll write my sbatch stuff, mamba activate an environment, then run the python script by adding ` /usr/bin/time -v ` `./pythonscript.py [argparse variables]`
    - `chmod bashscript.sh` and `chmod pythonscript.py`
    - `sbatch bashscript.sh`

After peer review, i need to:
- Answer those questions
- Add functions and unit tests to bioinfo.py

------

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

    I want to create two functions in the bioinfo.py script (and make unit tests for both in the script to make sure it's working).
    
        def reverse_complement(str):
            '''
            Return the reverse complement of a DNA sequence string.
            Example: 'GTAGCGTA' → 'TACGCTAC'
            '''
            return rev_com
        Input: AGCT
        Expected output: AGCT
        

        def convert_mean_phred(str):
            '''
            Convert a quality string into a mean Phred score.
            Assumes Phred+33 encoding.
            '''
            return mean_qscore
        Input: IIII
        Expected output: 40.0
