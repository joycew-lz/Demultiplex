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
    R1 and R4 are a read pair while R2 and R3 are their respective indexes, so their index pair. Each group of four lines from all four files represents one read-pair record with their indexes.
    We must keep in mind that R3 is the reverse compliment to R2. 
    We are also given 24 known/dual-matched indexes in a text file. These are the expected, correctly paired barcodes.
    ----
    We have to demultiplex the data by using the barcode information in order to know which sequences came from which samples after they were all sequenced together.
    To demultiplex the data, I need to go through every read (record) of the four files and demultiplex by index-pair (R2 and R3).
    * If R2 and R3 index sequences are not in the list of known indexes, or contain "N", they are categorized as "unknown".
    * If R2 and R3 index sequences are in the list of known indexes and match each other, they are categorized as "dual-matched".
    * If R2 and R3 index sequences do not match each other, they are categorized as "index-hopping". 


2. Describe output:

    FOR EACH RECORD, look at the two index reads (index-pairs) R2 and R3. SORT into one of the three categories below.

    * If one or both index sequences (the barcodes) includes "N's" OR do not match the 24 known indexes (unknown/low quality category):
        * Write out the record (as in the record of R1 and R4, the read pair), with modified headers, into one read1 FASTQ file and one read2 FASTQ file.
        * 2 FASTQ files
        * ALSO, if either of the index reads' mean quality score is lower than the threshold determined in part 1, we'll put it in the same(unknown/low quality category).
    ----
    * If the index-pair's sequences match one of the 24 known indexes and each other (consider R3 is reverse complement) (dual-matched):
        * Write out the record (as in the record of the read pair), with modified headers, into one read1 FASTQ file and one read2 FASTQ file **per matching index-pair**.
            * named like: A1_R1.fq, A1_R2.fq (with A1 being the known barcode "nickname")
        * 48 FASTQ files 
    * If the index-pair's sequences DO NOT match each other (index-hopping):
        * Write out the record (as in the record of the read pair), with modified headers, into one read1 FASTQ file and one read2 FASTQ file.
        * 2 FASTQ files

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
            
            - (*) modify the header for everything:
            look at index1 and index2's sequence line: line2[2] and line2[3]
            grab the sequences from index1 and index2 to get the index-pair:
                index1_seq = line2[2].strip()  # reminder to keep using strip
                index2_seq = line2[3].strip()
                rev_com_index2_seq = reverse_complement(index2_seq)
                index_pair = f"{index1_seq}-{rev_com_index2_seq}"
            append index_pair to the end of readpair1 and readpair 2's headers in line1: line1[0] and line1[1]
                line1[0] = line1[0].strip() + f" {index_pair}\n" 
                line1[1] = line1[1].strip() + f" {index_pair}\n"

            - (1) do we sort into the unknown/low quality category?:
            look at index1 and index2's sequence line: line2[2] and line2[3]
            if index1 seq line or index2 seq line NOT found in the 24 known matches OR contains "N":
                write this record (all four lines) from readpair1 to a read1 FASTQ file for all unknown/low-qual
                write this record (all four lines) from readpair2 to a read2 FASTQ file for all unknown/low-qual

            - (2) or do we sort into the unknown/low quality category given the threshold qc_cutoff?
            look at index1 and index2's quality score line: line4[2] and line4[3]
            convert quality score lines to mean phred scores using bioinfo.convert_mean_phred()
            if index1 OR index2's mean phred score is below the qc_cutoff:
                write this record (all four lines) from readpair1 to a read1 FASTQ file for all unknown/low-qual
                write this record (all four lines) from readpair2 to a read2 FASTQ file for all unknown/low-qual
            
            - (3) or do we sort into the dual-matched category?
            look at index1 and index2's sequence line: line2[2] and line2[3]
            if index1 seq line == reverse_complement(index2 seq line) AND index1 and index2 seq lines found in the 24 known indexes:
                write this record (all four lines) from readpair1 to a read1 FASTQ file, for each index-pair
                write this record (all four lines) from readpair2 to a read2 FASTQ file, for each index-pair
                (48 of these total for each index-pair)

            - (4) or do we sort the unmatched/index hopping category?
            look at index1 and index2's sequence line: line2[2] and line2[3]
            if index1 seq line != reverse_complement(index2 seq line) 
                write this record (all four lines) from readpair1 to a read1 FASTQ file for all unmatched
                write this record (all four lines) from readpair2 to a read2 FASTQ file for all unmatched
        
        empty the lists after looping through one record
```

```
At the end, write a summary:
    * the number of read-pairs with properly matched indexes (per index-pair),
    * the number of read pairs with index-hopping observed,
    * the numbers for each possible pair of indexes (both swapped and dual matched) (in a table),
    * the number of read-pairs with unknown index(es).
```
questions/more to come:
- is this method of storing each line of a record into a list, pulling what I need out the list through indexing, and emptying the list after looping through each four-line record feasible?
- how do i initialize and use a dictionary to help me write 52 files?
- this will be a python script, i can srun for the test files; sbatch for the actual files (should I create a new environment and activate this environment in my script)?

after peer review, i need to:
- answer those questions
- add functions to bioinfo.py and write unit tests for them

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
