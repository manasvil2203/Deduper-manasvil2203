#!/usr/bin/env python

import argparse
import re
import gzip


#Creating my varuables using Argparse
def get_args():
    parser = argparse.ArgumentParser(description= "Deduplicate a SAM (or SAM.gz) by (chromosome, 5' position, strand, UMI). INPUT MUST BE SORTED by RNAME then POS.")
    parser.add_argument("-f", "--file", help="Specify the filename of the sorted sam file", required=True)
    parser.add_argument("-o", "--output", help="Specify the filename of the output sam file", required=True)
    parser.add_argument("-u", "--umi", help="Specify the valid umi list", required=True)
    parser.add_argument("-s", "--summary", help="Specify the filename for the summary counts", required=True)
    

    return parser.parse_args()  

args = get_args()

#Making my variable global
f: str = args.file
o: str = args.output
u: str = args.umi
s: str  = args.summary



# initialize counters
#Counts the total number of header
header_count:int = 0
# Counts the total number of umis
wrong_umi_count:int = 0
# Count the number of pcr duplicate
duplicate_count:int = 0
# Counts the number of unique reads
unique_count:int = 0
# Counts the total number of reads
total_reads:int = 0

# dictionary to count reads per chromosome
chrom_counts:dict = {}

#Initialize an empty set 
umi_set :set = set()


# Open the umi list
with open(u, "r") as file:
    # For every line in the file
    for line in file:
        #Strip the line of the new line character
        item = line.strip()
        # Ensure it is not an empty line
        if item:  
            # Add the umi to the set
            umi_set.add(item)
       

#print(umi_set)

def determine_strand(flag:int) -> str:
    '''This function is taking in the flag and giving us the strand'''        

    #Determine the strand using the the flag
    if ((flag & 16) == 16):
        strand = "neg"
    else:
        strand = "pos"

    #print(strand)
    return strand
# Testing the function
#determine_strand(82)

def compute_five_prime(strand: str, pos: int, cigar: str) -> int:
    """
    This function figures out where the first base of the read (in the read's
    own 5' to 3' direction) maps on the reference genome. It accounts for the strand
    orientation, how far the alignment extends along the reference, and whether part of
    the read was soft-clipped off the start or end.
    """ 
    # Use regex to make a list of tuple of each number and character
    parts = re.findall(r"(\d+)([MIDNSHP])", cigar)
    #print(parts)

    # Initialize variable for leading S
    s_left = 0
    # Initialize variable for trailing S
    s_right = 0
    # Initialize variable for calculating refernce length 
    ref_len = 0

    # If the second element of the first tuple is S
    if parts[0][1] == "S":
        # Then s_left would be the number preceding S
        s_left = int(parts[0][0])
        #print(s_left)

    #If the second element of the last tuple is S
    if parts[-1][1] == "S":
        # Then s_right would be the number precding S
        s_right = int(parts[-1][0])
        #print(s_right)

    # For each tuple(length, op)
    for length, op in parts:
        # If the character is either M, D, N
        if op in ("M", "D", "N"):
            # Then add the preceding integer to ref_length
            ref_len += int(length)
            #print(ref_len)


    # If the strand is positive 
    if strand =="pos":
        # Then 5 prime is leftmost position - leading soft clip
        five_prime = pos - s_left
    # Else if it is negative
    else:
        # Then 5 prime is left most position + reference length - 1(technicaly right) + trailing soft clip
        five_prime = pos + ref_len - 1 + s_right

    #print(five_prime)
    return five_prime

#compute_five_prime("neg", 34, "6S5M6N7S")
def open_maybe_gzip(path: str, mode: str = "rt"):
    """Open a file, using gzip.open if it is compressed (.gz), otherwise open()."""
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

# Temp variable to store current chrom number
current_chrom = None

# Keeps track of UMI, strand and 5 prime
tracker_set: set = set()

#Open an output file
with open(o, "w") as out:
    #Parse through every line in file
    with open_maybe_gzip(f, "rt") as file:
        # For every line in the file
        for line in file: 
            # if the line starts with an @
            if line.startswith("@"):
                #Write the line into the output
                out.write(line)
                # Count the header line
                header_count += 1
                # Keep going
                continue

            #Count the totsal number of reads
            total_reads += 1
            # Grab the line that is not a header and strip and split
            bits = line.strip().split()
            #print(bits)
            # Get the umi from the first element of bits, strip it and grab 7th element
            umi = bits[0].split(":")[7]
            #print(umi)
                # If the umi is not in the valid umi set
            if umi not in umi_set:
                # increment counter
                wrong_umi_count += 1
                # then keep going
                continue

            # Grab the chromosome number from col 3
            chrom = bits[2]
            # print(chrom)
            # break

            # Grab from coloumn 2
            flag =  int(bits[1])
            # print(flag)
            # break

            # Use the determine_strand function to determine the strand
            strand = determine_strand(flag)
            # print(strand)
            # break

            # Grab coloumn 4 for left most start position 
            pos = int(bits[3])
            # print(pos)
            # break

            # Grab coloumn 6 to determine 5' position and ref length
            cigar = bits[5] 
            # print(cigar)
            # break

            # Use compute_five_prime function to determine 5' position
            five_prime = compute_five_prime(strand, pos, cigar)
            # print(five_prime)
            # break


            # If there is a new chromosome
            if chrom != current_chrom:
                #Empty the set
                tracker_set: set = set()
                #current_chrom will be now set to chromosome
                current_chrom = chrom
            

            #Create a variable to store a tuple of the 5', strand, UMI
            mini_key:tuple = (five_prime, strand, umi)

            # If the mini key is not in my set
            if  mini_key not in tracker_set:
                # Write out the line in the output file 
                out.write(line)
                #And add the mini_key into my set
                tracker_set.add(mini_key)
                # Increment counter of unique reads
                unique_count += 1
                # Check if this chromosome has been seen before
                if chrom in chrom_counts:
                    # If yes, increment its count by 1
                    chrom_counts[chrom] += 1
                else:
                    # If not, initialize its count to 1
                    chrom_counts[chrom] = 1
            # Else increment the counter of number of pcr duplicates
            else:
                duplicate_count += 1
            

with open(s, "w") as sum:
    sum.write("\n===== SUMMARY REPORT =====\n")
    sum.write(f"Header lines: {header_count}\n")
    sum.write(f"Total reads processed: {total_reads}\n")
    sum.write(f"Unique reads kept: {unique_count}\n")
    sum.write(f"Wrong UMIs skipped: {wrong_umi_count}\n")
    sum.write(f"PCR duplicates removed: {duplicate_count}\n")
    sum.write("\nReads per chromosome:\n")
    for chrom, count in chrom_counts.items():
        sum.write(f"{chrom}\t{count}\n")
    sum.write("==========================\n")