import argparse
import random

'''
Assignment 2, Part A
Generate random coding DNA given the start codon ATG and
1) the stop codon (coding cannot contain internal stop codons)
2) the sequence length
'''

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", type=int, help="length of the sequence")
    parser.add_argument("-s", help="stop codon")
    parser.add_argument("-l", help="display the sequence in lowercase")

    bases = ["T", "C", "A", "G"]
    seq = "ATG"
    # Initialize random length of sequence between 25-50 (incl)
    len = random.randint(25, 50)
    # Initialize random stop codon
    stop = ""
    for i in range(3):
        stop += bases[random.randint(0, 3)]

    # Parse command-line arguments
    args = parser.parse_args()
    if args.c:
        len = args.c
        if len < 6:
            print("Invalid length")
    if args.s:
        s = args.s.upper()
        if valid(s, bases):
            stop = s
        else:
            print("Invalid stop codon")

    # Keep track of last two generated bases in the sequence
    prev_bases = "TG"
    # If last two bases match the first two bases of the stop codon, remove last stop codon base from generation
    bases1 = ["T", "C", "A", "G"]
    bases1.remove(stop[2])

    for i in range(len - 6):
        base = bases[random.randint(0, 3)]
        if prev_bases == stop[0:2]:
            base = bases1[random.randint(0, 2)]
        seq += base
        prev_bases = prev_bases[1] + base
    # Add stop codon
    seq += stop
    if args.l:
        seq = seq.lower()
    print(seq)

# Returns true if a given uppercase stop codon is valid
def valid(s, bases):
    if len(s) == 3:
        for i in range(3):
            if not (s[i] in bases):
                return False
    else:
        return False
    return True

if __name__ == "__main__":
    main()