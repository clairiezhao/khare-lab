import argparse
import random

'''
Assignment 2, Part B
Translate a given DNA fragment according to a given reading frame (<0 for complement)
A negative frame corresponds to reading the complementary strand in the reverse direction
Possible reading frames: 1, 2, 3, -1, -2, -3
Translation includes * for STOP codon, X for invalid codons
'''

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", type=int, help="reading frame")
    parser.add_argument("-s", help="DNA fragment")

    # Parse command-line arguments, assume always given
    args = parser.parse_args()
    if args.f:
        if args.f == 0 or args.f > 3:
            f = 1
        elif args.f < -3:
            f = -1
        else:
            f = args.f
    else:
        print("No reading frame given")
        return
    if args.s:
        s = args.s
    else:
        print("No sequence given")
        return
    
    # dictionary
    table = init_translation_table()
    # translated sequence
    trans_s = ""

    # Frame on complementary strand
    if f < 0:
        s = complement(s)
        s = s[::-1]
    # Move start position according to given frame
    s = s[abs(f) - 1:]
    
    i = 0
    while (i+3) <= len(s):
        codon = s[i:i+3]
        if codon in table:
            trans_s += table[codon]
        else:
            trans_s += "X"
        i += 3

    print(trans_s)

    # Test: GACCGGCTCGAGTGCTACGCGCCACCCTCTCTACTACGACTAATT

# Create translation dictionary
def init_translation_table():
    table = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", # First base T
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
            "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
            "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", # First base C
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", # First base A
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", # First base G
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
    return table

# Return complement of a given DNA strand
def complement(s):
    complement = ""
    for char in s:
        upper_char = char.upper()
        # Match to bases
        if upper_char == "G":
            complement += "C"
        elif upper_char == "C":
            complement += "G"
        elif upper_char == "A":
            complement += "T"
        elif upper_char == "T":
            complement += "A"
        # Otherwise convert to X
        else:
            complement += "X"
    return complement

if __name__ == "__main__":
    main()