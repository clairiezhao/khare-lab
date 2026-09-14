import argparse
import sys

'''
Assignment 3
Convert FASTA file into tab delimited format
Output includes GI no., version no., gene name and symbol, and (opt) sequence
Tab = 4 spaces
'''
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="name of input file")
    parser.add_argument("-o", help="name of output file")
    parser.add_argument("-s", action="store_true", default=False, help="set sequence to be included in output")

    args = parser.parse_args()
    input_file = sys.stdin
    output_file = sys.stdout
    if args.i:
        try:
            input_file = open(args.i, 'r')
        except FileNotFoundError:
            input_file = sys.stdin
            print(f"{args.i} not found")
    if args.o:
        output_file = open(args.o, 'w')

    process_fasta_file(input_file, output_file, args.s)

# Given input file, output file, and boolean s = print sequence enabled
# Output input file in tab delimited format
def process_fasta_file(input_file, output_file, s):
    first_seq = True

    for line in input_file:
        line = line.strip()
        # Read record label
        if line[0] == ">":
            labels = line.split("|")
            if not first_seq:
                output_file.write("\n")
            output_file.write(f"{labels[1]}    {labels[3]}    {labels[4][1:]}")
            if s:
                output_file.write("    ")
            first_seq = False
        # Else read sequence line
        elif s:
            output_file.write(line)

if __name__ == "__main__":
    main()