'''
Assignment 1
Calculate and display 
1) the GC-content of an entered DNA fragment
2) its Watson-Crick complement
3) the reverse of the complement.
'''

def main():
    seq = input("Enter a sequence: ")
    complement = ""
    gc_count = 0
    at_count = 0

    for char in seq:
        # Convert char to uppercase
        upper_char = char.upper()
        # Match to bases
        if upper_char == "G":
            gc_count += 1
            complement += "C"
        elif upper_char == "C":
            gc_count += 1
            complement += "G"
        elif upper_char == "A":
            at_count += 1
            complement += "T"
        elif upper_char == "T":
            at_count += 1
            complement += "A"
        # Otherwise convert to X
        else:
            complement += "X"
    
    gc_content = float(gc_count) / (gc_count + at_count) * 100
    print("%GC content: " + f"{gc_content:.2f}%")
    print("Complement: " + complement)
    print("Reverse complement: " + complement[::-1])

    # Test: GACCGGCTCGAGTGCTACGCGCCACCCTCTCTACTACGACTAATT


if __name__ == '__main__':
    main()