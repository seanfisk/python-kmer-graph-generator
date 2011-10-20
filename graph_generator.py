import sys
import time

# generate the reverse complement of a sequence
def reverse_complement(sequence):
    rc = ""
    # iterate on string from beginning to end (I think?)
    for i in sequence[::-1]:
        if i == 'A':
            rc = rc + 'T'
        elif i == 'C':
            rc = rc + 'G'
        elif i == 'G':
            rc = rc + 'C'
        elif i == 'T':
            rc = rc + 'A'
        else:
            rc = rc + 'N'
            
    return rc

# check arguments
if len(sys.argv) < 5:
    print "Expected format: graph_generator.py <k> <input_file> <output_file> (<num_lines_to_read>)"
    exit
    
# open input file
input_file = open(sys.argv[2], "r")
# open output file
output_file = open(sys.argv[3], "w")

# initialize k from the first argument, convert to an integer
k = int(sys.argv[1])

# set num_lines_to_read to the argument if it'subsequence been given, otherwise default to -1
# num_lines_to_read denotes the amount of sequence lines to read from the file
if len(sys.argv) == 5:
    num_lines_to_read = int(sys.argv[4]);
else:
    num_lines_to_read = -1
    
# initialize kmer hash
k_mers = {}
# overlap_index = list() # seems unused

overlap_count = 0
k_mer_index = 0
read_index = 0
k_mer_count = 0

# assume all reads are same length as first read
line = input_file.readline()
# remove whitespace from beginning and end of string so length is accurate
line = line.strip()
read_length = len(line)

input_file.seek(0) # start reading at beginning of file again

# fields to indicate the length of things
# read_length - k + 1 : number of k-mer per read
output_file.write(str(num_lines_to_read) + "," + str(read_length) + "," + str(read_length - k + 1) + "\n")
# start timing
start_time = time.time()
for line in input_file:
    # remove whitespace from beginning and end of string
    line = line.strip()
    # get the reverse complement, see above
    rc = reverse_complement(line)
    
    # iterate over all k-mers in the sequence
    for i in range(read_length - k + 1):
        # grab k-mer substring from front of line
        subsequence = line[i : k + i]
        # grab k-mer substring from end_time of reverse complement
        reverse_subsequence = rc[read_length - k - i : read_length - i]
        
        # check if the k-mer overlaps with an existing k-mer
        if(subsequence[:-1] in k_mers):
            overlap_count = overlap_count + 1
            output_file.write(str(read_index) + "," + str(i) + ","  + subsequence + "," + reverse_subsequence +"," +  str(k_mers[subsequence[:-1]][0]) + "," + str(k_mers[subsequence[:-1]][1]) + "\n")
        else:
            output_file.write(str(read_index) + "," + str(i) + ","  + subsequence + "," + reverse_subsequence +","+ str(read_index) + "," + str(i) + "\n")
            
        # k-mer is not in the list
        if not (subsequence[1:] in k_mers):
            # add it to the list
            k_mers[subsequence[1:]] = [read_index, i]
            k_mers[reverse_subsequence[1:]] = [read_index, i]
            k_mer_count = k_mer_count + 2
            
    # advance indices
    read_index = read_index + 1
    num_lines_to_read = num_lines_to_read - 1
    # quit when num_lines_to_read is zero
    if(num_lines_to_read == 0):
        break
    
# stop timing
end_time = time.time()

# print timing
print "Determine overlap time:", end_time - start_time

# print results
print "Overlapping reads: ", overlap_count
print "Number of unique " +str(k) + "-mers: ", k_mer_count

# empty array
k_mers.clear()

# close files
input_file.close()
output_file.close()
