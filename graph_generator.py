import sys
import time

def reverse_comp(s):
    rc = ""
    for i in s[::-1]:
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

if len(sys.argv) < 5:
    print "Expected format: graph_generator.py <k> <input_file> <output_file> (<read_count>)"
    exit

fr = open(sys.argv[2], "r")
fw = open(sys.argv[3], "w")

k = int(sys.argv[1])

if len(sys.argv) == 5:
    read_count = int(sys.argv[4]);
else:
    read_count = -1

k_mers = {}
overlap_idx = list()

overlap_count = 0
k_mer_idx = 0
read_idx = 0
k_mer_count = 0

# assume all reads are same length as first read
line = fr.readline()
line = line.strip()
r_len = len(line)

fr.seek(0) # start reading at beginning of file again

fw.write(str(read_count) + "," + str(r_len) + "," + str(r_len - k + 1) + "\n")
start = time.time()
for line in fr:
    line = line.strip()
    rc = reverse_comp(line)

    for i in range(r_len - k + 1):
        s = line[i : k + i]
        rs = rc[r_len - k - i : r_len - i]
        
        if(s[:-1] in k_mers):
            overlap_count = overlap_count + 1
            fw.write(str(read_idx) + "," + str(i) + ","  + s + "," + rs +"," +  str(k_mers[s[:-1]][0]) + "," + str(k_mers[s[:-1]][1]) + "\n")
        else:
             fw.write(str(read_idx) + "," + str(i) + ","  + s + "," + rs +","+ str(read_idx) + "," + str(i) + "\n")

        if not (s[1:] in k_mers):
            k_mers[s[1:]] = [read_idx, i]
            k_mers[rs[1:]] = [read_idx, i]
            k_mer_count = k_mer_count + 2
           
            
    read_idx = read_idx + 1
    read_count = read_count - 1
    if(read_count == 0):
        break

end = time.time()

print "Determine overlap time:", end - start

print "Overlapping reads: ", overlap_count
print "Number of unique " +str(k) + "-mers: ", k_mer_count

k_mers.clear()

fr.close()
fw.close()
