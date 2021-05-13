# Assembled fastq to fasta
awk '{if(NR%4==2) print $0}' N41-I14_S14_I14.assembled.fastq > N41-I14_S14.assembled.seq

# Read counts
sort N41-I14_S14.assembled.seq | uniq -dc  | sort -nr > N41-I14_S14_read_count.txt