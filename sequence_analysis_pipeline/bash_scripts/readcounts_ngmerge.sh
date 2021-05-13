# Assembled fastq to fasta
awk '{if(NR%4==2) print $0}' N41-I14_S14_I14 > N41-I14_S14_I14.seq

# Read counts
sort N41-I14_S14_I14.seq | uniq -dc  | sort -nr > N41-I14_S14_I14_ngmerge_read_count.txt