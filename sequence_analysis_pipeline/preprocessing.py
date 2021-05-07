import subprocess
import sys
import re
import time

outputfilename = 'N41-I14_S14'
pear_directory = '/home/boydpeters/Software/pear-src-0.9.11/src/pear'
data_directory = '/sequence_analysis_pipeline/data/NGS/'

with open(data_directory + 'R1_files_processed.txt', 'r') as rf:
    lines = rf.readlines()
    R1_files_done = set()
    for line in lines:
        R1_files_done.add(line.strip())

with open(data_directory + 'R2_files_processed.txt', 'r') as rf:
    lines = rf.readlines()
    R2_files_done = set()
    for line in lines:
        R2_files_done.add(line.strip())


# retrieve all the files containing _R1 or _R2 in the file
R1_files = subprocess.check_output(
    'ls ' + data_directory + '/*_R1*', shell=True)
R2_files = subprocess.check_output(
    'ls ' + data_directory + '/*_R2*', shell=True)

R1_joined = R1_files.join('')
R2_joined = R2_files.join('')

R1_list = R1_joined.split('\n')
R2_list = R2_joined.split('\n')

# Sort the lists
R1_list.sort()
R2_list.sort()

# This is the regex for the pattern that finds the number placed next to the I,
# e.g. 'N41-I14_S15' -> 14
inum_pattern = re.compile('(?<=I)([0-9]+)')

with open('pearfrompy.sh', 'a') as bashfile:

    for i in range(len(R1_list)):
        inum_pattern
        cmdstr = pear_directory + ' -f ' + \
            R1_list[i] + ' -r ' + R2_list[i] + ' -o ' + outputfilename
