import os
import subprocess
import sys
import re
import time

outputfilename = 'N41-I14_S14'
pear_directory = '/home/boydpeters/pear-0.9.11-linux-x86_64/bin/pear'
data_directory = os.path.join('sequence_analysis_pipeline', 'data', 'NGS')


with open(os.path.join(data_directory, 'R1_files_processed.txt'), 'r') as rf:
    lines = rf.readlines()
    R1_files_processed = set()
    for line in lines:
        R1_files_processed.add(line.strip())

with open(os.path.join(data_directory, 'R2_files_processed.txt'), 'r') as rf:
    R2_files_processed = set()
    for line in lines:
        R2_files_processed.add(line.strip())

print('check')

# retrieve all the files containing _R1 or _R2 in the file and decode the bytes into a string
R1_files = subprocess.check_output(
    'ls ' + data_directory + '/*_R1*', shell=True).decode()
R2_files = subprocess.check_output(
    'ls ' + data_directory + '/*_R2*', shell=True).decode()

R1_set = set(R1_files.split('\n'))
R2_set = set(R2_files.split('\n'))

# Remove files already processed
R1_uniq = R1_set.difference(R1_files_processed)
R2_uniq = R2_set.difference(R2_files_processed)

# Change it back to lists
R1_list = list(R1_uniq)
R2_list = list(R2_uniq)

# Remove empty strings from the lists
R1_list = [R1_file for R1_file in R1_list if R1_file != '']
R2_list = [R2_file for R2_file in R2_list if R2_file != '']

# Sort the lists
R1_list.sort()
R2_list.sort()

# This is the regex for the pattern that finds the number placed next to the I,
# e.g. 'N41-I14_S15' -> 14
Inum_pattern = re.compile('(?<=I)([0-9]+)')

bashfile_directory = os.path.join(
    'sequence_analysis_pipeline', 'pearfrompy.sh')

with open(bashfile_directory, 'a') as bashfile:
    for i in range(len(R1_list)):

        if i < len(R2_list):
            # Pattern.search() returns a match object, so grab the first match using 'group(0)'
            R1_Inum = Inum_pattern.search(R1_list[i]).group(0)
            R2_Inum = Inum_pattern.search(R2_list[i]).group(0)

            if R1_Inum == R2_Inum:
                cmdstr = pear_directory + ' -f ' + \
                    R1_list[i] + ' -r ' + R2_list[i] + ' -o ' + outputfilename + '_I' + \
                    R1_Inum + ' > ' + outputfilename + '_I' + R1_Inum + '_pear.log'
                bashfile.write(cmdstr + '\n')
            else:
                print(
                    f'R1_Inum = {R1_Inum} and R2_Inum = {R2_Inum} -> R1_Inum != R2_Inum\nThe R1_filename is {R1_list[i]}')
        else:
            print(
                'There is probably a file missing, the length of R1_list is larger than the length of R2_list')

print('check 1')
# Make the bashfile an executable with 'chmod +x'
subprocess.call(
    'chmod +x ' + bashfile_directory, shell=True)
subprocess.call(bashfile_directory, shell=True)
print('check2')

# Maybe considering NGmerge instead of PEAR, since NGmerge is open-source software?
