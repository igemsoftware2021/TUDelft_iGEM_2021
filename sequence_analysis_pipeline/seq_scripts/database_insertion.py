
import sqlite3
import re
import numpy as np
import matplotlib.pyplot as plt
import csv


def read_ngs_references(path):
    ngs_references = {}
    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ngs_references[row["sequence"]] = row["name"]
    return ngs_references


# def is_reference(sequence, references):
#     patterns = []
#     for key in references:
#         patterns.append(re.compile(key))

#     for pat in patterns:
#         match = pat.search(sequence)
#         if bool(match):
#             return (True, sequence, match)
#     return (False,)


ngs_references = read_ngs_references(
    "sequence_analysis_pipeline/ngs_references.csv")

# print(ngs_references)

# inputfilename = snakemake.input[0]
inputfilename = "sequence_analysis_pipeline/data/NGS/processed/N41-I14_S14_read_count.txt"

clvd_prefix_pattern = re.compile("CTTTTCCGTATATCTCGCCAG")  # A
clvd_suffix_pattern = re.compile("AAAAAGAAACAGTC")
unclvd_prefix_pattern = re.compile("GGGAAACAAACAAA")  # W
unclvd_suffix_pattern = re.compile("AAAAAGAAACAGTC")

wrong_array = np.zeros(shape=(913709), dtype=np.int32)
count_array = np.arange(0, 913709, step=1)

with open(inputfilename) as rf:
    count = 0
    wrong = 0
    lines = rf.readlines()
    for line in lines:
        count += 1
        sequence = line.strip()
        # result = is_reference(sequence, ngs_references)
        # if result[0]:
        #     print(result[1], result[2])
        prefix_match = clvd_prefix_pattern.search(sequence)
        # print(int(bool(prefix_match)))
        if bool(prefix_match):
            clvd_prefix = 1
        else:
            prefix_match = unclvd_prefix_pattern.search(sequence)
            if bool(prefix_match):
                clvd_prefix = 0
            else:
                # clvd_prefix = 2
                # print("NOOOO CLVD OR UNCLVD PREFIX")
                # print(sequence)
                with open("sequence_analysis_pipeline/data/NGS/processed/bad_N41-I14_S14_read_count.txt", 'a+') as wf:
                    wf.write(sequence + '\n')
                wrong += 1
        # if count == 50000:
        #     break
        wrong_array[count] = wrong
        # print(clvd_prefix)
    print(wrong)


# conn = sqlite3.connect(':memory:')

# c = conn.cursor()

# # A database is created with the following columns:
# # reads: the number of reads, original_sequence: orignal sequence (TEXT)
# # cleaned_sequence: DNA sequence with barcode prefix and suffix removed (TEXT),
# # barcode: barcode of the sequence (TEXT), cleaved_prefix: yes(1)/no(0) (INTEGER),
# # cleaved_suffix: yes(1)/no(0) (INTEGER), reference: indicates whether sequence is a
# # reference sequence yes(1)/no(0) (INTEGER), round: round when the sequence was
# # sequenced(INTEGER), ligand: ligand present yes(1)/no(0) (INTEGER),
# # sensor: indicates whether sequence is a possible biosensor yes(1)/no(0) (INTEGER)

# c.execute("""CREATE TABLE sequences (
#             reads INTEGER,
#             original_sequence TEXT,
#             cleaned_sequence TEXT,
#             barcode TEXT,
#             cleaved_prefix INTEGER,
#             cleaved_suffix INTEGER,
#             reference INTEGER,
#             selection TEXT,
#             round INTEGER,
#             ligand_present INTEGER,
#             cleavage_fraction REAL,
#             fold_change REAL,
#             sensor INTEGER
#             )""")

fig, ax = plt.subplots()
ax.plot(wrong_array)
plt.show()
# plt.savefig('test')
