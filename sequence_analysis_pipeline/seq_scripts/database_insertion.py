
import sqlite3
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

print(matplotlib.get_backend())

# inputfilename = snakemake.input[0]
inputfilename = "sequence_analysis_pipeline/data/NGS/processed/N41-I14_S14_ngmerge_read_count.txt"

clvd_prefix_pattern = re.compile("CTTTTCCGTATATCTCGCCAG")
clvd_suffix_pattern = re.compile("AAAAAGAAACAGTC")
unclvd_prefix_pattern = re.compile("GGGAAACAAACAAA")
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
                wrong += 1
        # if count == 20000:
        #     break
        wrong_array[count] = wrong
        # print(clvd_prefix)
    print(wrong)


# conn = sqlite3.connect(':memory:')

# c = conn.cursor()

# # A database is created with the following columns:
# # sequence: DNA sequence with barcode prefix and suffix removed (TEXT),
# # barcode: barcode of the sequence (TEXT), cleaved_prefix: yes(1)/no(0) (INTEGER),
# # cleaved_suffix: yes(1)/no(0) (INTEGER), reference: indicates whether sequence is a
# # reference sequence yes(1)/no(0) (INTEGER), round: round when the sequence was
# # sequenced(INTEGER), ligand: ligand present yes(1)/no(0) (INTEGER),
# # sensor: indicates whether sequence is a possible biosensor yes(1)/no(0) (INTEGER)

# c.execute("""CREATE TABLE sequences (
#             sequence TEXT,
#             barcode TEXT,
#             cleaved_prefix INTEGER,
#             cleaved_suffix INTEGER,
#             reference INTEGER,
#             round INTEGER,
#             ligand INTEGER,
#             sensor INTEGER
#             )""")

fig, ax = plt.subplots()
ax.plot(wrong_array)
plt.show()
# plt.savefig('test')
