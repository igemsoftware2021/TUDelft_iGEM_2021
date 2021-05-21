
import sqlite3
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

print(matplotlib.get_backend())

# inputfilename = snakemake.input[0]
inputfilename = "/home/minke/iGEM_2021/TUDelft_iGEM_2021/sequence_analysis_pipeline/data/test_prefix/N38-I1_S1.txt"

clvd_prefix_pattern = re.compile("AGATCTTTTCCGTATATCTCGCCAG")
clvd_suffix_pattern = re.compile("AAAAAGAAA")
unclvd_prefix_pattern = re.compile("AGATGGGAAACAAACAAA")
unclvd_suffix_pattern = re.compile("AAAAAGAAA")

wrong_array = np.zeros(shape=(3277315), dtype=np.int32)
count_array = np.arange(0, 3277315, step=1)

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
                #print("NOOOO CLVD OR UNCLVD PREFIX")
                print(sequence)
                wrong += 1
        # if count == 20000:
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
