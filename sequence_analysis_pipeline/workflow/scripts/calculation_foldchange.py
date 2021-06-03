import sqlite3
from database_interface import DatabaseInterfaceCleanSequences
from tqdm import tqdm

# add right path for the table with the clean sequences
database_path = "sequence_analysis_pipeline/data/NGS/processed/S1_D80_database.db"
# database_path = snakemake.input[0]

TABLE_NAME = "clean_sequences"

# 1 - calculate the fold-change for k=1

# temporary storage for the calculated fold-change with k=1
fc_list = []

with DatabaseInterfaceCleanSequences(path=database_path) as db:
    # get the sequences with ligand present and uncleaved prefix (those are all unique sequences) and for biosensors, uncleaved is definitaly present in +round
    sequences_data = db.get_sequences(table=TABLE_NAME, cleaved_prefix=0)

    for i in tqdm(range(len(sequences_data))):
        k = 1
        seq = sequences_data[i][2]         # get the sequence
        cs_pos = sequences_data[i][9]      # get cleavage fraction

        # get info of the same sequence but then in the ligand NOT present round, with cleaved_prefix, cause for biosensors, cleaved is definitaly present in -round
        sequences_data_neg = db.get_sequence_negligand(table=TABLE_NAME, cleaned_sequence=seq)
        cs_neg = sequences_data_neg[0][9]  # get cleavage fraction

        # calculate fold-change with k=1
        fc_temp = k * ((1-cs_pos) / (1-cs_neg))
        fc_list.append(fc_temp)


# 2 - Find median value
# sort the list with temporary fold change values
fc_list.sort()

# take the median value
length = int(len(fc_list))
if (length % 2) == 0:  # even number of sequences
    ind = length / 2
    median = (fc_list[ind] + fc_list[ind+1]) / 2

else:
    ind = int(round(length/2))
    median = fc_list[ind]  


# 3 - Calculate new k with median fold-change value
k_new = 1 / median
with DatabaseInterfaceCleanSequences(path=database_path) as db:
    for i in tqdm(range(len(sequences_data))):
        seq = sequences_data[i][2]            # get the sequence
        cs_pos = sequences_data[i][9]         # get cleavage fraction
        seq_ID = sequences_data[i][0]

        # get info of the same sequence but then in the ligand NOT present round, with cleaved_prefix, cause for biosensors, cleaved is definitaly present in -round
        sequences_data_neg = db.get_sequence_negligand(table=TABLE_NAME, cleaned_sequence=seq)
        cs_neg = sequences_data_neg[0][9]     # get cleavage fraction
        seq_ID_neg = sequences_data_neg[0][0] # get sequence ID of negative round

        # calculate fold-change with k_new
        fc = k_new * ((1-cs_pos) / (1-cs_neg))
        
        # store the fold change in the database
        db.update_fold_change(table=TABLE_NAME, rowid=seq_ID, fold_change=fc)
        db.update_fold_change(table=TABLE_NAME, rowid=seq_ID_neg, fold_change=fc)

        
    