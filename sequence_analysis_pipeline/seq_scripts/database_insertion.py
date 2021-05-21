
import sqlite3
import re
import numpy as np
import matplotlib.pyplot as plt
import csv
from database_interface import DatabaseInterface
import seq_helper

# Filename is given by snakemake
inputfilename = 'fjasafakl'
#inputfiles = [snakemake.input[0], snakemake.input[1]]

# Set through config file
clvd_prefix_seq = "CTTTTCCGTATATCTCGCCAG"
clvd_prefix_name = "A"
clvd_suffix_seq = "AAAAAGAAA"
unclvd_prefix_seq = "GGGAAACAAACAAA"
unclvd_prefix_name = "W"
unclvd_suffix_seq = "AAAAAGAAA"

# Should be set through a config file
driver_round_pattern = "_R"
selection_pattern = "S_"
ligand_pattern = "L_"

driver_round = int(re.match(driver_round_pattern, inputfilename))
selection = re.match(selection_pattern, inputfilename)
ligand_present = int(re.match(ligand_pattern, inputfilename))

database_path = ":memory:"

ngs_references = seq_helper.read_ngs_references(
    "sequence_analysis_pipeline/ngs_references.csv")
ngs_references_pattern_dict = seq_helper.create_ngs_references_patterns(
    ngs_references)

# Create the database

with DatabaseInterface(database_path) as db:

    if not db.table_exists("sequences"):
        db.query("""CREATE TABLE sequences (
                    reads INTEGER,
                    original_sequence TEXT,
                    cleaned_sequence TEXT,
                    barcode TEXT,
                    prefix_name TEXT,
                    cleaved_prefix INTEGER,
                    reference_name TEXT,
                    selection TEXT,
                    round INTEGER,
                    ligand_present INTEGER,
                    cleavage_fraction REAL,
                    fold_change REAL,
                    sensor INTEGER
                    )""")
    else:
        print("Table already exists, check if files are already processed.")
        while True:
            user_input = input(
                "Do you want to continue filling the database? (Y/n)\t")
            if user_input.lower().startswith('y'):
                print("Database filling up...")
            elif user_input.lower().startswith('n'):
                exit()

with DatabaseInterface(database_path) as db:

    for inputfile in inputfiles:

        driver_round = int(re.match(driver_round_pattern, inputfilename))
        selection = re.match(selection_pattern, inputfilename)
        ligand_present = int(re.match(ligand_pattern, inputfilename))

        with open(inputfile) as rf:
            lines = rf.readlines()
            for line in lines:
                sequence = line.strip()

                # First check whether sequence is a reference sequence
                reference_name = seq_helper.is_reference_seq(
                    sequence, ngs_references_pattern_dict)

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

                # fig, ax = plt.subplots()
                # ax.plot(wrong_array)
                # plt.show()
                # plt.savefig('test')
