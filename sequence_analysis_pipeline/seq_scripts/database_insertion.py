import sqlite3
import re
import numpy as np
import matplotlib.pyplot as plt
import csv
from database_interface import DatabaseInterface
import seq_helper

# Filename is given by snakemake
inputfiles = [
    'sequence_analysis_pipeline/data/NGS/processed/N41-I14_S14_read_count.txt']
#inputfiles = [snakemake.input[0], snakemake.input[1]]

# Set through config file
clvd_prefix_seq = "CTTTTCCGTATATCTCGCCAG"
clvd_prefix_name = "A"
clvd_suffix_seq = "AAAAAGAAA"

unclvd_prefix_seq = "GGGAAACAAACAAA"
unclvd_prefix_name = "W"
unclvd_suffix_seq = "AAAAAGAAA"

clvd_prefix_info = {"seq": clvd_prefix_seq, "name": clvd_prefix_name}
unclvd_prefix_info = {"seq": unclvd_prefix_seq, "name": unclvd_prefix_name}


# Should be set through a config file
driver_round_pattern = "_R"
selection_pattern = "S_"
ligand_pattern = "L_"

# database_path = ":memory:"
database_path = "sequence_analysis_pipeline/data/NGS/processed/N41-I14_S14_database.db"

ngs_references = seq_helper.read_ngs_references(
    "sequence_analysis_pipeline/ngs_references.csv")
ngs_references_pattern_dict = seq_helper.create_ngs_references_patterns(
    ngs_references)

# Create the database
with DatabaseInterface(path=database_path) as db:

    if not db.table_exists("sequences"):
        db.query("""CREATE TABLE sequences (
                    read_count INTEGER,
                    original_sequence TEXT,
                    cleaned_sequence TEXT,
                    barcode TEXT,
                    cleaved_prefix INTEGER,
                    prefix_name TEXT,
                    reference_name TEXT,
                    selection TEXT,
                    driver_round INTEGER,
                    ligand_present INTEGER,
                    cleavage_fraction REAL,
                    fold_change REAL,
                    possible_sensor INTEGER
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

with DatabaseInterface(path=database_path) as db:
    for inputfile in inputfiles:

        # driver_round = int(re.match(driver_round_pattern, inputfilename))
        # selection = re.match(selection_pattern, inputfilename)
        # ligand_present = int(re.match(ligand_pattern, inputfilename))

        driver_round = 80
        selection = "S4"
        ligand_present = 0

        with open(inputfile) as rf:
            lines = rf.readlines()
            for line in lines:
                # Create a dictionary and store all general information for an unique sequence
                sequence_info = {"driver_round": driver_round, "selection": selection, "ligand_present": ligand_present,
                                 "cleavage_fraction": "NULL", "fold_change": "NULL", "possible_sensor": 0}

                read_count, sequence = line.strip().split()
                sequence_info["read_count"] = read_count
                sequence_info["original_sequence"] = sequence

                # Determine whether sequence is a reference sequence
                sequence_info["reference_name"] = seq_helper.reference_seq(
                    sequence, ngs_references_pattern_dict)

                clvd_prefix, prefix_name, prefix = seq_helper.determine_clvd_prefix(
                    sequence, clvd_prefix_info=clvd_prefix_info, unclvd_prefix_info=unclvd_prefix_info)

                sequence_info["cleaved_prefix"] = clvd_prefix
                sequence_info["prefix_name"] = prefix_name

                sequence_info["barcode"] = seq_helper.retrieve_barcode(
                    sequence, prefix)
                sequence_info["cleaned_sequence"] = seq_helper.cleanup_sequence(
                    sequence, prefix, clvd_suffix_seq)

                db.query("""INSERT INTO sequences VALUES (
                            :read_count, :original_sequence, :cleaned_sequence,
                            :barcode, :cleaved_prefix, :prefix_name,:reference_name,
                            :selection, :driver_round, :ligand_present, :cleavage_fraction,
                            :fold_change, :possible_sensor)""", parameters=sequence_info)

                # # A database is created with the following columns:
                # # reads: the number of reads, original_sequence: orignal sequence (TEXT)
                # # cleaned_sequence: DNA sequence with barcode prefix and suffix removed (TEXT),
                # # barcode: barcode of the sequence (TEXT), cleaved_prefix: yes(1)/no(0) (INTEGER),
                # # cleaved_suffix: yes(1)/no(0) (INTEGER), reference: indicates whether sequence is a
                # # reference sequence yes(1)/no(0) (INTEGER), round: round when the sequence was
                # # sequenced(INTEGER), ligand: ligand present yes(1)/no(0) (INTEGER),
                # # sensor: indicates whether sequence is a possible biosensor yes(1)/no(0) (INTEGER)


with DatabaseInterface(path=database_path) as db:
    results = db.get(
        "sequences", ["original_sequence", "cleaned_sequence"], limit=10)
    print(results)
