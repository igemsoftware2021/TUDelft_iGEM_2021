import yaml
import regex
from tqdm import tqdm
from database_interface import DatabaseInterfaceCleanSequences


database_path = "sequence_analysis_pipeline/data/NGS/processed/S1_D80_database.db"
# database_path = snakemake.output[0]


with DatabaseInterface(path=database_path) as db:
    # Fetchall rows only columns (cleaned_sequence, prefix_name, ligand_present)
    pre_data_seq = db.get(table="raw_sequences", columns=["cleaned_sequence", "cleaved_prefix", "ligand_present"])
    
    # Create a set of the list of tuples of the retrieved sequence data
    unique_rows = list(set(pre_data_seq))


with DatabaseInterfaceCleanSequences(path=database_path) as db:



# Go over all the tuples in the set, query the data, add the read counts and put all info
# in new table.
# Then do the cleavage fraction and stuff