configfile: "config.yaml"

# Example
# inputfiles:
# "S1_D80_l0_read_counts.txt"
# "S1_D80_l1_read_counts.txt"
# Then dataset = S1_D80 and sample = L0/L1

rule all:
    pass


# rule calculate_cleavage_fraction:
#     input:
#         "/data/processed/{dataset}_database.db"
#     output:
#         "/data/processed/{dataset}_database.db"

rule insert_counts_into_database:
    input:
        "/data/processed/{dataset}_L0_read_counts.txt"
        "/data/processed/{dataset}_L1_read_counts.txt"
    output:
        "/data/processed/{dataset}_database.db"
    script:
        "/seq_scripts/database_insertion.py"