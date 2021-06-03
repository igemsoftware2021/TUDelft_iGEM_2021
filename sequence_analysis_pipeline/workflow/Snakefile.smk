configfile: "config/config.yaml"

# Example
# inputfiles:
# "S1_D80_l0_read_counts.txt"
# "S1_D80_l1_read_counts.txt"
# Then dataset = S1_D80 and sample = L0/L1

rule all:
    input:


rule fastqc:
    input:
        "../../data/{sample}.fastq"
    output:
        "../../results/fastqc_reports/{sample}_fastqc.html"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        #!/bin/bash
        fastqc {input} --outdir={output}
        """

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
        "/result/databases/{dataset}_database.db"
    script:
        "/seq_scripts/database_insertion.py"