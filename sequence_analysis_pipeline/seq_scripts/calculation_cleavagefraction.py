import sqlite3
from database_interface import DatabaseInterfaceSequences

database_path = "sequence_analysis_pipeline/data/NGS/processed/S1_D80_database.db"
# database_path = snakemake.input[0]

with DatabaseInterfaceSequences(path=database_path) as db:

    # 1 - Get number of reads of the references with cleaved prefix in +ligand round
    clv_ref_info = db.get_ref_sequences()
    r_clv_ref = sum([s[1] for s in clv_ref_info])      #sums up second column (read_counts)

    # 2 - Get number of reads of the references with uncleaved prefix in +ligand round
    unclv_ref_info = db.get_ref_sequences(cleaved_prefix=0)
    r_unclv_ref = sum([p[1] for p in unclv_ref_info])  #sums up second column (read_counts)

    # 3 - Do the same with the references for the -ligand round
    # cleaved prefix
    clv_ref_info_neg = db.get_ref_sequences(ligand_present=0)
    r_clv_ref_neg = sum([k[1] for k in clv_ref_info_neg])      #sums up second column (read_counts)
    
    # uncleaved prefix
    unclv_ref_info_neg = db.get_ref_sequences(cleaved_prefix=0, ligand_present=0)
    r_unclv_ref_neg = sum([a[1] for a in unclv_ref_info])  #sums up second column (read_counts)


    # 4 - Select and count the reads of the cleaved sequences in the +ligand round
    clv_seq_info = db.get_sequences()                  #returns list of tuples with one value?
    
    for i in range(len(clv_seq_info)):
        one_seq_info = clv_seq_info[i]                 #get specific tuple with information of the sequence
        one_seq = one_seq_info[3]                      #get cleaned sequence 
        one_ID = one_seq_info[0]                       #get key ID 

        r_clv = one_seq_info[1]                        #get read count

        #find uncleaved reads of the same sequence
        unclv_seq_info = db.get_uncleaved_sequence(cleaned_sequence=one_seq)
        unclv_ID = unclv_seq_info[0]                   #get key ID of uncleaved sequence

        r_unclv = unclv_seq_info[1]                    #get read count


        #calculate the cleavage fraction in +ligand round
        clvg_frac = (r_clv / r_clv_ref)  /  ( (r_unclv / r_unclv_ref) + (r_clv / r_clv_ref) )
        
        #put it in database with ID's !!!!


        #-ligand round
        clv_seq_info_neg = db.get_sequence_negligand()                    #get info of cleaved sequence in negative round 
        unclv_seq_info_neg = db.get_sequence_negligand(cleaved_prefix=0)  #get info of uncleaved sequence in negative round

        #cleaved sequence 
        clv_ID_neg = clv_seq_info_neg[0]               #get ID of cleaved sequence in negative round 
        r_clv_neg = clv_seq_info_neg[1]                #get read count cleaved sequence in negative round

        #uncleaved sequence
        unclv_ID_neg = unclv_seq_info_neg[0]           #get ID of uncleaved sequence in negative round
        r_unclv_neg = unclv_seq_info_neg[1]            #get read count uncleaved sequence in negative round


        #calculate the cleavage fraction in -ligand round
        clvg_frac_neg = (r_clv_neg / r_clv_ref_neg)  /  ((r_unclv_neg / r_unclv_ref_neg) + (r_clv_neg / r_clv_ref_neg))
