import sqlite3

conn = sqlite3.connect(':memory:')

c = conn.cursor()

# A database is created with the following columns:
# sequence: DNA sequence with barcode prefix and suffix removed (TEXT),
# barcode: barcode of the sequence (TEXT), cleaved_prefix: yes(1)/no(0) (INTEGER),
# cleaved_suffix: yes(1)/no(0) (INTEGER), reference: indicates whether sequence is a
# reference sequence yes(1)/no(0) (INTEGER), round: round when the sequence was
# sequenced(INTEGER), ligand: ligand present yes(1)/no(0) (INTEGER),
# sensor: indicates whether sequence is a possible biosensor yes(1)/no(0) (INTEGER)

c.execute("""CREATE TABLE sequences (
            sequence TEXT,
            barcode TEXT,
            cleaved_prefix INTEGER,
            cleaved_suffix INTEGER,
            reference INTEGER,
            round INTEGER,
            ligand INTEGER,
            sensor INTEGER
            )""")
