import sqlite3


class DatabaseInterface:
    """Class to allow for easy interface with a database."""

    def __init__(self, path=None):
        self.path = path
        self.conn = None
        self.cursor = None

        if path is not None:
            self.open(path)

    def open(self, path: str):
        """Function tries to connect to a database."""
        try:
            self.conn = sqlite3.connect(path)
            self.cursor = self.conn.cursor()
        except sqlite3.Error:
            raise Exception("Error connecting to the database")

    def close(self):
        """Function closes the database connection if there is one."""
        if self.conn is not None:
            self.conn.commit()
            self.cursor.close()
            self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self.close()

    def is_open(self):
        """Function checks whether there is a database connection."""
        return self.conn is not None

    def table_exists(self, table: str) -> bool:
        """
        Function checks whether the table with 'table_name' exists in the sqlite database.\n
        args:\n
        table_name: (str) name of the table you want to check for existence.\n
        \n
        return values:\n
        boolean: tells whether the table exists
        """

        if not self.is_open():
            return False

        self.cursor.execute(
            f"SELECT name FROM sqlite_master WHERE type='table' and name='{table}'")
        tables = self.cursor.fetchall()
        if len(tables) > 0:
            return True
        # else
        return False

    def get(self, table: str, columns: list, limit: int = None, offset: int = 0) -> list:
        """
        Function to query data from a database.\n
        args:\n
        table: (str) name of the table you want to query.\n
        columns: (list) list of strings of the column names.\n
        limit: (int) optional keyword argument where you can limit the number of rows retrieved. The default value returns all rows.\n
        offset: (int) optional keyword argument that allows you to offset your search query. This is only used if a non-default limit is set.\n
        \n
        return values:\n
        list of tuples containing all information from the queried rows.
        """

        if not self.is_open():
            raise Exception("There is no database connection")

        if isinstance(columns, list):
            columns = ",".join(columns)
            if limit is None:
                self.cursor.execute(f"SELECT {columns} FROM {table} ")
            else:
                self.cursor.execute(
                    f"SELECT {columns} FROM {table} LIMIT {limit} OFFSET {offset}")
            return self.cursor.fetchall()
        else:
            raise TypeError(
                "The input variable 'columns' needs to be a list of column names as strings")

    def get_all(self, table: str, limit: int = None, offset: int = 0) -> list:
        """
        Function to query rows and all columns from a database.\n
        args:\n
        table: (str) name of the table you want to query.\n
        limit: (int) optional keyword argument where you can limit the number of rows retrieved. The default value returns all rows.\n
        offset: (int) optional keyword argument that allows you to offset your search query. This is only used if a non-default limit is set.\n
        \n
        return values:\n
        list of tuples containing all information from the queried rows.
        """

        if not self.is_open():
            raise Exception("There is no database connection")

        if limit is None:
            self.cursor.execute(f"SELECT * FROM {table} ")
        else:
            self.cursor.execute(
                f"SELECT * FROM {table} LIMIT {limit} OFFSET {offset}")
        return self.cursor.fetchall()

    def query(self, sql: str, parameters=None):
        """Function to query any SQL statement."""

        if not self.is_open():
            raise Exception("There is no database connection")

        if parameters is None:
            self.cursor.execute(sql)
        else:
            self.cursor.execute(sql, parameters)


class DatabaseInterfaceRawSequences(DatabaseInterface):

    def __init__(self, path=None):
        super().__init__(path)
        self.table = None

    def create_table(self, table: str):
        """
        Function creates a table in the database with the name given by the variable table.\n
        \n
        The created database has the following columns:\n
        id: (INTEGER PRIMARY KEY) unique integer for every row\n
        read_count: (INTEGER) number of reads\n
        original_sequence: (TEXT) sequence with barcode, prefix and suffix still attached\n
        cleaned_sequence: (TEXT) sequence with barcode, prefix and suffix removed\n
        barcode: (TEXT) barcode sequence\n
        cleaved_prefix: (INTEGER) indicates whether the prefix corresponds to the cleaved prefix. Yes(1)/No(0)/Don't know(2)\n
        prefix_name: (TEXT) name of the prefix\n
        reference_name: (TEXT) is name of the reference sequence, if not a reference sequence value is NULL\n
        selection: (TEXT) name of the selection\n
        driver_round: (INTEGER) round of driver at the moment of sequencing\n
        ligand_present: (INTEGER) indicates whether the ligand was present before sequencing. Yes(1)/No(0)\n
        cleavage_fraction: (REAL) value of the cleavage fraction for a sequence\n
        fold_change: (REAL) value of the fold change for a sequence\n
        possible_sensor: (INTEGER) indicates whether the sequence is a possible sensor. Yes(1)/No(1)\n
        mutated_prefix: (INTEGER) indicates whether the prefix had a mutation. Yes(1)/No(0)\n
        \n
        args:\n
        table: (str) name of the table to be created.
        """
        self.table = table
        self.query(f"""CREATE TABLE IF NOT EXISTS {table} (
                    id INTEGER PRIMARY KEY,
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
                    possible_sensor INTEGER,
                    mutated_prefix INTEGER
                    )""")

    def insert_sequence_info(self, table: str, sequence_info: dict):
        """
        Function inserts sequence info data into the database table 'table'.\n
        args:\n
        table: (str) name of the table.\n
        sequence_info: (dict) dictionary where the keyword is the name of the column in the table and the value is the value you want to insert in that column.
        """
        if not self.is_open():
            raise Exception("There is no database connection")

        self.query(f"""INSERT INTO {table}(read_count, original_sequence, cleaned_sequence,
                            barcode, cleaved_prefix, prefix_name, reference_name,
                            selection, driver_round, ligand_present, cleavage_fraction,
                            fold_change, possible_sensor, mutated_prefix) VALUES (
                            :read_count, :original_sequence, :cleaned_sequence,
                            :barcode, :cleaved_prefix, :prefix_name, :reference_name,
                            :selection, :driver_round, :ligand_present, :cleavage_fraction,
                            :fold_change, :possible_sensor, :mutated_prefix)""", parameters=sequence_info)
    
    

class DatabaseInterfaceCleanSequences(DatabaseInterface):

    def __init__(self, path=None):
        super().__init__(path)
        self.table = None

    # TODO remove unneeded columns
    def create_table(self, table: str):
        """
        Function creates a table in the database with the name given by the variable table.\n
        \n
        The created database has the following columns:\n
        id: (INTEGER PRIMARY KEY) unique integer for every row\n
        read_count: (INTEGER) number of reads\n
        # original_sequence: (TEXT) sequence with barcode, prefix and suffix still attached\n
        cleaned_sequence: (TEXT) sequence with barcode, prefix and suffix removed\n
        # barcode: (TEXT) barcode sequence\n
        cleaved_prefix: (INTEGER) indicates whether the prefix corresponds to the cleaved prefix. Yes(1)/No(0)/Don't know(2)\n
        prefix_name: (TEXT) name of the prefix\n
        reference_name: (TEXT) is name of the reference sequence, if not a reference sequence value is NULL\n
        selection: (TEXT) name of the selection\n
        driver_round: (INTEGER) round of driver at the moment of sequencing\n
        ligand_present: (INTEGER) indicates whether the ligand was present before sequencing. Yes(1)/No(0)\n
        cleavage_fraction: (REAL) value of the cleavage fraction for a sequence\n
        fold_change: (REAL) value of the fold change for a sequence\n
        possible_sensor: (INTEGER) indicates whether the sequence is a possible sensor. Yes(1)/No(1)\n
        mutated_prefix: (INTEGER) indicates whether the prefix had a mutation. Yes(1)/No(0)\n
        \n
        args:\n
        table: (str) name of the table to be created.
        """
        self.table = table
        self.query(f"""CREATE TABLE IF NOT EXISTS {table} (
                    id INTEGER PRIMARY KEY,
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
                    possible_sensor INTEGER,
                    mutated_prefix INTEGER
                    )""")

    # TODO remove unneeded columns
    def insert_sequence_info(self, table: str, sequence_info: dict):
        """
        Function inserts sequence info data into the database table 'table'.\n
        args:\n
        table: (str) name of the table.\n
        sequence_info: (dict) dictionary where the keyword is the name of the column in the table and the value is the value you want to insert in that column.
        """
        if not self.is_open():
            raise Exception("There is no database connection")

        self.query(f"""INSERT INTO {table}(read_count, original_sequence, cleaned_sequence,
                            barcode, cleaved_prefix, prefix_name, reference_name,
                            selection, driver_round, ligand_present, cleavage_fraction,
                            fold_change, possible_sensor, mutated_prefix) VALUES (
                            :read_count, :original_sequence, :cleaned_sequence,
                            :barcode, :cleaved_prefix, :prefix_name, :reference_name,
                            :selection, :driver_round, :ligand_present, :cleavage_fraction,
                            :fold_change, :possible_sensor, :mutated_prefix)""", parameters=sequence_info)

    

    def update_cleavage_fraction(self, table: str, rowid: int, cleavage_fraction: float):
        """
        Function updates the cleavage fraction at the row rowid.\n
        args:\n
        table: (str) name of the table\n
        rowid: (int) row to update\n
        cleavage_fraction: (float) value to insert in the cleavage_fraction column\n
        """
        self.query(
            f"UPDATE {table} SET cleavage_fraction={cleavage_fraction} WHERE id={rowid}")

    def update_fold_change(self, table: str, rowid: int, fold_change: float):
        """
        Function updates the cleavage fraction at the row rowid.\n
        args:\n
        table: (str) name of the table\n
        rowid: (int) row to update\n
        fold_change: (float) value to insert in the fold_change column\n
        """
        self.query(
            f"UPDATE {table} SET fold_change={fold_change} WHERE id={rowid}")


    def get_sequences(self, cleaved_prefix: int = 1, ligand_present: int = 1):
        """
        1 = yes
        0 = no
        """
        self.cursor.execute(
            f"SELECT * FROM sequences WHERE cleaved_prefix={cleaved_prefix} AND ligand_present={ligand_present}")
        return self.cursor.fetchall()

    def get_ref_sequences(self, cleaved_prefix: int = 1, ligand_present: int=1):
        """
        cleaved_prefix = 1 --> reference sequence with prefix corresponding to cleaved sequences
        cleaved_prefix = 0 --> reference sequence with prefix corresponding to uncleaved sequences
        """
        self.cursor.execute(
            f"SELECT * FROM sequences WHERE cleaved_prefix={cleaved_prefix} AND ligand_present={ligand_present} AND reference_name IS NOT NULL")
        return self.cursor.fetchall()
    
    def get_uncleaved_sequence(self, cleaned_sequence: str, cleaved_prefix: int=0, ligand_present: int=1):
        """
        Function finds the information of the uncleaved variant of a specific sequence
        cleaned_sequence = target sequence
        cleaved_prefix
        """

        self.cursor.execute(
            f"SELECT * FROM sequences WHERE cleaned_sequence={cleaned_sequence} AND cleaved_prefix={cleaved_prefix} AND ligand_present={ligand_present}")
        return self.cursor.fetchall()

    def get_sequence_negligand(self, cleaned_seqence: str, cleaved_prefix: int=1, ligand_present: int=0):
        """
        Function finds the information of specific sequence in the negative ligand rounds
        """
        self.cursor.execute(
            f"SELECT * FROM sequences WHERE cleaned_sequence={cleaned_sequence} AND ligand_present={ligand_present} AND cleaved_prefix={cleaved_prefix}")
        return self.cursor.fetchall()
