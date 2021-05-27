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

    def query(self, sql: str, parameters=None):
        """Function to query any other SQL statement."""

        if not self.is_open():
            raise Exception("There is no database connection")

        if parameters is None:
            self.cursor.execute(sql)
        else:
            self.cursor.execute(sql, parameters)


class DatabaseInterfaceSequences(DatabaseInterface):

    def __init__(self, path=None):
        super().__init__(path)

    def get_sequences(self, cleaved_prefix: int = 1, ligand_present: int = 1):
        """
        1 = yes
        0 = no
        """
        self.cursor.execute(
            f"SELECT * FROM sequences WHERE cleaved_prefix={cleaved_prefix} AND ligand_present={ligand_present}")
        return self.cursor.fetchall()
