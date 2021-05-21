import sqlite3


class Database:

    def __init__(self, path=None):
        self.conn = None
        self.cursor = None

        if path is not None:
            self.open(path)

    def open(self, path):
        """Function tries to connect to a database"""
        try:
            self.conn = sqlite3.connect(path)
            self.cursor = self.conn.cursor()
        except sqlite3.Error:
            print("Error connecting to the database")

    def close(self):
        """Function closes the database connection if there is one"""
        if self.conn is not None:
            self.conn.commit()
            self.cursor.close()
            self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self.close()

    def is_open(self):
        """Function checks whether there is a database connection"""
        return self.conn is not None

    def table_exists(self, table):
        """
        Function checks whether the table with 'table_name' exists in the sqlite database.
        args:
        table_name (str)

        return:
        boolean
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

    def get(self, table, columns, limit=None):

        if not self.is_open():
            raise Exception("There is no database connection")

        if isinstance(columns, list):
            columns = ",".join(columns)
            if limit is None:
                self.cursor.execute(f"SELECT {columns} FROM {table}")
            else:
                self.cursor.execute(
                    f"SELECT {columns} FROM {table} LIMIT {limit}")
            return self.cursor.fetchall()
        else:
            raise TypeError(
                "The input variable 'columns' needs to be a list of column names as strings")

    def query(self, sql):
        """Function to query any other SQL statement"""

        if not self.is_open():
            raise Exception("There is no database connection")

        self.cursor.execute(sql)
