import json
import sqlite3
from pathlib import Path


class SearchCache:
    """Own SQLite lifetime; None differs from a cached no-hit []."""
    def __init__(self,path):
        self.path=Path(path)
        self.connection=None

    def __enter__(self):
        self.path.parent.mkdir(parents=True,exist_ok=True)
        self.connection=sqlite3.connect(self.path)
        try:
            self.connection.execute('CREATE TABLE IF NOT EXISTS hits (key TEXT PRIMARY KEY, result TEXT NOT NULL)')
        except BaseException:
            self.connection.close()
            self.connection=None
            raise
        return self

    def __exit__(self,exc_type,exc,traceback):
        try:
            if exc_type is None:
                self.connection.commit()
            else:
                self.connection.rollback()
        finally:
            self.connection.close()
            self.connection=None

    def get(self,key):
        row=self.connection.execute('SELECT result FROM hits WHERE key=?',(key,)).fetchone()
        return None if row is None else json.loads(row[0])

    def put(self,key,result):
        self.connection.execute('INSERT OR REPLACE INTO hits VALUES (?,?)',(key,json.dumps(result)))
