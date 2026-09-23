"""Cross-platform durable documents and ordered events."""

import contextlib
import json
import sqlite3
import threading
import time
from pathlib import Path

_LOCKS = {}
_GUARD = threading.Lock()


class Store:
    def __init__(self, root):
        self.root = Path(root).expanduser().resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.path = self.root / "client.db"
        with self.connect() as db:
            db.executescript("""
            CREATE TABLE IF NOT EXISTS docs (kind TEXT, id TEXT, body TEXT NOT NULL,
              PRIMARY KEY(kind,id));
            CREATE TABLE IF NOT EXISTS events (seq INTEGER PRIMARY KEY AUTOINCREMENT,
              task TEXT NOT NULL, stamp REAL NOT NULL, body TEXT NOT NULL);
            """)
        with _GUARD:
            self.lock = _LOCKS.setdefault(str(self.path), threading.RLock())

    @contextlib.contextmanager
    def connect(self):
        db = sqlite3.connect(self.path, timeout=120)
        db.execute("PRAGMA journal_mode=WAL")
        db.execute("PRAGMA busy_timeout=120000")
        try:
            with db:
                yield db
        finally:
            db.close()

    def get(self, kind, key, default=None):
        with self.connect() as db:
            row = db.execute(
                "SELECT body FROM docs WHERE kind=? AND id=?", (kind, key)
            ).fetchone()
        return json.loads(row[0]) if row else default

    def put(self, kind, key, value):
        with self.connect() as db:
            db.execute(
                "INSERT OR REPLACE INTO docs VALUES (?,?,?)",
                (kind, key, json.dumps(value)),
            )
        return value

    def all(self, kind):
        with self.connect() as db:
            return [
                json.loads(r[0])
                for r in db.execute(
                    "SELECT body FROM docs WHERE kind=? ORDER BY id", (kind,)
                )
            ]

    def event(self, task, value):
        with self.connect() as db:
            stamp = time.time()
            cur = db.execute(
                "INSERT INTO events(task,stamp,body) VALUES (?,?,?)",
                (task, stamp, json.dumps(value)),
            )
            return {"seq": cur.lastrowid, "time": stamp, **value}

    def events(self, task, after=0):
        with self.connect() as db:
            return [
                {"seq": s, "time": t, **json.loads(v)}
                for s, t, v in db.execute(
                    "SELECT seq,stamp,body FROM events WHERE task=? AND seq>? ORDER BY seq LIMIT 500",
                    (task, after),
                )
            ]

    @contextlib.contextmanager
    def exclusive(self):
        # A separate DB supplies a cross-process lease without blocking document writes.
        with self.lock:
            db = sqlite3.connect(
                self.root / "lease.db", timeout=120, isolation_level=None
            )
            try:
                db.execute("BEGIN EXCLUSIVE")
                yield
                db.execute("COMMIT")
            finally:
                db.close()
