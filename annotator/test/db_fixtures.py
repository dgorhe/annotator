"""Build sqlite databases for committed GTF slices."""

from pathlib import Path

from annotator.build_gffutils_db import build_db

TEST_DATA = Path(__file__).parent / "data"
DB_CACHE = Path(__file__).parent / ".db_cache"


def ensure_db(gtf_name):
    DB_CACHE.mkdir(exist_ok=True)
    db_path = DB_CACHE / (Path(gtf_name).name + ".db")
    gtf_path = TEST_DATA / gtf_name
    if not db_path.exists():
        build_db(str(gtf_path), str(db_path), force=True)
    return str(db_path)
