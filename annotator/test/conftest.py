from pathlib import Path

import pytest

from annotator.test.db_fixtures import ensure_db


@pytest.fixture(scope="session")
def chr19_db():
    return ensure_db("gencode.v19.annotation.chr19.10000000.gtf")


@pytest.fixture(scope="session")
def chr21_db():
    return ensure_db("gencode.v19.annotation.chr21.46M-47M.gtf")


@pytest.fixture(scope="session")
def chr11_db():
    return ensure_db("gencode.v19.annotation.chr11.70M-71M.gtf")


@pytest.fixture(scope="session")
def proxdist_db():
    return ensure_db("proxdist_test.gtf")


@pytest.fixture(scope="session")
def celegans_chrIII_db():
    return ensure_db("c_elegans.PRJNA13758.WS257.canonical_geneset.chrIII.25000.gtf")
