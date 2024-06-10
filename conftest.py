import os
from pathlib import Path
import pytest

@pytest.fixture(scope="session")
def project_root():
    return Path(os.path.dirname(os.path.abspath(__file__)))
    