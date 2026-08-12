from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_manuscript_evidence_check() -> None:
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts/check_manuscript_evidence.py")],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert '"status": "valid"' in result.stdout
