import datetime as dt
from pathlib import Path

import numpy as np
import pytest

from gammaspec.rawio import DEFAULT_LSB_VOLTS

REAL_DATA = Path(__file__).resolve().parents[1] / "ADC_raw_output"


def write_capture(path: Path, codes: np.ndarray, created: dt.datetime,
                  lsb: float = DEFAULT_LSB_VOLTS, dt_ms: float = 5.555555555555556e-04) -> Path:
    """Write a capture file in the exact format of the 2021 acquisition software."""
    body = b"".join(b"%+.6e\n" % v for v in codes.astype(np.float64) * lsb)
    header = (f"Creation: {created.isoformat(' ')}\n"
              f"Signals:  {codes.size}\n"
              f"Dt:       {dt_ms!r} ms\n").encode()
    path.write_bytes(header + body)
    return path


@pytest.fixture
def capture_dir(tmp_path):
    """Five small captures whose file numbers are deliberately NOT in time order."""
    rng = np.random.default_rng(1)
    t0 = dt.datetime(2021, 8, 26, 13, 34, 48, 617564)
    order = [3, 0, 4, 1, 2]  # file number -> position in time
    data = {}
    for k, num in enumerate(order):
        codes = rng.integers(-65000, 3800, size=4096).astype(np.int32)
        codes[:3] = [0, 1, -1]
        write_capture(tmp_path / f"raw_data_{num}.txt", codes, t0 + dt.timedelta(seconds=2.3 * k))
        data[num] = codes
    return tmp_path, order, data
