import numpy as np
import pytest

from gamma_spectrometry.cli import main
from gamma_spectrometry.rawio import DEFAULT_LSB_VOLTS, build_index
from gamma_spectrometry.store import CaptureStore, convert


@pytest.mark.parametrize("workers", [1, 2])
def test_convert_roundtrip(capture_dir, tmp_path, workers):
    d, order, data = capture_dir
    out = convert(build_index(d), tmp_path / "out" / "s.h5", workers=workers)
    with CaptureStore(out) as s:
        assert len(s) == 5
        assert list(s.file_number) == order
        assert np.all(np.diff(s.time_s) > 0)
        assert s.lsb_volts == DEFAULT_LSB_VOLTS
        assert s.live_time_s == pytest.approx(5 * 4096 / 1.8e6)
        for i, num in enumerate(order):
            np.testing.assert_array_equal(s.codes(i), data[num])
            np.testing.assert_allclose(s.volts(i), data[num] * DEFAULT_LSB_VOLTS)


def test_convert_refuses_to_overwrite(capture_dir, tmp_path):
    d, _, _ = capture_dir
    out = convert(build_index(d), tmp_path / "s.h5", workers=1)
    with pytest.raises(FileExistsError):
        convert(build_index(d), out, workers=1)
    convert(build_index(d), out, workers=1, overwrite=True)


def test_cli_convert_limit(capture_dir, tmp_path, capsys):
    d, order, _ = capture_dir
    out = tmp_path / "s.h5"
    assert main(["convert", str(d), str(out), "--limit", "3", "--workers", "1"]) == 0
    with CaptureStore(out) as s:
        assert list(s.file_number) == order[:3]
    assert "live time" in capsys.readouterr().out
