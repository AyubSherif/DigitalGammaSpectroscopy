import numpy as np
import pytest

from gammaspec.rawio import (DEFAULT_LSB_VOLTS, RawFormatError, build_index, estimate_lsb,
                             parse_samples, read_header, read_raw_file, volts_to_codes)

from .conftest import REAL_DATA


def test_header(capture_dir):
    d, order, _ = capture_dir
    h = read_header(d / "raw_data_3.txt")
    assert h.file_number == 3
    assert h.n_samples == 4096
    assert h.dt_s == pytest.approx(1 / 1.8e6)
    assert h.duration_s == pytest.approx(4096 / 1.8e6)


def test_codes_roundtrip_exact(capture_dir):
    d, _, data = capture_dir
    for num, codes in data.items():
        _, got = read_raw_file(d / f"raw_data_{num}.txt")
        assert got.dtype == np.int32
        np.testing.assert_array_equal(got, codes)


def test_fast_parser_matches_generic_float_parsing(capture_dir):
    d, _, _ = capture_dir
    p = d / "raw_data_0.txt"
    h = read_header(p)
    body = p.read_bytes()[h.header_bytes:]
    np.testing.assert_allclose(parse_samples(body), np.array(body.split(), dtype=np.float64), rtol=0, atol=1e-15)


def test_generic_fallback_for_other_formats():
    body = b"0.5\r\n-1.25e-3\r\n2\r\n"
    np.testing.assert_allclose(parse_samples(body), [0.5, -1.25e-3, 2.0])


def test_rejects_values_off_the_adc_grid():
    with pytest.raises(RawFormatError):
        volts_to_codes(np.array([0.0, 0.4 * DEFAULT_LSB_VOLTS]), DEFAULT_LSB_VOLTS)


def test_sample_count_mismatch(capture_dir):
    d, _, _ = capture_dir
    p = d / "raw_data_1.txt"
    p.write_bytes(p.read_bytes()[:-14])  # drop the last sample line
    with pytest.raises(RawFormatError, match="samples"):
        read_raw_file(p)


def test_index_is_chronological(capture_dir):
    d, order, _ = capture_dir
    headers = build_index(d)
    assert [h.file_number for h in headers] == order
    assert all(a.created < b.created for a, b in zip(headers, headers[1:]))


def test_estimate_lsb(capture_dir):
    d, _, _ = capture_dir
    _, volts = read_raw_file(d / "raw_data_2.txt", lsb=None)
    assert estimate_lsb(volts) == pytest.approx(DEFAULT_LSB_VOLTS, rel=1e-7)


@pytest.mark.skipif(not (REAL_DATA / "raw_data_0.txt").exists(), reason="real data not present")
def test_real_file_codes_reproduce_text():
    """Every printed value of a real capture is reproduced by code*LSB to 7 significant digits,
    and each code maps to exactly one printed string (i.e. the conversion is lossless)."""
    p = REAL_DATA / "raw_data_0.txt"
    h, codes = read_raw_file(p)
    lines = np.array(p.read_bytes()[h.header_bytes:].split())
    volts = lines.astype(np.float64)
    np.testing.assert_allclose(codes * DEFAULT_LSB_VOLTS, volts, rtol=1e-6, atol=1e-12)
    order = np.argsort(codes, kind="stable")
    c, s = codes[order], lines[order]
    same_code = c[1:] == c[:-1]
    assert np.all(s[1:][same_code] == s[:-1][same_code])
