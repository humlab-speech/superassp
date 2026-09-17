import tempfile
import unittest
from pathlib import Path

from snack import Sound


def _find_fixture() -> Path:
    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        candidate = parent / "demos" / "tcl" / "ex1.wav"
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Unable to locate demos/tcl/ex1.wav")


FIXTURE = _find_fixture()
RATE_CONVERT_SAMPLE_DELTA = 2.0
# Reference outputs from the legacy Tcl tests on ex1.wav.
REF_POWER_0 = 57.668
REF_POWER_1 = 58.916
REF_SPECTRUM_0 = -27.449
REF_SPECTRUM_1 = -29.165


class SnackSoundFromTclTests(unittest.TestCase):
    def test_sample_and_length_from_fixture(self):
        s = Sound()
        s.read(str(FIXTURE))
        self.assertEqual(s.length(), 15820)
        self.assertEqual(s.get_sample(0, 0), 1177.0)

    def test_crop_and_cut(self):
        s = Sound()
        s.length(300000)
        s.set_sample(0, 1000, 17)
        s.set_sample(0, 299000, 17)
        s.crop(1000, 299000)
        self.assertEqual(s.get_sample(0, 0), 17.0)
        self.assertEqual(s.get_sample(0, 298000), 17.0)

        s = Sound()
        s.length(300000)
        s.set_sample(0, 999, 17)
        s.set_sample(0, 299001, 17)
        s.cut(1000, 299000)
        self.assertEqual(s.get_sample(0, 999), 17.0)
        self.assertEqual(s.get_sample(0, 1000), 17.0)

    def test_concatenate_copy_insert(self):
        source = Sound()
        source.read(str(FIXTURE))

        dst = Sound()
        dst.length(1000)
        dst.concatenate(source)
        self.assertEqual(dst.get_sample(0, 999), 0.0)
        self.assertEqual(dst.get_sample(0, 6000), 7443.0)
        self.assertEqual(dst.get_sample(0, 11000), 779.0)

        copied = Sound()
        copied.copy(source, start=5000, end=10000)
        self.assertEqual(copied.get_sample(0, 0), 7443.0)
        self.assertEqual(copied.get_sample(0, 5000), 779.0)

        inserted = Sound()
        inserted.length(1000)
        inserted.insert(copied, 500)
        self.assertEqual(inserted.get_sample(0, 499), 0.0)
        self.assertEqual(inserted.get_sample(0, 500), 7443.0)
        self.assertEqual(inserted.get_sample(0, 5500), 779.0)
        self.assertEqual(inserted.get_sample(0, 5501), 0.0)

    def test_convert(self):
        s = Sound()
        s.read(str(FIXTURE))
        s.convert(channels=2)
        self.assertEqual(s.get_sample(0, 0), 1177.0)
        self.assertEqual(s.get_sample(1, 0), 1177.0)

        s = Sound()
        s.read(str(FIXTURE))
        s.convert(encoding="Mulaw")
        self.assertEqual(s.get_sample(0, 0), 1180.0)

        s = Sound()
        s.read(str(FIXTURE))
        s.convert(rate=32000)
        self.assertAlmostEqual(s.get_sample(0, 5000), 1688.0, delta=RATE_CONVERT_SAMPLE_DELTA)

    def test_data_and_append_data(self):
        s = Sound()
        s.data(b"\x00\x00\xff\x00\x00\x00")
        info = s.info()
        self.assertEqual(info["length"], 3)
        self.assertEqual(info["max"], 255)
        self.assertEqual(info["min"], 0)

        s = Sound()
        s.length(1000)
        with open(FIXTURE, "rb") as f:
            s.append_data(f.read())
        self.assertEqual(s.get_sample(0, 999), 0.0)
        self.assertEqual(s.get_sample(0, 6000), 7443.0)
        self.assertEqual(s.get_sample(0, 11000), 779.0)

    def test_max_min(self):
        s = Sound()
        s.read(str(FIXTURE))
        self.assertEqual(s.max_sample(), 14264.0)
        self.assertEqual(s.min_sample(), -8288.0)
        self.assertEqual(s.max_sample(start=0, end=4), 1201.0)
        self.assertEqual(s.min_sample(start=0, end=4), 1169.0)

    def test_reverse(self):
        s = Sound(channels=2)
        s.length(4)
        s.set_sample(0, 0, 1)
        s.set_sample(1, 0, 2)
        s.set_sample(0, 1, 3)
        s.set_sample(1, 1, 4)
        s.set_sample(0, 2, 5)
        s.set_sample(1, 2, 6)
        s.set_sample(0, 3, 7)
        s.set_sample(1, 3, 8)
        s.reverse()
        self.assertEqual((s.get_sample(0, 0), s.get_sample(1, 0)), (7.0, 8.0))
        self.assertEqual((s.get_sample(0, 1), s.get_sample(1, 1)), (5.0, 6.0))
        self.assertEqual((s.get_sample(0, 2), s.get_sample(1, 2)), (3.0, 4.0))
        self.assertEqual((s.get_sample(0, 3), s.get_sample(1, 3)), (1.0, 2.0))

    def test_fileio(self):
        with tempfile.TemporaryDirectory() as td:
            wav = Path(td) / "snackTest.wav"
            au = Path(td) / "snackTest.au"
            aiff = Path(td) / "snackTest.aiff"

            s = Sound()
            s.length(300)

            s.write(str(wav))
            s.read(str(wav))
            self.assertEqual(s.info()["file_format"], "WAV")

            s.write(str(au))
            s.read(str(au))
            self.assertEqual(s.info()["file_format"], "AU")

            s.write(str(aiff))
            s.read(str(aiff))
            self.assertEqual(s.info()["file_format"], "AIFF")

    def test_power_and_spectrum(self):
        s = Sound()
        s.read(str(FIXTURE))
        power = s.power()
        spectrum = s.power_spectrum()
        self.assertGreater(len(power), 0)
        self.assertEqual(len(spectrum), 256)
        self.assertAlmostEqual(power[0], REF_POWER_0, places=2)
        self.assertAlmostEqual(power[1], REF_POWER_1, places=2)
        self.assertAlmostEqual(spectrum[0], REF_SPECTRUM_0, places=2)
        self.assertAlmostEqual(spectrum[1], REF_SPECTRUM_1, places=2)


if __name__ == "__main__":
    unittest.main()
