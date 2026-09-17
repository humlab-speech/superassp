# snack-sound

A Python C extension for the [Snack Sound Toolkit](http://www.speech.kth.se/snack/)
(v2.2) developed at KTH — Royal Institute of Technology, Sweden.

Snack provides high-level audio objects for playing, recording, and analysing
sound. This package exposes Snack's non-graphical functionality directly to
Python without requiring Tcl/Tk or Tkinter at runtime.

## Requirements

* Python ≥ 3.8
* A Tcl installation at build time (`tcl-dev` on Debian/Ubuntu, `tcl-tk` via
  Homebrew on macOS)
* ALSA (`libasound2-dev`) on Linux for audio I/O

## Installation

```
pip install snack-sound
```

Or build from source:

```
cd python/
pip install .
```

## Quick start

```python
from snack import Sound

s = Sound()
s.read("hello.wav")
print(s.info())       # (length, rate, maxval, channels, encoding, type)
s.play(blocking=True)
s.write("copy.wav")
```

### Tone generation

```python
from snack import Sound, Filter

s = Sound()
filt = Filter("generator", 440, 30000, 0.0, "sine", 8000)
s.play(filter=filt, blocking=True)
```

### Signal analysis

```python
from snack import Sound

s = Sound()
s.read("speech.wav")
pitch = s.pitch()          # list of F0 values per frame
power = s.power()          # RMS power per frame
formants = s.formant()     # list of (F1, F2, F3, …) per frame
spectrum = s.dBPowerSpectrum()
```

## License

BSD-style Snack Toolkit license; see `../LICENSE.txt`.
