"""
snack - Python bindings for the Snack Sound Toolkit.

This package provides a direct Python API over the Snack C library without
any Tcl/Tk dependency visible to callers.  The C extension (_snack) embeds
a Tcl interpreter internally for resource management, but you never need
to install or interact with Tcl or Tk yourself.

Quick start::

    from snack import Sound

    s = Sound()
    s.read("speech.wav")
    print(s.length, s.sample_rate, s.channels)
    s.play()

    pitch = s.pitch()
    formants = s.formant()          # list of [f1, f2, f3, ...] per frame
    spectrum = s.power_spectrum(fft_length=1024)

Filters::

    from snack import Sound, Filter

    filt = Filter("iir", ...)
    s = Sound()
    s.read("speech.wav")
    s.apply_filter(filt)
"""

from ._snack import (
    Sound,
    Filter,
    get_output_devices,
    get_input_devices,
    audio_select_output,
    audio_select_input,
    audio_frequencies,
    audio_encodings,
    audio_play,
    audio_stop,
    audio_pause,
    audio_elapsed_time,
    audio_play_gain,
    audio_record_gain,
    audio_play_latency,
    LIN16,
    ALAW,
    MULAW,
    LIN8OFFSET,
    LIN8,
    LIN24,
    LIN32,
    SNACK_FLOAT,
    SNACK_DOUBLE,
)

__all__ = [
    "Sound",
    "Filter",
    "get_output_devices",
    "get_input_devices",
    "audio_select_output",
    "audio_select_input",
    "audio_frequencies",
    "audio_encodings",
    "audio_play",
    "audio_stop",
    "audio_pause",
    "audio_elapsed_time",
    "audio_play_gain",
    "audio_record_gain",
    "audio_play_latency",
    "LIN16",
    "ALAW",
    "MULAW",
    "LIN8OFFSET",
    "LIN8",
    "LIN24",
    "LIN32",
    "SNACK_FLOAT",
    "SNACK_DOUBLE",
]

__version__ = "2.2.10"
