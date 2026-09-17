"""
tkSnack
An interface to Kare Sjolander's Snack Tcl extension via Tkinter.
http://www.speech.kth.se/snack/index.html

by Kevin Russell and Kare Sjolander
Ported to Python 3 by contributors

Note: for new code, prefer the `snack` package (snack/_snack.c) which
provides direct C bindings and does not require a Tk window.
"""

import tkinter

Tkroot = None
audio = None
mixer = None

def initializeSnack(newroot):
    global Tkroot, audio, mixer
    Tkroot = newroot
    Tkroot.tk.call('eval', 'package require snack')
    Tkroot.tk.call('snack::createIcons')
    Tkroot.tk.call('snack::setUseOldObjAPI')
    audio = AudioControllerSingleton()
    mixer = MixerControllerSingleton()


def _cast(astring):
    """Convert a string returned by a Tcl call to int, float, or str.
    Returns None for empty strings."""
    try:
        return int(astring)
    except (ValueError, TypeError):
        try:
            return float(astring)
        except (ValueError, TypeError):
            return astring if astring else None


class TkObject:
    """Mixin providing Python/Tk communication helpers for Snack objects.
    Modelled on tkinter.Misc."""

    def _getboolean(self, astring):
        if astring:
            return self.tk.getboolean(astring)

    def _getints(self, astring):
        if astring:
            return tuple(map(int, self.tk.splitlist(astring)))

    def _getdoubles(self, astring):
        if astring:
            return tuple(map(float, self.tk.splitlist(astring)))

    def _options(self, cnf, kw=None):
        if kw:
            cnf = tkinter._cnfmerge((cnf, kw))
        else:
            cnf = tkinter._cnfmerge(cnf)
        res = ()
        for k, v in cnf.items():
            if v is not None:
                if k[-1] == '_':
                    k = k[:-1]
                res = res + ('-' + k, v)
        return res

    def configure(self, cnf=None, **kw):
        self._configure(cnf, kw)

    def _configure(self, cnf=None, kw=None):
        if kw is None:
            kw = {}
        if kw:
            cnf = tkinter._cnfmerge((cnf, kw))
        elif cnf:
            cnf = tkinter._cnfmerge(cnf)
        if cnf is None:
            cnf = {}
            for x in self.tk.split(
                    self.tk.call(self.name, 'configure')):
                cnf[x[0][1:]] = (x[0][1:],) + x[1:]
            return cnf
        if isinstance(cnf, str):
            x = self.tk.split(self.tk.call(self.name, 'configure', '-' + cnf))
            return (x[0][1:],) + x[1:]
        self.tk.call((self.name, 'configure') + self._options(cnf))

    config = configure

    def cget(self, key):
        return _cast(self.tk.call(self.name, 'cget', '-' + key))

    __getitem__ = cget

    def __setitem__(self, key, value):
        self.configure({key: value})

    def keys(self):
        return [x[0][1:] for x in
                self.tk.split(self.tk.call(self.name, 'configure'))]

    def __str__(self):
        return self.name


class Sound(TkObject):

    def __init__(self, name=None, master=None, **kw):
        self.name = None
        if not master:
            if Tkroot:
                master = Tkroot
            else:
                raise RuntimeError(
                    'Tk not initialised or not registered with Snack')
        self.tk = master.tk
        if not name:
            self.name = self.tk.call(('sound',) + self._options(kw))
        else:
            self.name = self.tk.call(('sound', name) + self._options(kw))

    def append(self, binarydata, **kw):
        """Append binary string data to the end of the sound."""
        self.tk.call((self.name, 'append', binarydata) + self._options(kw))

    def concatenate(self, othersound):
        """Concatenate sample data from othersound onto the end of this sound."""
        self.tk.call(self.name, 'concatenate', othersound.name)

    def configure(self, **kw):
        """Set options for a sound."""
        self.tk.call((self.name, 'configure') + self._options(kw))

    def copy(self, sound, **kw):
        """Copy sample data from another sound into self."""
        self.tk.call((self.name, 'copy', sound.name) + self._options(kw))

    def changed(self, flag):
        """Inform Snack that the sound object has been modified externally."""
        self.tk.call((self.name, 'changed', flag))

    def convert(self, **kw):
        """Convert sample encoding, sample rate, or channel count."""
        self.tk.call((self.name, 'convert') + self._options(kw))

    def crop(self, start=1, end=None, **kw):
        """Remove all samples outside [start..end]."""
        if end is None:
            end = self.length()
        self.tk.call((self.name, 'crop', start, end) + self._options(kw))

    def cut(self, start=1, end=None, **kw):
        """Remove all samples inside [start..end]."""
        if end is None:
            end = self.length()
        self.tk.call((self.name, 'cut', start, end) + self._options(kw))

    def data(self, binarydata=None, **kw):
        """Load sound data from, or write to, a binary string."""
        if binarydata:
            self.tk.call((self.name, 'data', binarydata) + self._options(kw))
        else:
            return self.tk.call((self.name, 'data') + self._options(kw))

    def destroy(self):
        """Remove the Tcl command and free the storage for this sound."""
        self.tk.call(self.name, 'destroy')

    def dBPowerSpectrum(self, **kw):
        """Compute the log FFT power spectrum; return a list of dB values."""
        result = self.tk.call((self.name, 'dBPowerSpectrum') + self._options(kw))
        return self._getdoubles(result)

    def powerSpectrum(self, **kw):
        """Compute the FFT power spectrum; return a list of magnitude values."""
        result = self.tk.call((self.name, 'powerSpectrum') + self._options(kw))
        return self._getdoubles(result)

    def filter(self, filt, **kw):
        """Apply the given filter to the sound."""
        return self.tk.call((self.name, 'filter', filt.name) + self._options(kw))

    def formant(self, **kw):
        """Return a list of formant trajectories."""
        result = self.tk.call((self.name, 'formant') + self._options(kw))
        return [self._getdoubles(row) for row in self.tk.splitlist(result)]

    def flush(self):
        """Remove all audio data from the sound."""
        self.tk.call(self.name, 'flush')

    def info(self, fmt='string'):
        """Return information about the sound.

        Fields: length, rate, max, min, encoding, channels, fileFormat, headerSize.
        """
        result = self.tk.call(self.name, 'info')
        if fmt == 'list':
            return [_cast(x) for x in result.split()]
        return result

    def insert(self, sound, position, **kw):
        """Insert sound at position."""
        self.tk.call((self.name, 'insert', sound.name, position) + self._options(kw))

    def length(self, n=None, **kw):
        """Get/set the length of the sound in samples (or seconds via -unit)."""
        if n is not None:
            result = self.tk.call((self.name, 'length', n) + self._options(kw))
        else:
            result = self.tk.call((self.name, 'length') + self._options(kw))
        return _cast(result)

    def load(self, filename, **kw):
        """Read new sound data from a file (synonym for read)."""
        self.tk.call((self.name, 'read', filename) + self._options(kw))

    def max(self, **kw):
        """Return the largest positive sample value."""
        return _cast(self.tk.call((self.name, 'max') + self._options(kw)))

    def min(self, **kw):
        """Return the largest negative sample value."""
        return _cast(self.tk.call((self.name, 'min') + self._options(kw)))

    def mix(self, sound, **kw):
        """Mix sample data from another sound into self."""
        self.tk.call((self.name, 'mix', sound.name) + self._options(kw))

    def pause(self):
        """Toggle pause for current play/record operation."""
        self.tk.call(self.name, 'pause')

    def pitch(self, method=None, **kw):
        """Return a list of pitch values."""
        if method is None or method in ('amdf', 'AMDF'):
            result = self.tk.call((self.name, 'pitch') + self._options(kw))
            return self._getdoubles(result)
        result = self.tk.call(
            (self.name, 'pitch', '-method', method) + self._options(kw))
        return [self._getdoubles(row) for row in self.tk.splitlist(result)]

    def play(self, **kw):
        """Play the sound."""
        self.tk.call((self.name, 'play') + self._options(kw))

    def power(self, **kw):
        """Return a list of power values."""
        result = self.tk.call((self.name, 'power') + self._options(kw))
        return self._getdoubles(result)

    def read(self, filename, **kw):
        """Read new sound data from a file."""
        self.tk.call((self.name, 'read', filename) + self._options(kw))

    def record(self, **kw):
        """Start recording from the audio device."""
        self.tk.call((self.name, 'record') + self._options(kw))

    def reverse(self, **kw):
        """Reverse the sample order."""
        self.tk.call((self.name, 'reverse') + self._options(kw))

    def sample(self, index, left=None, right=None):
        """Get or set the sample value at index."""
        if right is not None:
            if left is None:
                left = '?'
            opts = (left, right)
        elif left is not None:
            opts = (left,)
        else:
            opts = ()
        return _cast(self.tk.call((self.name, 'sample', index) + opts))

    def stop(self):
        """Stop the current play or record operation."""
        self.tk.call(self.name, 'stop')

    def stretch(self, **kw):
        self.tk.call((self.name, 'stretch') + self._options(kw))

    def write(self, filename, **kw):
        """Write sound data to a file."""
        self.tk.call((self.name, 'write', filename) + self._options(kw))


class AudioControllerSingleton(TkObject):
    """Controls audio device properties.  One instance is created by
    initializeSnack() and stored in the module-level ``audio`` variable."""

    def __init__(self):
        self.tk = Tkroot.tk

    def encodings(self):
        """Return supported sample encoding formats for the current device."""
        return self.tk.splitlist(self.tk.call('snack::audio', 'encodings'))

    def rates(self):
        """Return supported sample rates for the current device."""
        return self._getints(self.tk.call('snack::audio', 'frequencies'))

    def frequencies(self):
        """Return supported sample rates for the current device."""
        return self._getints(self.tk.call('snack::audio', 'frequencies'))

    def inputDevices(self):
        """Return available audio input devices."""
        return self.tk.splitlist(self.tk.call('snack::audio', 'inputDevices'))

    def outputDevices(self):
        """Return available audio output devices."""
        return self.tk.splitlist(self.tk.call('snack::audio', 'outputDevices'))

    def playLatency(self, latency=None):
        """Get/set (in ms) how much audio is queued to the device."""
        if latency is not None:
            return _cast(self.tk.call('snack::audio', 'playLatency', latency))
        return _cast(self.tk.call('snack::audio', 'playLatency'))

    def pause(self):
        """Toggle play/pause for all playback on the audio device."""
        self.tk.call('snack::audio', 'pause')

    def play(self):
        """Resume paused playback."""
        self.tk.call('snack::audio', 'play')

    def play_gain(self, gain=None):
        """Get/set the current play gain (0–100)."""
        if gain is not None:
            return _cast(self.tk.call('snack::audio', 'play_gain', gain))
        return _cast(self.tk.call('snack::audio', 'play_gain'))

    def selectOutput(self, device):
        """Select a default audio output device."""
        self.tk.call('snack::audio', 'selectOutput', device)

    def selectInput(self, device):
        """Select a default audio input device."""
        self.tk.call('snack::audio', 'selectInput', device)

    def stop(self):
        """Stop all playback on the audio device."""
        self.tk.call('snack::audio', 'stop')

    def elapsedTime(self):
        """Return elapsed time since the audio device started playback."""
        return self.tk.getdouble(self.tk.call('snack::audio', 'elapsedTime'))


class Filter(TkObject):

    def __init__(self, name, *args, **kw):
        global Tkroot
        self.name = None
        if not Tkroot:
            raise RuntimeError('Tk not initialised or not registered with Snack')
        self.tk = Tkroot.tk
        self.name = self.tk.call(('snack::filter', name) + args +
                                 self._options(kw))

    def configure(self, *args):
        """Configure the filter."""
        self.tk.call((self.name, 'configure') + args)

    def destroy(self):
        """Remove the Tcl command for the filter and free its storage."""
        self.tk.call(self.name, 'destroy')


class MixerControllerSingleton(TkObject):
    """Controls mixer device properties.  One instance is created by
    initializeSnack() and stored in the module-level ``mixer`` variable."""

    def __init__(self):
        self.tk = Tkroot.tk

    def channels(self, line):
        """Return channel names for the specified mixer line."""
        return self.tk.splitlist(self.tk.call('snack::mixer', 'channels', line))

    def devices(self):
        """Return available mixer devices."""
        return self.tk.splitlist(self.tk.call('snack::mixer', 'devices'))

    def input(self, jack=None, tclVar=None):
        """Get/set the current input jack."""
        opts = ()
        if jack is not None:
            opts = opts + (jack,)
        if tclVar is not None:
            opts = opts + (tclVar,)
        return self.tk.call(('snack::mixer', 'input') + opts)

    def inputs(self):
        """Return available input ports."""
        return self.tk.splitlist(self.tk.call('snack::mixer', 'inputs'))

    def lines(self):
        """Return mixer line names."""
        return self.tk.splitlist(self.tk.call('snack::mixer', 'lines'))

    def output(self, jack=None, tclVar=None):
        """Get/set the current output jack."""
        opts = ()
        if jack is not None:
            opts = opts + (jack,)
        if tclVar is not None:
            opts = opts + (tclVar,)
        return self.tk.call(('snack::mixer', 'output') + opts)

    def outputs(self):
        """Return available output ports."""
        return self.tk.splitlist(self.tk.call('snack::mixer', 'outputs'))

    def update(self):
        """Update all linked variables to reflect mixer device status."""
        self.tk.call('snack::mixer', 'update')

    def volume(self, line, leftVar=None, rightVar=None):
        if self.channels(line)[0] == 'Mono':
            return self.tk.call('snack::mixer', 'volume', line, rightVar)
        return self.tk.call('snack::mixer', 'volume', line, leftVar, rightVar)

    def select(self, device):
        """Select a device as the default mixer."""
        self.tk.call('snack::mixer', 'select', device)


class SnackCanvas(tkinter.Canvas):

    def __init__(self, master=None, cnf=None, **kw):
        super().__init__(master, **(cnf or {}), **kw)

    def create_spectrogram(self, *args, **kw):
        """Draw a spectrogram of a sound on the canvas."""
        return self._create('spectrogram', args, kw)

    def create_section(self, *args, **kw):
        """Draw an FFT log power spectrum section."""
        return self._create('section', args, kw)

    def create_waveform(self, *args, **kw):
        """Draw a waveform."""
        return self._create('waveform', args, kw)


def createSpectrogram(canvas, *args, **kw):
    """Draw a spectrogram of a sound on canvas."""
    return canvas._create('spectrogram', args, kw)

def createSection(canvas, *args, **kw):
    """Draw an FFT log power spectrum section on canvas."""
    return canvas._create('section', args, kw)

def createWaveform(canvas, *args, **kw):
    """Draw a waveform on canvas."""
    return canvas._create('waveform', args, kw)


class SoundFrame(tkinter.Frame):
    """A simple tape-recorder widget."""

    def __init__(self, parent=None, sound=None, *args, **kw):
        super().__init__(parent, *args, **kw)
        self.sound = sound or Sound()
        self.canvas = SnackCanvas(self, height=100)
        self.canvas.create_waveform(0, 0, sound=self.sound.name)
        self.canvas.pack(side='top')
        bbar = tkinter.Frame(self)
        bbar.pack(side='left')
        tkinter.Button(bbar, bitmap='snackPlay',   command=self.play  ).pack(side='left')
        tkinter.Button(bbar, bitmap='snackRecord', command=self.record).pack(side='left')
        tkinter.Button(bbar, bitmap='snackStop',   command=self.stop  ).pack(side='left')
        tkinter.Button(bbar, text='Info', command=self.info).pack(side='left')

    def load(self):
        f = Tkroot.tk.call('eval', 'snack::getOpenFile')
        self.sound.read(f, progress='snack::progressCallback')

    def play(self):
        self.sound.play()

    def stop(self):
        self.sound.stop()

    def record(self):
        self.sound.record()

    def info(self):
        print(self.sound.info())


if __name__ == '__main__':
    root = tkinter.Tk()
    initializeSnack(root)
    frame = SoundFrame(root)
    frame.pack(expand=0)
    root.mainloop()
