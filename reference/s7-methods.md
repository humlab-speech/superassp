# S7 Method System for DSP Functions

This file sets up S7 method dispatch for all lst\_\* and trk\_\*
functions, allowing them to accept both character vectors (file paths)
and AVAudio objects.

## Details

The system works by:

1.  Converting each existing function to an S7 generic

2.  Registering the original implementation as the character method

3.  Adding an AVAudio method that converts to temp file and calls
    original

This preserves full backward compatibility while adding AVAudio support.

Note: listOfFiles is now a mandatory parameter in all DSP functions (no
default value).
