###############################################
# Convert PraatPic to PDF using Praat
# For Praat 6.0+ (no external tools needed)
# Place this script in same directory as .prapic files
#
# Use with Parselmouth workflow:
# 1. Python/Parselmouth saves pictures as .prapic
# 2. Run this Praat script to convert to PDF
###############################################

# Get current directory
current_dir$ = "./"

# Create list of PraatPic files
pic_list = Create Strings as file list: "pic_files", current_dir$ + "*.prapic"
n_files = Get number of strings

if n_files = 0
    exitScript: "No .prapic files found in current directory"
endif

writeInfoLine: "Found ", n_files, " Praat picture file(s)"
appendInfoLine: "Converting to PDF..."
appendInfoLine: ""

# Process each file
for i to n_files
    selectObject: pic_list
    filename$ = Get string: i
    basename$ = filename$ - ".prapic"

    pic_file$ = current_dir$ + filename$
    pdf_file$ = current_dir$ + basename$ + ".pdf"

    # Clear picture window
    Erase all

    # Read PraatPic file into Picture window
    Read from praat picture file: pic_file$

    # Save as PDF
    Save as PDF file: pdf_file$

    appendInfoLine: "✓ ", filename$, " → ", basename$, ".pdf"
endfor

removeObject: pic_list

appendInfoLine: ""
appendInfoLine: "Done! Converted ", n_files, " image file(s) to PDF"
