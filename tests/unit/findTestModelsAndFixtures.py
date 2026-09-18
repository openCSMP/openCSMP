# SKM query: I have a bunch of test files that read input files from the test directory, 
# but the latter contains test artifacts and other stuff. 
# So my idea is to "grep" the source code for the file names to obtain a complete list 
# of the input files required.
# Most tests construct models using constructors with the syntax:
# 
# ( "BoxHalfs2D", "CSMP-variables.txt" );
#
# so the files I am looking for are indicated by the constructor arguments i quotation marks. 
# Test fixtures include "-variables.txt" text files that must also be present as well config files
# with the added string "-configuration.txt" and "-regions.txt" files.
# the following script that I can run on osx from the command line in the tests directory 
# to get a list of these files?

#!/usr/bin/env python3

# HOW TO RUN: python3 findTestModelsAndFixtures.py

#!/usr/bin/env python3

import os
import re

# ---------------------------------------------------------
# regex patterns
# ---------------------------------------------------------

constructor_pattern = re.compile(
    r'\b\w+\s*\(\s*"([^"]+)"\s*,\s*"([^"]+\.txt)"\s*\)',
    re.DOTALL
)

assignment_pattern = re.compile(
    r'\b\w*name\w*\s*=\s*"([^"]+)"',
    re.IGNORECASE
)

string_ctor_pattern = re.compile(
    r'\b\w*name\w*\s*\(\s*"([^"]+)"\s*\)',
    re.IGNORECASE
)

# remove comments
comment_pattern = re.compile(
    r'//.*?$|/\*.*?\*/',
    re.DOTALL | re.MULTILINE
)

extensions = {".cpp", ".cc", ".cxx", ".hpp", ".h"}

required = set()

# ---------------------------------------------------------
# scan source files
# ---------------------------------------------------------

for root, dirs, files in os.walk("."):
    for file in files:

        if os.path.splitext(file)[1] not in extensions:
            continue

        path = os.path.join(root, file)

        try:
            text = open(path, encoding="utf-8").read()
        except Exception:
            continue

        # remove comments
        text = comment_pattern.sub("", text)

        # constructor calls
        for name, txtfile in constructor_pattern.findall(text):

            required.add(name)
            required.add(txtfile)

            base = txtfile[:-4]

            required.add(base + "-variables.txt")
            required.add(base + "-configuration.txt")
            required.add(base + "-regions.txt")

        # assignments
        for name in assignment_pattern.findall(text):
            required.add(name)

        # string constructors
        for name in string_ctor_pattern.findall(text):
            required.add(name)

# ---------------------------------------------------------
# detect existing files in directory
# ---------------------------------------------------------

existing = set()

for f in os.listdir("."):
    if os.path.isfile(f):
        existing.add(f)

missing = sorted(required - existing)
unused = sorted(existing - required)

# ---------------------------------------------------------
# output
# ---------------------------------------------------------

print("\nREQUIRED INPUTS\n----------------")
for f in sorted(required):
    print(f)

#print("\nMISSING FILES\n-------------")
#for f in missing:
#    print(f)

#print("\nUNUSED FILES IN DIRECTORY\n-------------------------")
#for f in unused:
#    print(f)        
    
# missing file check (run from command line)
# ./list_test_inputs.py | while read f; do
#     [ -e "$f" ] || echo "MISSING: $f"
# done