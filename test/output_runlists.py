import os
import re

# Directory containing the files
directory = '/sphenix/tg/tg01/jets/egm2153/JetValOutput'

# Regular expression to match the numbers after 'output'
pattern = re.compile(r'output_(\d+)_')

# Set to store the extracted numbers
numbers_set = set()

# Iterate through all files in the directory
for filename in os.listdir(directory):
    match = pattern.search(filename)
    if match:
        # Add the extracted number to the set
        numbers_set.add(int(match.group(1)))

print(sorted(numbers_set))
