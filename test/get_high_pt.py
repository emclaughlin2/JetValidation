import os
import ROOT
import fnmatch

def count_instances(filename):
    # Open the root file
    file = ROOT.TFile(filename)
    # Access the tree
    tree = file.Get("T")
    
    # Define the criteria
    criteria = "pt > 60"
    
    # Use the Draw method to apply the selection and count entries
    #tree.Draw(">>event_list", criteria, "entrylist")
    #event_list = ROOT.gDirectory.Get("event_list")
    #count = event_list.GetN()

    # Use the GetEntries method with a selection string to count entries
    count = tree.GetEntries(criteria)
    
    return count

def process_root_files(directory):
    results = {}
    
    for filename in os.listdir(directory):
        if fnmatch.fnmatch(filename, 'output_47289_*.root'):
            print(f'{filename}')
            filepath = os.path.join(directory, filename)
            count = count_instances(filepath)
            results[filename] = count
            print(f'{filename}: {count} instances')
    
    return results

# Example usage:
# Assuming root files are in the 'root_files' directory
directory = '/sphenix/tg/tg01/jets/egm2153/JetValOutput'
results = process_root_files(directory)

# Print the results
for filename, count in results.items():
    print(f"{filename}: {count} instances")
