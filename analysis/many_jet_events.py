import os
import ROOT

def count_instances(filename):
    # Open the root file
    file = ROOT.TFile(filename)
    # Access the tree
    tree = file.Get("T")
    
    # Define the criteria
    criteria = "nJet > 7"
    
    # Use the Draw method to apply the selection and count entries
    tree.Draw(">>event_list", criteria, "entrylist")
    event_list = ROOT.gDirectory.Get("event_list")
    count = event_list.GetN()
    
    return count

def process_root_files(directory):
    results = {}
    
    for filename in os.listdir(directory):
        if filename.startswith("output_47360") and filename.endswith('.root'):
            filepath = os.path.join(directory, filename)
            count = count_instances(filepath)
            print(f"{filename}: {count} instances")
            results[filename] = count
    
    return results

# Example usage:
# Assuming root files are in the 'root_files' directory
directory = '/sphenix/tg/tg01/jets/egm2153/JetValOutput'
results = process_root_files(directory)

totalcount = 0
n = 0

# Print the results
for filename, count in results.items():
    print(f"{filename}: {count} instances")
    totalcount += count
    n += 1

print(totalcount/n)
