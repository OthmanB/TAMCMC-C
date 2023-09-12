import os
import re
from termcolor import colored

def detect_finished_files(main_directory, file_prefix, file_extension):
    finished_files = []
    not_finished_files = []
    not_finished_info=[]
    # Traverse through all the directories and subdirectories
    for root, dirs, files in os.walk(main_directory):
        for filename in files:
            # Check if the file has the correct prefix and extension
            if filename.startswith(file_prefix) and filename.endswith(file_extension):
                file_path = os.path.join(root, filename)
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    # Check the last 4 lines of the file
                    last_lines = lines[-4:]
                    if any("Calculation finished in:" in line for line in last_lines):
                        # File has finished
                        finished_files.append((filename, file_path))
                    else:
                        # File has not finished
                        not_finished_files.append((filename, file_path))
                        info={}
                        i=0
                        for line in lines:
                            #if "swaping rate" in line:
                            #    info["swapping_rate"] =line.strip()
                            if "Total time so far" in line:
                                # Only store the last occurrence of the term
                                info["total_time"] = line.strip()
                                i_last_samples=i-1
                            if "Processing rate:" in line:
                                # Only store the last occurrence of the term
                                info["processing_rate"] = line.strip()
                            i=i+1
                        info["processed_samples"]=lines[i_last_samples].strip()
                        info["last_line"]=lines[-1].strip()
                        not_finished_info.append(info)
    return finished_files, not_finished_files, not_finished_info

 
def display_file_paths(files_list, info=None):
    # Create a dictionary to store the branches
    branches = {}
    max_shift_spaces = 0  # Variable to track the maximum shift spaces needed
    for i, file_tuple in enumerate(files_list):
        filename, file_path = file_tuple
        # Normalize file path for consistency
        file_path = os.path.normpath(file_path)
        # Split file_path into directory components
        directories = file_path.split(os.sep)
        # Get the key for the branch
        branch_key = os.sep.join(directories[:-1])
        # Get the branch from the dictionary, or create a new one if it doesn't exist
        if branch_key in branches:
            branch = branches[branch_key]
        else:
            branch = []
            branches[branch_key] = branch
        # Add the file to the branch
        if info is not None:
            branch.append((directories[-1], info[i]))
            shift_spaces = len(directories[-1]) - len(filename) - 4
            max_shift_spaces = max(max_shift_spaces, shift_spaces)
        else:
            branch.append((directories[-1]))
    # Traverse the branches and display the tree
    for branch_key, branch in branches.items():
        # Split the branch_key into directories
        directories = branch_key.split(os.sep)
        num_directories = len(directories)
        # Print the graphical representation of the branch
        for i, directory in enumerate(directories):
            if i == num_directories - 1:
                # Last directory
                print(f"{directory}/")
            elif i == 0:
                # First directory
                print(f"+--{directory}/")
            else:
                # Intermediate directories
                print(f"|  " * i + f"+--{directory}/")
        # Print the files in the branch
        if info is not None:
            shift_spaces=num_directories
            for file_name, info_name in branch:
                print(f"|  " * num_directories + f"+--{file_name}")
                shift_spaces = len(file_name) - len("slurm") + max_shift_spaces
                try:
                    line = info_name["processed_samples"]
                    value = re.findall(r'\[(\d+)\]', line)[0]
                    colored_line = colored(f'Number of samples processed : {value}', 'yellow')
                    print(f"|  " * (num_directories + 1) + " " * shift_spaces + colored_line)
                except:
                    colored_line = colored("Info on swapping rate not found", "red")
                    print(f"|  " * (num_directories + 1) + " " * shift_spaces + colored_line)
                #
                try:
                    line = info_name["total_time"]
                    colored_line = colored(line, 'yellow')
                    print(f"|  " * (num_directories + 1) + " " * shift_spaces + colored_line)
                except:
                    colored_line = colored("Info on total time not found", "red")
                    print(f"|  " * (num_directories + 1) + " " * shift_spaces + colored_line)
                #
                try:
                    line = info_name["processing_rate"]
                    value = re.findall(r'\((.*?)\)', line)[0]
                    colored_line = line.replace(f'({value})', colored(f'({value})', 'yellow'))
                    print(f"|  " * (num_directories + 1) + " " * shift_spaces + colored_line)
                except:
                    colored_line = colored("Info on processing rate not found", "red")
                    print(f"|  " * (num_directories + 1) + " " * shift_spaces + colored_line)
                line = info_name["last_line"]
                colored_line = colored("last line: " + line, 'blue')
                print(f"|  " * (num_directories + 1) + " " * shift_spaces + colored_line)
        else:
            for file_name in branch:
                print(f"|  " * num_directories + f"+--{file_name}")

'''
  Set root_dir as the highest level directory that contain all of the slurm files
  you want to scan
'''
root_dir="/scratch/ob19/ajAlm_project/"
file_prefix="slurm-"
finished_files, not_finished_files, not_finished_info= detect_finished_files(root_dir, file_prefix, ".out")

print(colored(" --------------------------------------------- ","red"))
print(colored("                 FINISHED FILES                ","red"))
print(colored(" --------------------------------------------- ","red"))
print(" ")
display_file_paths(finished_files)

print(colored(" --------------------------------------------- ","red"))
print(colored("               NOT FINISHED FILES              ","red"))
print(colored(" --------------------------------------------- ","red"))
print(" ")
display_file_paths(not_finished_files, info=not_finished_info)
