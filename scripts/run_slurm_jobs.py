import os

def run_jobs_in_subdirs(main_directory, filename, command):
    '''
        function that scan a main_directory. If it finds the 
        filename, it runs the command
    '''
    for subdir in os.listdir(main_directory):
        full_path = os.path.join(main_directory, subdir)
        if os.path.isdir(full_path) and filename in os.listdir(full_path):
            os.chdir(full_path)
            os.system(command)

main_directory = "/scratch/ob19/ajAlm_project/"
filename = "job_array.sh"
command = "sbatch job_array.sh"

run_jobs_in_subdirs(main_directory, filename, command)
