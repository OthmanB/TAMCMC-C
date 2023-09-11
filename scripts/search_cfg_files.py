import os
from termcolor import colored

def search_cfg_files(main_directory, cfg_file, sentences, cfg_dir="Config"):
    for subdir in os.listdir(main_directory):
        full_path = os.path.join(main_directory, subdir)
        if cfg_dir != "" :
            cfg_dir0 = os.path.join(full_path, cfg_dir)
        else:
            cfg_dir0=full_path
        if os.path.isdir(cfg_dir0):
            cfg_file_path = os.path.join(cfg_dir0, cfg_file)
            if os.path.isfile(cfg_file_path):
                print("File:", colored(subdir + "/" + cfg_file, 'red'), "\n")
                try:   
                    with open(cfg_file_path, "r") as file:
                        for line in file:
                            line = line.lstrip()
                            if not line.startswith("#"):
                                if ";" in line:
                                    line = line[:line.index(";")]
                                for sentence in sentences:
                                    if sentence in line:
                                        try:
                                            equal_index = line.index("=")
                                            before_equal = line[:equal_index + 1]
                                            after_equal = line[equal_index + 1:]
                                            colored_before_equal = colored(before_equal, 'blue')
                                            colored_after_equal = colored(after_equal, 'yellow')
                                            print("      " + colored_before_equal + colored_after_equal.rstrip())
                                        except ValueError:
                                            print("      " + colored(line, "cyan"))
                                        break
                except Exception as e:
                    print("Error while reading file:", str(e))
                print(colored("="*30 + "\n", "red"))