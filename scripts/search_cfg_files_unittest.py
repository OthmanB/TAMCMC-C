import unittest
import os

from search_cfg_files import search_cfg_files

class SearchCfgFilesTestCase(unittest.TestCase):
    def test_search_cfg_files(self):
        main_directory = "tests/"
        cfg_dir = "Config"
        cfg_file = "test_cfg_file.txt"
        sentences = ["kplr001435467_kasoc-wpsd_slc_v1_21000_SK", "cfg_models_dir", "cfg_out_dir", "processing", "Nsamples",
                     "c0", "restore", "core_out", "core_in", 
                     "start_index_processing", "last_index_processing"]
        search_cfg_files(main_directory, cfg_file, sentences, cfg_dir=cfg_dir)

        # Add your assertions here

if __name__ == '__main__':
    unittest.main()

