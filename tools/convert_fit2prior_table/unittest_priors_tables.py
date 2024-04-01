import unittest
import time
import os
import sys
#from io import StringIO
from PIL import Image, ImageChops
#from contextlib import redirect_stdout

from make_priors_tables import write_table, make_prior_table

class TestClassMain(unittest.TestCase):
    def assertImagesEqual(self, image_path1, image_path2):
        image1 = Image.open(image_path1)
        image2 = Image.open(image_path2)
        diff = ImageChops.difference(image1, image2)
        self.assertEqual(diff.getbbox(), None, "The generated plot is incorrect")

    def testWriteTable1D(self):
        xvals = [1, 2, 3, 4, 5]
        yvals = [0.1, 0.2, 0.3, 0.4, 0.5]
        label1 = "x"
        unit1 = "m"
        output_file = "test_data/Outputs/test_output1D.txt"
        process_name = "Process1"
        #
        expected_data = ["1.000000000000\t0.100000000000\n", 
                         "2.000000000000\t0.200000000000\n",
                         "3.000000000000\t0.300000000000\n",
                         "4.000000000000\t0.400000000000\n",
                         "5.000000000000\t0.500000000000\n"]

        err = write_table(xvals, yvals, label1, unit1, output_file, process_name)
        with open(output_file, "r") as f:
            output_lines = f.readlines()
        self.assertEqual(output_lines[0],"# Tabulated priors for Process1\n")
        self.assertEqual(output_lines[1],"# Table created by TAMCMC/tools/convert_fit2prior_table/make_priors_tables.py\n")
        self.assertEqual(output_lines[2],"# dimensions: 1\n")
        self.assertEqual(output_lines[3], "!\t"+label1+"\t PDF\n")
        self.assertEqual(output_lines[4], "*\t"+unit1+"\t (no_unit)\n")
        self.assertEqual(output_lines[5:], expected_data)

    def testWriteTable2D(self):
        xvals = [1, 2, 3, 4, 5]
        yvals = [0., 45, 90]
        zvals = [[1.1, 1.2, 1.3, 1.4, 1.5], [2.1, 2.2, 2.3, 2.4, 2.5], [3.1, 3.2, 3.3, 3.4, 3.5]]
        label1 = "x"
        unit1 = "(Hz)"
        label2 = "y"
        unit2 = "(deg)"
        output_file = "test_data/Outputs/test_output2D.txt"
        process_name = "Process2"
        #
        expected_data = ["{:12}\t1.000000000000\t2.000000000000\t3.000000000000\t4.000000000000\t5.000000000000\t\n".format("NA"),
                         "0.000000000000\t1.100000000000\t1.200000000000\t1.300000000000\t1.400000000000\t1.500000000000\t\n",
                         "45.000000000000\t2.100000000000\t2.200000000000\t2.300000000000\t2.400000000000\t2.500000000000\t\n",
                         "90.000000000000\t3.100000000000\t3.200000000000\t3.300000000000\t3.400000000000\t3.500000000000\t\n"]
        write_table(xvals, yvals, label1, unit1, output_file, process_name, zvals=zvals, label2=label2, unit2=unit2)
        with open(output_file, "r") as f:
            output_lines = f.readlines()
        self.assertEqual(output_lines[0],"# Tabulated priors for Process2\n")
        self.assertEqual(output_lines[1],"# Table created by TAMCMC/tools/convert_fit2prior_table/make_priors_tables.py\n")
        self.assertEqual(output_lines[2],"# dimensions: 2\n")
        self.assertEqual(output_lines[3], "!\t"+label1+"\t"+label2+"\t PDF\n")
        self.assertEqual(output_lines[4], "*\t"+unit1+"\t"+unit2+"\t (no_unit)\n")
        self.assertEqual(output_lines[5:], expected_data)

    def testMakePriorTable1D(self):
        binning = 30
        work_dir = os.getcwd()
        process_name = "10280410_Gaussfit"
        dir_tamcmc_outputs = work_dir + "/test_data/" #+ process_name + "/outputs/"
        output_dir = work_dir + "/test_data/Outputs/"
        index1 = 7
        period = 5
        #cpp_prg = "../../bin/"
        cpp_prg = "../bin/"
        #
        make_prior_table(index1, binning, output_dir, dir_tamcmc_outputs, process_name, phase='A', chain=0,
                         first_index=0, last_index=-1, period=period, index2=None, cpp_prg=cpp_prg)
        #
        expected_image_path = work_dir + "/test_data/Expectations/10280410_Gaussfit1d_EXPECTED.jpg"
        generated_image_path = output_dir + "10280410_Gaussfit_0.priors.jpg"
        self.assertImagesEqual(expected_image_path, generated_image_path)

    def testMakePriorTable2D(self):
        binning = 30
        work_dir = os.getcwd()
        process_name = "10280410_Gaussfit"
        dir_tamcmc_outputs = work_dir + "/test_data/" #+ process_name + "/outputs/"
        output_dir = work_dir + "/test_data/Outputs/"
        index1 = 8
        index2 = 7
        period = 5
        #cpp_prg = "../../bin/"
        cpp_prg = "../bin/"
        #
        # Add a check to wait for the file to be created
        generated_image_path = output_dir + "10280410_Gaussfit_0.priors.jpg"
        while not os.path.isfile(generated_image_path):
            time.sleep(1)  # Wait for 1 second before checking again

        make_prior_table(index1, binning, output_dir, dir_tamcmc_outputs, process_name, phase='A', chain=0,
                         first_index=0, last_index=-1, period=period, index2=index2, cpp_prg=cpp_prg)
        #
        expected_image_path = work_dir + "/test_data/Expectations/10280410_Gaussfit2d_EXPECTED.jpg"
        generated_image_path = output_dir + "10280410_Gaussfit_1.priors.jpg"
        self.assertImagesEqual(expected_image_path, generated_image_path)

if __name__ == "__main__":
    work_dir = os.getcwd()
    output_dir = work_dir + "/test_data/Outputs/"
    print("Clearing the output test directory...")
    # Iterate over all the files and subdirectories in the given directory
    for root, dirs, files in os.walk(output_dir):
        # Delete each file
        for file in files:
            file_path = os.path.join(root, file)
            os.remove(file_path)
        # Delete each subdirectory
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            os.rmdir(dir_path)
    time.sleep(1)  # Wait for 1 second before checking again
    print("Running tests...")
    unittest.main()