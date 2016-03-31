/*
 * Config.cpp
 *
 * Contains all kind of functions
 * used to process and/or encapsulate data
 * 
 *  Created on: 20 Mar 2016
 *      Author: obenomar
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "config.h"

Config::Config(string data_file_in, string MALA_cfg_file_in, string output_cfg_file_in){ // The constructor

	data.data_file=data_file_in;
	MALA.MALA_cfg_file=MALA_cfg_file_in;
	outputs.output_cfg_file=output_cfg_file_in;

	bool verbose_data=0;
	string delimiter=" ";
	int x_ind=0; // column index of the x values
        int y_ind=1; // column index of the y values
	int ysig_ind=-1; // negative index means that it is not provided
	MALA.proposal_type="Random"; // Exact_on_grid does not yet exist

	Data_Nd data_in=read_data_ascii_Ncols(data.data_file, delimiter, verbose_data);
	if(x_ind >=0){
		data.data.x=data_in.data.col(x_ind);
		data.data.xlabel=data_in.labels[x_ind];
		data.data.xunit=data_in.units[x_ind];
	}
	if(y_ind >=0){
		data.data.y=data_in.data.col(y_ind);
		data.data.ylabel=data_in.labels[y_ind];
		data.data.yunit=data_in.units[y_ind];
	}
	if(ysig_ind >= 0){
		data.data.sigma_y=data_in.data.col(ysig_ind);
	} else {
		std::cout << "Warning: Config::Config.data.data.sigma_y was not specified in the data file" << std::endl;
		std::cout << "         In case of a chi_square fit, this variable is required" << std::endl;
		std::cout << "         To avoid code failure, this value will be fixed to 1" << std::endl;
		std::cout << "         If you wish to use sigma_y != 1, then add a column in you data file: " << data_file_in << std::endl;
		VectorXd tmp(data.data.x.size());
		data.data.sigma_y=tmp.setConstant(1);
	}
	data.data.header=data_in.header;
	data.data.Nx=data.data.x.size();

}


Data_Nd Config::read_data_ascii_Ncols(const string file_in_name, const string delimiter, const bool verbose_data){
/*
 * This function read an input file (file_in_name) which may or may not contain a header, labels, units indicator.
 * The number of column can be as small as 1
*/

    int cpt, Nrows;
    string token, line, line0, subline0;
    vector<string> header, labels, units, data_str;
    long double tmp_val;
    Eigen::MatrixXd data;
    Data_Nd all_data_out; // The structure that encapsulate all the data, the header, the labels and units
    int data_Maxsize=1000000;
 
    ifstream file_in;
    std::cout << "  Assumption of this program: " << std::endl;
    std::cout << "       - header lines appear on the top of the file and are indicated by a # for first character. The header is however optional (then no # are found)" << std::endl;
    std::cout << "       - labels appear in one single line that is just after the  header and is indicated by a ! for first character. Labels are optional (then no ! are found)" << std::endl;
    std::cout << "       - units appear in one single line that is just after the  labels or the header (if the labels are missing) and is indicated by a * for first character. Units are optional (then no * are found)" << std::endl;
    std::cout << "       - Maximum number of lines for the data: " << data_Maxsize << std::endl;
    file_in.open(file_in_name.c_str());
    if (file_in.is_open()) {
	std::cout << "File opened... processing lines" << std::endl;

		// [1] Get the header
		cpt=0;
		std::getline(file_in, line0);
		line0=strtrim(line0); // remove any white space at the begining/end of the string
		subline0=strtrim(line0.substr(0, 1)); // pick the first character
		subline0=subline0.c_str();
		if (subline0 == "#"){
			while(subline0 == "#"){		
				header.push_back(strtrim(line0.substr(1, std::string::npos))); // add all characters except the first one (and any space at the begining)

				std::getline(file_in, line0);
				line0=strtrim(line0); // remove any white space at the begining/end of the string
				subline0=strtrim(line0.substr(0, 1)); // pick the first character
				subline0=subline0.c_str();
				
				cpt=cpt+1;
			}
			std::cout << "   [1] " << cpt << " header lines found..." << std::endl;
		} else{
			header.push_back("");
			std::cout << "   [1] Header not found. Header vector set to a blank vector<string> of size 1. Pursuing operations..." << std::endl;
		}

		// [2] Read the labels... these are expected just after the header... If not found, then labels is left blank of size 1
		// No need to read line... already read by the previous loop
		line0=strtrim(line0); // remove any white space at the begining/end of the string
		subline0=strtrim(line0.substr(0, 1)); // pick the first character
		subline0=subline0.c_str();
		if(subline0 == "!"){
			
			labels=strsplit(strtrim(line0.substr(1,std::string::npos)), " \t"); // remove either when you found a white space or a tabulation
			std::cout << "   [2] " << labels.size() << " labels found..." << std::endl;

			for(int i=0; i<labels.size();i++){
				labels[i]=strtrim(labels[i]); // remove spaces at begining/end of each substring
				std::cout << "         - " << labels[i] << std::endl;
			}

		} else {
			labels.push_back("");
			std::cout << "   [2] No labels found. Label vector set to a blank vector<string> of size 1. Pursuing operations..." << std::endl;
		}
		
		// [3] Read the units... these are expected just after the labels... If not found, then units is left blank of size 1
		// Read a new line only if we found a label. Otherwise, it is not necessary...
		if (labels[0] != ""){std::getline(file_in, line0);}
		line0=strtrim(line0); 
		subline0=strtrim(line0.substr(0, 1)); 
		subline0=subline0.c_str();
		if(subline0 == "*"){
			units=strsplit(strtrim(line0.substr(1,std::string::npos)), " \t"); 
			std::cout << "   [3] " << units.size() << " units found..." << std::endl;
			for(int i=0; i<units.size();i++){
				units[i]=strtrim(units[i]);
				std::cout << "         - " << units[i] << std::endl;
			}
		} else {
			units.push_back("");
			std::cout << "   [3] No units found. Units vector set to a blank vector<string> of size 1." << std::endl;
		}

		// [4] Read the data...
		std::cout <<  "   [4] Now processing the data..." << std::endl;
		if (labels[0] != "" || units[0] != ""){  // case where we need to read a new line before looping
			std::getline(file_in, line0);
		} 
		Nrows=0;
	   while(!file_in.eof()){
		data_str=strsplit(strtrim(line0), " \t"); 
		if (Nrows == 0) {
			data.resize(data_Maxsize, data_str.size());
		}
		for(int i=0; i<data_str.size();i++){
			if ( ! (istringstream(data_str[i]) >> tmp_val) ){tmp_val = nan("");} // If the number can be converted, then tmp_val=value. Otherwise tmp_val = NaN
			data(Nrows, i)=tmp_val;		
		}
		if (verbose_data == 1) {std::cout << data.row(Nrows) << std::endl;} // Show all entries only if requested
		std::getline(file_in, line0);
		Nrows=Nrows+1;
	    }
	file_in.close();
	data.conservativeResize(Nrows, data_str.size());
	std::cout << "         - Number of lines found: " << Nrows << std::endl;
	std::cout << "         - Number of columns found: " << data_str.size() << std::endl;
     } else {
	std::cout << "Could not open the file!" << std::endl;
	std::cout << "use the full path of the file in order to read it without issue " << std::endl;
	exit(EXIT_FAILURE);
     }

     all_data_out.data=data;
     all_data_out.header=header;
     all_data_out.labels=labels;
     all_data_out.units=units;

return all_data_out;
}


std::string Config::strtrim(const std::string& str){
/*
 * Small function that remove white space at the end and at
 * the begining of a string. 
 * The original program was taken from http://stackoverflow.com/questions/1798112/removing-leading-and-trailing-spaces-from-a-string
 * Modified in order to not use the C++11 standard:
 *  	- const auto --> const string
 *      - optional argument for the separator is now hardcoded
*/
    string whitespace = " \t";
    size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    size_t strEnd = str.find_last_not_of(whitespace);
    size_t strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

vector<std::string> Config::strsplit(const std::string str, const std::string delimiters){
/*
 * Take a string and split it each time one of the listed delimiters is detected
*/
	string str0=strtrim(str);
	size_t pos=0;
	vector<string> str_splitted;

	while ((pos = str0.find_first_of(delimiters)) != std::string::npos) {
		    
		str_splitted.push_back(str0.substr(0, pos)); // get the substring
		str0.erase(0, pos + delimiters.length());
		str0=strtrim(str0); // remove any extra white space at the begining before iterating
	}
	str_splitted.push_back(str0); // do not forget to add the end of the string

return str_splitted;
}

inline bool Config::file_exists(const std::string& name) {
/* 
 * Test whether a file exists
*/
    return ( access( name.c_str(), F_OK ) != -1 );
}

