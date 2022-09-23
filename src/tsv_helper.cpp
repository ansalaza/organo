#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "antimestamp.hpp"
#include "tsv_helper.hpp"

void load_ha_tagged_reads(const std::string& file, std::unordered_map<std::string, int>& map){
  std::ifstream inputfile(file.c_str());
  std::string line;
  //iterate through end of file
  while (std::getline(inputfile, line)){
  	int index = 0;
  	std::string read;
  	std::string hap;
		std::string value;
		std::istringstream stream(line);
		while(std::getline(stream, value, '\t')){
			switch(index){
		      case 0: {
		        read = value;
		        break;
		      }
		      case 1: {
		        hap = value;
		        break;
		      }
		      default: {
		        break;
		      }
		    }
		    ++index;
		}
		if(index < 2) std::cerr << "(" << antimestamp() << "): WARNING: unexpected number of columns in Haplotype-tagged line: " << line << std::endl;  
		else if(hap != "NA"){
			auto it = map.find(read);
			if(it != map.end()) std::cerr << "(" << antimestamp() << "): WARNING: read " << read << " already tagged. Overwriting... " << std::endl;
			map.insert({read, std::stoi(hap)});
		}
  }
  inputfile.close();
}