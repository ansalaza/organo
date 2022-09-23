#include <iostream>
#include <unistd.h>
#include <string>
#include "commands.hpp"


/**
 *  Function to pring out general (help) message of gvt
 */
void printHelp(){
   std::cout << ("Usage:\n organo [command]") << std::endl;
   std::cout << ("      genotype      Cluster target-sequences in given FASTA file and genotype.\n") << std::endl;
 }

/**
 * Main method. Parse input commands
 */
int main(int argc, char **argv){
  if(argc == 1) printHelp();
  else{
    if(std::string(argv[1]) == "generate") command_generate(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "realign") command_realign(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "cluster") command_cluster(argc - 1, &argv[1]);
    else printHelp();
  }

  return 0;
}
