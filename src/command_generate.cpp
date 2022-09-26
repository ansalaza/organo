#include <vector>
#include "commands.hpp"
#include "cxxopts.hpp"
#include "generate.hpp"
#include "organo_opts.hpp"

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_generate_parser(int argc, char** argv){
  try{
    cxxopts::Options options(argv[0]);
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <BAM>")
      .add_options()
      ("b, bed", "BED-formatted file of target regions.", cxxopts::value<std::string>())
      ("hap", "Two-column TSV file specifying haplotype-tagged reads.", cxxopts::value<std::string>()->default_value(""))
      ("fasta", "Output in FASTA format.", cxxopts::value<bool>()->default_value("false"))
      ("m, mapq", "Minimum mapping quality.", cxxopts::value<uint8_t>()->default_value("20"))
      ("t, threads", "Total number of threads.", cxxopts::value<uint32_t>()->default_value("4"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    //no BAM files provided, output help message
    if(inputs.empty()) std::cout << options.help();
    else {
      organo_opts params;
      const std::string bed = result["bed"].as<std::string>();
      const std::string hap = result["hap"].as<std::string>();
      params.fasta = result["fasta"].as<bool>();

      //process mapq accordingly
      uint8_t mapq = result["mapq"].as<uint8_t>();
      if(mapq >= 0 && mapq <= 60) params.min_mapq = mapq;
      else{
        std::cout << "Invalid minimum mapq: " << (int)mapq <<  options.help() <<  std::endl;
        exit(0);
      }

      //process threads accordingly
      int threads = result["threads"].as<uint32_t>();
      if(threads < 0){
        std::cout << "No. of threads must be >0 " << '\n' << options.help() <<  std::endl;
        exit(0);
      }
      else if(threads == 1) params.threads = 1; 
      else params.threads = threads - 1;
      
      generate(inputs.front(), bed, hap, params);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "Error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}

int command_generate(int argc, char **argv){
  //parse CLI arguments
  command_generate_parser(argc, argv);
  return 0;
}
