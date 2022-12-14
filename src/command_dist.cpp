#include <vector>
#include "commands.hpp"
#include "cxxopts.hpp"
#include "organo_opts.hpp"
#include "dist.hpp"

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_dist_parser(int argc, char** argv){
  try{
    cxxopts::Options options(argv[0]);
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <FASTA>")
      .add_options()
      ("c, min-cov", "Only consider alleles with at least this coverage.", cxxopts::value<int>()->default_value("0"))
      ("f, min-af", "Only consider alleles with at least this allele-frequency.", cxxopts::value<double>()->default_value("0"))
      ("h, hp", "Compress homopolymers.", cxxopts::value<bool>()->default_value("false"))
      ("d, hpd", "Compress homopolymers during distance calculations only.", cxxopts::value<bool>()->default_value("false"))
      ("t, threads", "Total number of threads.", cxxopts::value<int>()->default_value("4"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    //no BAM files provided, output help message
    if(inputs.size() < 2) std::cout << options.help();
    else {
      organo_opts params;
      params.hp = result["hp"].as<bool>();
      params.hpd = result["hpd"].as<bool>();

      //process mincov accordingly
      int mincov = result["min-cov"].as<int>();
      if(mincov >= 0) params.mincov = mincov;
      else{
        std::cout << "Min-coverage must be >= 0: " << mincov << '\n' << options.help() <<  std::endl;
        exit(1);
      }

      //process mincov accordingly
      double minaf = result["min-af"].as<double>();
      if(minaf >= 0) params.minaf = minaf;
      else{
        std::cout << "Min-af must be >= 0: " << minaf << '\n' << options.help() <<  std::endl;
        exit(1);
      }
      
      //process threads accordingly
      int threads = result["threads"].as<int>();
      if(threads < 0){
        std::cout << "No. of threads must be >0 " << '\n' << options.help() <<  std::endl;
        exit(0);
      }
      else if(threads == 1) params.threads = 1; 
      else params.threads = threads - 1;

      dist(inputs, params);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "Error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}

int command_dist(int argc, char **argv){
  //parse CLI arguments
  command_dist_parser(argc, argv);
  return 0;
}
