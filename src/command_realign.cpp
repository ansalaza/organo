#include <vector>
#include "commands.hpp"
#include "cxxopts.hpp"
#include "realign.hpp"
#include "organo_opts.hpp"

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_realign_parser(int argc, char** argv){
  try{
    cxxopts::Options options(argv[0]);
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <FASTA>")
      .add_options()
      ("e, error", "Maximum sequence error.", cxxopts::value<double>()->default_value("0.05"))
      ("s, minsim", "Minimum sequence similarity for non-spanning assignment.", cxxopts::value<double>()->default_value("0.99"))
      ("k, kmer-size", "Matching kmer-size for realignments.", cxxopts::value<int>()->default_value("5"))

      ("m, max-cov", "Ignore target regions with more than this number of sequences", cxxopts::value<int>()->default_value("50"))
      ("h, hp", "Compress homopolymers.", cxxopts::value<bool>()->default_value("false"))
      ("d, hpd", "Compress homopolymers during distance calculations only.", cxxopts::value<bool>()->default_value("false"))
      
      ("t, threads", "Total number of threads.", cxxopts::value<uint32_t>()->default_value("4"))
      ("block", "Number of target regions to load in memory (used for multithreading).", cxxopts::value<uint32_t>()->default_value("50000"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    //no BAM files provided, output help message
    if(inputs.empty()) std::cout << options.help();
    else {
      organo_opts params;
      params.hp = result["hp"].as<bool>();
      params.hpd = result["hpd"].as<bool>();
      params.block_size = result["block"].as<uint32_t>();
      

      //process minsim accordingly
      double minsim = result["minsim"].as<double>();
      if(minsim >= 0.0 || minsim <= 1.0) params.minsim = minsim;
      else{
        std::cout << "Min-similarity must be [0.0, 1.0]: " << minsim << '\n' << options.help() <<  std::endl;
        exit(1);
      }

      //process maxerror accordingly
      double maxerror = result["error"].as<double>();
      if(maxerror >= 0.0 || maxerror <= 1.0) params.maxerror = maxerror;
      else{
        std::cout << "Max error must be [0.0, 1.0]: " << maxerror << '\n' << options.help() <<  std::endl;
        exit(1);
      }

      //process kmer-size accordingly
      int k = result["kmer-size"].as<int>();
      if(k > 0) params.k = k;
      else{
        std::cout << "Kmer-size must be >0: " << k << '\n' << options.help() <<  std::endl;
        exit(1);
      }

      //process max seqs accordingly
      int max = result["max-cov"].as<int>();
      if(max > 0) params.maxcov = max;
      else{
        std::cout << "Max coverage must be >0" << max << '\n' << options.help() <<  std::endl;
        exit(1);
      }
    
      
      //process threads accordingly
      int threads = result["threads"].as<uint32_t>();
      if(threads < 0){
        std::cout << "No. of threads must be >0 " << '\n' << options.help() <<  std::endl;
        exit(0);
      }
      else if(threads == 1) params.threads = 1; 
      else params.threads = threads - 1;

      realign(inputs.front(), params);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "Error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}

int command_realign(int argc, char **argv){
  //parse CLI arguments
  command_realign_parser(argc, argv);
  return 0;
}
