#include <vector>
#include "commands.hpp"
#include "cxxopts.hpp"
#include "organo_opts.hpp"
#include "cluster.hpp"

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_cluster_parser(int argc, char** argv){
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
      ("a, max-alleles", "Maximum number of alleles.", cxxopts::value<int>()->default_value("0"))
      ("s, min-sim", "Minimum sequence similarity for non-spanning assignment.", cxxopts::value<double>()->default_value("0.99"))

      ("c, min-cov", "Minimum coverage support (used for guiding cluster-merging).", cxxopts::value<int>()->default_value("2"))
      ("l, low-cov", "Automatically switch to 'low-coverage' mode when per-region cov is <= X.", cxxopts::value<int>()->default_value("6"))
      ("m, max-cov", "Ignore region above this coverage threshold", cxxopts::value<int>()->default_value("50"))

      ("g, group", "Group sequences based on region ('TG'-tag).", cxxopts::value<bool>()->default_value("false"))
      ("hap", "Group sequences based on haplotype ('HP'-tag).", cxxopts::value<bool>()->default_value("false"))
      ("h, hp", "Compress homopolymers.", cxxopts::value<bool>()->default_value("false"))
      ("d, hpd", "Compress homopolymers during distance calculations only.", cxxopts::value<bool>()->default_value("false"))
      
      ("tsv", "Output  pairwise sequence distance information in TSV format.", cxxopts::value<std::string>()->default_value(""))
      ("matrix", "Output  pairwise sequence distance information in matrix format.", cxxopts::value<std::string>()->default_value(""))
      
      ("t, threads", "Total number of threads.", cxxopts::value<int>()->default_value("4"))
      ("block", "Number of target regions to load in memory (used for multithreading).", cxxopts::value<uint32_t>()->default_value("50000"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    //no BAM files provided, output help message
    if(inputs.empty()) std::cout << options.help();
    else {
      organo_opts params;
      params.group = result["group"].as<bool>();
      params.haplotype = result["hap"].as<bool>();
      params.hp = result["hp"].as<bool>();
      params.hpd = result["hpd"].as<bool>();
      params.block_size = result["block"].as<uint32_t>();
      params.tsv = result["tsv"].as<std::string>();
      params.matrix = result["matrix"].as<std::string>();
      params.maxalleles = result["max-alleles"].as<int>();

      //process mincov accordingly
      int mincov = result["min-cov"].as<int>();
      if(mincov >= 0) params.mincov = mincov;
      else{
        std::cout << "Min-coverage must be >= 0: " << mincov << '\n' << options.help() <<  std::endl;
        exit(1);
      }

      //process low-cov accordingly
      int lowcov = result["low-cov"].as<int>();
      if(lowcov >= 0) params.lowcov = lowcov;
      else{
        std::cout << "Low-coverage must be >= 0: " << lowcov << '\n' << options.help() <<  std::endl;
        exit(1);
      }

      //process minsim accordingly
      double minsim = result["min-sim"].as<double>();
      if(minsim >= 0.0 || minsim <= 1.0) params.minsim = minsim;
      else{
        std::cout << "Min-similarity must be [0.0, 1.0]: " << minsim << '\n' << options.help() <<  std::endl;
        exit(1);
      }

      //process maxerror accordingly
      double maxerror = result["error"].as<double>();
      if(maxerror >= 0.0 || maxerror <= 1.0) params.maxerror = maxerror;
      else{
        std::cout << "Min-similarity must be [0.0, 1.0]: " << maxerror << '\n' << options.help() <<  std::endl;
        exit(1);
      }

      //process max seqs accordingly
      int max = result["max-cov"].as<int>();
      if(max > 0) params.maxcov = max;
      else{
        std::cout << "Max number of sequences must be >0" << max << '\n' << options.help() <<  std::endl;
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

      cluster(inputs.front(), params);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "Error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}

int command_cluster(int argc, char **argv){
  //parse CLI arguments
  command_cluster_parser(argc, argv);
  return 0;
}
