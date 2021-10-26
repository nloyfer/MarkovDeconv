//
// Created by nloyfer on 5/27/19.
//

#include "pat2rlen.h"


/**
 *  Main
 *
 */

/*
 * Input Arguments Parsering class
 */
class InputParser {
public:
    InputParser(int &argc, char **argv) {
        for (int i = 1; i < argc; ++i)
            this->tokens.emplace_back(std::string(argv[i]));
    }

    const std::string &getCmdOption(const std::string &option) const {
        std::vector<std::string>::const_iterator itr;
        itr = std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
            return *itr;
        }
        static const std::string empty_string;
        return empty_string;
    }

    bool cmdOptionExists(const std::string &option) const {
        return std::find(this->tokens.begin(), this->tokens.end(), option)
               != this->tokens.end();
    }

private:
    std::vector <std::string> tokens;
};


std::string get_param_str(InputParser &input, std::string name, std::string defval) {
    std::string param_str = input.getCmdOption(name);
    if (!param_str.empty()) {
        return param_str;
    }
    return defval;
}


//void print_help() {
//    std::cerr << "" << std::endl;
//}


int main(int argc, char *argv[]) {

    if ((argc < 3)) {
        std::cerr << "Usage: EXEC PAT_PATH OUTPUT_PREFIX [-k MARKOV_ORDER] [-d] [-b BLOCKS_PATH]\n"
                  << "-k      markov order. Generate m0, m1,...,mk files. Default is k=5.\n"
                  << "-b      Blocks file. No header. Columns are:\n"
                  << "        <chr, start, end, startCpG, endCpG>. Default is nps20.\n"
                  << "        File may be gzipped. Note sometimes bgzip causes problems.\n"
                  << "--gzip  Compress output files.\n"
                  << "-d      Debug mode. Only use the first 2,500 blocks in the blocks file."
                  << std::endl;
        return -1;
    }

    InputParser input(argc, argv);

//    bool help = input.cmdOptionExists("-h");
//    if (help) {
//        print_help();
//        return 1;
//    }
    std::string pat_path = std::string(argv[1]);
    std::string output_path = std::string(argv[2]);
    std::string blocks_path = get_param_str(input, "-b", DEFAULT_BLOCKS);
    int markov_order = std::stoi(get_param_str(input, "-k", DEFAULT_MARKOV_ORDER));

    bool debug = input.cmdOptionExists("-d");
    bool compress_outputs = input.cmdOptionExists("--gzip");

    try {
        Pat2rlen(output_path, blocks_path, markov_order, compress_outputs, debug).parse(pat_path);
    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
    return 0;
}

