//
// Created by nloyfer on 5/27/19.
//

#ifndef PATS_READS_LENS_PAT2RLEN_H
#define PATS_READS_LENS_PAT2RLEN_H

#include <boost/algorithm/string/predicate.hpp>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <string>
#include <stdlib.h>

#define MAX_LEN 11

#define SEP "\t"
#define DEFAULT_BLOCKS "/cs/cbio/netanel/blocks/outputs/nps20_genome.tsv.gz"
#define DEFAULT_MARKOV_ORDER "5"

class Pat2rlen {
    int morder = 0;
    int32_t *counts;
    int32_t **marr;

    bool debug;
    bool compress;
    std::string blocks_path;
    std::string output_prefix;
    std::vector<int> borders_starts;
    std::vector<int> borders_ends;
    int nr_blocks = 0;
    int cur_block_ind = 0;

    std::vector<int> nr_states;

    int read_blocks();

    void dump(int *data, int width, std::string out_path);

    int proc_line(std::vector <std::string> tokens);

    void update_block(int ind, std::string pat, int count);

    int32_t* init_array(int len);

    void debug_print(int ind, std::vector<std::string> tokens);

    int blocks_helper(std::istream &instream);

public:
    Pat2rlen(std::string rlen_path, std::string blocks_path, int m_order, bool compress, bool deb);

    ~Pat2rlen();

    void parse(std::string pat_path);
};

bool hasEnding (std::string const &fullString, std::string const &ending);
#endif //PATS_READS_LENS_PAT2RLEN_H
