

#include "pat2rlen.h"


//  g++ -std=c++11 pat2rlen.cpp -o pat2rlen
//  g++ -std=c++11 pat2rlen.cpp -o pat2rlen -lz -lboost_iostreams


std::vector <std::string> line2tokens(std::string &line) {
    /**
     * Break string line to tokens, return it as a vector of strings
     */
    std::vector <std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while (getline(lineStream, cell, '\t'))
        result.push_back(cell);
    return result;
}

void print_vec(std::vector <std::string> &vec) {
    /** print a vector to stderr, tab separated */
    for (auto &j: vec)
        std::cerr << j << "\t";
    std::cerr << std::endl;
}


int32_t *Pat2rlen::init_array(int len) {
    int *arr = new int32_t[nr_blocks * len];
    std::fill_n(arr, nr_blocks * len, 0);
    return arr;
}

// Constructor
Pat2rlen::Pat2rlen(std::string in_output_prefix, std::string in_blocks_path, int m_order, bool comprs, bool deb) {
    morder = m_order;
    output_prefix = in_output_prefix;
    blocks_path = in_blocks_path;
    debug = deb;
    compress = comprs;

    int k = 1;
    for (int i = 0; i < morder + 1; i++) {
        nr_states.push_back(k *= 2);
    }

    // load blocks file
    int r = read_blocks();
    nr_blocks = borders_starts.size();

    // Init arrays to zeros
    counts = init_array(MAX_LEN);

    marr = new int32_t *[morder + 1];
    for (int i = 0; i < morder + 1; i++) {
        marr[i] = init_array(nr_states[i]);
    }
}

Pat2rlen::~Pat2rlen() {
    for (int i = 0; i < morder + 1; i++) {
        delete[] marr[i];
    }
    delete[] marr;
    delete[] counts;
}


int Pat2rlen::blocks_helper(std::istream &instream) {
    //Iterate lines
    std::vector <std::string> tokens;
    std::string line;
    int cur_start = 0, cur_end = 0;
    int abc = 0;
    while (std::getline(instream, line)) {
//        std::cerr << ++abc << std::endl;
        // skip empty lines and comments
        if (line.empty() || (!(line.rfind("#", 0)))) { continue; }

        tokens = line2tokens(line);
        if (tokens.size() < 5) {
            std::cerr << "Invalid blocks file format. ";
            std::cerr << "Should be: chr start end startCpG endCpG\n";
            abort();
        }

        // validate blocks:

        cur_start = std::stoi(tokens[3]);
        //  negative block start, unsorted blocks file
        if ((cur_start <= 0) || (cur_start < cur_end)) {
            std::cerr << "Invalid block start: " << cur_start << std::endl;
            std::cerr << "Make sure blocks are sorted by CpG-Index and monotonic (blockEnd > blockStart).\n";
            throw std::invalid_argument("Invalid blocks");
        }

        // If block end is smaller than its start:
        cur_end = std::stoi(tokens[4]);
        if ((cur_end <= cur_start)) {
            std::cerr << "Invalid block: " << cur_start << "\t" << cur_end << std::endl;
            throw std::invalid_argument("Invalid block");
        }
        borders_starts.push_back(cur_start);
        borders_ends.push_back(cur_end);
        if (debug && (borders_starts.size() >= 2500)) {
//        if (debug && (borders_starts.size() >= 20)) {
            break;
        }
    }
    return 0;
}

int Pat2rlen::read_blocks() {
    /**
     * Load blocks gzipped file into vector<int> borders_starts, borders_ends.
     */


    // std::cout << "loading blocks..." << std::endl;

    if (hasEnding(blocks_path, ".gz")) {
        // Open the gzipped file:
        std::ifstream file(blocks_path, std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_streambuf <boost::iostreams::input> inbuf;
        inbuf.push(boost::iostreams::gzip_decompressor());
        inbuf.push(file);
        std::istream instream(&inbuf);
        blocks_helper(instream);
        //Cleanup
        file.close();

    } else {
        std::ifstream instream(blocks_path);
        blocks_helper(instream);
        instream.close();
    }

    //std::cerr << "nr_blocks: " << borders_starts.size() << ", last block: ";
    //std::cerr << borders_starts[borders_starts.size() - 1] << "-" << borders_ends[borders_ends.size() - 1] << std::endl;
    return 0;
}


void Pat2rlen::dump(int32_t *data, int width, std::string out_path) {
    /**
     */
    std::ofstream bofs;
    bofs.open(out_path);
    if (!(bofs.is_open())) {
        std::cerr << "could not open output file " << out_path << std::endl;
        return;
    }

    std::cerr << "dumping to " << out_path << "... " << std::endl;
    for (int i = 0; i < nr_blocks; i++) {
        bofs << borders_starts[i] << SEP << borders_ends[i] << SEP;
        for (int j = 0; j < width; j++) {
            bofs << data[i * width + j];
            if (j < width - 1)
                bofs << SEP;
        }
        bofs << std::endl;
    }
    bofs.close();

    // gzip:
    if (compress) {
        system(("gzip " + out_path).c_str());
    }
}


void update_m2(int32_t *arr, int block_ind, std::string pat, int count, int nr_states, int log_ns) {

    int block_ind_times_nr_states = block_ind * nr_states;

    for (int i = log_ns - 1; i < pat.length(); i++) {
        bool valid_kmer = true;
        for (int k = 0; k < log_ns; k++) {
            if (pat[i - k] == '.') {
                valid_kmer = false;
                break;
            }
        }
        if (valid_kmer) {   // if current k-mer contains no dots
            arr[block_ind_times_nr_states + std::stoi(pat.substr(i - log_ns + 1, log_ns), nullptr, 2)] += count;
        }
    }
}

void Pat2rlen::update_block(int block_ind, std::string pat, int32_t count) {

//    std::cerr << block_ind << ") updating: " << pat << std::endl;

    // update lengths counts array:
    int len = pat.length();
    if (len == 0)
        return;
    len = (len > MAX_LEN) ? MAX_LEN : len;
    counts[block_ind * MAX_LEN + len - 1] += count;

    std::replace(pat.begin(), pat.end(), 'C', '1');
    std::replace(pat.begin(), pat.end(), 'T', '0');

    for (int i = 0; i < morder + 1; i++) {
        update_m2(marr[i], block_ind, pat, count, nr_states[i], i + 1);
    }
}

int Pat2rlen::proc_line(std::vector <std::string> tokens) {
    /**
     * Given one line of the form "chr1 1245345 CCT 3", update counts array -
     * Increase in the corresponding cells.
     */
    if (tokens.size() < 4) {
        throw std::invalid_argument("Invalid site in input file. too few columns");
    }

    int read_start = std::stoi(tokens[1]);
    std::string pattern = tokens[2];
    int count = std::stoi(tokens[3]);
    int read_end = read_start + (int) pattern.length() - 1;

//    print_vec(tokens);


    // read starts later than the last border - finish
    //                                          CTCTCT
    // |---|   |-----| |--|  ...  |-----| |-|
    if (read_start > borders_ends[borders_ends.size() - 1]) {
        //std::cerr << "Reads reached the end of blocks: " << read_start << " > "
        //          << borders_ends[borders_ends.size() - 1] << ". Breaking.\n";
        return 1;
    }

    // read starts after current block ends - move on to next relevant block
    //                   CTCTCT
    //         |-----|
    while (read_start >= borders_ends[cur_block_ind])
        cur_block_ind++;
    // todo: dump block line on the fly

    // read ends before current block starts: skip read.
    //  CTCTCT
    //         |-----|
    if (read_end < borders_starts[cur_block_ind])
        return 0;

    // read starts before current block, but continues to the block (and possibly beyond) - clip beginning of read
    //      CTCTCT
    //         |-----|
    if (read_start < borders_starts[cur_block_ind]) {
        pattern = pattern.substr(borders_starts[cur_block_ind] - read_start);
        read_start = borders_starts[cur_block_ind];
    }

    // If we reached here, current reads starts in cur_block_ind
    //        CTCTCT..CCTTTCTCT
    //     |-----|   |--|    |---|
    // tmp_block_ind index is running for the current read, from cur_block_ind until the read ends.
    int tmp_block_ind = cur_block_ind;
    while ((!pattern.empty()) && (tmp_block_ind < nr_blocks)) {
        if ((read_start >= borders_starts[tmp_block_ind]) && (read_start < borders_ends[tmp_block_ind])) {
            int head_size = std::min(borders_ends[tmp_block_ind] - read_start, (int) pattern.length());
            update_block(tmp_block_ind, pattern.substr(0, head_size), (int32_t) count);
            pattern = pattern.substr(head_size);
            read_start += head_size;
            tmp_block_ind++;
        }
        else if (read_end < borders_starts[tmp_block_ind]) {
            break;
        }
        else {
            pattern = pattern.substr(borders_starts[tmp_block_ind] - read_start);
            read_start = borders_starts[tmp_block_ind];
        }
    }
    return 0;
}

void Pat2rlen::debug_print(int ind, std::vector<std::string> tokens){
    std::cerr << ind << ") " << borders_starts[ind] << "-" << borders_ends[ind] << std::endl;
    print_vec(tokens);
}


bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

void Pat2rlen::parse(std::string pat_path) {

    try {

        // Open the gzipped pat file:
        //std::ifstream file(pat_path, std::ios_base::in | std::ios_base::binary);
        //boost::iostreams::filtering_streambuf <boost::iostreams::input> inbuf;
        //inbuf.push(boost::iostreams::gzip_decompressor());
        //inbuf.push(file);
        //std::istream instream(&inbuf);
        bool empty_pat = true;
        int line_i = 0;
        for (std::string line_str; std::getline(std::cin, line_str);) {
            if (line_str.empty()) { continue; } // skip empty lines
            empty_pat = false;

            std::vector <std::string> tokens = line2tokens(line_str);
            if (tokens.size() < 4) {
                print_vec(tokens);
                throw std::invalid_argument("too few columns in file, line " + std::to_string(line_i));
            } else if (!(tokens.empty())) {
                int r = proc_line(tokens);
                if (r == 1) {
                    break;
                }
            } else {
                std::cerr << "something went wrong... tokens is empty" << std::endl;
            }
            line_i++;
//            std::cerr << line_i << std::endl;
            if (line_i % 10000000 == 0) {
                std::cerr << line_i / 1000000 << "M" << std::endl;
            }
        }
        //file.close();
        // if pat file is empty, don't produce m files:
        if (empty_pat) {
            std::cerr << "Warning: empty pat file: " << output_prefix << std::endl;
            //return;
        }
        // dump markov counts:
        for (int i = 0; i < morder + 1; i++) {
            dump(marr[i], nr_states[i], output_prefix + ".m" + std::to_string(i));
        }

        // dump reads lengths file:
        //dump(counts, MAX_LEN, output_prefix + ".reads_lengths");
    }
    catch (std::exception &e) {
        std::cerr << "failed calculating rlen" << std::endl;
        std::cerr << e.what() << std::endl;
        return;
    }
}


