// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <chrono> // timing functions

#include <sharg/all.hpp> // parser

#include <vector>
#include <unordered_map>
#include <string>
#include <ranges>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/all.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

#include <seqan3/io/sequence_file/all.hpp>
#include <filesystem>

#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/search/all.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/unordered_map.hpp>
#include <fstream>
#include <iomanip>


constexpr size_t Sigma_unmin = 6;
constexpr size_t Sigma_min = 6;

std::string dir = "/storage/mi/eflecks/fMINdex/";


void write_stats_file(std::string filename, int const window_size, int const kmer_length, size_t const hits_unminimized, size_t const hits_minimized, auto const time_unminimized, auto const time_minimized, unsigned int const db_size_minimized, unsigned int const query_size_minimized) {
    std::ofstream outfile;
    
    std::ifstream f(dir + filename);
    if (f.good()) { // if file exists - append
        outfile.open(dir + filename, std::ios_base::app); // append instead of overwrite
        outfile << window_size << "," << kmer_length << ","
                << hits_unminimized << "," << hits_minimized << "," 
                << time_unminimized << "," << time_minimized << ","
                << db_size_minimized << "," << query_size_minimized << "\n";
    } else { // if file doesn't exist - create, make header
        outfile.open(dir + filename);
        outfile << "window_size" << "," << "kmer_length" << ","
                << "hits_unminimized" << "," << "hits_minimized" << "," 
                << "time_unminimized" << "," << "time_minimized" << "," 
                << "db_size_minimized" << "query_size_minimized" << "\n";
        outfile << window_size << "," << kmer_length << ","
                << hits_unminimized << "," << hits_minimized << "," 
                << time_unminimized << "," << time_minimized << ","
                << db_size_minimized << "," << query_size_minimized << "\n";

    }
    
}





// read fasta into dna5 vector
std::tuple<std::vector<seqan3::dna5_vector>, unsigned int> read_fasta(std::string const filename){

    seqan3::debug_stream << std::left << std::setw(30) << "loading" << filename;

    unsigned int sum = 0;
    
    seqan3::sequence_file_input file_input{dir + filename}; // contains multiple sequences
    std::vector<seqan3::dna5_vector> dna_vector{}; // vector of sequences
    for (auto & record : file_input) {
        dna_vector.push_back(record.sequence());
        sum += record.sequence().size();
    }


    seqan3::debug_stream << "-" << std::right << std::setw(10) << "LOADED\n";
    seqan3::debug_stream << "containing " << dna_vector.size() << " sequences, " << sum << " bases\n";
    
    if (dna_vector.size()<10) {
        seqan3::debug_stream << "sequences:\n";
        for (auto seq : dna_vector)
            seqan3::debug_stream << seq << "\n";
    }
    return std::tuple<std::vector<seqan3::dna5_vector>, unsigned int>(dna_vector, sum);
}




std::vector<uint8_t> to_base_5(auto input) {

    std::vector<uint8_t> output;
    
    while (input > 0) {
        output.push_back(input % 5 + 1);
        input = input / 5;
    }
    auto result = output | std::views::reverse | std::views::common;

    return std::vector(result.begin(), result.end());
}





// minimize text, fill hash map
std::vector<uint8_t> min_make_map(seqan3::dna5_vector const &input, std::unordered_map<size_t, uint8_t> &map, uint8_t const kmer_length, unsigned int const window_size) {
    
    auto minimized = input | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{kmer_length}}, seqan3::window_size{window_size});

    // Get hash values
    uint64_t seed = 0x8F3F73B5CF1C9ADE; // The default seed from minimiser_hash
    // Use XOR on all minimiser values
    auto hash_values = minimized | std::views::transform(
                           [seed](uint64_t i)
                           {
                               return to_base_5(i ^ seed);
                           });

    std::vector<uint8_t> result{};
    for(auto && item : hash_values){
        result.insert(result.end(), item.begin(), item.end());
    }



    /*std::vector<uint8_t> result{};

    for (uint64_t item : minimized) {
        //auto [find, success] = map.try_emplace(item, map.size());
        //result.push_back(find->second);

        auto find = map.find(item);
        if (find==map.end()) {
            if (map.size() >= Sigma_min) {
                throw std::runtime_error("too many different minimizers, adjust k, w or Sigma");
            }
            map.insert({item, map.size()+1});
            result.push_back(map.size());
        } else { 
            result.push_back(find->second);
        }
    }*/
    
    return result; 
}




// minimize query, use given hashmap
std::vector<uint8_t> min_use_map(seqan3::dna5_vector const &query, std::unordered_map<size_t, uint8_t> const &map, uint8_t const kmer_length, unsigned int const window_size) {
        
    auto minimized = query | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{kmer_length}}, seqan3::window_size{window_size});
    
    /*std::vector<uint8_t> result{};

    for (uint64_t item : minimized) {
        auto find = map.find(item);
        if (find==map.end()) {
            return {}; // no instances of query exist in text
        } else { 
            result.push_back(find->second);
        }
    }*/

    //
    // Get hash values
    uint64_t seed = 0x8F3F73B5CF1C9ADE; // The default seed from minimiser_hash
    // Use XOR on all minimiser values
    auto hash_values = minimized | std::views::transform(
                           [seed](uint64_t i)
                           {
                               return to_base_5(i ^ seed);
                           });

    std::vector<uint8_t> result{};
    for(auto && item : hash_values){
        result.insert(result.end(), item.begin(), item.end());
    }

    
    return result; 
}


// minimize entire databank (fasta file) with multiple sequences
std::tuple<std::vector<std::vector<uint8_t>>, unsigned int> minimise_db(std::string const filename, int const kmer_length, int const window_size){
    
    auto [input, db_size_unminimized] = read_fasta(filename);
    std::unordered_map<size_t, uint8_t> map{};
    std::vector<std::vector<uint8_t>> result{};
    result.reserve(input.size());


    seqan3::debug_stream << std::left << std::setw(30) << "minimizing database";


    unsigned int sum = 0;
    std::vector<uint8_t> minseq{};

    for (seqan3::dna5_vector sequence : input){
        //minq.clear();
        minseq = min_make_map(sequence, map, kmer_length, window_size);
        sum += minseq.size();
        result.push_back(std::move(minseq));
        
    }
    

    // save map to file
    auto ofs     = std::ofstream(dir + "map", std::ios::binary);
    auto archive = cereal::BinaryOutputArchive{ofs};
    archive(map);
    

    seqan3::debug_stream << "-" << std::right << std::setw(10) << "DONE\n";
    seqan3::debug_stream << "minimized: " << result.size() << " sequences, " << sum << " bases\n";
    if (result.size()<10) {
        seqan3::debug_stream << "sequences:\n";
        for (auto seq : result)
            seqan3::debug_stream << seq << "\n";
    }
    return std::tuple<std::vector<std::vector<uint8_t>>, unsigned int>(result, sum);
}





// minimize all queries (from fasta file)
std::tuple<std::vector<std::vector<uint8_t>>, unsigned int> minimise_queries(std::string const filename, int const kmer_length, int const window_size){
    
    // load map from file
    auto ifs     = std::ifstream(dir + "map", std::ios::binary);
    auto archive = cereal::BinaryInputArchive{ifs};
    auto map = std::unordered_map<size_t, uint8_t>{};
    archive(map);

    auto [queries, query_size_unminimized] = read_fasta(filename);
    std::vector<std::vector<uint8_t>> result{};
    result.reserve(queries.size());

    seqan3::debug_stream << std::left << std::setw(30) << "minimizing queries";

    unsigned int sum = 0;
    std::vector<uint8_t> minq{};
    for (seqan3::dna5_vector query : queries){
        //minq.clear();
        minq = min_use_map(query, map, kmer_length, window_size);
        sum += minq.size();
        result.push_back(std::move(minq));
        
    }
    



    seqan3::debug_stream << "-" << std::right << std::setw(10) << "DONE\n";
    seqan3::debug_stream << "minimized: " << sum << " bases\n";

    if (result.size()<10) {
        seqan3::debug_stream << "sequences:\n";
        for (auto seq : result)
            seqan3::debug_stream << seq << "\n";
    }
    return std::tuple<std::vector<std::vector<uint8_t>>, unsigned int>(result, sum);
}





// convert vector of dna5 sequences to vector of vector of uint8_t 
std::vector<std::vector<uint8_t>> to_uint8_t(std::vector<seqan3::dna5_vector> const& dna_vector) {

    std::vector<std::vector<uint8_t>> result{};

    for (seqan3::dna5_vector sequence : dna_vector) {

        std::vector<uint8_t> converted_sequence{};
        converted_sequence.reserve(sequence.size());

        for (seqan3::dna5 base : sequence) {
            converted_sequence.push_back(base.to_rank()+1);
        }
        result.push_back(std::move(converted_sequence));
    }

    if (result.size()<10) {
        seqan3::debug_stream << "sequences:\n";
        for (auto seq : result)
            seqan3::debug_stream << seq << "\n";
    }
    return result;
}






template <size_t TSigma>
void get_fmindex(std::vector<std::vector<uint8_t>> const& sequences, std::string const filename){
    using Table = fmindex_collection::occtable::Interleaved_16<6>;

    unsigned int sum{0};
    for (auto seq : sequences) sum += seq.size();

    seqan3::debug_stream << sequences.size() << " sequences   " << sum << " bases\n";
    
    if (sequences.size()<10) {
        seqan3::debug_stream << "sequences:\n";
        for (auto seq : sequences)
            seqan3::debug_stream << seq << "\n";
    }

    // get fmindex
    seqan3::debug_stream << std::left << std::setw(30) << "generating fm-index";
    auto index = fmindex_collection::FMIndex<Table>{sequences, /*samplingRate*/16, /*threadNbr*/20};
    seqan3::debug_stream << "-" << std::right << std::setw(10) << "DONE\n";
    
    seqan3::debug_stream << "FMindex size: " << index.size() << "\n";
    
    // save to file
    seqan3::debug_stream << std::left << std::setw(30) << "saving index to file";
    auto ofs     = std::ofstream(dir + filename, std::ios::binary);
    auto archive = cereal::BinaryOutputArchive{ofs};
    archive(index);
    seqan3::debug_stream << "-" << std::right << std::setw(10) << "SAVED\n";

}





template<size_t Sigma>
std::tuple<std::vector<int>, std::chrono::duration<long int, std::ratio<1,1000>>> fmindex_search( std::vector<std::vector<uint8_t>> const& queries, std::string const filename){
    using Table = fmindex_collection::occtable::Interleaved_16<6>;

    // load index from file
    seqan3::debug_stream << std::left << std::setw(30) << "loading index from file";
    
    auto ifs     = std::ifstream(dir + filename, std::ios::binary);
    auto archive = cereal::BinaryInputArchive{ifs};
    auto index = fmindex_collection::FMIndex<Table>{};
    archive(index);
    
    seqan3::debug_stream << "-" << std::right << std::setw(10) << "LOADED\n";
    seqan3::debug_stream << "FMindex size: " << index.size() << "\n";
    

    // search
    std::vector<std::chrono::milliseconds> search_times{};
    std::vector<int> counts{};
    for(int i=0; i<5; i++) {
        auto before_minimized_search = std::chrono::high_resolution_clock::now();

        counts.clear();
        for (size_t queryId=0; queryId<queries.size(); queryId++){
            if( queries[queryId].size() == 0) {
                //seqan3::debug_stream << "query " << queryId << " - # hits: 0\n";
                counts.push_back(0);
                continue;
            }

            auto cursor = fmindex_collection::search_no_errors::search(index, queries[queryId]);
            //seqan3::debug_stream << "found something " << queryId << " " << cursor.count() << "\n";
            /*for (auto i : cursor) {
                auto [chr, pos] = index.locate(i);
                //seqan3::debug_stream << "chr/pos: "<< chr << " " << pos << "\n";
            }*/ 
            // note: bitvector w/ select capabilities needed to get back to original positions from minimizer,
            // select call would count towards time eval
            

            //seqan3::debug_stream << "query " << queryId << " - # hits: " << cursor.len << "\n";
            counts.push_back(cursor.len);

        }
        
        auto after_minimized_search = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(after_minimized_search - before_minimized_search);
        search_times.push_back(time);
    }
    
        

    std::chrono::milliseconds min_search_time = search_times[0];
    for (std::chrono::milliseconds time : search_times) 
        if (time<min_search_time) min_search_time = time;
    

    return std::tuple(counts, min_search_time);
}







int main(int argc, char ** argv)
{   
    // parser
    sharg::parser parser{"minimize", argc, argv};
    
    // flag for minimizing
    bool minimize{false};
    parser.add_flag(minimize, sharg::config{.short_id = 'm', .long_id = "minimize", .description = "minimize before constructing fmindex"});
    
    // flag to run full test
    bool test{false};
    parser.add_flag(test, sharg::config{.short_id = 't', .long_id = "test", .description = "run full test: generate fm-index on minimizers and unminimized database and compare"});
    
    // flag to generate fm-index instead of reading it from file
    bool generate{false};
    parser.add_flag(minimize, sharg::config{.short_id = 'g', .long_id = "generate", .description = "generate fm-index new instead of reading from file, should it exist. If the index file does not exist the index will be generated in any case."});

    // file inputs
    std::string db_filename{};
    parser.add_positional_option(db_filename, sharg::config{.description = "database file (FASTA)"});
    //parser.add_option(db_filename, sharg::config{.short_id = 'd', .long_id = "database", .description = "filename of database, don't specify to not construct fm-Index"});
    std::string search_filename{};
    parser.add_positional_option(search_filename, sharg::config{.description = "query file (FASTA)"});
    //parser.add_option(search_filename, sharg::config{.short_id = 'q', .long_id = "queries", .description = "filename for queries, leave empty/don't specify for no search"});

    // window size and kmer 
    unsigned int window_size{};
    parser.add_positional_option(window_size, sharg::config{.description = "window size the minimizer uses"});
    uint8_t kmer_length{};
    parser.add_positional_option(kmer_length, sharg::config{.description = "kmer length for minimizer"});



    try
    {
        parser.parse(); // Trigger command line parsing.
    }
    catch (sharg::parser_error const & ext) // Catch user errors.
    {
        std::cerr << "Parsing error. " << ext.what() << '\n'; // Give error message.
        return -1;
    }








    if (test==1) {
        if(db_filename.length()==0 || search_filename.length()==0) throw std::runtime_error("missing database or query file in call");


        // run without minimizer
        seqan3::debug_stream << "\n########     without minimizers     ########\n";
        // get fm index
        std::ifstream f(dir + db_filename);
        if (!f.good() || generate) { // if index isn't saved or generate flag is set to true
            auto [sequences, db_size_unminimized] = read_fasta(db_filename);
            auto seqs = to_uint8_t(sequences);
            get_fmindex<Sigma_unmin>(seqs, db_filename + "_index");
        } else {
            unsigned int db_size_unminimized = 0;
        }
        

        // search queries
        auto [queries, query_size_unminimized] = read_fasta(search_filename);
        auto qs = to_uint8_t(queries);

        //auto before_unminimized_search = std::chrono::high_resolution_clock::now();
        auto [counts, time] = fmindex_search<Sigma_unmin>(qs, db_filename + "_index");
        //auto after_unminimized_search = std::chrono::high_resolution_clock::now();




        // run with minimizer
        seqan3::debug_stream << "\n########      with  minimizers      ########\n";
        // minimize db, create fm index
        auto [minseqs, db_size_minimized] = minimise_db(db_filename, kmer_length, window_size);
        get_fmindex<Sigma_min>(minseqs, db_filename + "_mindex");
        
        // minimize & search queries
        auto [minqs, query_size_minimized] = minimise_queries(search_filename, kmer_length, window_size);
        //seqan3::debug_stream << minqs << "\n\n";
        
        //auto before_minimized_search = std::chrono::high_resolution_clock::now();
        auto [mincounts, mintime] = fmindex_search<Sigma_min>(minqs, db_filename + "_mindex");
        //auto after_minimized_search = std::chrono::high_resolution_clock::now();
        


        // compare search time
        /*auto ms_unminimized_search = std::chrono::duration_cast<std::chrono::milliseconds>(after_unminimized_search - before_unminimized_search);
        seqan3::debug_stream << "time for unminimized search: " << ms_unminimized_search.count() << "ms\n";
        auto ms_minimized_search = std::chrono::duration_cast<std::chrono::milliseconds>(after_minimized_search - before_minimized_search);
        seqan3::debug_stream << "time for minimized search:   " << ms_minimized_search.count() << "ms\n";
        */

        seqan3::debug_stream << "time for unminimized search: " << time << "\n";
        seqan3::debug_stream << "time for minimized search:   " << mintime << "\n";
        


        // compare # of hits (FP)
        size_t count_unmin = 0;
        size_t count_min = 0;
        for (long unsigned int i=0; i<counts.size(); i++) {
            count_unmin += counts[i];
            count_min += mincounts[i];
            //seqan3::debug_stream << "query " << i << " - " << counts[i] << " vs " << mincounts[i] << "\n";
        }
        seqan3::debug_stream << "# hits unminimized: " << count_unmin << "\n";
        seqan3::debug_stream << "# hits minimized:   " << count_min << "\n";
        seqan3::debug_stream << "# FP:               " << count_min - count_unmin << "\n";


        write_stats_file("stats.csv",
                        window_size, kmer_length, 
                        count_unmin, count_min, 
                        //ms_unminimized_search.count(), ms_minimized_search.count());
                        time, mintime,
                        db_size_minimized, query_size_minimized);




    } else if (minimize == 0) {

        // run without minimizer
        seqan3::debug_stream << "\n########     without minimizers     ########\n";
        // get fm index
        std::ifstream f(dir + db_filename);
        if (!f.good() || generate) { // if index isn't saved or generate flag is set to true
            auto [sequences, db_size_unminimized] = read_fasta(db_filename);
            auto seqs = to_uint8_t(sequences);
            get_fmindex<Sigma_unmin>(seqs, db_filename + "_index");
        }
        

        // search queries
        auto [queries, query_size_unminimized] = read_fasta(search_filename);
        auto qs = to_uint8_t(queries);

        auto [counts, time] = fmindex_search<Sigma_unmin>(qs, db_filename + "_index");
        seqan3::debug_stream << "time: " << time << "ms\n";

        int count_unmin{0};
        for (long unsigned int i=0; i<counts.size(); i++) {
            count_unmin += counts[i];
            }
        seqan3::debug_stream << "hits: " << count_unmin << "\n";



    } else { // with minimizing
        seqan3::debug_stream << "running WITH minimizing\n";

        // get fm index
        auto [minseqs, db_size_minimized]  = minimise_db(db_filename, kmer_length, window_size);
        get_fmindex<Sigma_min>(minseqs, db_filename + "_mindex");


    
        // search
        auto [minqs, query_size_minimized]  = minimise_queries(search_filename, kmer_length, window_size);
        seqan3::debug_stream << "searching queries\n";
            
        //auto before_minimized_search = std::chrono::high_resolution_clock::now();
        auto [mincounts, mintime] = fmindex_search<Sigma_min>(minqs, db_filename + "_mindex");
        //auto after_minimized_search = std::chrono::high_resolution_clock::now();
        //auto ms_minimized_search = std::chrono::duration_cast<std::chrono::milliseconds>(after_minimized_search - before_minimized_search);
        seqan3::debug_stream << "time for minimized search:   " << mintime << "ms\n";


    }
    
    
    
    
    return 0;
}
