// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <chrono> // timing functions

#include <sharg/all.hpp> // parser

#include <vector>
#include <unordered_map>
#include <string>

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

int const window_size = 20;
int const kmer_length = 4;
constexpr size_t Sigma_unmin = 6;
constexpr size_t Sigma_min = 150;




// read fasta into dna5 vector
std::vector<seqan3::dna5_vector> read_fasta(std::string const filename){

    seqan3::debug_stream << std::left << std::setw(30) << "loading fasta file ";
    
    seqan3::sequence_file_input file_input{"/storage/mi/eflecks/fMINdex/" + filename}; // contains multiple sequences
    std::vector<seqan3::dna5_vector> dna_vector{}; // vector of sequences
    for (auto & record : file_input) dna_vector.push_back(record.sequence());


    seqan3::debug_stream << "-" << std::right << std::setw(10) << "LOADED\n";
    return dna_vector;
}





// minimize text, fill hash map
std::vector<uint8_t> min_make_map(seqan3::dna5_vector const &input, std::unordered_map<size_t, uint8_t> &map) {
    
    auto minimized = input | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{kmer_length}}, seqan3::window_size{window_size});

    std::vector<uint8_t> result{};

    for (uint64_t item : minimized) {
        //auto [find, success] = map.try_emplace(item, map.size());
        //result.push_back(find->second);

        auto find = map.find(item);
        if (find==map.end()) {
            if (map.size() >= 255) {
                throw std::runtime_error("too many different minimizers");
            }
            map.insert({item, map.size()+1});
            result.push_back(map.size());
        } else { 
            result.push_back(find->second);
        }
    }

    return result; 
}




// minimize query, use given hashmap
std::vector<uint8_t> min_use_map(seqan3::dna5_vector const &query, std::unordered_map<size_t, uint8_t> const &map) {
        
    auto minimized = query | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{kmer_length}}, seqan3::window_size{window_size});
    
    std::vector<uint8_t> result{};

    for (uint64_t item : minimized) {
        auto find = map.find(item);
        if (find==map.end()) {
            return {}; // no instances of query exist in text
        } else { 
            result.push_back(find->second);
        }
    }
    
    return result; 
}


// minimize entire databank (fasta file) with multiple sequences
std::vector<std::vector<uint8_t>> minimise_db(std::string const filename){
    
    std::vector<seqan3::dna5_vector> input = read_fasta(filename);
    std::unordered_map<size_t, uint8_t> map{};
    std::vector<std::vector<uint8_t>> result{};

    seqan3::debug_stream << std::left << std::setw(30) << "minimizing database";
    for (seqan3::dna5_vector sequence : input){
        std::vector<uint8_t> minseq = min_make_map(sequence, map);
        result.push_back(minseq);
    }

    // save map to file
    auto ofs     = std::ofstream("/storage/mi/eflecks/fMINdex/map", std::ios::binary);
    auto archive = cereal::BinaryOutputArchive{ofs};
    archive(map);
    

    seqan3::debug_stream << "-" << std::right << std::setw(10) << "DONE\n";
    seqan3::debug_stream << "map size: " << map.size() << "\n";
    return result;
}





// minimize all queries (from fasta file)
std::vector<std::vector<uint8_t>> minimise_queries(std::string const filename){
    
    // load map from file
    auto ifs     = std::ifstream("/storage/mi/eflecks/fMINdex/map", std::ios::binary);
    auto archive = cereal::BinaryInputArchive{ifs};
    auto map = std::unordered_map<size_t, uint8_t>{};
    archive(map);
    seqan3::debug_stream << "map   " << map << "\n";

    std::vector<seqan3::dna5_vector> queries = read_fasta(filename);
    std::vector<std::vector<uint8_t>> result{};

    seqan3::debug_stream << std::left << std::setw(30) << "minimizing queries";
    for (seqan3::dna5_vector query : queries){
        std::vector<uint8_t> minq = min_use_map(query, map);
        result.push_back(minq);
    }


    seqan3::debug_stream << "-" << std::right << std::setw(10) << "DONE\n";
    return result;
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

    return result;
}

template <size_t TSigma>
void get_fmindex(std::vector<std::vector<uint8_t>> const& sequences, std::string const filename){
    using Table = fmindex_collection::occtable::Interleaved_16<TSigma>;
    
    // get fmindex
    seqan3::debug_stream << std::left << std::setw(30) << "generating fm-index";
    auto index = fmindex_collection::FMIndex<Table>{sequences, /*samplingRate*/16, /*threadNbr*/16};
    seqan3::debug_stream << "-" << std::right << std::setw(10) << "DONE\n";
    

    // save to file
    seqan3::debug_stream << std::left << std::setw(30) << "saving index to file";
    auto ofs     = std::ofstream("/storage/mi/eflecks/fMINdex/" + filename, std::ios::binary);
    auto archive = cereal::BinaryOutputArchive{ofs};
    archive(index);


    seqan3::debug_stream << "-" << std::right << std::setw(10) << "SAVED\n";

}





template<size_t Sigma>
std::vector<int> fmindex_search( std::vector<std::vector<uint8_t>> const& queries, std::string const filename){
    
    using Table = fmindex_collection::occtable::Interleaved_16<Sigma>;

    // load index from file
    seqan3::debug_stream << std::left << std::setw(30) << "loading index from file";
    
    auto ifs     = std::ifstream("/storage/mi/eflecks/fMINdex/" + filename, std::ios::binary);
    auto archive = cereal::BinaryInputArchive{ifs};
    auto index = fmindex_collection::FMIndex<Table>{};
    archive(index);
    
    seqan3::debug_stream << "-" << std::right << std::setw(10) << "LOADED\n";


    // search
    std::vector<int> counts{};
    for (size_t queryId=0; queryId<queries.size(); queryId++){
        if( queries[queryId].size() == 0) {
            //seqan3::debug_stream << "query " << queryId << " - # hits: 0\n";
            counts.push_back(0);
            continue;
        }

        auto cursor = fmindex_collection::search_no_errors::search(index, queries[queryId]);
        seqan3::debug_stream << "found something " << queryId << " " << cursor.count() << "\n";
        for (auto i : cursor) {
            auto [chr, pos] = index.locate(i);
            //seqan3::debug_stream << "chr/pos: "<< chr << " " << pos << "\n";
        } // note: bitvector w/ select capabilities needed to get back to original positions from minimizer,
        // select call would count towards time eval
        

        //seqan3::debug_stream << "query " << queryId << " - # hits: " << cursor.len << "\n";
        counts.push_back(cursor.len);

    }
    return counts;
}


int main(int argc, char ** argv)
{   
    // parser
    sharg::parser parser{"minimize", argc, argv};
    // flags for minimizing/not minimizing
    bool minimize{false};
    parser.add_flag(minimize, sharg::config{.short_id = 'm', .long_id = "min", .description = "minimize before constructing fmindex"});
    
    bool test{false};
    parser.add_flag(test, sharg::config{.short_id = 't', .long_id = "test", .description = "run comparison test"});

    std::string db_filename{};
    parser.add_option(db_filename, sharg::config{.short_id = 'd', .long_id = "database", .description = "filename of database, don't specify to not construct fm-Index"});
    std::string search_filename{};
    parser.add_option(search_filename, sharg::config{.short_id = 'q', .long_id = "queries", .description = "filename for queries, leave empty/don't specify for no search"});

    //parser.add_positional_option();

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
        auto sequences = read_fasta(db_filename);
        auto seqs = to_uint8_t(sequences);
        get_fmindex<Sigma_unmin>(seqs, db_filename + "_index");

        // search queries
        auto queries = read_fasta(search_filename);
        auto qs = to_uint8_t(queries);

        auto before_unminimized_search = std::chrono::high_resolution_clock::now();
        std::vector<int> counts = fmindex_search<Sigma_unmin>(qs, db_filename + "_index");
        auto after_unminimized_search = std::chrono::high_resolution_clock::now();




        // run with minimizer
        seqan3::debug_stream << "\n########      with  minimizers      ########\n";
        // minimize db, create fm index
        auto minseqs = minimise_db(db_filename);
        get_fmindex<Sigma_min>(minseqs, db_filename + "_mindex");
        
        // minimize & search queries
        auto minqs = minimise_queries(search_filename);
        seqan3::debug_stream << minqs << "\n\n";
        
        auto before_minimized_search = std::chrono::high_resolution_clock::now();
        std::vector<int> mincounts = fmindex_search<Sigma_min>(minqs, db_filename + "_mindex");
        auto after_minimized_search = std::chrono::high_resolution_clock::now();
        




        // compare search time
        auto ms_unminimized_search = std::chrono::duration_cast<std::chrono::microseconds>(after_unminimized_search - before_unminimized_search);
        seqan3::debug_stream << "time for unminimized search: " << ms_unminimized_search.count() << "microseconds\n";
        auto ms_minimized_search = std::chrono::duration_cast<std::chrono::microseconds>(after_minimized_search - before_minimized_search);
        seqan3::debug_stream << "time for minimized search:   " << ms_minimized_search.count() << "microseconds\n";


        // compare # of hits (FP)
        size_t count_unmin = 0;
        size_t count_min = 0;
        for (int i=0; i<counts.size(); i++) {
            count_unmin += counts[i];
            count_min += mincounts[i];
            seqan3::debug_stream << "query " << i << " - " << counts[i] << " vs " << mincounts[i] << "\n";
        }
        seqan3::debug_stream << "# hits unminimized: " << count_unmin << "\n";
        seqan3::debug_stream << "# hits minimized:   " << count_min << "\n";
        seqan3::debug_stream << "# FP:               " << count_min - count_unmin << "\n";







    } else if (minimize == 0) {

        seqan3::debug_stream << "running WITHOUT minimizing\n";

        // get fm index
        if (db_filename.length()>0) {
            auto sequences = read_fasta(db_filename);
            auto seqs = to_uint8_t(sequences);

            long unsigned int sum = 0;
            for (auto seq : seqs) sum += seq.size();
            seqan3::debug_stream << "db: # seqs =  " << seqs.size() << "  # bases = " << sum << "\n";
        
            get_fmindex<Sigma_unmin>(seqs, "test_index");
        }
        

        // search
        if (search_filename.length()>0) {
            auto queries = read_fasta(search_filename);
            auto qs = to_uint8_t(queries);

            seqan3::debug_stream << "searching queries\n";
            fmindex_search<Sigma_unmin>(qs, "test_index");
        }
        



    } else { // with minimizing
        seqan3::debug_stream << "running WITH minimizing\n";

        // get fm index
        if (db_filename.length()>0) {
            auto minseqs = minimise_db(db_filename);
            
            

            int sum = 0;
            for (auto seq : minseqs) sum += seq.size();
            seqan3::debug_stream << "db: # seqs =  " << minseqs.size() << "  # bases = " << sum << "\n";

            get_fmindex<Sigma_min>(minseqs, "test_mindex");

            
        }
    
        // search
        if (search_filename.length()>0) {
            auto minqs = minimise_queries(search_filename);
            seqan3::debug_stream << "searching queries\n";
            fmindex_search<Sigma_min>(minqs, "test_mindex");
        }
    }
    
    
    
    
    return 0;
}
