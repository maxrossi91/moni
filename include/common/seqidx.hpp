/* seqidx - an index fo the sequence names in a fasta file
    Copyright (C) 2020 Massimiliano Rossi
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file seqidx.cpp
   \brief seqidx.cpp an index fo the sequence names in a fasta file.
   \author Massimiliano Rossi
   \date 07/08/2021
*/
#ifndef _SEQIDX_HH
#define _SEQIDX_HH

#include <common.hpp>

#include <sdsl/sd_vector.hpp>

#include <kseq.h>
#include <zlib.h>

// KSEQ_DECLARE(gzFile);

class seqidx
{
public:
    seqidx()
    {
        u = 0;
    }
    /**
     * @brief Construct a new seqidx object
     * 
     * @param filename filepath of the fasta/q file
     */
    seqidx(std::string filename)
    {
        gzFile fp(gzopen(filename.c_str(), "r"));
        if (fp == nullptr)
            error("gzopen() file " + std::string(filename) + " failed");

        kseq_t *seq = kseq_init(fp);

        std::vector<size_t> onset(1,0);
        u = 0;

        while (kseq_read(seq) >= 0)
        {
            u += seq->seq.l;
            names.push_back(std::string(seq->name.s));
            onset.push_back(u);
        }

        kseq_destroy(seq);
        gzclose(fp);

        sdsl::sd_vector_builder builder(u, onset.size());
        for (auto idx : onset)
            builder.set(idx);

        starts = sdsl::sd_vector<>(builder);
        rank1 = sdsl::sd_vector<>::rank_1_type(&starts);
        select1 = sdsl::sd_vector<>::select_1_type(&starts);
    }

    /**
     * @brief Construct a new seqidx object from onset, list of sequence names and total length
     * 
     * @param onset the popsitions
     * @param names_ 
     * @param l 
     */
    seqidx(const std::vector<size_t>& onset, const std::vector<std::string>& names_, const size_t l)
    {
        assert(onset.size() == names_.size());
        assert(onset[0] == 0);
        assert(onset.back() < l);
        assert(std::is_sorted(onset.begin(), onset.end()));

        u = l;
        names = std::vector<std::string>(names_);


        sdsl::sd_vector_builder builder(u, onset.size());
        for (auto idx : onset)
            builder.set(idx);
        
        builder.set(u);

        starts = sdsl::sd_vector<>(builder);
        rank1 = sdsl::sd_vector<>::rank_1_type(&starts);
        select1 = sdsl::sd_vector<>::select_1_type(&starts);
    }

    

    /**
     * @brief Return the length of the i-th sequence
     * 
     * @param i 
     * @return size_t 
     */
    inline size_t length(const size_t i)
    {
        assert(i < names.size());
        return select1(i+1) - select1(i);
    }

    /**
     * @brief return the name of the sequence pos belongs.
     * 
     * @param pos the position in the set of sequences.
     * @return std::string the name of the sequence pos belongs.
     */
    inline std::string operator[](const size_t pos)
    {
        return names[rank1(pos)-1];
    }

    /**
     * @brief Check if the substring [pos.pos+len-1] does not span two sequences.
     * 
     * @param pos the position of the substring.
     * @param len the length of the substring.
     * @return true if the substring does not span two sequences.
     * @return false if the substring spans two sequences.
     */
    inline bool valid(size_t pos, size_t len)
    {
        return (pos + len <= select1(rank1(pos)+1));
    }

    /**
     * @brief return the SAM header description of the reference file
     * 
     * @return std::string 
     */
    std::string to_sam()
    {
        std::string res = "";
        for (size_t i = 0; i < names.size(); ++i)
            res += "@SQ\tSN:" + names[i] + "\tLN:" + std::to_string(length(i)) + "\n";
        return res;    
    }

    size_t serialize(std::ostream &out)
    {

        size_t w_bytes = 0;

        out.write((char *)&u, sizeof(u));

        w_bytes += sizeof(u);

        if (u == 0)
            return w_bytes;

        w_bytes += starts.serialize(out);
        w_bytes += sdsl::serialize(names.size(), out);
        for(size_t i = 0; i < names.size(); ++i)
        {
            w_bytes += sdsl::serialize(names[i].size(), out);
            w_bytes = my_serialize_array<char, std::string::size_type>(names[i].data(), names[i].size(), out);
        }
        return w_bytes;
    }

    void load(std::istream &in)
    {

        in.read((char *)&u, sizeof(u));

        if (u == 0)
            return;

        starts.load(in);
        rank1 = sdsl::sd_vector<>::rank_1_type(&starts);
        select1 = sdsl::sd_vector<>::select_1_type(&starts);

        std::vector<std::string>::size_type names_size;
        sdsl::load(names_size, in);
        names.resize(names_size);
        for (size_t i = 0; i < names.size(); ++i)
        {
            std::string::size_type string_size;
            sdsl::load(string_size, in);
            names[i].resize(string_size);
            my_load_array<char, std::string::size_type>(&names[i][0], names[i].size(), in);
        }
    }

    std::string get_file_extension() const
    {
        return  ".idx";
    }

protected:
    size_t u;
    
    sdsl::sd_vector<> starts;
    sdsl::sd_vector<>::rank_1_type rank1;
    sdsl::sd_vector<>::select_1_type select1;
    
    std::vector<std::string> names;

};

#endif /* end of include guard: _SEQIDX_HH */
