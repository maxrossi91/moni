/* extender_klib - Extend the MEMs of the reads to the reference using the klib library for SW
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
   \file extender_klib.cpp
   \brief extender_klib.cpp Extend the MEMs of the reads to the reference using the klib library for SW
   \author Massimiliano Rossi
   \date 13/07/2020
*/

#ifndef _EXTENDER_KLIB_HH
#define _EXTENDER_KLIB_HH

#include <common.hpp>

#include <sdsl/io.hpp>

#include <ms_pointers.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <PlainSlp.hpp>
#include <FixedBitLenCode.hpp>

#include <ksw.h>
#include <ssw.h>

#include <libgen.h>
#include <seqidx.hpp>

////////////////////////////////////////////////////////////////////////////////
/// SLP definitions
////////////////////////////////////////////////////////////////////////////////

using SelSd = SelectSdvec<>;
using DagcSd = DirectAccessibleGammaCode<SelSd>;
using Fblc = FixedBitLenCode<>;

using shaped_slp_t = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;
using plain_slp_t = PlainSlp<uint32_t, Fblc, Fblc>;

template <typename slp_t>
std::string get_slp_file_extension()
{
    return std::string(".slp");
}

template <>
std::string get_slp_file_extension<shaped_slp_t>()
{
    return std::string(".slp");
}

template <>
std::string get_slp_file_extension<plain_slp_t>()
{
    return std::string(".plain.slp");
}
////////////////////////////////////////////////////////////////////////////////

template <typename slp_t>
class extender
{
public:
    extender(std::string filename,
            size_t min_len_ = 50,
            bool forward_only_ = true) : min_len(min_len_),
                                         forward_only(forward_only_)
    {
        verbose("Loading the matching statistics index");
        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string filename_ms = filename + ms.get_file_extension();

        ifstream fs_ms(filename_ms);
        ms.load(fs_ms);
        fs_ms.close();

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Matching statistics index construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

        verbose("Loading random access");
        t_insert_start = std::chrono::high_resolution_clock::now();

        std::string filename_slp = filename + get_slp_file_extension<slp_t>();

        ifstream fs(filename_slp);
        ra.load(fs);
        fs.close();

        n = ra.getLen();

        t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Matching statistics index loading complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

        std::string filename_idx = filename + idx.get_file_extension();
        verbose("Loading fasta index file: " + filename_idx);
        t_insert_start = std::chrono::high_resolution_clock::now();


        ifstream fs_idx(filename_idx);
        idx.load(fs_idx);
        fs_idx.close();

        t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Fasta index loading complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());


        verbose("Initialize the local aligner");
        t_insert_start = std::chrono::high_resolution_clock::now();

        if (minsc > 0xffff)
            minsc = 0xffff;
        xtra |= KSW_XSUBO | minsc;
        // initialize scoring matrix
        for (i = k = 0; i < 4; ++i)
        {
            for (j = 0; j < 4; ++j)
                mat[k++] = i == j ? sa : -sb;
            mat[k++] = 0; // ambiguous base
        }
        for (j = 0; j < 5; ++j)
            mat[k++] = 0;

        t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Local aligner initialization complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

        verbose("Minimum MEM length: ", min_len);
    }

    bool extend(kseq_t *read, FILE *out, uint8_t strand)
    {
        size_t mem_pos = 0;
        size_t mem_len = 0;
        size_t mem_idx = 0;

        bool extended = false;

        auto pointers = ms.query(read->seq.s, read->seq.l);
        std::vector<size_t> lengths(pointers.size());
        size_t l = 0;
        for (size_t i = 0; i < pointers.size(); ++i)
        {
            size_t pos = pointers[i];
            while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
                ++l;

            lengths[i] = l;
            l = (l == 0 ? 0 : (l - 1));

            // Update MEM
            if (lengths[i] > mem_len)
            {
                mem_len = lengths[i];
                mem_pos = pointers[i];
                mem_idx = i;
            }
        }

        // Align the read
        if (mem_len >= min_len)
        {
            char *str = (char *)malloc(400);

            int32_t maskLen = read->seq.l / 2;
            maskLen = maskLen < 15 ? 15 : maskLen;

            // Extract the context from the reference
            size_t left_occ = (mem_pos > 100 ? mem_pos - 100 : 0);
            size_t len = mem_len + 100 + (mem_pos > 100 ? 100 : 100 - mem_pos);
            ra.expandSubstr(left_occ, len, str);

            size_t min_score = 20 + 8 * log(read->seq.l);

            uint8_t *seq = (uint8_t *)malloc(read->seq.l);
            // Convert A,C,G,T,N into 0,1,2,3,4
            for (i = 0; i < (int)read->seq.l; ++i)
                seq[i] = seq_nt4_table[(int)read->seq.s[i]];
            // for (i = 0; i < (int)read->seq.l; ++i)
            //   read->seq.s[i] = seq_nt4_table[(int)read->seq.s[i]];

            for (i = 0; i < (int)len; ++i)
                str[i] = seq_nt4_table[(int)str[i]];

            int score;

            kswq_t *q = 0;
            kswr_t r;

            r = ksw_align(read->seq.l, (uint8_t *)seq, len, (uint8_t *)str, 5, mat, gapo, gape, xtra, &q);
            // score = ksw_global(read->seq.l, (uint8_t *)read->seq.s, len, (uint8_t *)str, 5, mat, gapo, gape, w, &n_cigar, &cigar);

            int n_cigar;
            uint32_t *cigar;

            size_t new_seq_len = r.qe - r.qb;
            size_t new_ref_len = r.te - r.tb;
            uint8_t *new_seq = (uint8_t *)(seq + r.qb);
            // uint8_t *new_seq = (uint8_t *)(read->seq.s + r.qb);
            uint8_t *new_ref = (uint8_t *)(str + r.tb);

            score = ksw_global(new_seq_len, (uint8_t *)new_seq, new_ref_len, new_ref, 5, mat, gapo, gape, w, &n_cigar, &cigar);

            std::string cig;

            // for(size_t i = 0; i < n_cigar; ++i)
            // {
            //   // for (i = 0; i < ez->n_cigar; ++i)
            //   //   printf("%d%c", ez->cigar[i] >> 4, "MID"[ez->cigar[i] & 0xf]);
            //   cig += std::to_string(cigar[i] >> 4) + "MID"[cigar[i] & 0xf];
            // }

            size_t mismatch = mark_mismatch(r.tb, r.qb, r.qe, (int8_t *)str, (int8_t *)seq, read->seq.l, &cigar, &n_cigar);
            for (c = 0; c < (n_cigar); ++c)
            {
                char letter = cigar_int_to_op(cigar[c]);
                uint32_t length = cigar_int_to_len(cigar[c]);
                // fprintf(out, "%lu%c", (unsigned long)length, letter);
                cig += std::to_string((unsigned long)length) + letter;
            }

            // if(r.score > 0)
            //   printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n", "human", r.tb, r.te + 1, read->name.s, r.qb, r.qe + 1, r.score, r.score2, r.te2);
            //   // std::cout << "\rCurrent score... "<< r.score << std::flush;

            // // Declares a default Aligner
            // StripedSmithWaterman::Aligner aligner;
            // // Declares a default filter
            // StripedSmithWaterman::Filter filter;
            // // StripedSmithWaterman::Filter filter(true, true, min_score, 32767);
            // // Declares an alignment that stores the result
            // StripedSmithWaterman::Alignment alignment;
            // // Aligns the query to the ref
            // aligner.Align(read->seq.s, str, len, filter, &alignment, maskLen);

            // // Update alignment method
            r.tb += left_occ;
            r.te += left_occ;
            r.te2 += left_occ;

            if (r.score >= min_score)
            {
                ssw_write_sam(r, idx[r.tb].c_str(), read, strand, out, cig, mismatch);
                extended = true;
            }

            // extended_reads++;
            free(cigar);
            free(q);
            delete str;
            delete seq;
        }
        return extended;
    }

    size_t get_extended_reads()
    {
        return extended_reads;
    }

    // Adapted from SSW
    static void ssw_write_sam(kswr_t &a,
                              const char *ref_seq_name,
                              const kseq_t *read,
                              int8_t strand,
                              FILE *out,
                              std::string cigar,
                              size_t mismatches) // 0: forward aligned ; 1: reverse complement aligned
    {
        // Sam format output
        fprintf(out, "%s\t", read->name.s);
        if (a.score == 0)
            fprintf(out, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
        else
        {
            int32_t c, p;
            uint32_t mapq = -4.343 * log(1 - (double)abs(a.score - a.score2) / (double)a.score);
            mapq = (uint32_t)(mapq + 4.99);
            mapq = mapq < 254 ? mapq : 254;
            if (strand)
                fprintf(out, "16\t");
            else
                fprintf(out, "0\t");
            // TODO: Find the correct reference name.
            fprintf(out, "%s\t%d\t%d\t", ref_seq_name, a.tb + 1, mapq);
            // size_t mismatch = mark_mismatch(a.tb, a.qb, a.qe, (int8_t*)ref, (int8_t*)read_, read->seq.l, cigar, cigarLen);
            // for (c = 0; c < (*cigarLen); ++c)
            // {
            //   char letter = cigar_int_to_op((*cigar)[c]);
            //   uint32_t length = cigar_int_to_len((*cigar)[c]);
            //   fprintf(out, "%lu%c", (unsigned long)length, letter);
            // }
            // fprintf(out, "\t*\t");
            // fprintf(out, "%s", a.cigar_string.c_str());
            fprintf(out, "%s", cigar.c_str());
            fprintf(out, "\t*\t0\t0\t");
            fprintf(out, "%s", read->seq.s);
            fprintf(out, "\t");
            if (read->qual.s && strand)
            {
                for (p = read->qual.l - 1; p >= 0; --p)
                    fprintf(out, "%c", read->qual.s[p]);
            }
            else if (read->qual.s)
                fprintf(out, "%s", read->qual.s);
            else
                fprintf(out, "*");
            fprintf(out, "\tAS:i:%d", a.score);
            fprintf(out, "\tNM:i:%d\t", mismatches);
            // fprintf(out, "\tNM:i:%d\t", a.mismatches);
            if (a.score2 > 0)
                fprintf(out, "ZS:i:%d\n", a.score2);
            else
                fprintf(out, "\n");
        }
    }

    std::string to_sam()
    {
        std::string res = "@HD VN:1.6 SO:unknown\n";
        res += idx.to_sam();
        res += "@PG ID:moni PN:moni VN:0.1.0\n";
        return res; 
    }

protected:
    ms_pointers<> ms;
    slp_t ra;
    seqidx idx;

    size_t min_len = 0;
    size_t extended_reads = 0;
    size_t n = 0;

    unsigned char seq_nt4_table[256] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

    int c, sa = 2, sb = 2, i, j, k, max_rseq = 0;
    int w = 4000;
    int8_t mat[25];
    int gapo = 5, gape = 2, minsc = 0, xtra = KSW_XSTART;
    uint8_t *rseq = 0;

    bool forward_only;
};

#endif /* end of include guard: _EXTENDER_KLIB_HH */
