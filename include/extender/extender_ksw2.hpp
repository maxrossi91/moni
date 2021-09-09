/* extender_ksw2 - Extend the MEMs of the reads to the reference using the ksw2 library for SW
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
   \file extender_ksw2.cpp
   \brief extender_ksw2.cpp Extend the MEMs of the reads to the reference using the ksw2 library for SW
   \author Massimiliano Rossi
   \date 13/07/2020
*/

#ifndef _EXTENDER_KSW2_HH
#define _EXTENDER_KSW2_HH

#include <common.hpp>

#include <sdsl/io.hpp>

#include <ms_pointers.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <PlainSlp.hpp>
#include <FixedBitLenCode.hpp>

#include <ksw2.h>

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

template <typename slp_t,
          typename ms_t>
class extender
{
public:
    // using SelSd = SelectSdvec<>;
    // using DagcSd = DirectAccessibleGammaCode<SelSd>;
    // using SlpT = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;
    using ll = long long int;

    typedef struct{

        size_t min_len = 25;    // Minimum MEM length
        size_t ext_len = 100;   // Extension length
        size_t top_k = 1;       // Report the top_k alignments

        // ksw2 parameters
        int8_t smatch = 2;      // Match score default
        int8_t smismatch = 4;   // Mismatch score default
        int8_t gapo = 4;        // Gap open penalty
        int8_t gapo2 = 13;      // Gap open penalty
        int8_t gape = 2;        // Gap extension penalty
        int8_t gape2 = 1;       // Gap extension penalty
        int end_bonus = 400;    // Bonus to add at the extension score to declare the alignment

        int w = -1;             // Band width
        int zdrop = -1;         // Zdrop enable

        bool forward_only = true;      // Align only 

    } config_t;

    // extender(std::string filename,
    //         size_t min_len_ = 50,
    //         bool forward_only_ = true) : min_len(min_len_),
    //                                      forward_only(forward_only_)

    extender(std::string filename,
            config_t config = config_t()) : 
                min_len(config.min_len),        // Minimum MEM length
                ext_len(config.ext_len),        // Extension length
                top_k(config.top_k),            // Report the top_k alignments
                smatch(config.smatch),          // Match score default
                smismatch(config.smismatch),    // Mismatch score default
                gapo(config.gapo),              // Gap open penalty
                gapo2(config.gapo2),            // Gap open penalty
                gape(config.gape),              // Gap extension penalty
                gape2(config.gape2),            // Gap extension penalty
                end_bonus(config.end_bonus),    // Bonus to add at the extension score to declare the alignment
                w(config.w),                    // Band width
                zdrop(config.zdrop),            // Zdrop enable
                forward_only(config.forward_only)
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

        std::string filename_slp = filename + get_slp_file_extension<slp_t>();
        verbose("Loading random access file: " + filename_slp);
        t_insert_start = std::chrono::high_resolution_clock::now();

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

        ksw_gen_simple_mat(m, mat, smatch, -smismatch);

        t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Local aligner initialization complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

        verbose("Minimum MEM length: ", min_len);

    }

    // Destructor
    ~extender()
    {
        // NtD
    }

    bool extend(kseq_t *read, FILE *out, uint8_t strand)
    {

        bool extended = false;

        mem_t mem = find_longest_mem(read);

        // Extend the read
        if (mem.len >= min_len)
        {

            // Extractin left and right context of the read
            // lcs: left context sequence
            size_t lcs_len = mem.idx;
            uint8_t *lcs = (uint8_t *)malloc(lcs_len);
            // verbose("lcs: " + std::string(read->seq.s).substr(0,lcs_len));
            // Convert A,C,G,T,N into 0,1,2,3,4
            // The left context is reversed
            for (size_t i = 0; i < lcs_len; ++i)
                lcs[lcs_len - i - 1] = seq_nt4_table[(int)read->seq.s[i]];

            // rcs: right context sequence
            size_t rcs_occ = (mem.idx + mem.len); // The first character of the right context
            size_t rcs_len = read->seq.l - rcs_occ;
            uint8_t *rcs = (uint8_t *)malloc(rcs_len);
            // verbose("rcs: " + std::string(read->seq.s).substr(rcs_occ,rcs_len));
            // Convert A,C,G,T,N into 0,1,2,3,4
            for (size_t i = 0; i < rcs_len; ++i)
                rcs[i] = seq_nt4_table[(int)read->seq.s[rcs_occ + i]];

            int32_t min_score = 20 + 8 * log(read->seq.l);
            // verbose("Number of occurrences: " + std::to_string(occs.size()));

            int32_t score = extend(
                mem.pos,
                mem.len,
                lcs,     // Left context of the read
                lcs_len, // Left context of the read lngth
                rcs,     // Right context of the read
                rcs_len  // Right context of the read length
            );

            if (score > min_score)
            {
                extend(mem.pos, mem.len, lcs, lcs_len, rcs, rcs_len, false, 0, min_score, read, strand, out);
                extended = true;
            }
            
            delete lcs;
            delete rcs;
        }
        return extended;
    }

    typedef struct mem_t
    {
        size_t pos = 0;           // Position in the reference
        size_t len = 0;           // Length
        size_t idx = 0;           // Position in the pattern

        mem_t(size_t p, size_t l, size_t i)
        {
            pos = p; // Position in the reference
            len = l; // Length of the MEM
            idx = i; // Position in the read
        }

    } mem_t;

    inline mem_t find_longest_mem(kseq_t *read)
    {
        size_t mem_pos = 0;
        size_t mem_len = 0;
        size_t mem_idx = 0;

        auto pointers = ms.query(read->seq.s, read->seq.l);
        std::vector<size_t> lengths(pointers.size());
        size_t l = 0;
        size_t n_Ns = 0;
        for (size_t i = 0; i < pointers.size(); ++i)
        {
            size_t pos = pointers[i];
            while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
            {
                if (read->seq.s[i + l] == 'N')
                    n_Ns++;
                else
                    n_Ns = 0;
                ++l;
            }

            lengths[i] = l;
            l = (l == 0 ? 0 : (l - 1));

            // Update MEM
            if (lengths[i] > mem_len and n_Ns < lengths[i])
            {
                mem_len = lengths[i];
                mem_pos = pointers[i];
                mem_idx = i;
            }
        }

        return mem_t(mem_pos, mem_len, mem_idx);
    }

    size_t get_extended_reads()
    {
        return extended_reads;
    }

    // If score_only is true we compute the score of the alignment.
    // If score_only is false, we extend again the read and we write the result
    // in the SAM file, so we need to give the second best score.
    int32_t extend(
        const size_t mem_pos,
        const size_t mem_len,
        const uint8_t *lcs,           // Left context of the read
        const size_t lcs_len,         // Left context of the read lngth
        const uint8_t *rcs,           // Right context of the read
        const size_t rcs_len,         // Right context of the read length
        const bool score_only = true, // Report only the score
        const int32_t score2 = 0,     // The score of the second best alignment
        const int32_t min_score = 0,  // The minimum score to call an alignment
        const kseq_t *read = nullptr, // The read that has been aligned
        int8_t strand = 0,            // 0: forward aligned ; 1: reverse complement aligned
        FILE *out = nullptr,          // The SAM file pointer
        const bool realign = false    // Realign globally the read
    )
    {
        int flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT;

        if (score_only)
            flag = KSW_EZ_SCORE_ONLY;

        int score_lc = 0;
        int score_rc = 0;

        ksw_extz_t ez_lc;
        ksw_extz_t ez_rc;
        ksw_extz_t ez;
        memset(&ez_lc, 0, sizeof(ksw_extz_t));
        memset(&ez_rc, 0, sizeof(ksw_extz_t));
        memset(&ez, 0, sizeof(ksw_extz_t));
        // TODO: Update end_bonus according to the MEM contribution to the score

        // Extract the context from the reference
        // lc: left context
        ksw_reset_extz(&ez_lc);
        if (lcs_len > 0)
        {
            size_t lc_occ = (mem_pos > ext_len ? mem_pos - ext_len : 0);
            size_t lc_len = (mem_pos > ext_len ? ext_len : ext_len - mem_pos);
            char *tmp_lc = (char *)malloc(ext_len);
            ra.expandSubstr(lc_occ, lc_len, tmp_lc);
            // verbose("lc: " + std::string(lc));
            // Convert A,C,G,T,N into 0,1,2,3,4
            // The left context is reversed
            uint8_t *lc = (uint8_t *)malloc(ext_len);
            for (size_t i = 0; i < lc_len; ++i)
                lc[lc_len - i - 1] = seq_nt4_table[(int)tmp_lc[i]];

            delete tmp_lc;

            // Query: lcs
            // Target: lc
            // verbose("aligning lc and lcs");
            ksw_extz2_sse(km, lcs_len, (uint8_t *)lcs, lc_len, (uint8_t *)lc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_lc);
            score_lc = ez_lc.mqe;
            // verbose("lc score: " + std::to_string(score_lc));
            // Check if the extension reached the end or the query
            assert(score_only or ez_lc.reach_end);

            // std::string blc = print_BLAST_like((uint8_t*)lc,(uint8_t*)lcs,ez_lc.cigar,ez_lc.n_cigar);
            // std::cout<<blc;

            delete lc;
        }

        // rc: right context
        ksw_reset_extz(&ez_rc);
        if (rcs_len > 0)
        {
            size_t rc_occ = mem_pos + mem_len;
            size_t rc_len = (rc_occ < n - ext_len ? ext_len : n - rc_occ);
            char *rc = (char *)malloc(ext_len);
            ra.expandSubstr(rc_occ, rc_len, rc);
            // verbose("rc: " + std::string(rc));
            // Convert A,C,G,T,N into 0,1,2,3,4
            for (size_t i = 0; i < rc_len; ++i)
                rc[i] = seq_nt4_table[(int)rc[i]];

            // Query: rcs
            // Target: rc
            // verbose("aligning rc and rcs");
            ksw_extz2_sse(km, rcs_len, (uint8_t *)rcs, rc_len, (uint8_t *)rc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_rc);
            score_rc = ez_rc.mqe;
            // verbose("rc score: " + std::to_string(score_rc));
            // Check if the extension reached the end or the query
            assert(score_only or ez_rc.reach_end);

            // std::string brc = print_BLAST_like((uint8_t*)rc,(uint8_t*)rcs,ez_rc.cigar,ez_rc.n_cigar);
            // std::cout<<brc;
            delete rc;
        }

        // Compute the final score
        int32_t score = mem_len * smatch + score_lc + score_rc;

        if (not score_only)
        {
            // Compute starting position in reference
            size_t ref_pos = mem_pos - (lcs_len > 0 ? ez_lc.mqe_t + 1 : 0);
            size_t ref_len = (lcs_len > 0 ? ez_lc.mqe_t + 1 : 0) + mem_len + (rcs_len > 0 ? ez_rc.mqe_t + 1 : 0);
            char *ref = (char *)malloc(ref_len);
            ra.expandSubstr(ref_pos, ref_len, ref);
            // Convert A,C,G,T,N into 0,1,2,3,4
            for (size_t i = 0; i < ref_len; ++i)
                ref[i] = seq_nt4_table[(int)ref[i]];

            // Convert the read
            size_t seq_len = read->seq.l;
            uint8_t *seq = (uint8_t *)malloc(seq_len);
            for (size_t i = 0; i < seq_len; ++i)
                seq[i] = seq_nt4_table[(int)read->seq.s[i]];

            char *tmp = (char *)calloc(max(ref_len, seq_len), 1);

            if (realign)
            {
                // Realign the whole sequence globally
                flag = KSW_EZ_RIGHT;
                ksw_reset_extz(&ez);
                ksw_extz2_sse(km, seq_len, (uint8_t *)seq, ref_len, (uint8_t *)ref, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez);

                // std::string bfull = print_BLAST_like((uint8_t*)ref,seq,ez.cigar,ez.n_cigar);
                // std::cout << bfull;

                // Example were ez.score is lrger than score:
                // Left context alignment
                // 22333022233022233302223302223
                // ||||  ||||||||| ||||||*|||||*
                // 2233  222330222 3302220302220
                // Right context alignment
                // 33022233022233022233022233022233      0222334
                // *||||||||||||||||*||||||||||||||      ||||||*
                // 130222330222330220330222330222332222330222330
                // [INFO] 16:26:16 - Message: old score:  130  new score:  140
                // Global alignment
                // 2233    3022233022233302223  30222330222330222330222330222  330222  330222330222330222330222330222330222334
                // ||||    |||||||||||*| ||||*  |||||||||||||||||||||||||||||  *|||||  |||||||||||*||||||||||||||*|||||||||||*
                // 223322233022233022203 02220  30222330222330222330222330222  130222  330222330220330222330222332222330222330
                // The original occurrence of the MEM has been shifted to the left by 6 positions,
                // reducing the gap in the right context, and moving in to the left context.

                assert(ez.score >= score);

                // Concatenate the CIGAR strings

                std::string cigar_s;
                for (size_t i = 0; i < ez.n_cigar; ++i)
                    cigar_s += std::to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];

                // Compute the MD:Z field and thenumber of mismatches
                std::pair<std::string, size_t> md_nm = write_MD_core((uint8_t *)ref, seq, ez.cigar, ez.n_cigar, tmp, 0);
                std::pair<std::string,size_t> pos = idx.index(ref_pos);
                write_sam(ez.score, score2, min_score, pos.second, pos.first.c_str(), read, strand, out, cigar_s, md_nm.first, md_nm.second);
            }
            else
            {
                // Concatenate the CIGAR strings
                size_t n_cigar = ez_lc.n_cigar + ez_rc.n_cigar + 1;
                uint32_t *cigar = (uint32_t *)calloc(n_cigar, sizeof(uint32_t));
                size_t i = 0;

                for (size_t j = 0; j < ez_lc.n_cigar; ++j)
                    cigar[i++] = ez_lc.cigar[ez_lc.n_cigar - j - 1];

                if (ez_lc.n_cigar > 0 and ((cigar[i - 1] & 0xf) == 0))
                { // If the previous operation is also an M then merge the two operations
                    cigar[i - 1] += (((uint32_t)mem_len) << 4);
                    --n_cigar;
                }
                else
                    cigar[i++] = (((uint32_t)mem_len) << 4);

                if (ez_rc.n_cigar > 0)
                {
                    if ((ez_rc.cigar[0] & 0xf) == 0)
                    { // If the next operation is also an M then merge the two operations
                        cigar[i - 1] += ez_rc.cigar[0];
                        --n_cigar;
                    }
                    else
                        cigar[i++] = ez_rc.cigar[0];
                }

                for (size_t j = 1; j < ez_rc.n_cigar; ++j)
                    cigar[i++] = ez_rc.cigar[j];

                assert(i <= n_cigar);

                // std::string bfull = print_BLAST_like((uint8_t*)ref,seq,cigar,n_cigar);
                // std::cout << bfull;

                std::string cigar_s;
                for (size_t i = 0; i < n_cigar; ++i)
                    cigar_s += std::to_string(cigar[i] >> 4) + "MID"[cigar[i] & 0xf];

                // Compute the MD:Z field and thenumber of mismatches
                std::pair<std::string, size_t> md_nm = write_MD_core((uint8_t *)ref, seq, cigar, n_cigar, tmp, 0);
                std::pair<std::string,size_t> pos = idx.index(ref_pos);
                write_sam(score, score2, min_score, pos.second, pos.first.c_str(), read, strand, out, cigar_s, md_nm.first, md_nm.second);

                delete cigar;
            }
            delete tmp;
            delete ref;
            delete seq;
        }

        if (ez_lc.m_cigar > 0)
            delete ez_lc.cigar;
        if (ez_rc.m_cigar > 0)
            delete ez_rc.cigar;
        if (ez.m_cigar > 0)
            delete ez.cigar;

        return score;
    }


    // Readapted from https://github.com/lh3/minimap2/blob/c9874e2dc50e32bbff4ded01cf5ec0e9be0a53dd/format.c
    // tmp is a string of length max(reference length, query length)
    static std::pair<std::string, size_t> write_MD_core(const uint8_t *tseq, const uint8_t *qseq, const uint32_t *cigar, const size_t n_cigar, char *tmp, int write_tag)
    {
        std::string mdz;
        int i, q_off, t_off, l_MD = 0, NM = 0;
        if (write_tag)
            mdz += "MD:Z:"; //printf("MD:Z:");
        for (i = q_off = t_off = 0; i < (int)n_cigar; ++i)
        {
            int j, op = cigar[i] & 0xf, len = cigar[i] >> 4;
            assert((op >= 0 && op <= 3) || op == 7 || op == 8);
            if (op == 0 || op == 7 || op == 8)
            { // match
                for (j = 0; j < len; ++j)
                {
                    if (qseq[q_off + j] != tseq[t_off + j])
                    {
                        mdz += std::to_string(l_MD) + "ACGTN"[tseq[t_off + j]];
                        // printf("%d%c", l_MD, "ACGTN"[tseq[t_off + j]]);
                        l_MD = 0;
                        ++NM;
                    }
                    else
                        ++l_MD;
                }
                q_off += len, t_off += len;
            }
            else if (op == 1)
            { // insertion to ref
                q_off += len;
                NM += len;
            }
            else if (op == 2)
            { // deletion from ref
                for (j = 0, tmp[len] = 0; j < len; ++j)
                    tmp[j] = "ACGTN"[tseq[t_off + j]];
                mdz += std::to_string(l_MD) + "^" + std::string(tmp);
                // printf("%d^%s", l_MD, tmp);
                l_MD = 0;
                t_off += len;
                NM += len;
            }
            else if (op == 3)
            { // reference skip
                t_off += len;
            }
        }
        if (l_MD > 0)
            mdz += std::to_string(l_MD); //printf("%d", l_MD);
        // assert(t_off == r->re - r->rs && q_off == r->qe - r->qs);
        return make_pair(mdz, NM);
    }

    // From https://github.com/lh3/ksw2/blob/master/cli.c
    static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
    {
        int i, j;
        a = a < 0 ? -a : a;
        b = b > 0 ? -b : b;
        for (i = 0; i < m - 1; ++i)
        {
            for (j = 0; j < m - 1; ++j)
                mat[i * m + j] = i == j ? a : b;
            mat[i * m + m - 1] = 0;
        }
        for (j = 0; j < m; ++j)
            mat[(m - 1) * m + j] = 0;
    }

    // Adapted from https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/master/src/main.c
    void write_sam(const int32_t score,
                   const int32_t score2,
                   const int32_t min_score,
                   size_t ref_pos,
                   const char *ref_seq_name,
                   const kseq_t *read,
                   int8_t strand, // 0: forward aligned ; 1: reverse complement aligned
                   FILE *out,
                   std::string &cigar,
                   std::string &md,
                   size_t mismatches)
    {
        // Sam format output
        fprintf(out, "%s\t", read->name.s);
        if (score == 0)
            fprintf(out, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
        else
        {
            int32_t c, p;
            // uint32_t mapq = -4.343 * log(1 - (double)abs(score - score2) / (double)score);
            // mapq = (uint32_t)(mapq + 4.99);
            // mapq = mapq < 254 ? mapq : 254;
            uint32_t mapq = compute_mapq(score, score2, min_score, read->seq.l);
            if (strand)
                fprintf(out, "16\t");
            else
                fprintf(out, "0\t");
            // TODO: Find the correct reference name.
            fprintf(out, "%s\t%d\t%d\t", ref_seq_name, ref_pos + 1, mapq);
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
            fprintf(out, "\tAS:i:%d", score);
            fprintf(out, "\tNM:i:%d", mismatches);
            if (score2 > 0)
                fprintf(out, "\tZS:i:%d", score2);
            fprintf(out, "\tMD:Z:%s\n", md.c_str());
        }
    }

    /*!
    Compute the mapping quality of the alignment
    Inspired from https://github.com/BenLangmead/bowtie2/blob/4512b199768e562e8627ffdfd9253affc96f6fc6/unique.h
  */
    size_t compute_mapq(
        const int32_t score,     // Best alignment score
        const int32_t score2,    // Second best alignemt score
        const int32_t min_score, // Minimum alignemt score
        const size_t read_l      // Read length
    )
    {
        int32_t max_score = read_l * smatch;
        int32_t best = max_score - score;
        size_t best_bin = (size_t)((double)best * (10.0 / (double)(max_score - min_score)) + 0.5);
        if (score2 > min_score)
        {
            int32_t diff = score - score2;
            size_t diff_bin = (size_t)((double)diff * (10.0 / (double)(max_score - min_score)) + 0.5);
            if (best == max_score)
                return unp_sec_perf[best_bin];
            else
                return unp_sec[diff_bin][best_bin];
        }
        else
        {
            if (best == max_score)
                return unp_nosec_perf;
            else
                return unp_nosec[best_bin];
        }
    }

    static std::string print_BLAST_like(const uint8_t *tseq, const uint8_t *qseq, const uint32_t *cigar, const size_t n_cigar)
    {
        std::string target_o;
        std::string bars_o;
        std::string seq_o;

        int i, q_off, t_off, l_MD = 0;
        for (i = q_off = t_off = 0; i < (int)n_cigar; ++i)
        {
            int j, op = cigar[i] & 0xf, len = cigar[i] >> 4;
            assert((op >= 0 && op <= 3) || op == 7 || op == 8);
            if (op == 0 || op == 7 || op == 8)
            { // match
                for (j = 0; j < len; ++j)
                {
                    if (qseq[q_off + j] != tseq[t_off + j])
                    {
                        bars_o += "*";
                    }
                    else
                    {
                        bars_o += "|";
                    }
                    target_o += std::to_string(tseq[t_off + j]);
                    seq_o += std::to_string(qseq[q_off + j]);
                }
                q_off += len, t_off += len;
            }
            else if (op == 1)
            { // insertion to ref
                for (j = 0; j < len; ++j)
                {
                    target_o += " ";
                    bars_o += " ";
                    seq_o += std::to_string(qseq[q_off + j]);
                }
                q_off += len;
            }
            else if (op == 2)
            { // deletion from ref
                for (j = 0; j < len; ++j)
                {
                    seq_o += " ";
                    bars_o += " ";
                    target_o += std::to_string(tseq[t_off + j]);
                }
                t_off += len;
            }
            else if (op == 3)
            { // reference skip
                for (j = 0; j < len; ++j)
                {
                    seq_o += " ";
                    bars_o += " ";
                    target_o += std::to_string(tseq[t_off + j]);
                }
                t_off += len;
            }
        }
        return target_o + "\n" + bars_o + "\n" + seq_o + "\n";
    }

    std::string to_sam()
    {
        std::string res = "@HD VN:1.6 SO:unknown\n";
        res += idx.to_sam();
        res += "@PG\tID:moni\tPN:moni\tVN:0.1.0\n";
        return res; 
    }

protected:
    ms_t ms;
    slp_t ra;
    seqidx idx;
    // SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;

    const size_t min_len = 0;
    const size_t ext_len = 100;   // Extension length
    size_t extended_reads = 0;
    size_t n = 0;
    const size_t top_k = 1; // report the top_k alignments

    const unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
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

    const int8_t smatch = 2;    // Match score default
    const int8_t smismatch = 4; // Mismatch score default
    const int8_t gapo = 4;      // Gap open penalty
    const int8_t gapo2 = 13;    // Gap open penalty
    const int8_t gape = 2;      // Gap extension penalty
    const int8_t gape2 = 1;     // Gap extension penalty
    const int end_bonus = 400;  // Bonus to add at the extension score to declare the alignment

    const int w = -1; // Band width
    const int zdrop = -1;

    void *km = 0; // Kalloc

    // int8_t max_rseq = 0;

    const int m = 5;
    int8_t mat[25];
    // int minsc = 0, xtra = KSW_XSTART;
    // uint8_t *rseq = 0;

    const bool forward_only;

    // From https://github.com/BenLangmead/bowtie2/blob/4512b199768e562e8627ffdfd9253affc96f6fc6/unique.cpp
    // There is no valid second-best alignment and the best alignment has a
    // perfect score.
    const uint32_t unp_nosec_perf = 44;

    // There is no valid second-best alignment.  We stratify the alignment
    // score of the best alignment into 10 bins.
    const uint32_t unp_nosec[11] = {
        43, 42, 41, 36, 32, 27, 20, 11, 4, 1, 0};

    // The best alignment has a perfect score, and we stratify the distance
    // between best and second-best alignment scores into 10 bins.
    const uint32_t unp_sec_perf[11] = {
        2, 16, 23, 30, 31, 32, 34, 36, 38, 40, 42};

    // The best alignment has a non-perfect score, and we stratify both by best
    // alignment score (specifically, the maximum score minus the best "best")
    // and by the distance between the best and second-best alignment scores
    // ("difference").  Each is stratified into 10 bins.  Each row is a
    // difference (smaller elts = smaller differences) and each column is a
    // best score (smaller elts = higher best alignment scores).
    const uint32_t unp_sec[11][11] = {
        {2, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0},
        {20, 14, 7, 3, 2, 1, 0, 0, 0, 0, 0},
        {20, 16, 10, 6, 3, 1, 0, 0, 0, 0, 0},
        {20, 17, 13, 9, 3, 1, 1, 0, 0, 0, 0},
        {21, 19, 15, 9, 5, 2, 2, 0, 0, 0, 0},
        {22, 21, 16, 11, 10, 5, 0, 0, 0, 0, 0},
        {23, 22, 19, 16, 11, 0, 0, 0, 0, 0, 0},
        {24, 25, 21, 30, 0, 0, 0, 0, 0, 0, 0},
        {30, 26, 29, 0, 0, 0, 0, 0, 0, 0, 0},
        {30, 27, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    };

    //
    // Paired mapping quality:
    //

    // There is no valid second-best alignment and the best alignment has a
    // perfect score.
    const uint32_t pair_nosec_perf = 44;
};

#endif /* end of include guard: _EXTENDER_KSW2_HH */
