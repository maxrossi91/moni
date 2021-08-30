/* extender_reads_dispatcher - Dispatches the reads in single and multithread.
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
   \file extender_reads_dispatcher.cpp
   \brief extender_reads_dispatcher.cpp Dispatches the reads in single and multithread.
   \author Massimiliano Rossi
   \date 29/04/2021
*/

#ifndef _READS_DISPATCHER_HH
#define _READS_DISPATCHER_HH

extern "C"{
#include <xerrors.h>
}

#include <common.hpp>
#include <kseq.h>
#include <zlib.h>

////////////////////////////////////////////////////////////////////////////////
/// kseq extra
////////////////////////////////////////////////////////////////////////////////

static inline size_t ks_tell(kseq_t *seq)
{
    return gztell(seq->f->f) - seq->f->end + seq->f->begin;
}

void copy_kstring_t(kstring_t &l, kstring_t &r)
{
    l.l = r.l;
    l.m = r.m;
    l.s = (char *)malloc(l.m);
    for (size_t i = 0; i < r.m; ++i)
        l.s[i] = r.s[i];
}
void copy_kseq_t(kseq_t *l, kseq_t *r)
{
    copy_kstring_t(l->name, r->name);
    copy_kstring_t(l->comment, r->comment);
    copy_kstring_t(l->seq, r->seq);
    copy_kstring_t(l->qual, r->qual);
    l->last_char = r->last_char;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// xerror extra (conditions)
////////////////////////////////////////////////////////////////////////////////

#ifndef Thread_error_wait
    #define Thread_error_wait 5
#endif

// cond
int xpthread_cond_init(pthread_cond_t *cond, const pthread_condattr_t *attr, int linea, const char *file)
{
    int e = pthread_cond_init(cond, attr);
    if (e != 0)
    {
        xperror(e, "Error in pthread_cond_init");
        fprintf(stderr, "== %d == Line: %d, File: %s\n", getpid(), linea, file);
        sleep(Thread_error_wait); // do not kill immediately other threads
        exit(1);
    }
    return e;
}

int xpthread_cond_destroy(pthread_cond_t *cond, int linea, const char *file)
{
    int e = pthread_cond_destroy(cond);
    if (e != 0)
    {
        xperror(e, "Error in pthread_cond_destroy");
        fprintf(stderr, "== %d == Line: %d, File: %s\n", getpid(), linea, file);
        sleep(Thread_error_wait); // do not kill immediately other threads
        exit(1);
    }
    return e;
}

int xpthread_cond_wait(pthread_cond_t *cond, pthread_mutex_t *mutex, int linea, const char *file)
{
    int e = pthread_cond_wait(cond, mutex);
    if (e != 0)
    {
        xperror(e, "Error in pthread_cond_lock");
        fprintf(stderr, "== %d == Line: %d, File: %s\n", getpid(), linea, file);
        sleep(Thread_error_wait); // do not kill immediately other threads
        exit(1);
    }
    return e;
}

int xpthread_cond_signal(pthread_cond_t *cond, int linea, const char *file)
{
    int e = pthread_cond_signal(cond);
    if (e != 0)
    {
        xperror(e, "Error in pthread_cond_unlock");
        fprintf(stderr, "== %d == Line: %d, File: %s\n", getpid(), linea, file);
        sleep(Thread_error_wait); // do not kill immediately other threads
        exit(1);
    }
    return e;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Parallel computation
////////////////////////////////////////////////////////////////////////////////

// This should be done using buffering.
size_t next_start_fastq(gzFile fp)
{
    int c;
    // Special case when we arr at the beginning of the file.
    if ((gztell(fp) == 0) && ((c = gzgetc(fp)) != EOF) && c == '@')
        return 0;

    // Strart from the previous character
    gzseek(fp, -1, SEEK_CUR);

    std::vector<std::pair<int, size_t>> window;
    // Find the first new line
    for (size_t i = 0; i < 4; ++i)
    {
        while (((c = gzgetc(fp)) != EOF) && (c != (int)'\n'))
        {
        }
        if (c == EOF)
            return gztell(fp);
        if ((c = gzgetc(fp)) == EOF)
            return gztell(fp);
        window.push_back(std::make_pair(c, gztell(fp) - 1));
    }

    for (size_t i = 0; i < 2; ++i)
    {
        if (window[i].first == '@' && window[i + 2].first == '+')
            return window[i].second;
        if (window[i].first == '+' && window[i + 2].first == '@')
            return window[i + 2].second;
    }

    return gztell(fp);
}

// test if the file is gzipped
static inline bool is_gzipped(std::string filename)
{
    FILE *fp = fopen(filename.c_str(), "rb");
    if (fp == NULL)
        error("Opening file " + filename);
    int byte1 = 0, byte2 = 0;
    fread(&byte1, sizeof(char), 1, fp);
    fread(&byte2, sizeof(char), 1, fp);
    fclose(fp);
    return (byte1 == 0x1f && byte2 == 0x8b);
}

// Return the length of the file
// Assumes that the file is not compressed
static inline size_t get_file_size(std::string filename)
{
    if (is_gzipped(filename))
    {
        std::cerr << "The input is gzipped!" << std::endl;
        return -1;
    }
    FILE *fp = fopen(filename.c_str(), "r");
    fseek(fp, 0L, SEEK_END);
    size_t size = ftell(fp);
    fclose(fp);
    return size;
}

std::vector<size_t> split_fastq(std::string filename, size_t n_threads)
{
    //Precondition: the file is not gzipped
    // scan file for start positions and execute threads
    size_t size = get_file_size(filename);

    gzFile fp = gzopen(filename.c_str(), "r");
    if (fp == Z_NULL)
    {
        throw new std::runtime_error("Cannot open input file " + filename);
    }

    std::vector<size_t> starts(n_threads + 1);
    for (int i = 0; i < n_threads; ++i)
    {
        size_t start = (size_t)((size * i) / n_threads);
        gzseek(fp, start, SEEK_SET);
        starts[i] = next_start_fastq(fp);
    }
    starts[n_threads] = size;
    gzclose(fp);
    return starts;
}

inline char complement(const char n)
{
    switch (n)
    {
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    default:
        return n;
    }
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Merge SAMs
////////////////////////////////////////////////////////////////////////////////


// Merges te file in filename in the file pointed by fp
void append_file(const std::string filename, FILE *fp){
    const size_t buff_size = 16384;

    uint8_t buff[buff_size];
    size_t size = 0;

    struct stat filestat;
    FILE *fd;

    if ((fd = fopen(filename.c_str(), "r")) == nullptr)
        error("open() file " + std::string(filename) + " failed");

    // int fn = fileno(fd);
    // if (fstat(fn, &filestat) < 0)
    //     error("stat() file " + std::string(filename) + " failed");

    // size_t length = filestat.st_size;
    size_t length = 0;

    while((length = fread(buff, sizeof(uint8_t), buff_size, fd)) == buff_size)
        if ((fwrite(buff, sizeof(uint8_t), buff_size, fp)) != buff_size)
            error("fwrite() file " + std::string(filename) + " failed");
    
    assert(length < buff_size);
    if(length > 0)
        if ((fwrite(buff, sizeof(uint8_t), length, fp)) != length)
            error("fwrite() file " + std::string(filename) + " failed");


    fclose(fd);
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Multithreads workers
////////////////////////////////////////////////////////////////////////////////

pthread_mutex_t mutex_reads_dispatcher;
pthread_cond_t cond_reads_dispatcher;
// Critical variables
size_t n_active_threads = 0;
// std::vector<bool> active_threads;

template <typename extender_t>
struct mt_param_t
{
    // Parameters
    extender_t *extender;
    std::string pattern_filename;
    std::string sam_filename;
    size_t start;
    size_t end;
    size_t wk_id;
    // Return values
    size_t n_reads;
    size_t n_extended_reads;
};

template <typename extender_t>
void *mt_extend_worker(void *param)
{
    mt_param_t<extender_t> *p = (mt_param_t<extender_t> *)param;
    size_t n_reads = 0;
    size_t n_extended_reads = 0;

    FILE *sam_fd;
    gzFile fp;

    if ((sam_fd = fopen(p->sam_filename.c_str(), "w")) == nullptr)
        error("open() file " + p->sam_filename + " failed");

    if ((fp = gzopen(p->pattern_filename.c_str(), "r")) == Z_NULL)
        error("open() file " + p->pattern_filename + " failed");

    gzseek(fp, p->start, SEEK_SET);

    kseq_t rev;
    int l;

    kseq_t *seq = kseq_init(fp);
    while ((ks_tell(seq) < p->end) && ((l = kseq_read(seq)) >= 0))
    {

        bool fwd_extend = p->extender->extend(seq, sam_fd, 0);

        //copy seq
        copy_kseq_t(&rev, seq);

        for (size_t i = 0; i < seq->seq.l; ++i)
            rev.seq.s[i] = complement(seq->seq.s[seq->seq.l - i - 1]);

        if (rev.seq.m > rev.seq.l)
            rev.seq.s[rev.seq.l] = 0;

        bool rev_extend = p->extender->extend(&rev, sam_fd, 1);

        if (fwd_extend or rev_extend)
            n_extended_reads++;
        n_reads++;

        free(rev.name.s);
        free(rev.comment.s);
        free(rev.seq.s);
        free(rev.qual.s);
    }

    verbose("Number of extended reads block ", p->wk_id, " : ", n_extended_reads, "/", n_reads);
    p->n_reads = n_reads;
    p->n_extended_reads = n_extended_reads;
    kseq_destroy(seq);
    gzclose(fp);
    fclose(sam_fd);

    // Update the number of active threads
    xpthread_mutex_lock(&mutex_reads_dispatcher, __LINE__, __FILE__);
    {
        --n_active_threads;
        xpthread_cond_signal(&cond_reads_dispatcher, __LINE__, __FILE__);
    }
    xpthread_mutex_unlock(&mutex_reads_dispatcher, __LINE__, __FILE__);

    return NULL;
}

template <typename extender_t>
size_t mt_extend(extender_t *extender, std::string pattern_filename, std::string sam_filename, size_t n_threads, size_t k)
{
    xpthread_mutex_init(&mutex_reads_dispatcher, NULL, __LINE__, __FILE__);
    xpthread_cond_init(&cond_reads_dispatcher, NULL, __LINE__, __FILE__);

    // active_threads = std::vector<bool>(n_threads, false);
    pthread_t t[k * n_threads] = {0};
    mt_param_t<extender_t> params[k * n_threads];
    std::vector<size_t> starts = split_fastq(pattern_filename, k * n_threads);
    for (size_t i = 0; i < k * n_threads; ++i)
    {
        // Get the number of active threads
        xpthread_mutex_lock(&mutex_reads_dispatcher, __LINE__, __FILE__);
        {
            while(n_active_threads >= n_threads)
                xpthread_cond_wait(&cond_reads_dispatcher, &mutex_reads_dispatcher, __LINE__, __FILE__);
            assert(n_active_threads < n_threads);
            // Create a new thread
            params[i].extender = extender;
            params[i].pattern_filename = pattern_filename;
            params[i].sam_filename = sam_filename + "_" + std::to_string(i) + ".sam";
            params[i].start = starts[i];
            params[i].end = starts[i + 1];
            params[i].wk_id = i;
            xpthread_create(&t[i], NULL, &mt_extend_worker<extender_t>, &params[i], __LINE__, __FILE__);
            // Update the number of active threads
            ++n_active_threads;
        }
        xpthread_mutex_unlock(&mutex_reads_dispatcher, __LINE__, __FILE__);
    }

    size_t tot_reads = 0;
    size_t tot_extended_reads = 0;

    for (size_t i = 0; i < k * n_threads; ++i)
    {
        xpthread_join(t[i], NULL, __LINE__, __FILE__);
    }

    // sleep(5);
    verbose("Merging temporary SAM files");

    FILE *fd;

    if ((fd = fopen(std::string(sam_filename + ".sam").c_str(), "w")) == nullptr)
        error("open() file " + std::string(sam_filename + ".sam") + " failed");

    fprintf(fd, "%s", extender->to_sam().c_str());

    for (size_t i = 0; i < k * n_threads; ++i)
    {
        tot_reads += params[i].n_reads;
        tot_extended_reads += params[i].n_extended_reads;

        append_file(params[i].sam_filename, fd);
        if (std::remove(params[i].sam_filename.c_str()) != 0)
            error("remove() file " + params[i].sam_filename + " failed");
    }
    
    xpthread_mutex_destroy(&mutex_reads_dispatcher, __LINE__, __FILE__);
    xpthread_cond_destroy(&cond_reads_dispatcher, __LINE__, __FILE__);

    verbose("Number of extended reads: ", tot_extended_reads, "/", tot_reads);
    return tot_extended_reads;
}

////////////////////////////////////////////////////////////////////////////////
/// Single Thread
////////////////////////////////////////////////////////////////////////////////
template <typename extender_t>
size_t st_extend(extender_t *extender, std::string pattern_filename, std::string sam_filename)
{
    size_t n_reads = 0;
    size_t n_extended_reads = 0;
    kseq_t rev;
    int l;
    FILE *sam_fd;

    sam_filename += ".sam";

    if ((sam_fd = fopen(sam_filename.c_str(), "w")) == nullptr)
        error("open() file " + sam_filename + " failed");

    fprintf(sam_fd, "%s", extender->to_sam().c_str());

    gzFile fp = gzopen(pattern_filename.c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0)
    {

        bool fwd_extend = extender->extend(seq, sam_fd, 0);

        //copy seq
        copy_kseq_t(&rev, seq);

        for (size_t i = 0; i < seq->seq.l; ++i)
            rev.seq.s[i] = complement(seq->seq.s[seq->seq.l - i - 1]);

        if (rev.seq.m > rev.seq.l)
            rev.seq.s[rev.seq.l] = 0;

        bool rev_extend = extender->extend(&rev, sam_fd, 1);

        if (fwd_extend or rev_extend)
            n_extended_reads++;
        n_reads++;

        free(rev.name.s);
        free(rev.comment.s);
        free(rev.seq.s);
        free(rev.qual.s);
    }

    verbose("Number of extended reads: ", n_extended_reads, "/", n_reads);
    kseq_destroy(seq);
    gzclose(fp);
    fclose(sam_fd);

    // sleep(5);

    return n_extended_reads;
}

#endif /* end of include guard: _READS_DISPATCHER_HH */
