/* mems - Computes the MEMs from BWT and Thresholds
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
   \file mems.cpp
   \brief mems.cpp Computes the MEMs from BWT and Thresholds.
   \author Massimiliano Rossi
   \date 13/07/2020
*/

extern "C" {
#include <xerrors.h>
}

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/io.hpp>

#include <ms_pointers.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <PlainSlp.hpp>
#include <FixedBitLenCode.hpp>

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
  if(fp == NULL) error("Opening file " + filename);
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
  for (int i = 0; i < n_threads + 1; ++i)
  {
    size_t start = (size_t)((size * i) / n_threads);
    gzseek(fp, start, SEEK_SET);
    starts[i] = next_start_fastq(fp);
  }
  gzclose(fp);
  return starts;
}

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
class mems_c
{
public:

  mems_c(std::string filename)
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
  }

  // Destructor
  ~mems_c() 
  {
      // NtD
  }

  // The outfile has the following format. The first size_t integer store the
  // length l of the name. Then the following l characters stores the name of
  // the read, the next size_t integer stores the number m of MEMs, and the 
  // following m size_t pairs of integers stores the positions and lengths of 
  // the MEMs.
  void maxrimal_exact_matches(kseq_t *read, FILE* out)
  {
    auto pointers = ms.query(read->seq.s, read->seq.l);
    std::vector<size_t> lengths(pointers.size());
    std::vector<std::pair<size_t,size_t>> mems;

    size_t l = 0;
    for (size_t i = 0; i < pointers.size(); ++i)
    {
      size_t pos = pointers[i];
      while ((i + l) < read->seq.l && (pos + l) < n && (i < 1 || pos != (pointers[i-1] + 1) ) && read->seq.s[i + l] == ra.charAt(pos + l))
        ++l;

      lengths[i] = l;
      l = (l == 0 ? 0 : (l - 1));
 
      if((i == 0) or (lengths[i] >= lengths[i-1]))
        mems.push_back(make_pair(i,lengths[i]));
    }

    // Original MS computation
    // for (size_t i = 0; i < pointers.size(); ++i)
    // {
    //   size_t pos = pointers[i];
    //   while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
    //     ++l;

    //   lengths[i] = l;
    //   l = (l == 0 ? 0 : (l - 1));
    // }

    assert(lengths.size() == pointers.size());

    size_t h_length = read->name.l;
    fwrite(&h_length, sizeof(size_t), 1,out);
    fwrite(read->name.s, sizeof(char),h_length,out);
    size_t q_length = mems.size();
    fwrite(&q_length, sizeof(size_t), 1,out);
    fwrite(mems.data(), sizeof(std::pair<size_t,size_t>),q_length,out);
  }

protected:
  ms_pointers<> ms;
  slp_t ra;
  size_t n = 0;
};



char complement(char n)
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

template <typename ms_t>
struct mt_param_t
{
  // Parameters
  ms_t *ms;
  std::string pattern_filename;
  std::string out_filename;
  size_t start;
  size_t end;
  size_t wk_id;
};

template <typename ms_t>
void *mt_ms_worker(void *param)
{
  mt_param_t<ms_t> *p = (mt_param_t<ms_t>*) param;
  size_t n_reads = 0;
  size_t n_processed_reads = 0;

  FILE *out_fd;
  gzFile fp;

  if ((out_fd = fopen(p->out_filename.c_str(), "w")) == nullptr)
    error("open() file " + p->out_filename + " failed");

  if ((fp = gzopen(p->pattern_filename.c_str(), "r")) == Z_NULL)
    error("open() file " + p->pattern_filename + " failed");

  gzseek(fp, p->start, SEEK_SET);

  kseq_t rev;
  int l;

  kseq_t *seq = kseq_init(fp);
  while ((ks_tell(seq) < p->end) && ((l = kseq_read(seq)) >= 0))
  {

    p->ms->maxrimal_exact_matches(seq,out_fd);

  }

  kseq_destroy(seq);
  gzclose(fp);
  fclose(out_fd);

  return NULL;
}

template <typename ms_t>
void mt_ms(ms_t *ms, std::string pattern_filename, std::string out_filename, size_t n_threads)
{
  pthread_t t[n_threads] = {0};
  mt_param_t<ms_t> params[n_threads];
  std::vector<size_t> starts = split_fastq(pattern_filename, n_threads);
  for(size_t i = 0; i < n_threads; ++i)
  {
    params[i].ms = ms;
    params[i].pattern_filename = pattern_filename;
    params[i].out_filename = out_filename + "_" + std::to_string(i) + ".mems.tmp.out";
    params[i].start = starts[i];
    params[i].end = starts[i+1];
    params[i].wk_id = i;
    xpthread_create(&t[i], NULL, &mt_ms_worker<ms_t>, &params[i], __LINE__, __FILE__);
  }

  for(size_t i = 0; i < n_threads; ++i)
  {
    xpthread_join(t[i],NULL,__LINE__,__FILE__);
  }

  // sleep(5);


  return;
}


////////////////////////////////////////////////////////////////////////////////
/// Single Thread
////////////////////////////////////////////////////////////////////////////////
template <typename ms_t>
size_t st_ms(ms_t *ms, std::string pattern_filename, std::string out_filename)
{
  size_t n_reads = 0;
  size_t n_processed_reads = 0;
  kseq_t rev;
  int l;
  FILE *out_fd;

  out_filename += "_0.mems.tmp.out";

  if ((out_fd = fopen(out_filename.c_str(), "w")) == nullptr)
    error("open() file " + out_filename + " failed");

  gzFile fp = gzopen(pattern_filename.c_str(), "r");
  kseq_t* seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0)
  {

    ms->maxrimal_exact_matches(seq, out_fd);

  }

  kseq_destroy(seq);
  gzclose(fp);
  fclose(out_fd);

  // sleep(5);

  return n_processed_reads;
}


typedef std::pair<std::string, std::vector<uint8_t>> pattern_t;

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  std::string patterns = ""; // path to patterns file
  std::string output   = ""; // output file prefix
  size_t l = 25;             // minumum MEM length
  size_t th = 1;             // number of threads
  bool shaped_slp = false;   // use shaped slp
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-p patterns] [-o output] [-t threads] [-l len] [-q shaped_slp] [-b batch]\n\n" +
                    "Copmputes the matching statistics of the reads in the pattern against the reference index in infile.\n" +
                    "shaped_slp: [boolean] - use shaped slp. (def. false)\n" +
                    "   pattens: [string]  - path to patterns file.\n" +
                    "    output: [string]  - output file prefix.\n" +
                    "       len: [integer] - minimum MEM lengt (def. 25)\n" +
                    "    thread: [integer] - number of threads (def. 1)\n");

  std::string sarg;
  while ((c = getopt(argc, argv, "l:hp:o:t:")) != -1)
  {
    switch (c)
    {
    case 'p':
      arg.patterns.assign(optarg);
      break;
    case 'o':
      arg.output.assign(optarg);
      break;
    case 'l':
      sarg.assign(optarg);
      arg.l = stoi(sarg);
      break;
    case 't':
      sarg.assign(optarg);
      arg.th = stoi(sarg);
      break;
    case 'q':
      arg.shaped_slp = true;
      break;
    case 'h':
      error(usage);
    case '?':
      error("Unknown option.\n", usage);
      exit(1);
    }
  }
  // the only input parameter is the file name
  if (argc == optind + 1)
  {
    arg.filename.assign(argv[optind]);
  }
  else
  {
    error("Invalid number of arguments\n", usage);
  }
}

//********** end argument options ********************

template <typename ms_t>
void dispatcher(Args &args)
{
  verbose("Construction of the matching statistics data structure");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  ms_t ms(args.filename);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::string base_name = basename(args.filename.data());
  std::string out_filename = args.patterns + "_" + base_name;
  if(args.output != "")
    out_filename = args.output;

  if (is_gzipped(args.patterns))
  {
    verbose("The input is gzipped - forcing single thread matching statistics.");
    args.th = 1;
  }

  if (args.th == 1)
    st_ms<ms_t>(&ms, args.patterns, out_filename);
  else
    mt_ms<ms_t>(&ms, args.patterns, out_filename, args.th);

  // TODO: Merge the SAM files.

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  verbose("Printing plain output");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::ofstream f_mems(out_filename + ".mems");

  if (!f_mems.is_open())
    error("open() file " + std::string(out_filename) + ".mems failed");

  size_t n_seq = 0;
  for (size_t i = 0; i < args.th; ++i)
  {
    std::string tmp_filename = out_filename + "_" + std::to_string(i) + ".mems.tmp.out";
    FILE *in_fd;

    if ((in_fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
      error("open() file " + tmp_filename + " failed");

    size_t length = 0;
    size_t m = 100; // Reserved size for pointers and lengths
    std::vector<std::pair<size_t,size_t>> mem(m);
    size_t s = 100; // Reserved size for read name
    char* rname = (char *)malloc(s * sizeof(char));
    while (!feof(in_fd) and fread(&length, sizeof(size_t), 1, in_fd) > 0)
    {
      // Reading read name
      if (s < length)
      {
        // Resize lengths and pointers
        s = length;
        rname = (char *)realloc(rname, m * sizeof(char));
      }

      if ((fread(rname, sizeof(char), length, in_fd)) != length)
        error("fread() file " + std::string(tmp_filename) + " failed");

      // TODO: Store the fasta headers somewhere
      f_mems << ">" + std::string(rname,length) << endl;

      // Reading MEMs
      if ((fread(&length, sizeof(size_t), 1, in_fd)) != 1)
        error("fread() file " + std::string(tmp_filename) + " failed");

      if (m < length)
      {
        // Resize lengths and pointers
        m = length;
        mem.resize(m);
      }

      if ((fread(mem.data(), sizeof(std::pair<size_t,size_t>), length, in_fd)) != length)
        error("fread() file " + std::string(tmp_filename) + " failed");

      // TODO: Store the fasta headers somewhere
      // f_mems << ">" + std::to_string(n_seq) << endl;
      for (size_t i = 0; i < length; ++i)
        f_mems << "(" << mem[i].first << "," << mem[i].second << ") ";
      f_mems << endl;

      n_seq++;
    }
    fclose(in_fd);
    if (std::remove(tmp_filename.c_str()) != 0)
      error("remove() file " + tmp_filename + " failed");
  }

  f_mems.close();

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());
}

int main(int argc, char *const argv[])
{
  Args args;
  parseArgs(argc, argv, args);

  if (args.shaped_slp)
  {
    dispatcher<mems_c<shaped_slp_t>>(args);
  }
  else
  {
    dispatcher<mems_c<plain_slp_t>>(args);
  }
  return 0;
}