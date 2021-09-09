/* extend_ksw2 - Extend the MEMs of the reads to the reference
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
   \file extend_ksw2.cpp
   \brief extend_ksw2.cpp Extend the MEMs of the reads to the reference.
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
#include <extender_ksw2.hpp>
#include <extend_reads_dispatcher.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <PlainSlp.hpp>
#include <FixedBitLenCode.hpp>

#include <libgen.h>

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  std::string patterns = ""; // path to patterns file
  std::string output   = ""; // output file prefix
  size_t l = 25;             // minumum MEM length
  size_t th = 1;             // number of threads
  size_t b = 1;              // number of batches per thread pool
  bool shaped_slp = false;   // use shaped slp
  size_t ext_len = 100;      // Extension length
  // size_t top_k = 1;       // Report the top_k alignments

  // ksw2 parameters
  int8_t smatch = 2;      // Match score default
  int8_t smismatch = 4;   // Mismatch score default
  int8_t gapo = 4;        // Gap open penalty
  int8_t gapo2 = 13;      // Gap open penalty
  int8_t gape = 2;        // Gap extension penalty
  int8_t gape2 = 1;       // Gap extension penalty
  // int end_bonus = 400;    // Bonus to add at the extension score to declare the alignment

  // int w = -1;             // Band width
  // int zdrop = -1;         // Zdrop enable
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-p patterns] [-t threads] [-l len] [-q shaped_slp] [-b batch] [-L ext_l] [-A smatch] [-B smismatc] [-O gapo] [-E gape]\n\n" +
                    "Extends the MEMs of the reads in the pattern against the reference index in infile.\n" +
                    "shaped_slp: [boolean] - use shaped slp. (def. false)\n" +
                    "   pattens: [string]  - path to patterns file.\n" +
                    "    output: [string]  - output file prefix.\n" +
                    "       len: [integer] - minimum MEM lengt (def. 25)\n" +
                    "    thread: [integer] - number of threads (def. 1)\n" +
                    "     ext_l: [integer] - length of reference substring for extension (def. " + std::to_string(arg.ext_len) + ")\n" +
                    "    smatch: [integer] - match score value (def. " + std::to_string(arg.smatch) + ")\n" +
                    " smismatch: [integer] - mismatch penalty value (def. " + std::to_string(arg.smismatch) + ")\n" +
                    "      gapo: [integer] - gap open penalty value (def. " + std::to_string(arg.gapo) + "," + std::to_string(arg.gapo2) + ")\n" +
                    "      gape: [integer] - gap extension penalty value (def. " + std::to_string(arg.gape) + "," + std::to_string(arg.gape2) + ")\n" +
                    "     batch: [integer] - number of batches per therad pool (def. 1)\n");

  std::string sarg;
  char* s;
  while ((c = getopt(argc, argv, "l:hp:o:b:t:qA:B:O:E:L:")) != -1)
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
    case 'b':
      sarg.assign(optarg);
      arg.b = stoi(sarg);
      break;
    case 'L':
      sarg.assign(optarg);
      arg.ext_len = stoi(sarg);
      break;
    case 'A':
      sarg.assign(optarg);
      arg.smatch = stoi(sarg);
      break;
    case 'B':
      sarg.assign(optarg);
      arg.smismatch = stoi(sarg);
      break;
    case 'O':
      arg.gapo = arg.gapo2 = strtol(optarg, &s, 10);
      if (*s == ',') arg.gapo2 = strtol(s+1, &s, 10);
      break;
    case 'E':
      arg.gape = arg.gape2 = strtol(optarg, &s, 10);
      if (*s == ',') arg.gape2 = strtol(s+1, &s, 10);
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


template<typename extender_t>
typename extender_t::config_t configurer(Args &args){
  typename extender_t::config_t config;
  
  config.min_len    = args.l;           // Minimum MEM length
  config.ext_len    = args.ext_len;     // Extension length

  // ksw2 parameters
  config.smatch     = args.smatch;      // Match score default
  config.smismatch  = args.smismatch;   // Mismatch score default
  config.gapo       = args.gapo;        // Gap open penalty
  config.gapo2      = args.gapo2;       // Gap open penalty
  config.gape       = args.gape;        // Gap extension penalty
  config.gape2      = args.gape2;       // Gap extension penalty

  return config;
}

template<typename extender_t>
void dispatcher(Args &args){
  verbose("Construction of the extender");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();


  extender_t extender(args.filename, configurer<extender_t>(args));

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::string base_name = basename(args.filename.data());
  std::string sam_filename = args.patterns + "_" + base_name + "_" + std::to_string(args.l);
  if(args.output != "")
    sam_filename = args.output;

  if (is_gzipped(args.patterns))
  {
    verbose("The input is gzipped - forcing single thread extension.");
    args.th = 1;
  }

  if (args.th == 1)
    st_extend<extender_t>(&extender, args.patterns, sam_filename);
  else
    mt_extend<extender_t>(&extender, args.patterns, sam_filename, args.th, args.b);

  // TODO: Merge the SAM files.

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());
}

int main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

  if (args.shaped_slp)
  {
    dispatcher<extender<shaped_slp_t, ms_pointers<>>>(args);
  }
  else
  {
    dispatcher<extender<plain_slp_t, ms_pointers<>>>(args);
  }

  return 0;
}