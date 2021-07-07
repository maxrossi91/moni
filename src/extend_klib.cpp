/* extend_klib - Extend the MEMs of the reads to the reference
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
   \file extend_klib.cpp
   \brief extend_klib.cpp Extend the MEMs of the reads to the reference.
   \author Massimiliano Rossi
   \date 30/04/2021
*/

extern "C" {
#include <xerrors.h>
}

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/io.hpp>

#include <ms_pointers.hpp>
#include <extender_klib.hpp>
#include <extend_reads_dispatcher.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <PlainSlp.hpp>
#include <FixedBitLenCode.hpp>

#include <ksw.h>
#include <ssw.h>

#include <libgen.h>

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  std::string patterns = ""; // path to patterns file
  size_t l = 25;             // minumum MEM length
  size_t th = 1;             // number of threads
  size_t b = 1;              // number of batches per thread pool
  bool shaped_slp = false;   // use shaped slp
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-p patterns] [-t threads] [-l len] [-q shaped_slp] [-b batch]\n\n" +
                    "Extends the MEMs of the reads in the pattern against the reference index in infile.\n" +
                    "shaped_slp: [boolean] - use shaped slp. (def. false)\n" +
                    "   pattens: [string]  - path to patterns file.\n" +
                    "       len: [integer] - minimum MEM lengt (def. 25)\n" +
                    "    thread: [integer] - number of threads (def. 1)\n" +
                    "     batch: [integer] - number of batches per therad pool (def. 1)\n");

  std::string sarg;
  while ((c = getopt(argc, argv, "l:hp:b:t:")) != -1)
  {
    switch (c)
    {
    case 'p':
      arg.patterns.assign(optarg);
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
void dispatcher(Args &args){
  verbose("Construction of the extender");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  extender_t extender(args.filename, args.l);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::string base_name = basename(args.filename.data());
  std::string sam_filename = args.patterns + "_" + base_name + "_" + std::to_string(args.l);

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
    dispatcher<extender<shaped_slp_t>>(args);
  }
  else
  {
    dispatcher<extender<plain_slp_t>>(args);
  }

  return 0;
}