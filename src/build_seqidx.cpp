/* build_seqidx - Builds the sequence index for the reference
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
   \file build_seqidx.cpp
   \brief build_seqidx.cpp Builds the sequence index for the reference.
   \author Massimiliano Rossi
   \date 07/08/2021
*/

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <malloc_count.h>
#include <kseq.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread);

#include <seqidx.hpp>

#include <filesystem>
namespace fs = std::filesystem;

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  std::string outpath = ""; // path where to output the file
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-o outpath]\n\n" +
                    "Computes the .idx file storing the sequence names and starting positions.\n" +
                    "outpath: [string]  - path to where to output the file.\n");

  std::string sarg;
  while ((c = getopt(argc, argv, "o:")) != -1)
  {
    switch (c)
    {
    case 'o':
      arg.outpath.assign(optarg);
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

int main(int argc, char *const argv[])
{
  Args args;
  parseArgs(argc, argv, args);

  // Building the sequence idx

  verbose("Building the sequence index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  seqidx idx(args.filename);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Sequence index construction complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());


  std::string outfile = "";
  if(args.outpath == "") outfile = args.filename;
  else outfile = args.outpath + fs::path(args.filename).filename().string();
  outfile += idx.get_file_extension();

  std::ofstream out(outfile);
  idx.serialize(out);

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Sequence index serialzation complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());
  return 0;
}