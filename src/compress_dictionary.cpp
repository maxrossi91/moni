/* compress_dictionary - Computes the compressed dictionary from prefix-free parse dictionary
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
   \file compress_dictionary.cpp
   \brief compress_dictionary.cpp Computes the compressed dictionary from prefix-free parse dictionary.
   \author Massimiliano Rossi
   \date 16/09/2020
*/

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <malloc_count.h>

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  size_t w = 10;             // sliding window size and its default
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-w wsize]\n\n" +
                    "Computes the pfp data structures of infile, provided that infile.parse, infile.dict, and infile.occ exists.\n" +
                    "  wsize: [integer] - sliding window size (def. 10)\n");

  std::string sarg;
  while ((c = getopt(argc, argv, "w:smcfl:rhp:t:")) != -1)
  {
    switch (c)
    {
    case 'w':
      sarg.assign(optarg);
      arg.w = stoi(sarg);
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

std::string execute_cmd(const char* cmd) {
    std::array<char, 256> buffer{};
    std::string output = "";

    std::string cmd_plus_stderr = std::string(cmd) + " 2>&1";
    FILE* pipe = popen(cmd_plus_stderr.data(), "r"); // Extract stderr as well
    if (!pipe) {throw std::runtime_error("popen() failed!");}

    try {
        std::size_t bytes;
        while ((bytes = fread(buffer.data(), sizeof(char), sizeof(buffer), pipe))) {
            output += std::string(buffer.data(), bytes);
        }
    } catch (...) {
        pclose(pipe);
        throw std::runtime_error("Error occurred while reading popen() stream.");
    }
    pclose(pipe);
    return output;
}

int main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

  // Building the r-index

  verbose("Compressing the dictionary");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  // Open output files
  std::string dicz_filename = args.filename + ".dicz";
  std::string dicz_len_filename = args.filename + ".dicz.len";

  FILE *dicz;
  FILE *dicz_len;

  if ((dicz = fopen(dicz_filename.c_str(), "w")) == nullptr)
    error("open() file " + std::string(dicz_filename) + " failed");

  if ((dicz_len = fopen(dicz_len_filename.c_str(), "w")) == nullptr)
    error("open() file " + std::string(dicz_len_filename) + " failed");

  // Open the dictionary
  std::string dict_filename = args.filename + ".dict";
  std::vector<uint8_t> dict;
  read_file(dict_filename.c_str(), dict);

  // Start processing

  
  // Generating phrase lengths
  verbose("Generating phrase lengths");
  std::vector<size_t> lengths(1,0);
  
  // Counting the number of Dollars at the beginning
  size_t i = 0, j = 0;
  while(dict[i++] == Dollar)
    j++;
  dict.erase(dict.begin(), dict.begin() + j);

  for(auto chr: dict)
  {
    // Skip the Dollars
    if(chr == EndOfDict)
      continue;

    // Hit end of phrase
    if(chr == EndOfWord)
      lengths.push_back(0);
    else
      lengths.back()++;
  }

  if (lengths.back()==0)
    lengths.pop_back();

  verbose("Found", lengths.size(), " phrases ");

  verbose("Generating phrases");
  uint8_t* ptr = dict.data(); // Beginning of the current phrase

  bool empty_first_phrase = false;
  for(size_t i = 0; i < lengths.size(); i++)
  {
    size_t compressed_length = lengths[i] - args.w;

    // special case: starts with a trigger string
    if (i==0 && compressed_length == 0) {
      ptr += lengths[i] + 1; 
      empty_first_phrase = true;
      continue;
    } else if (i > 0 && compressed_length == 0) {
      error("encountered a length=0 phrase after removing trigger string, which should not occur.");
    }

    if ((fwrite(&compressed_length, 4, 1, dicz_len)) != 1)
      error("fwrite() file " + std::string(dicz_len_filename) + " failed");

    if ((fwrite(ptr, sizeof(uint8_t), compressed_length, dicz)) != compressed_length)
      error("fwrite() file " + std::string(dicz_filename) + " failed");

    ptr += lengths[i] + 1;
  }
  fclose(dicz);
  fclose(dicz_len);

  // re-writes parse file to shift down all the phrase ids by 1 
  // since we removed the empty beginning phrase
  if (empty_first_phrase) {
    verbose("alert: found that the first phrase length is 0"
            " so we will rewrite *.parse file to generated correct SLP.");

    // read in all the phrase ids in parse
    std::string parse_filename = args.filename + ".parse";
    std::vector<uint32_t> parse_arr;
    read_file(parse_filename.c_str(), parse_arr);

    // make sure first phrase is lowest lexicographically and then remove it
    if (parse_arr[0] != 1)
      error("parse should being with lowest lexicographic phrase.");
    parse_arr.erase(parse_arr.begin());

    // rename the old parse file as *.parse_with_empty_phrase
    std::ostringstream command_stream;
    command_stream << "mv " << parse_filename << " " << (args.filename + ".parse_with_empty_phrase");
    auto mv_log = execute_cmd(command_stream.str().c_str());

    verbose("executed this command: " + command_stream.str());

    // open new parse file for writing
    FILE* new_parse_file;
    if ((new_parse_file = fopen((args.filename + ".parse").c_str(), "w")) == nullptr)
      verbose("open() file " + std::string(args.filename + ".parse" + " failed"));

    // iterate through each element of parse and decrement by 1
    for (size_t i = 0; i < parse_arr.size(); i++) {
      if (parse_arr[i] == 1)
        error("issue occurred when creating new parse file.");
      parse_arr[i]--;
      
      // write it out 
      if ((fwrite(&parse_arr[i], 4, 1, new_parse_file)) != 1)
        verbose("fwrite() file " + std::string(args.filename + ".parse") + " failed"); 
    }
    fclose(new_parse_file);
  }







  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  return 0;
}