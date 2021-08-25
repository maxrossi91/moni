/* pfp-ds - prefix free parsing data structures
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
   \file common.hpp
   \brief common.hpp contains common features.
   \author Massimiliano Rossi
   \date 12/03/2020
*/

#ifndef _COMMON_HH
#define _COMMON_HH

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <assert.h>

#include <sys/time.h>

#include <sys/mman.h> // for mmap
#include <unistd.h>
#include <sys/stat.h>
 #include <fcntl.h>

#include <sstream>      // std::stringstream

#include <vector>      // std::vector

#include <chrono>       // high_resolution_clock

#include <sdsl/io.hpp>  // serialize and load
#include <type_traits>  // enable_if_t and is_fundamental

//**************************** From  Big-BWT ***********************************
// special symbols used by the construction algorithm:
//   they cannot appear in the input file
//   the 0 symbol is used in the final BWT file as the EOF char

#define Dollar 2     // special char for the parsing algorithm, must be the highest special char
#define EndOfWord 1  // word delimiter for the plain dictionary file
#define EndOfDict 0  // end of dictionary delimiter
//******************************************************************************

#define THRBYTES 5 // The number of bytes for the thresholds
#define SSABYTES 5 // The number of bytes for the thresholds

std::string NowTime();
void _internal_messageInfo(const std::string message);
void _internal_messageWarning( const std::string file, const unsigned int line, const std::string message);
void _internal_messageError( const std::string file, const unsigned int line,const std::string message);


std::string NowTime()
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    char buffer[100];
    tm r;
    strftime(buffer, sizeof(buffer), "%X", localtime_r(&tv.tv_sec, &r));
    char result[100];
    snprintf(result, 100, "%s"/*.%06ld"*/, buffer/*, (long)tv.tv_usec*/);
    return result;
}


template<typename T>
inline void _internal_message_helper(std::stringstream &ss, T const &first) { ss << first; }
template<typename T, typename... Args>
inline void _internal_message_helper(std::stringstream &ss, T const &first, const Args&... args) { ss << first << " "; _internal_message_helper(ss,args...); }
template<typename T, typename... Args>
inline std::string _internal_message(T const &first, const Args&... args) { std::stringstream ss; _internal_message_helper(ss,first,args...); return ss.str(); }


void _internal_messageInfo(const std::string message)
{
  std::cout << "[INFO] " << NowTime() << " - " << "Message: " << message << std::endl;
}

void _internal_messageWarning( const std::string file, const unsigned int line,
  const std::string message)
{
  std::cout << "[WARNING] " << NowTime() << " - "
  << "File: " << file << '\n'
  << "Line: " << line << '\n'
  << "Message: " << message << std::endl;
}

void _internal_messageError( const std::string file, const unsigned int line,
  const std::string message)
{
  std::cerr << "[ERROR] " << NowTime() << " - "
  << "File: " << file << '\n'
  << "Line: " << line << '\n'
  << "Message: " << message << std::endl;
  assert( false );
  exit( 1 );
}



#define info( args... ) \
    _internal_messageInfo( _internal_message(args) )

#ifdef VERBOSE
  #define verbose( args... ) \
      _internal_messageInfo( _internal_message(args) )
#else
  #define verbose( args... )
#endif

#define warning( args... ) \
    _internal_messageWarning( __FILE__, __LINE__, _internal_message(args) )

#define error( args... ) \
    _internal_messageError( __FILE__, __LINE__, _internal_message(args) )


// converts elemens in csv format
template <typename T>
inline void csv_helper(std::stringstream &ss, T const &first){ss << first;}
template <typename T, typename... Args>
inline void csv_helper(std::stringstream &ss, T const &first, const Args &... args){ ss << first << ", "; csv_helper(ss, args...);}
template <typename T, typename... Args>
inline std::string csv(T const &first, const Args &... args){std::stringstream ss;csv_helper(ss, first, args...); return ss.str();}

//*********************** File I/O *********************************************
template<typename T>
void map_file(const char *filename, T*& ptr, size_t& length){
    struct stat filestat;
    int fd;

    if ((fd = open(filename, O_RDONLY)) < 0)
        error("open() file " + std::string(filename) + " failed" );

    if (fstat(fd, &filestat) < 0)
        error("stat() file " + std::string(filename) + " failed" );

    if(filestat.st_size % sizeof(T) != 0)
        error("invilid file " + std::string(filename));

    length = filestat.st_size / sizeof(T);

    if ((ptr = mmap(NULL, filestat.st_size, PROT_READ, MAP_SHARED, fd, 0)) == MAP_FAILED)
        error("mmap() file " + std::string(filename) + " failed");
}

template<typename T>
void read_file(const char *filename, T*& ptr, size_t& length){
    struct stat filestat;
    FILE* fd;

    if ((fd = fopen(filename, "r")) == nullptr)
        error("open() file " + std::string(filename) + " failed" );

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        error("stat() file " + std::string(filename) + " failed" );

    if(filestat.st_size % sizeof(T) != 0)
        error("invilid file " + std::string(filename));

    length = filestat.st_size / sizeof(T);
    ptr = new T[length];

    if ((fread(ptr, sizeof(T), length, fd)) != length)
        error("fread() file " + std::string(filename) + " failed");

    fclose(fd);
}

template<typename T>
void read_file(const char *filename, std::vector<T>& ptr){
    struct stat filestat;
    FILE* fd;

    if ((fd = fopen(filename, "r")) == nullptr)
        error("open() file " + std::string(filename) + " failed" );

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        error("stat() file " + std::string(filename) + " failed" );

    if(filestat.st_size % sizeof(T) != 0)
        error("invilid file " + std::string(filename));

    size_t length = filestat.st_size / sizeof(T);
    ptr.resize(length);

    if ((fread(&ptr[0], sizeof(T), length, fd)) != length)
        error("fread() file " + std::string(filename) + " failed");

    fclose(fd);
}

void read_file(const char *filename, std::string &ptr)
{
  struct stat filestat;
  FILE *fd;

  if ((fd = fopen(filename, "r")) == nullptr)
    error("open() file " + std::string(filename) + " failed");

  int fn = fileno(fd);
  if (fstat(fn, &filestat) < 0)
    error("stat() file " + std::string(filename) + " failed");

  if (filestat.st_size % sizeof(char) != 0)
    error("invilid file " + std::string(filename));

  size_t length = filestat.st_size / sizeof(char);
  ptr.resize(length);

  if ((fread(&ptr[0], sizeof(char), length, fd)) != length)
    error("fread() file " + std::string(filename) + " failed");

  fclose(fd);
}

template<typename T>
void read_fasta_file(const char *filename, std::vector<T>& v){
    FILE* fd;

    if ((fd = fopen(filename, "r")) == nullptr)
        error("open() file " + std::string(filename) + " failed" );

    v.clear();

    char c;
    while (fread( &c, sizeof(char), 1,fd) == 1) {
      if(c == '>'){
        while(fread( &c, sizeof(char), 1,fd) == 1 && c != '\n');
      }else{
        v.push_back(c);
        while(fread( &c, sizeof(char), 1,fd) == 1 && c!= '\n') v.push_back(c);
      }
  	}
  	fclose(fd);
}

template <typename T>
void write_file(const char *filename, std::vector<T> &ptr)
{
  struct stat filestat;
  FILE *fd;

  if ((fd = fopen(filename, "w")) == nullptr)
    error("open() file " + std::string(filename) + " failed");

  size_t length = ptr.size(); 
  if ((fwrite(&ptr[0], sizeof(T), length, fd)) != length)
    error("fwrite() file " + std::string(filename) + " failed");

  fclose(fd);
}

//*********************** Time resources ***************************************

/*!
 * op the operation that we want to measure
 */
#define _elapsed_time(op)                                                                                               \
  ({                                                                                                                    \
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();          \
    op;                                                                                                                 \
    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();            \
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count()); \
    std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count();                                \
  })



//********** begin my serialize edit from sdsl ********************
// Those are wrapper around most of the serialization functions of sdsl


template <class T, typename size_type>
uint64_t
my_serialize_array(const T* p, const size_type size, std::ostream &out, typename std::enable_if<std::is_fundamental<T>::value>::type * = 0)
{
  size_t written_bytes = 0;
  if (size > 0)
  {

    size_type idx = 0;
    while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (size))
    {
      out.write((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
      written_bytes += sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T);
      p += sdsl::conf::SDSL_BLOCK_SIZE;
      idx += sdsl::conf::SDSL_BLOCK_SIZE;
    }
    out.write((char *)p, ((size) - idx) * sizeof(T));
    written_bytes += ((size) - idx) * sizeof(T);

  }
  return written_bytes;
}

//! Serialize each element of an std::vector
/*!
 * \param vec The vector which should be serialized.
 * \param out Output stream to which should be written.
 * \param v   Structure tree node. Note: If all elements have the same
 *            structure, then it is tried to combine all elements (i.e.
 *            make one node w with size set to the cumulative sum of all
 *           sizes of the children)
 */
// specialization for fundamental types
template <class T>
uint64_t
my_serialize_vector(const std::vector<T> &vec, std::ostream &out, sdsl::structure_tree_node *v, std::string name, typename std::enable_if<std::is_fundamental<T>::value>::type * = 0)
{
  if (vec.size() > 0)
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, "std::vector<" + sdsl::util::class_name(vec[0]) + ">");
    size_t written_bytes = 0;

    const T *p = &vec[0];
    typename std::vector<T>::size_type idx = 0;
    while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (vec.size()))
    {
      out.write((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
      written_bytes += sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T);
      p += sdsl::conf::SDSL_BLOCK_SIZE;
      idx += sdsl::conf::SDSL_BLOCK_SIZE;
    }
    out.write((char *)p, ((vec.size()) - idx) * sizeof(T));
    written_bytes += ((vec.size()) - idx) * sizeof(T);

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }
  else
  {
    return 0;
  }
}

template <typename X>
uint64_t
my_serialize(const std::vector<X> &x,
             std::ostream &out, sdsl::structure_tree_node *v = nullptr,
             std::string name = "", typename std::enable_if<std::is_fundamental<X>::value>::type * = 0)
{
  return sdsl::serialize(x.size(), out, v, name) + my_serialize_vector(x, out, v, name);
}


/**
 * @brief Load an array of size elements into p. p should be preallocated.
 * 
 * \tparam T 
 * \tparam size_type 
 * @param p 
 * @param size 
 * @param in 
 */
template <class T, typename size_type>
void my_load_array(T *p, const size_type size, std::istream &in, typename std::enable_if<std::is_fundamental<T>::value>::type * = 0)
{
  size_type idx = 0;
  while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (size))
  {
    in.read((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
    p += sdsl::conf::SDSL_BLOCK_SIZE;
    idx += sdsl::conf::SDSL_BLOCK_SIZE;
  }
  in.read((char *)p, ((size) - idx) * sizeof(T));
}

//! Load all elements of a vector from a input stream
/*! \param vec  Vector whose elements should be loaded.
 *  \param in   Input stream.
 *  \par Note
 *   The vector has to be resized prior the loading
 *   of its elements.
 */
template <class T>
void my_load_vector(std::vector<T> &vec, std::istream &in, typename std::enable_if<std::is_fundamental<T>::value>::type * = 0)
{
  T *p = &vec[0];
  typename std::vector<T>::size_type idx = 0;
  while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (vec.size()))
  {
    in.read((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
    p += sdsl::conf::SDSL_BLOCK_SIZE;
    idx += sdsl::conf::SDSL_BLOCK_SIZE;
  }
  in.read((char *)p, ((vec.size()) - idx) * sizeof(T));
}

template <typename X>
void my_load(std::vector<X> &x, std::istream &in, typename std::enable_if<std::is_fundamental<X>::value>::type * = 0)
{
  typename std::vector<X>::size_type size;
  sdsl::load(size, in);
  x.resize(size);
  my_load_vector(x, in);
}




#endif /* end of include guard: _COMMON_HH */
