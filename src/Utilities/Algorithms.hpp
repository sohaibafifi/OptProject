// Part of OptProject: A basis project for RO problems
// By Sohaib Afifi <me [at] sohaibafifi.com>
// Copyright (c) 2012-2013 University of Technology of Compiegne, France


#pragma once
#include "../libs.hpp"
#include "mtrand/mtrand.h"
namespace Utilities{
class Algorithms
{
public:
    static string itos(int i);
    static string dtos(double i, unsigned precision = 3);
    static void debug(string message) ;

    /**
     * @brief strip white space from a string
     * @param s the input string
     * @param ws characters to be stripped, default are " \t\n\r"
     * @return the string stripped
     */
    static string stripWhiteSpace(const string & s, const string & ws = " \t\r\n");

    /**
     * @brief split a string with a delimiter character
     */
    static void splitString(string line, vector<string> & split, char sep=';');

    /**
     * @brief get euclidean norm of a vector
     */
    static double getVectorLength(vector<double> v);

    /**
     * @brief get average from a set of value
     */
    static double getAverage(vector<double> v);

    /**
     * @brief generate random order indices
     */
    static void generateRandomIndices(vector<unsigned int> & indices, const unsigned int & n);

    /**
     *@brief initialize random seed
     */
    static unsigned long initializeRandomSeed(unsigned long _seed = 0);
    static unsigned long seed;
    static MTRand_int32 irand;
    static MTRand drand;

    /**
     *@brief create n distinct random seeds
     */
    static void generateRandomSeed(vector<unsigned int> & seed, unsigned int n);

    /**
      * @briefget standard deviation from a set of value
      */
    static double getStandardDeviation(vector<double> v);

    /**
      * @brief get relative deviation (sdv/avg) from a set of value
      */
    static double getRelativeDeviation(vector<double> v);

    /**
      * @brief get the directory from file path
      */
    static std::string directory(string path);

    /**
     * @brief get the file name from the file path
    **/
    static std::string filename(string path);

    /**
     * @brief check if a file exists
     */
    static bool fileExists(string filepath);

    };
} // end  namespace
// define a shortcut
typedef Utilities::Algorithms __;
// some useful macros
#define lambda(return_type, ...) \
  __extension__ \
  ({ \
    return_type __fn__ __VA_ARGS__ \
    __fn__; \
  })
#define QUOTE(name) #name
#define STR(name) QUOTE(name)
#define ISNAN(x) ((x) != (x))

#ifdef DEBUG
#define dout cout
#else
#define dout 0 && cout
#endif
