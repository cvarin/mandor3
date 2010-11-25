/** \file misc_cfgReader.h
  * \brief Reads parameters from ascii configuration files (one parameter in a row plus comment).
  * Enhances functionality of the usual printf by the skipping of the rest of line (usually some
  * comments sit there) untill next signal character \b @ will appear at the start of the line.
  * Initially was done for setup.out, now is used everywhere. It alse contains routines for safe
  * opening of file (fopen + check against NULL and error messaging).
  *
  * \todo: Use lexx and yacc for this stuff - much more powerfull and also code will be kept clean.
  */

#ifndef cfgReader_header
#define cfgReader_header

#include <stdio.h>

void   cfg_nextLine (FILE *fp);
double cfg_readDouble (FILE *fp, const char *caller, int line);
int    cfg_readInt (FILE *fp, const char *caller, int line);
double cfg_readOptDouble (FILE *fp, const char *caller, int line);
int    cfg_readOptInt (FILE *fp, const char *caller, int line);
int    cfg_isParameter (FILE *fp);
int    cfg_isOption (FILE *fp);
char*  cfg_readLine (FILE *fp);
const char*  cfg_readWord (FILE *fp, const char *caller, int line);
const char*  cfg_readOptWord (FILE *fp, const char *caller, int line);
int    cfg_identifyWord (const char *word, ... );

FILE*  cfg_open (const char *name, const char *attr, const char *parentFunc);

/// Terminator of the list for cfg_identifyWord() routine.
#define mc_cfgTermGuesses		(char*) NULL, -1
/*
 * Here I use ISO C99 compiler standart - name of the function is a predefined
 * macro __func__ of the type 'static const char *', name of the file - __FILE__, etc.
 */
#define cfg_readDouble(file)    cfg_readDouble(file, __FILE__, __LINE__)		///< Adds two default arguments to the call.
#define cfg_readInt(file)       cfg_readInt(file, __FILE__, __LINE__)			///< Adds two default arguments to the call.
#define cfg_readOptDouble(file) cfg_readOptDouble(file, __FILE__, __LINE__)	///< Adds two default arguments to the call.
#define cfg_readOptInt(file)    cfg_readOptInt(file, __FILE__, __LINE__)		///< Adds two default arguments to the call.
#define cfg_readWord(file)      cfg_readWord(file, __FILE__, __LINE__)		///< Adds two default arguments to the call.
#define cfg_readOptWord(file)   cfg_readOptWord(file, __FILE__, __LINE__)		///< Adds two default arguments to the call.

#endif
