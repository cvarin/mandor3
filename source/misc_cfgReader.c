/** \file misc_cfgReader.c
  * \brief Reader of the parameters from ascii configuration files (one parameter in a row plus comment).
  */

#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "log.h"

#define cfgReader_header	///< Blocks the including of own header and corresponding define extension.

#define mc_bufferLineSize	(1024)		///< Size of the line buffer.

static char cfgLine[mc_bufferLineSize];		///< Local line buffer.

// ---------------------------------------------------------------------------
/// Checks if the first symbol to read is \b @.
// ---------------------------------------------------------------------------
int
cfg_isParameter (FILE *fp)
{
  int c = fgetc (fp);
  if (feof (fp))
    return 0;
  ungetc (c, fp);
  return (c == '@');
}

// ---------------------------------------------------------------------------
/// Checks if the first symbol to read is \b >.
// ---------------------------------------------------------------------------
int
cfg_isOption (FILE *fp)
{
  int c = fgetc (fp);
  if (feof (fp))
    return 0;
  ungetc (c, fp);
  return (c == '>');
}

// ---------------------------------------------------------------------------
/// Skips the rest of the line.
// ---------------------------------------------------------------------------
void
cfg_nextLine (FILE *fp)
{
  cfgLine[mc_bufferLineSize-2] = 0;						// Sets 'too big line' detector.
  fgets (cfgLine, mc_bufferLineSize, fp);
  ENSURE (cfgLine[mc_bufferLineSize-2] == 0,
          "lines longer that %d characters in input file - increase buffer",
          mc_bufferLineSize - 2);
}

// ---------------------------------------------------------------------------
/// \brief Reads double. If signal mark is missed then caller is identified using \b caller and \b line parameters (added automatically by macros).
// ---------------------------------------------------------------------------
double
cfg_readDouble (FILE *fp, const char *caller, int line)
{
  double tmp;
  int loaded_items = fscanf (fp, "@ %le %*[^\n] ", &tmp);
  ENSURE(loaded_items == 1, "%s/%d cannot get double number", caller, line);
  return tmp;
}

// ---------------------------------------------------------------------------
/// \brief Reads int. If signal mark is missed then caller is identified using \b caller and \b line coordinates, added automatically by macros.
// ---------------------------------------------------------------------------
int
cfg_readInt (FILE *fp, const char *caller, int line)
{
  int tmp;
  // XXX in new lib 'fscanf ("@ %d",) == 1' works better!
  ENSURE (fgetc (fp) == '@',
          "missing required int parameter (%s/line %d)", caller, line);
  fscanf (fp, "%d", &tmp);
  cfg_nextLine (fp);
  return tmp;
}

// ---------------------------------------------------------------------------
/// \brief Reads double as option. If signal mark is missed then caller is
/// identified using \b caller and \b line parameters (added automatically by
/// macros).
// ---------------------------------------------------------------------------
double
cfg_readOptDouble (FILE *fp, const char *caller, int line)
{
  double tmp;
  ENSURE (fgetc (fp) == '>',
          "missing optional double parameter (%s/line %d)", caller, line);
  fscanf (fp, "%le", &tmp);
  cfg_nextLine (fp);
  return tmp;
}

// ---------------------------------------------------------------------------
/// \brief Reads int as option. If signal mark is missed then caller is identified using \b caller and \b line coordinates, added automatically by macros.
// ---------------------------------------------------------------------------
int
cfg_readOptInt (FILE *fp, const char *caller, int line)
{
  int tmp;
  ENSURE (fgetc (fp) == '>',
          "missing optional int parameter (%s/line %d)", caller, line);
  fscanf (fp, "%d", &tmp);
  cfg_nextLine (fp);
  return tmp;
}

// ---------------------------------------------------------------------------
/// Opens file and checks that operation is done successfully.
// ---------------------------------------------------------------------------
FILE*
cfg_open (const char *name, const char *attr, const char *parentFunc)
{
  FILE *fp = fopen (name, attr);
  ENSURE (fp, "cannot open file '%s' with attributes '%s' in '%s'",
                                                      name, attr, parentFunc);
  return fp;
}

// ---------------------------------------------------------------------------
/// Reads line from the file and returns pointer to the line.
// ---------------------------------------------------------------------------
char*
cfg_readLine (FILE *fp)
{
  cfg_nextLine (fp);		// Gets line to the buffer.
  return cfgLine;
}

// ---------------------------------------------------------------------------
/// \brief Reads line from config file, analyzes content and returns pointer to
/// the first word. If signal mark is missed then caller is identified using
/// \b caller and \b line coordinates, added automatically by macros.
// ---------------------------------------------------------------------------
const char*
cfg_readWord (FILE *fp, const char *caller, int line)
{
    // Gets line to the buffer.
    char *word = cfg_readLine (fp);

    ENSURE (*word == '@', "parameter is abscent [%s/line %d]", caller, line);

    while (*word && (*word == ' ' || *word == '@' || *word == '>'))
        ++word;

    // Terminates string after the first word.
    char *end = word;
    while (*end > ' ')
        ++end;
    *end = 0;

    static char words[8][200];
    static int  slot = 0;
    slot = (slot + 1) & 7;
    strcpy (words[slot], word);
    return words[slot];
}

// ---------------------------------------------------------------------------
/// \brief Reads line from config file, analyzes content and returns pointer to the first word. If signal mark is missed then
/// caller is identified using \b caller and \b line coordinates, added automatically by macros.
// ---------------------------------------------------------------------------
const char *
cfg_readOptWord (FILE *fp, const char *caller, int line)
{
  char *word = cfg_readLine (fp);						// Gets line to the buffer.

  ENSURE (*word == '>', "parameter is abscent [%s/line %d]", caller, line);

  while (*word && (*word == ' ' || *word == '@' || *word == '>'))
    ++word;

  char *end = word;								// Terminates string after the first word.
  while (*end > ' ')
    ++end;
  *end = 0;

  return word;
}

// ---------------------------------------------------------------------------
/// \brief Gets word and list of choises (in form <b>"variant word", variantID</b>) ended by
/// <b>, NULL, -1</b> (use macros \b mc_cfgTermGuesses). Returns ID for the guess if one of the
/// variants is correct or -1.
// ---------------------------------------------------------------------------
int
cfg_identifyWord (const char *word, ... )
{
  va_list  argptr;								// Starts argument list dispatching.
  va_start (argptr, word);

  int retValue = -1;
  while (1)
  {
    char *guess = va_arg(argptr, char*);
    int   guessID = va_arg(argptr, int);

    if (!guess && guessID == -1)						// Terminator is detected.
      break;

    ENSURE (guessID != -1,
            "ID '-1' (suggested for '%s') is reserved for 'not found'", guess);

    if (! strcmp (word, guess))
    {
      retValue = guessID;
      break;
    }
  }
  va_end (argptr);								// Ends argument list dispatching.

  return retValue;
}
