/** \file commandLine.c
  * Command line parser tools.
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "log.h"
#include "type_vector.h"
#include "misc_cfgReader.h"

#define MC_MAX_KEY_LENGTH 200		///< Storage size for the 'key' string.

// ---------------------------------------------------------------------------
/// Key-value pair imported from command line using processing script.
// ---------------------------------------------------------------------------
typedef struct
{
  char key[MC_MAX_KEY_LENGTH];		///< 'key' string.
  int  type;				///< Type of the element stored in the value addressed memory.
  union
  {
    int     u_bool;
    int     u_int;
    double  u_double;
    int     u_vec3i[3];
    double  u_vec3d[3];
    char*   u_string;
  } value;				///< Union to keep it readable w/o typecasts.
} clPair_t;

static int        clPoolN = 0;		///< Number of pairs imported.
static clPair_t *clPool = NULL;	///< Elements themself.

/// Type aliases (used in clPair_t::type field).
enum {eBool = 0, eInt, eDouble, eVec3i, eVec3d, eString, eSequence};

// ---------------------------------------------------------------------------
/// Prints all 'key->value' pairs.
// ---------------------------------------------------------------------------
void
cl_dumpConfig (void)
{
  for (clPair_t* pair = clPool ; pair < clPool + clPoolN ; ++pair)
  {
    switch (pair->type)
    {
      case eBool:	SAY_DEBUG ("%s: %s", pair->key, (pair->value.u_bool) ? "true" : "false");	break;
      case eInt:	SAY_DEBUG ("%s: %d", pair->key, pair->value.u_int);				break;
      case eDouble:	SAY_DEBUG ("%s: %e", pair->key, pair->value.u_double);				break;
      case eString:
      case eSequence:	SAY_DEBUG ("%s: %s", pair->key, pair->value.u_string);				break;
      case eVec3i:
        SAY_DEBUG ("%s: (%d, %d, %d)", pair->key, pair->value.u_vec3i[0], pair->value.u_vec3i[1], pair->value.u_vec3i[2]);
      break;

      case eVec3d:
        SAY_DEBUG ("%s: (%e, %e, %e)", pair->key, pair->value.u_vec3d[0], pair->value.u_vec3d[1], pair->value.u_vec3d[2]);
      break;

      default:
        DIE ("unknown type '%d'", pair->type);
      break;
    }
  }
}

// ---------------------------------------------------------------------------
/// Imports command line options and resource file(s).
// ---------------------------------------------------------------------------
void
cl_import (const char *importPipe, int argc, char **argv, int skip)
{
  char *_ = (char*) calloc (30001, sizeof (char));									// Command line and 'magik' string.
  assert (_ && skip >= 0);

  strcpy (_, importPipe);												// Forms command line with all arguments.
  for (int i = 1 + skip ; i < argc ; ++i)
  {
    strncat (_, " \"", 30001 - strlen (_));
    strncat (_, argv[i], 30001 - strlen (_));
    strncat (_, "\"", 30001 - strlen (_));
    ENSURE (_[30000] == 0 && _[29999] == 0, "too long command line");
  }

  // Opens session with Perl parser and uses '_' as big input buffer.
  FILE *rc = popen (_, "r");												// Opens pipe to read parsed data.
  assert (rc);
  while (1)
  {
    /// Reads line of the code, strips new-line character and checks if there are no overflows.
    #define MF_GETLINE 										\
    assert (!feof (rc));									\
    fscanf (rc, "%30000[^\n]\n", _);								\
    assert (_[30000] == 0 && _[29999] == 0);

    MF_GETLINE;														// Reads token ('>>> DATA BEGIN' or 'BYE').
    if (!strcmp (_, "BYE"))
      break;

    ENSURE (!strcmp (_, ">>> DATA BEGIN"),
            "expected '>>> DATA BEGIN' instead of '%s'", _);

    clPool = (clPair_t *) realloc (clPool, (++clPoolN)*sizeof (clPair_t));						// Allocates new key.
    assert (clPool);

    clPair_t *pair = clPool + clPoolN - 1;										// Alias for data being imported.

    MF_GETLINE;														// Reads 'key' field.
    assert (strlen (_) < MC_MAX_KEY_LENGTH);
    strcpy (pair->key, _);

    MF_GETLINE;														// Reads 'type' field.
    pair->type = cfg_identifyWord (_, "bool", eBool, "int", eInt, "double", eDouble, 					// Recognizes type.
      "vec3i_t", eVec3i, "vec3d_t", eVec3d, "string", eString, "sequence", eSequence, mc_cfgTermGuesses);
    switch (pair->type)
    {
      case eBool:
      case eInt:
        MF_GETLINE;													// Reads one integer value.
        pair->value.u_int = atoi (_);
      break;

      case eDouble:
        MF_GETLINE;													// Reads one double value.
        pair->value.u_double = atof (_);
      break;

      case eVec3i:
        for (int i = 0 ; i < 3 ; ++i)											// Reads three integer values.
        {
          MF_GETLINE;
          pair->value.u_vec3i[i] = atoi (_);
        }
      break;

      case eVec3d:
        for (int i = 0 ; i < 3 ; ++i)											// Reads three double values.
        {
          MF_GETLINE;
          pair->value.u_vec3d[i] = atof (_);
        }
      break;

      case eString:
        MF_GETLINE;													// Reads string.
        pair->value.u_string = (char*) malloc (strlen (_) + 1);
        assert (pair->value.u_string);
        strcpy (pair->value.u_string, _);
      break;

      case eSequence:
      {
        MF_GETLINE;													// Reads length of packed sequence.
        int L = atoi (_), l = L + 2;											// +2: '\0' and '\n' after last element.
        assert (L > 0);

        pair->value.u_string = (char*) malloc (l*sizeof (char*));							// Allocates sequence container.
        assert (pair->value.u_string);

        MF_GETLINE;													// Reads number of elements in the sequence.
        int N = atoi (_);
        strcpy (pair->value.u_string, _);										// Copies 'number of elements' string.
        strcat (pair->value.u_string, "\n");
        l -= strlen (_) + 1;
        for (int i = 0 ; i < N ; ++i)
        {
          MF_GETLINE;													// Reads string.
          assert (strcmp (_, "<<< DATA END") && l > 0);									// Checks that terminator is not swallowed.
          strncat (pair->value.u_string, _, l);
          strcat (pair->value.u_string, "\n");
          l -= strlen (_) + 1;
        }
        assert (l == 1);												// Only room for NULL terminator is left.
        pair->value.u_string[L] = 0;											// Turns last '\n' into NULL terminator.
      }
      break;

      default:
        DIE ("unrecognized type '%s'", _);
      break;
    }

    MF_GETLINE;														// Checks that chunk is terminated.
    ENSURE (!strcmp (_, "<<< DATA END"),
            "expected '<<< DATA END' instead of '%s'", _);
    #undef MF_GETLINE
  }
  pclose (rc);

  cl_dumpConfig ();
}

// ---------------------------------------------------------------------------
/// Searches for the key-value pair.
// ---------------------------------------------------------------------------
static clPair_t *
cl_find (const char *key)
{
  for (clPair_t* pair = clPool ; pair < clPool + clPoolN ; ++pair)
    if (!strcmp (key, pair->key))
      return pair;
  return NULL;
}

/// Upper half of the common search routine (finds value for key and checks that type is correct.
#define MF_PRE(_type_, _Type_)								\
int											\
cl_find ## _Type_ (const char *key, _type_ *value)					\
{											\
  clPair_t *pair;									\
  if ((pair = cl_find (key)))								\
  {											\
    assert (pair->type == e ## _Type_);

/// Lower half of the common search routine (exports value through pointer in argumnents).
#define MF_POST										\
    return 1;										\
  }											\
  return 0;										\
}

/// Searches for the optional key-value pair.
MF_PRE(int, Int)
    *value = pair->value.u_int;
MF_POST

/// Searches for the optional key-value pair.
MF_PRE(double, Double)
    *value = pair->value.u_double;
MF_POST

/// Searches for the optional key-value pair.
MF_PRE(int, Vec3i)
    for (int i = 0 ; i < 3 ; ++i)
      value[i] = pair->value.u_vec3i[i];
MF_POST

/// Searches for the optional key-value pair.
MF_PRE(double, Vec3d)
    for (int i = 0 ; i < 3 ; ++i)
      value[i] = pair->value.u_vec3d[i];
MF_POST

/// Searches for the optional key-value pair.
MF_PRE(char*, String)
    *value = (char*) malloc (strlen (pair->value.u_string) + 1);
    assert (*value);
    strcpy (*value, pair->value.u_string);
MF_POST

/// Searches for the optional key-value pair.
MF_PRE(char*, Sequence)
    *value = (char*) malloc (strlen (pair->value.u_string) + 1);
    assert (*value);
    strcpy (*value, pair->value.u_string);
MF_POST
