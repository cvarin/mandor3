/** \file commandLine.h
  * Command line analysis tool: reads command line arguments and helps to extract options given.
  *
  * Options may have long name started with '--' and short name started with '-'. Arguments of the
  * option (if there are any) is separated by ':'.
  */

#ifndef MC_COMMANDLINE_HEADER
#define MC_COMMANDLINE_HEADER				///< \internal Multiple include flag.

void cl_import (const char *importPipe, int argc, char **argv, int skip);
void cl_dumpConfig (void);

int cl_findInt (char *key, int *value);
int cl_findDouble (char *key, double *value);
int cl_findVec3i (char *key, int *value);
int cl_findVec3d (char *key, double *value);
int cl_findString (char *key, char **value);
int cl_findSequence (char *key, char **value);

#endif
