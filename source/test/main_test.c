/** \file main_test.c
  * \brief Test of the system reg-exp engine used to parse options.
  */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "commandLine.h"

int
main (int argc, char *argv[])
{
  cl_import ("./source/import.pl", argc, argv, 0);
  cl_dumpConfig ();

  return 0;
}
