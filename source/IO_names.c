/** \file IO_names.c
 * Wrapper for \b sprintf used by IO routines to convert number of the record and/or cpu number into filename.
 */

#include <stdio.h>

#include "log.h"

// XXX - remove quit, replace all wrappers with '_'.

// ---------------------------------------------------------------------------
/// Engraves cpu number into the format string and returns the result.
// ---------------------------------------------------------------------------
const char*
IO_plasmaName (int record, int cpu)
{
  static char name[300];
  sprintf (name, "binData/plasma_%06d_%03d.bin", record, cpu);
  return name;
}

// ---------------------------------------------------------------------------
/// Engraves cpu number into the format string and returns the result.
// ---------------------------------------------------------------------------
const char *
IO_nameCpu (const char *format, int nodeNum)
{
  static char name[100];
  ENSURE (nodeNum >= 0 && nodeNum < 1000,
          "node number '%d' is longer than 3 digits", nodeNum);
  sprintf (name, format, nodeNum);
  return name;
}

// ---------------------------------------------------------------------------
/// Engraves cpu number and record number into the format string and returns the result.
// ---------------------------------------------------------------------------
const char *
IO_nameCpuRec (const char *format, int cpu_num, int recordNum)
{
  static char name[80];
  if (cpu_num < 0 || recordNum < 0 || cpu_num > 999 || recordNum > 999999)
    DIE ("bad parameters for format string cpuNum (%d > 3 digits) or recordNum (%d > 6 digits)", cpu_num, recordNum);
  sprintf (name, format, recordNum, cpu_num);
  return name;
}

// ---------------------------------------------------------------------------
/// Engraves record number into the format string and returns the result.
// ---------------------------------------------------------------------------
const char *
IO_nameRec (const char *format, int recordNum)
{
  static char name[80];
  ENSURE (recordNum >= 0, "bad recordNum (%d)", recordNum);
  sprintf (name, format, recordNum);
  return name;
}
