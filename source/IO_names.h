/** \file IO_names.h
 * Wrapper for \b sprintf used by IO routines to convert number of the record and/or cpu number into filename.
 */

#ifndef MC_IO_NAMES_HEADER
#define MC_IO_NAMES_HEADER				///< \internal Guard.

const char* IO_plasmaName (int record, int cpu);

const char * IO_nameCpu (const char *format, int nodeNum);
const char * IO_nameCpuRec (const char *format, int cpu_num, int recordNum);
const char * IO_nameRec (const char *format, int recordNum);

static const char sysName_info[50] = "binData/sys_%06d.txt";			// Format string for global parameters of the record files' manifold.
static const char sysName_prtn_full[50] = "binData/sys_prtn_%06d.bin";		// Format string for name of partition file (recNum).

static const char sysName_E_full[50] = "binData/sys_E_%06d(%03d).bin";		// Format string for name of electric field file (recNum, node).
static const char sysName_H_full[50] = "binData/sys_H_%06d(%03d).bin";		// Format string for name of magnetic field file (recNum, node).
static const char sysName_DF_full[50] = "binData/sys_DF_%06d(%03d).bin";	// Format string for name of distribution function file (recNum, node).

static const char sysName_E_map[50] = "binData/sys_E_%06d(%%03d).bin";		// Format string for generator of name of electric field file (recNum, node).
static const char sysName_H_map[50] = "binData/sys_H_%06d(%%03d).bin";		// Format string for generator of name of magnetic field file (recNum, node).
static const char sysName_DF_map[50] = "binData/sys_DF_%06d(%%03d).bin";	// Format string for generator of name of distribution function file (recNum, node).


#endif
