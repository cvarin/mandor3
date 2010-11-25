/** \file IO_fileMap.h
  * \brief Mapping of the nodes on files in the ./binData (see IO_fileMap.c).
  */

#ifndef IO_fileMap_header
#define IO_fileMap_header

#include "misc_partition.h"

/**
  * Parameters of the file-map grouped together.
  */
typedef struct
{
  regList_t  map;				///< Partitioning of the data sub-domains.
  int       *visMatrix;				///< Visibility matrix for all nodes for all pieces (main goal - clusters with local filesystems).
  int        cpu_total;				///< Used as another dimension for visibility matrix.
} fileMap_t;

/// Null file-map (for default initialization).
#define mc_fileMap_init { mc_regList_init, 0, 0 }

void fileMap_buildVisMatrix (fileMap_t *map, const char *format);
int  fileMap_readable (const fileMap_t *map, int node, int reader);

int  fileMap_save (void);
void fileMap_load (fileMap_t *map, int id);
void fileMap_free (fileMap_t *map);
void fileMap_init (int mode);

#define mc_fileMap_setup	529546		/**< Start-up initialization code for \b setup (record number starts from 0).			*/
#define mc_fileMap_start	163756		/**< Start-up initialization code for \b simulation-start (record number starts from 1).	*/
#define mc_fileMap_contn	359416		/**< Start-up initialization code for \b simulation-continue (from the last one).		*/

#endif
