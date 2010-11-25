/** \file parr_regLists.h
  * \brief Parallel syncronisation of the regList_t arrays over cluster (see parr_regLists.c, type_reg.c).
  */

#ifndef mc_parr_regLists_h
#define mc_parr_regLists_h					///< Guards against multiple include.

#include <stdio.h>
#include <stdlib.h>

#include "type_reg.h"

regList_t* syncReg_createArray (void);
regList_t* syncReg_sync (regList_t *local);
void       syncReg_addClaims (regList_t* list, const regList_t *claims);

#endif
