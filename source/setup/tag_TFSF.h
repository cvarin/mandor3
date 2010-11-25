/** \file tag_TFSF.h
  * Header and explanation for tag_TFSF.c.
  *
  * To see what <b>T</b>otal <b>F</b>ield / <b>S</b>cattered <b>F</b>ield" interface is please consult documentation for main calculation
  * module (see em_TFSF.c file documentation). Speaking shortly it is domain limited by 6 faces. Initially there are no field in the domain so
  * any non-zero solution is a wave coming through the boundaries. With some manipulation this solution on the boundary may be recorded and
  * added to any other computation as a specific pre-recorded source. Advantage of this source is that interface boundary is numerically
  * (means \b absolutely) transparent and can be used together with absorbing boundary condition as compound input-output boundary condition
  * (for laser pulses, for examples).
  *
  * This tag indicates request for TFSF module support and the size of the TFSF domain.
  *
  * Example of the config file entry is
    <pre>
    [TFSF]
    @ 7             TF/SF xMin [mesh nodes]
    @ 42            TF/SF xMax [mesh nodes]
    @ 7             TF/SF yMin [mesh nodes]
    @ 63            TF/SF yMax [mesh nodes]
    @ 35            TF/SF zMin [mesh nodes]
    @ 63            TF/SF zMax [mesh nodes]
    @ record        TF/SF interface regime (record, play, playBackward)

      TF/SF stands for Total Field/Scattered Field interface.
    </pre>
  *
  */

#ifndef mc_tag_TFSF_header
#define mc_tag_TFSF_header					///< \internal Multiple include guard.

#include <stdio.h>

void tag_TFSF (FILE *fp);

#endif
