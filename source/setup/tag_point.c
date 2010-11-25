/** \file tag_point.c
  * Adds point perturbation to magnetic and electric field at given point (direct usage is a testing of the ABC, TF/SF interface and so on).
  *
  * Example of the config file entry is
    <pre>
    [point]
    @ 10            Ex
    @ 10            Ey
    @ 10            Ez
    @ 0             Hx
    @ 0             Hy
    @ 0             Hz
    @ 5.0           x   [micron]
    @ 5.0           y   [micron]
    @ 4.1           z   [micron]

       Places point perturbation into the domain.
    </pre>
  *
  * Point perturbation corresponds to the smallest packet in the real space and to the widest packet in the Fourier space. Because of this
  * feature this tag is usefull primary for testing of the numerical properties of the schemes used (they stability and dispersion properties).
  */

#include "type_mesh.h"

#include "misc_units.h"
#include "log.h"
#include "misc_cfgReader.h"

// ---------------------------------------------------------------------------
/// Adds point perturbation to the EM-field.
// ---------------------------------------------------------------------------
void
tag_point (FILE *fp, meshVec_p E, meshVec_p H)
{
  int i, j, k, clear = 1;
  double pointEx, pointEy, pointEz, pointHx, pointHy, pointHz, pointX, pointY, pointZ;

  pointEx = cfg_readDouble (fp);
  pointEy = cfg_readDouble (fp);
  pointEz = cfg_readDouble (fp);
  pointHx = cfg_readDouble (fp);
  pointHy = cfg_readDouble (fp);
  pointHz = cfg_readDouble (fp);
  pointX  = cfg_readDouble (fp)*units (mc_micron);
  pointY  = cfg_readDouble (fp)*units (mc_micron);
  pointZ  = cfg_readDouble (fp)*units (mc_micron);

  i = pointX/h1;
  j = pointY/h2;
  k = pointZ/h3;

  say ("tag_point: ");

  if (!mf_mesh_pointIsOutside (E, i, j, k))
  {
    mv_fx(E, i, j, k) += pointEx;
    mv_fy(E, i, j, k) += pointEy;
    mv_fz(E, i, j, k) += pointEz;
    say ("  - added point perturbation to \\vec E (%.2e, %.2e, %.2e) to node [%d, %d, %d]", pointEx, pointEy, pointEz, i, j, k);
    clear = 0;
  }

  if (!mf_mesh_pointIsOutside (H, i, j, k))
  {
    mv_fx(H, i, j, k) += pointHx;
    mv_fy(H, i, j, k) += pointHy;
    mv_fz(H, i, j, k) += pointHz;
    say ("  - added point perturbation to \\vec H (%.2e, %.2e, %.2e) to node [%d, %d, %d]", pointHx, pointHy, pointHz, i, j, k);
    clear = 0;
  }

  if (clear)
    say ("  - no perturbation is added on this node.");
}
