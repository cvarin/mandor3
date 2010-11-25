/** \file main.c
  * Fourier transformer of any fourier data.
  */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "frame.h"
#include "type_reg.h"
#include "type_mesh.h"
#include "type_CFile.h"
#include "type_vector.h"

#include "commandLine.h"

#include "fft.h"

#include "reg_unroller.h"

#include "log.h"
#include "misc_cfgReader.h"

static const char  *fft_outFileName = NULL;		///< Name of the output file.
static int          fft_inFileNamesN = 0;		///< Number of input files.
static const char **fft_inFileNames = NULL;		///< Names of the input files.

// ---------------------------------------------------------------------------
/// Reads input/output files from command line.
// ---------------------------------------------------------------------------
static void
main_getFileNames (int argc, char **argv)
{
  assert (argc > 2);
  if (argv[1][0] == '-' || argv[1][0] == '+')										// Gets output file name.
    error ("%s: '%s' - looks like option instead of output file.", argv[0], argv[1]);

  int a = 2;														// Counts input files (second arg, ..).
  while (a < argc && argv[a][0] != '-' && argv[a][0] != '+')   a += 1;
  fft_inFileNamesN = a - 2;

  fft_outFileName = argv[1];												// Aliases name of the output files.
  fft_inFileNames = (const char **) argv + 2;										// Aliases name(s) of the input file(s).
  if (!fft_inFileNamesN)												// Checks input files' presence.
    error ("%s: '%s %s ...' - looks like option instead of input file(s).", argv[0], argv[1], argv[2]);

  say ("Output file:\n  %s", fft_outFileName);
  say ("Input file(s):");
  for (int i = 0 ; i < fft_inFileNamesN ; ++i)
    say ("  %s", fft_inFileNames[i]);
}

// ---------------------------------------------------------------------------
/// Reads mesh storage mapping.
// ---------------------------------------------------------------------------
static regList_t
main_storageDecomposition (void)
{
  assert (fft_inFileNamesN);

  regList_t map = mc_regList_init;
  for (int f = 0 ; f < fft_inFileNamesN ; ++f)										// Collects region sizes from files.
  {
    reg_t reg;
    FILE *fp = cfg_open (fft_inFileNames[f], "rb", __func__);								// Loads mesh sizes.
    fread (reg.min, sizeof (int), 3, fp);
    fread (reg.max, sizeof (int), 3, fp);
    fclose (fp);

    regList_add (&map, &reg);
  }

  return map;
}

// ---------------------------------------------------------------------------
/// Analizes domain decomposition and returns subdomain size.
// ---------------------------------------------------------------------------
static reg_t
main_configuredReg (const regList_t *map)
{
  reg_t dmn, subdmn = {{INT_MIN, INT_MIN, INT_MIN}, {INT_MAX, INT_MAX, INT_MAX}};

  if (cl_findVec3i ("fft->min r", subdmn.min))										// Gets space subdomain sizes.
    say ("%s: subdomain min is imported from command line.", __func__);

  if (cl_findVec3i ("fft->max r", subdmn.max))
    say ("%s: subdomain max is imported from command line.", __func__);

  MF_VEC_COPY (dmn.min, map->list[0].min);										// Initial guess for min/max search.
  MF_VEC_COPY (dmn.max, map->list[0].max);
  for (reg_t *r = map->list ; r < map->list + map->N ; ++r)								// Computes domain size.
    for (int axis = 0 ; axis < 3 ; ++axis)
    {
      dmn.min[axis] = (dmn.min[axis] <= r->min[axis]) ? dmn.min[axis] : r->min[axis];
      dmn.max[axis] = (dmn.max[axis] >= r->max[axis]) ? dmn.max[axis] : r->max[axis];
    }

  mf_reg_collapse (&dmn);												// Accounts for dimensionality.
  reg_overlap (&subdmn, &dmn, NULL);											// Fits subdomain to domain.

  say ("Domain: '%s'", reg_printRanges (&dmn));
  say ("Subdomain: '%s'", reg_printRanges (&subdmn));

  // Check that 'dmn' is covered by 'map' completely.
  regList_t pieces = mc_regList_init;											// Puts dmn in the list.
  regList_add (&pieces, &dmn);
  for (reg_t *r = map->list ; r < map->list + map->N ; ++r)								// Projects dmn on map.
    for (int axis = 0 ; axis < 3 ; ++axis)
    {
      regList_slice (&pieces, axis, r->min[axis]);
      regList_slice (&pieces, axis, r->max[axis]);
    }

  for (reg_t *r = pieces.list ; r < pieces.list + pieces.N ; ++r)							// Checks that all pieces have host.
  {
    int noHost = 1;
    for (reg_t *h = map->list ; h < map->list + map->N && noHost ; ++h)							// Checks all hosts.
      if (reg_isInside (r, h))
        noHost = 0;

    if (noHost)
      error ("%s: domain is not completely covered:\n  o domain '%s'\n  o homeless piece '%s'", __func__,
             reg_printRanges (&dmn), reg_printRanges (r));
  }

  return subdmn;
}

// ---------------------------------------------------------------------------
/// Reads parameters of the config file.
// ---------------------------------------------------------------------------
static double*
main_loadData (const reg_t *subdmn)
{
  char *type;
  if (!cl_findString ("fft->type", &type))
    error ("%s: type is not described (use '--fft:vector{x|y|z}' or '--fft:scalar' options).");

  const int skip = cfg_identifyWord (type, "scalar", 1, "vector", 3, mc_cfgTermGuesses);				// Gets spacing of the components.
  if (skip == -1)
    error ("%s: bad type (%s), should be 'scalar' or 'vector'.", __func__, type);
  free (type);

  int offset = 0;
  if (skip == 3)													// Gets component of the vector.
  {
    int typeKeyExists = cl_findString ("fft->component", &type);
    assert (typeKeyExists);
    offset = cfg_identifyWord (type, "x", 0, "y", 1, "z", 2, mc_cfgTermGuesses);					// Gets components' shift.
    assert (offset >= 0 && offset < 3);
    free (type);
  }

  double *data = (double *) malloc (reg_volume (subdmn)*sizeof (double));						// Allocates and NaNifies storage.
  assert (data);
  memset (data, -1, reg_volume (subdmn)*sizeof (double));

  for (int f = 0 ; f < fft_inFileNamesN ; ++f)										// Loads mesh subregion.
  {
    reg_t mesh;
    FILE *fp = cfg_open (fft_inFileNames[f], "rb", __func__);								// Loads mesh sizes.
    fread (mesh.min, sizeof (int), 3, fp);
    fread (mesh.max, sizeof (int), 3, fp);
    fseek (fp, sizeof (mesh_t), SEEK_CUR);										// Skips mesh_t structure.
    const long long int start = ftell (fp) + sizeof (double)*offset;							// Saves data segment position.

    reg_t load = *subdmn;
    if (reg_overlap (&load, &mesh, NULL))										// Loads overlapped part.
    {
      fclose (fp);
      continue;
    }

    long long int _s_[3], _o_, _S_[3], _O_;										// Unrolling coefficients for domain and mesh.
    MF_UNROLL(_S_, _O_, subdmn);
    MF_UNROLL(_s_, _o_, &mesh);

    double tmp[3];
    for (int i = load.min[0] ; i <= load.max[0] ; ++i)
      for (int j = load.min[1] ; j <= load.max[1] ; ++j)
      {
        fseek (fp, start + skip*sizeof (double)*MF_AIM(i, j, load.min[2], _s_, _o_), SEEK_SET);				// Positions to the start of span.
        for (int k = load.min[2] ; k <= load.max[2] ; ++k)								// Loads span.
        {
          fread (data + MF_AIM(i, j, k, _S_, _O_), sizeof (double), 1, fp);
          fread (tmp, sizeof (double), skip - 1, fp);
        }
      }

    fclose (fp);
  }

  return data;
}

// ---------------------------------------------------------------------------
/// Writes data in tecplot format.
// ---------------------------------------------------------------------------
static void
main_exportToTecplot (const reg_t *dmn, const double *data, FILE *fp, const double *scale)
{
  fprintf (fp, "variables = \"mx\", \"my\", \"mz\", \"f<sub>k</sub>\"\n");
  fprintf (fp, "zone t=\"fft\", i = %d, j = %d, k = %d, f = block\n",
            dmn->max[0] - dmn->min[0] + 1, dmn->max[1] - dmn->min[1] + 1, dmn->max[2] - dmn->min[2] + 1);

  long long int _s_[3], _o_;
  MF_UNROLL(_s_, _o_, dmn);

  #define MF_SCAN 										\
  for (int k = dmn->min[2] ; k <= dmn->max[2] ; ++k)						\
    for (int j = dmn->min[1] ; j <= dmn->max[1] ; ++j)						\
      for (int i = dmn->min[0] ; i <= dmn->max[0] ; ++i)

  MF_SCAN  fprintf (fp, "%.3e\n", i*scale[0]);
  MF_SCAN  fprintf (fp, "%.3e\n", j*scale[1]);
  MF_SCAN  fprintf (fp, "%.3e\n", k*scale[2]);
  MF_SCAN  fprintf (fp, "%.3e\n", data[MF_AIM(i, j, k, _s_, _o_)]);

  #undef MF_SCAN
}

// ---------------------------------------------------------------------------
/// Entry point for spectr.out diagnostic.
// ---------------------------------------------------------------------------
int
main (int argc, char *argv[])
{
  if (argc < 2)
    error ("Usage: %s <output file name> <input file(s)> [options].", argv[0]);

  msg_openLogFile ("output/logs/fft.log", 1);
  cl_import ("./source/import.pl", argc, argv, 2);									// Imports data to the command line parser.

  main_getFileNames (argc, argv);											// Inits local filename vars.
  regList_t map = main_storageDecomposition ();										// Extracts mesh partitioning.
  reg_t spaceReg = main_configuredReg (&map);										// Extracts loadable subdomain.
  double *data = main_loadData (&spaceReg);										// Loads data.

  reg_t spectrReg = {{INT_MIN, INT_MIN, INT_MIN}, {INT_MAX, INT_MAX, INT_MAX}};						// Default is the biggest domain possible.
  fft (&spaceReg, data, &spectrReg);

  double lambda[3] = {1, 1, 1};												// Computes scaling coefficients.
  cl_findVec3d ("fft->lambda", lambda);											// Overrides default by command line.
  say ("lambda_0 = [%.3f, %.3f, %.3f] mesh steps.", lambda[0], lambda[1], lambda[2]);
  for (int axis = 0 ; axis < 3 ; ++axis)
    lambda[axis] = fabs (lambda[axis])/(spaceReg.max[axis] - spaceReg.min[axis] + 1 - ACTIVATOR[axis]);

  FILE *fp = cfg_open (fft_outFileName, "wt", argv[0]);
  main_exportToTecplot (&spectrReg, data, fp, lambda);
  fclose (fp);

  free (data);
  regList_clean (&map);

  say ("Done.");

  return EXIT_SUCCESS;
}
