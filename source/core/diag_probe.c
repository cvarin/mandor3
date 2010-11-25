/** \file diag_probe.c
  * \brief Probe diagnostic.
  *
  * Probe is point all field in which are recorded with high resolution in time.
  * This diagnostic is useful to study high frequency part of the information
  * and also it may be used to study frequency spectrums. More details are
  * provided in the diag_probe.h file.
  *
  * \attention The most sutable time to call this diagnostic is after the first
  *            half-step for magnetic field - at this moment all field
  *            components are defined on the same time layer.
  */

#include <stdlib.h>
#include <dirent.h>
#include <unistd.h>

#include <mpi.h>

#include "diag_probe.h"

#include "misc_PIC.h"
#include "misc_units.h"
#include "log.h"
#include "misc_MPItags.h"
#include "misc_cfgReader.h"

// ---------------------------------------------------------------------------
/// Parameters of the probe (local compound, for MPI exchange/run but \b not
/// for file IO).
// ---------------------------------------------------------------------------
typedef struct
{
    FILE *file;
    int   i;
    int   j;
    int   k;
} probeParam_t;

static int           probesN = 0;	///< Number of probes in the system.
static probeParam_t *probes  = NULL;	///< Array of probes.

// ---------------------------------------------------------------------------
/// Opens output files for active probes.
// ---------------------------------------------------------------------------
void
probe_touch (void)
{
    for (int p = 0 ; p < probesN ; ++p)											// Closes all opened files.
        if (probes[p].file) {
            fclose (probes[p].file);
            probes[p].file = NULL;
        }

    MPI_Barrier (MPI_COMM_WORLD);

    // Master saves global parameters.
    if (!cpu_here) {
        FILE *fp = cfg_open ("binData/probe_globals.bin", "wb", __func__);
        probeGlobals_t header = {
            tau,
            h1*mc_have_x,
            h2*mc_have_y,
            h3*mc_have_z,
            { dmn_min[0]*h1*mc_have_x,
              dmn_min[1]*h2*mc_have_y,
              dmn_min[2]*h3*mc_have_z },
            { dmn_max[0]*h1*mc_have_x,
              dmn_max[1]*h2*mc_have_y,
              dmn_max[2]*h3*mc_have_z}
        };
        fwrite (&header, sizeof (probeGlobals_t), 1, fp);
        fclose (fp);
    }

    // Opens files for append.
    for (int p = 0 ; p < probesN ; ++p) {
        if (((probes[p].i < cpu_min[0] || probes[p].i >= cpu_max[0]) && mc_have_x) ||
            ((probes[p].j < cpu_min[1] || probes[p].j >= cpu_max[1]) && mc_have_y) ||
            ((probes[p].k < cpu_min[2] || probes[p].k >= cpu_max[2]) && mc_have_z))
            continue;

        char name[200];
        sprintf (name, "binData/probe_%d_%d_%d.bin",
                       probes[p].i, probes[p].j, probes[p].k);

        // Checks if file exists.
        if (!(probes[p].file = fopen (name, "rb"))) {
            probeHeader_t header = {
                "probes", probes[p].i, probes[p].j, probes[p].k
            };
            probes[p].file = cfg_open (name, "wb", __func__);									// Adds header to just created file.
            fwrite (&header, sizeof (probeHeader_t), 1, probes[p].file);
        } else {
            fclose (probes[p].file);												// Closes file opened by test.
            probes[p].file = cfg_open (name, "ab", __func__);									// Opens to continue output.
        }
    }
}

// ---------------------------------------------------------------------------
/// Deallocates probes.
// ---------------------------------------------------------------------------
static void
probe_deallocate (void)
{
    for (int p = 0 ; p < probesN ; ++p)											// Closes all files on all cpus.
        if (probes[p].file) {
            fclose (probes[p].file);
            probes[p].file = NULL;
        }

    probesN = 0;														// Deallocates memory.
    if (probes) {
        free (probes);
        probes = NULL;
    }
}

// ---------------------------------------------------------------------------
/// \brief Removes all probe data from binData folder.
// ---------------------------------------------------------------------------
void
probe_removeOldFiles (void)
{
    struct dirent *item   = NULL;
    DIR           *folder = opendir ("./binData");											// Opens filelist.
    ENSURE (folder, "cannot open 'binData' folder");

    int N = 0;
    SAY_DEBUG ("probe_removeOldFiles:");
    // Searches for the probe file.
    while ((item = readdir (folder))) {
        if (strncmp (item->d_name, "probe_", 6))										// Skips files without proper prefix.
            continue;

        if (strncmp (item->d_name + strlen (item->d_name) - 4, ".bin", 4))							// Skips files without proper extension.
            continue;

        char name[300];
        sprintf (name, "binData/%s", item->d_name);
        ++N;
        SAY_DEBUG ("  removing %s ...", name);
        remove (name);
    }
    SAY_DEBUG ("  %d files removed.", N);

    closedir (folder);													// Closes folder.
}

// ---------------------------------------------------------------------------
/// Turns <b>(x, y, z)</b> into <b>[i, j, k]</b> and stores new probe without
/// repetition test.
// ---------------------------------------------------------------------------
static void
probe_addDraft (double x, double y, double z)
{
    const double micron = units (mc_micron);
    probes = (probeParam_t *) realloc (probes,
                                       (++probesN)*sizeof (probeParam_t));
    probes[probesN-1].file = NULL;
    probes[probesN-1].i = mc_have_x*MF_DBL_TO_INT (x*micron/h1);
    probes[probesN-1].j = mc_have_y*MF_DBL_TO_INT (y*micron/h2);
    probes[probesN-1].k = mc_have_z*MF_DBL_TO_INT (z*micron/h3);
}

// ---------------------------------------------------------------------------
/// Allocates probes and assigns them to the cpus.
// ---------------------------------------------------------------------------
void
probe_allocate (int cont)
{
    probe_deallocate ();

    if (!cont)
        probe_removeOldFiles ();

    const double micron = units (mc_micron);

    // Master reads parameters.
    if (cpu_here == 0) {
        FILE *fp = fopen ("run_probes.cfg", "rt");
        if (!fp)
            SAY_DEBUG ("probe_allocate: file 'run_probes.cfg' is absent, "
                       "no probe diagnostic at all.");

        // Terminator is used to cut file.
        while (fp) {
            const char *word = cfg_readWord (fp);											// Reads type of the probe set.
            enum { point = 0, line, rectangle, volume, end };
            int probeSet = cfg_identifyWord (word,
                                             "point",     point,
                                             "line",      line,
                                             "rectangle", rectangle,
                                             "volume",    volume,
                                             "end",       end,
                                             mc_cfgTermGuesses);

            switch (probeSet)
            {
                default:
                case -1:
                    ENSURE (0, "bad probe group '%s' (use 'point', 'line', "
                                  "'rectangle', or 'volume')", word);
                break;

                case point: {
                    double x = cfg_readOptDouble (fp),
                           y = cfg_readOptDouble (fp),
                           z = cfg_readOptDouble (fp);
                    probe_addDraft (x, y, z);
                }
                break;

                case line: {
                    double x1 = cfg_readOptDouble (fp),
                           y1 = cfg_readOptDouble (fp),
                           z1 = cfg_readOptDouble (fp),
                           x2 = cfg_readOptDouble (fp),
                           y2 = cfg_readOptDouble (fp),
                           z2 = cfg_readOptDouble (fp);
                    int N = cfg_readOptInt (fp);											// Reads number of additional points in between.

                    N = 2 + abs (N);												// Abs guards against negative N.

                    x2 = (x2 - x1)/(N - 1);											// Turns x2 -> dx, .., z2 -> dz.
                    y2 = (y2 - y1)/(N - 1);
                    z2 = (z2 - z1)/(N - 1);
                    for (int i = 0 ; i < N ; ++i,
                                             x1 += x2, y1 += y2, z1 += z2)
                        probe_addDraft (x1, y1, z1);
                }
                break;

                case rectangle:	 {
                double x = cfg_readOptDouble (fp),	// Rectangular origin.
                       y = cfg_readOptDouble (fp),
                       z = cfg_readOptDouble (fp),
                       px = cfg_readOptDouble (fp),	// Vertex shifted along one side.
                       py = cfg_readOptDouble (fp),
                       pz = cfg_readOptDouble (fp),
                       qx = cfg_readOptDouble (fp),	// Vertex shifted along another side.
                       qy = cfg_readOptDouble (fp),
                       qz = cfg_readOptDouble (fp);
                int Np = cfg_readOptInt (fp),		// Number of additional points in between.
                    Nq = cfg_readOptInt (fp);

                Np = 2 + abs (Np);			// Guard against negative Np/Nq.
                Nq = 2 + abs (Nq);

                px = (px - x)/(Np - 1);			// Turns px to dx and so on.
                py = (py - y)/(Np - 1);
                pz = (pz - z)/(Np - 1);
                qx = (qx - x)/(Nq - 1);											// Turns qx to dx and so on.
                qy = (qy - y)/(Nq - 1);
                qz = (qz - z)/(Nq - 1);

                for (int p = 0 ; p < Np ; ++p)										// Allocates probes.
                    for (int q = 0 ; q < Nq ; ++q)
                        probe_addDraft (x + p*px + q*qx,
                                        y + p*py + q*qy,
                                        z + p*pz + q*qz);
                }
                break;

                case volume: {
                    double x = cfg_readOptDouble (fp),	// Volume origin.
                           y = cfg_readOptDouble (fp),
                           z = cfg_readOptDouble (fp),
                           px = cfg_readOptDouble (fp),	// Vertex shifted along side 1.
                           py = cfg_readOptDouble (fp),
                           pz = cfg_readOptDouble (fp),
                           qx = cfg_readOptDouble (fp),	// Vertex shifted along side 2.
                           qy = cfg_readOptDouble (fp),
                           qz = cfg_readOptDouble (fp),
                           rx = cfg_readOptDouble (fp),	// Vertex shifted along side 3.
                           ry = cfg_readOptDouble (fp),
                           rz = cfg_readOptDouble (fp);
                    int Np = cfg_readOptInt (fp),	// Reads number of additional points in between.
                        Nq = cfg_readOptInt (fp),
                        Nr = cfg_readOptInt (fp);

                    Np = 2 + abs (Np);			// Guard against negative Np/Nq.
                    Nq = 2 + abs (Nq);
                    Nr = 2 + abs (Nr);

                    px = (px - x)/(Np - 1);		// Turns px to dx and so on.
                    py = (py - y)/(Np - 1);
                    pz = (pz - z)/(Np - 1);
                    qx = (qx - x)/(Nq - 1);		// Turns qx to dx and so on.
                    qy = (qy - y)/(Nq - 1);
                    qz = (qz - z)/(Nq - 1);
                    rx = (rx - x)/(Nr - 1);		// Turns rx to dx and so on.
                    ry = (ry - y)/(Nr - 1);
                    rz = (rz - z)/(Nr - 1);

                    for (int p = 0 ; p < Np ; ++p)										// Allocates probes.
                        for (int q = 0 ; q < Nq ; ++q)
                            for (int r = 0 ; r < Nr ; ++r)
                                probe_addDraft (x + p*px + q*qx + r*rx,
                                                y + p*py + q*qy + r*ry,
                                                z + p*pz + q*qz + r*rz);
                }
                break;

                case end:													// Terminator is received.
                break;
            }

            // Terminates loop.
            if (probeSet == end) {
                fclose (fp);
                fp = NULL;
            }
        }

        for (int p1 = 0 ; p1 < probesN ; ++p1)										// Avoids repetitions.
            for (int p2 = p1 + 1 ; p2 < probesN ; ++p2)
                // Replaces repeated probe and steps back
                if ( probes[p1].i == probes[p2].i &&
                     probes[p1].j == probes[p2].j &&
                     probes[p1].k == probes[p2].k )
                    probes[p2--] = probes[--probesN];				//   in p2 to process replacing probe.

        fp = cfg_open ("output/probe_locations.dat", "wt", "probe_allocate");	// File with probe positions marked.
        fprintf (fp, "variables = \"x [<greek>m</greek>m]\", "
                                 "\"y [<greek>m</greek>m]\", "
                                 "\"z [<greek>m</greek>m]\"\n"
                                 "zone t=\"domain\", f= point\n");
        for (int i = 0 ; i < 2 ; ++i)											// Saves corners of the domain.
            for (int j = 0 ; j < 2 ; ++j)
                for (int k = 0 ; k < 2 ; ++k)
                    fprintf (fp, "%e %e %e\n", (dmn_min[0]*i + dmn_max[0]*(1 - i))*h1*mc_have_x/micron,
                                               (dmn_min[1]*j + dmn_max[1]*(1 - j))*h2*mc_have_y/micron,
                                               (dmn_min[2]*k + dmn_max[2]*(1 - k))*h3*mc_have_z/micron);

        fprintf (fp, "zone t=\"probes\", f= point\n");
        for (int p = 0 ; p < probesN ; ++p)
            fprintf (fp, "%e %e %e\n", probes[p].i*h1*mc_have_x/micron,
                                       probes[p].j*h2*mc_have_y/micron,
                                       probes[p].k*h3*mc_have_z/micron);
        fclose (fp);
    }

    MPI_Bcast (&probesN, 1, MPI_INT, 0, MPI_COMM_WORLD);									// Broadcasts parameters.
    if (cpu_here && probesN > 0)												// Allocates probes array on other cpus.
        probes = (probeParam_t *) malloc (probesN*sizeof (probeParam_t));

    if (probesN > 0)													// Gets data.
        MPI_Bcast (probes, probesN*sizeof (probeParam_t), MPI_BYTE, 0, MPI_COMM_WORLD);

    probe_touch ();													// Opens files for active probes.

    ENSURE (!atexit (probe_deallocate),"cannot submit memory deallocator");
}

// ---------------------------------------------------------------------------
/// Writes data.
// ---------------------------------------------------------------------------
void
probe_postData (double time, meshVec_RO_p E, meshVec_RO_p H)
{
  probeParam_t * p = probes;
  for (const probeParam_t * const end = p + probesN ; p < end ; ++p)
  {
    if (!p->file)
      continue;

    probeSample_t data;

    data.time = Time;
    data.Ex = 0.5*(mv_fx (E, p->i, p->j, p->k) + mv_fx (E, p->i + 1, p->j, p->k));
    data.Ey = 0.5*(mv_fy (E, p->i, p->j, p->k) + mv_fy (E, p->i, p->j + 1, p->k));
    data.Ez = 0.5*(mv_fz (E, p->i, p->j, p->k) + mv_fz (E, p->i, p->j, p->k + 1));
    data.Hx = 0.25*(mv_fx (H, p->i, p->j, p->k) + mv_fx (H, p->i, p->j + 1, p->k) + mv_fx (H, p->i, p->j, p->k + 1) + mv_fx (H, p->i, p->j + 1, p->k + 1));
    data.Hy = 0.25*(mv_fy (H, p->i, p->j, p->k) + mv_fy (H, p->i + 1, p->j, p->k) + mv_fy (H, p->i, p->j, p->k + 1) + mv_fy (H, p->i + 1, p->j, p->k + 1));
    data.Hz = 0.25*(mv_fz (H, p->i, p->j, p->k) + mv_fz (H, p->i + 1, p->j, p->k) + mv_fz (H, p->i, p->j + 1, p->k) + mv_fz (H, p->i + 1, p->j + 1, p->k));

    fwrite (&data, sizeof (probeSample_t), 1, p->file);
  }
}

// ---------------------------------------------------------------------------
/// \brief Loads data from saved files. \b data is a pointer to array which is allocated/deallocated here but belongs to caller.
// ---------------------------------------------------------------------------
int
probe_loadData (probeHeader_t *header, probeSample_t **data)
{
  static DIR           *folder = NULL;
  static struct dirent *item = NULL;

  if (! *data)
  {
    folder = opendir ("./binData");											// Opens filelist.
    item = readdir (folder);												// Gets entry.
  }

  ENSURE (folder, "provide NULL '(*data)' pointer to start loading");

  while ((item = readdir (folder)))											// Searches for the probe file.
  {
    if (strncmp (item->d_name, "probe_", 6))										// Skips files without proper prefix.
      continue;

    if (strncmp (item->d_name + strlen (item->d_name) - 4, ".bin", 4))							// Skips files without proper extension.
      continue;

    char name[300];
    sprintf (name, "binData/%s", item->d_name);

    FILE *fp = cfg_open (name, "rb", "probe_loadData");									// Opens file for reading.
    fseek (fp, 0, SEEK_END);
    ENSURE (ftell (fp) >= sizeof (probeHeader_t), "damaged file (too short)");

    fseek (fp, 0, SEEK_SET);												// Checks magic.
    fread (header, sizeof (probeHeader_t), 1, fp);
    if (!mf_probeMagic(header))												// If file has magic it is a good file.
    {
      fclose (fp);
      continue;
    }

    long long startPos = ftell (fp);											// Gets number of records.
    fseek (fp, 0, SEEK_END);
    const long long int N = (ftell (fp) - startPos)/sizeof(probeSample_t);
    fseek (fp, startPos, SEEK_SET);

    *data = (probeSample_t *) realloc (*data, sizeof(probeSample_t)*N);
    ENSURE (*data, "out of memory (%.3f Kb)", sizeof(probeSample_t)*N/1024.0);

    fread (*data, sizeof(probeSample_t), N, fp);									// Reads the data.
    fclose (fp);

    return N;
  }

  if (!item)
  {
    closedir (folder);													// Closes folder.
    free (*data);
    *data = NULL;
    return -1;
  }

  return -1;
}
