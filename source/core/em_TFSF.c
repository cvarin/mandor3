/** \file em_TFSF.c
  * \brief 'Total Field/Scattered Field' interface support (see em_TFSF.h).
  *
  * Data structure:
  *
  * - Interface is a parallelepiped with 6 faces around region which is
  *   initially field-free (to support initial fields the field in the volume
  *   should be loaded too).
  *
  * - All 6 faces are grouped in a list of nodeFace_s structures.
  *
  * - Each face has two subgroups to deal with two different polarizations; only
  *   tangential electric field nodes are located on the face planes, so normal
  *   component of the field is not processed.
  *
  * - Face is splitted on lines (nodeLine_s arrays). Line structure hosts
  *   detailed information about how data is unrolled in 1D buffers.
  *
  * How it works:
  *
  * - Parallel \b play is done by simply reading source for RHS directly from
  *   prerecorded file.
  *
  * - Parallel <b>backward play</b> is done by reading source for RHS directly
  *   from prerecorded file, changing time \b t into <b>T - t</b> (\b T is the
  *   time duration of the recording session) and changing sign of magnetic
  *   field.
  *
  * - Parallel \b recording is two pass process. Main pass is parallel with
  *   output from each node sended into separate file stream (to remove messy
  *   parallel syncronization in order to construct more solid output). That
  *   means that one more pass is required to merge all file streams into solid
  *   one (to avoid multiple partitioning issues). Separate merging pass is done
  *   in scalar mode because of small resources required: merging has the same
  *   speed like direct file<->file copy, memory pressure is small
  *   (dimensionality is smaller by one) and parallelization is not required.
  *
  * \note Electric and magnetic field, used for Lorentz force interpolation, are
  *       not corrected before applying to particles. The reason is simple: if
  *       particle is crossing an interface the solution is not physical anyway
  *       (or field is so small already that correction doesn't really change
  *       much).
  *
  * \warning For big simulations the size of data file may be bigger than few Gb
  *          so additional support for 64-bits positioning in file streams
  *          should be used (on SGI \b ftell64 and \b fseek64 functions are
  *          used).
  *
  * \warning Visualization artifact (see em_caps.h for explanation and example).
  */

#include <math.h>
#include <unistd.h>
#include <stdlib.h>

#include <mpi.h>

#include "misc_units.h"
#include "log.h"
#include "misc_cfgReader.h"
#include "misc_definedKeys.h"

#include "em_TFSF.h"

/// All defragmented recorded streams are assigned to the stream with
/// unrealistic cpu number.
#define mc_TFSF_joinedStreamsCPU	(99999)

/// Group of lines which forms face of the interface (not necessarily
/// rectangular).
typedef struct
{
  reg_t		face;		///< Total face region (one before any partitioning).
  reg_t		toCache;	///< Local part of the face to IO/cache.
  double       *buffer;		///< Buffer for unrolled (on \b tfsfFace_t::toCache region) data for solid IO.
  double        sourceFactor;	///< Weight of the initial field in the source term (accounted on the record pass after caching).
  int           component;	///< Component of field (one of two possible tangent components, -1 <=> face is deactivated).
  int           ep;		///< Tangential vector number 1. Used in the 2D unrolling.
  int           eq;		///< Tangential vector number 2. Used in the 2D unrolling.
  const char   *name;		///< Name of the face (used as part of the file name).
  FILE         *fp;		///< File pointer (to keep file stream and cached data grouped together).
} tfsfFace_t;

/// Storage order index generator for tfsfE, tfsfH: "axis/top/pq"; index = "axis*4 + top*2 + {\b 0 for \b p / \b 1 for \b q}".
#define mf_face(axis, top, comp) (((axis) << 2) + ((top) << 1) + (comp))

// ---------------------------------------------------------------------------
/// \brief Structures to process \b electric field sources on the TFSF faces.
/// Storage order: 'axis/top/pq', in other word index of array is
/// 'axis*4 + top*2 + {\b 0 for \b p / \b 1 for \b q}'.
// ---------------------------------------------------------------------------
static tfsfFace_t tfsfE[12] = {
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ey_xMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ez_xMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ey_xMax", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ez_xMax", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ez_yMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ex_yMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ez_yMax", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ex_yMax", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ex_zMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ey_zMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ex_zMax", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Ey_zMax", NULL }
};

// ---------------------------------------------------------------------------
/// Structures to process \b magnetic field sources on the TFSF faces (storage order like for tfsfE).
// ---------------------------------------------------------------------------
static tfsfFace_t tfsfH[12] = {
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hz_xMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hy_xMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hz_xMax", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hy_xMax", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hx_yMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hz_yMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hx_yMax", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hz_yMax", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hy_zMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hx_zMin", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hy_zMax", NULL },
    { mc_reg_init, mc_reg_init, NULL, 0, -1, 0, 0, "Hx_zMax", NULL }
};

static int   tfsfLastFrame   = 0;	///< Number of frames written, updated and written in \b record mode.
static int   tfsfFrame       = 0;	///< Current frame number, used in \b play modes.
static int   tfsfCachedFrame = -1;	///< Current cached frame number, used in \b play modes.
static reg_t tfsfInterface   = {	///< Interfaced region (barcode holds mode).
    .min     = {0, 0, 0},
    .max     = {0, 0, 0},
    .cpu     = 0,
    .barcode = mc_TFSF_uninitialized
};

// Declaration of the function to submit to atexit.
static void TFSF_shutDown (void);


// ---------------------------------------------------------------------------
/// Wrapper around \b sprintf for convinience.
// ---------------------------------------------------------------------------
static char*
TFSF_filename (tfsfFace_t *face, int cpu)
{
   static char name[200];
   sprintf (name, ".EM_sources/TFSF/%s_%03d.bin", face->name, cpu);									// Assembles file name.
   return name;
}

// ---------------------------------------------------------------------------
/// Reads parameters from the config file and distributes it over cluster.
// ---------------------------------------------------------------------------
static void
TFSF_readParameters (void)
{
    if (!cpu_here) {
        FILE *fp = fopen (".EM_sources/TFSF/TFSF.cfg", "rt");
        // Reads basic parameters.
        if (fp) {
            tfsfInterface.min[0] = cfg_readInt (fp);
            tfsfInterface.max[0] = cfg_readInt (fp);
            tfsfInterface.min[1] = cfg_readInt (fp);
            tfsfInterface.max[1] = cfg_readInt (fp);
            tfsfInterface.min[2] = cfg_readInt (fp);
            tfsfInterface.max[2] = cfg_readInt (fp);

            // Barcode holds mode (rec or play <-|->).
            tfsfInterface.barcode = cfg_readInt (fp);
            fclose (fp);
            mf_reg_collapse (&tfsfInterface);
        } else
            tfsfInterface.barcode = mc_TFSF_uninitialized;
    }

    MPI_Bcast (&tfsfInterface, sizeof (reg_t), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Checks if mode is one of the active.
    if (tfsfInterface.barcode == mc_TFSF_uninitialized) {
        SAY_DEBUG ("TFSF_readParameters: file 'TFSF/TFSF.cfg' doesn't exist, "
                   "no interface activated.");
        tfsfInterface.barcode = mc_TFSF_inactive;
        return;
    }

    say ("TFSF:\n  - interface mesh region is %s;",
              reg_printRanges (&tfsfInterface));
    say ("  - interface space region is [%.3f,%.3f]x[%.3f,%.3f]x[%.3f,%.3f] [r0];",
              tfsfInterface.min[0]*h1*mc_have_x, tfsfInterface.max[0]*h1*mc_have_x,
              tfsfInterface.min[1]*h2*mc_have_y, tfsfInterface.max[1]*h2*mc_have_y,
              tfsfInterface.min[2]*h3*mc_have_z, tfsfInterface.max[2]*h3*mc_have_z);
    say ("  - interface space region is [%.3f,%.3f]x[%.3f,%.3f]x[%.3f,%.3f] [micron];",
              tfsfInterface.min[0]*h1*mc_have_x/units (mc_micron), tfsfInterface.max[0]*h1*mc_have_x/units (mc_micron),
              tfsfInterface.min[1]*h2*mc_have_y/units (mc_micron), tfsfInterface.max[1]*h2*mc_have_y/units (mc_micron),
              tfsfInterface.min[2]*h3*mc_have_z/units (mc_micron), tfsfInterface.max[2]*h3*mc_have_z/units (mc_micron));

    ENSURE (reg_volume (&tfsfInterface) > 0, "bad interface region");
}

// ---------------------------------------------------------------------------
/// Truncates all faces to fit to the given domain.
// ---------------------------------------------------------------------------
static void
TFSF_shapeFaces (tfsfFace_t *face, const reg_t *domain, const reg_t *ghost)
{
    SAY_DEBUG ("  - Fitting geometry to domain %s with ghost %s...",
               reg_printRanges (domain),
               (ghost) ? reg_printRanges (ghost)
                       : "NULL");

    // Scans all faces.
    for (const tfsfFace_t * const end = face + 12 ; face < end ; ++face) {
        // Skips deactivated faces.
        if (face->component < 0)
            continue;

        // Gets full size to fit to domain.
        face->toCache = face->face;

        // Reduces face to fit into domain.
        if (reg_overlap (&face->toCache, domain, ghost)) {
            face->component = -1;
            SAY_DEBUG ("    face '%s' is excluded completely", face->name);
        } else {
            int size = reg_volume (&face->toCache)*sizeof (double);								// Allocates memory.
            SAY_DEBUG ("    face '%s' stays as %s (%f Kb)", face->name,
                       reg_printRanges (&face->toCache), size/1024.0);
            face->buffer = (double*) malloc (size);
            ENSURE (face->buffer, "out of memory (face %s, %.3f Kb)",
                       face->name, size/1024.0);
        }
    }
}

// ---------------------------------------------------------------------------
/// Prepares information for defragmentation of file streams.
// ---------------------------------------------------------------------------
static void
TFSF_defragInit (tfsfFace_t *face, const char *fileName)
{
    int  faceFlag[12],
        *flags = faceFlag;
    for (const tfsfFace_t * const end = face + 12 ; face < end ; ++face)
        *(flags++) = (face->component >= 0);
    face -= 12;														// Gets original face pointer.

    flags = (int*) calloc (cpu_total*12, sizeof (int));
    MPI_Gather (faceFlag, 12, MPI_INT, flags, 12, MPI_INT, 0, MPI_COMM_WORLD);						// Assembles map on master node.

    // Writes consistency check parameters.
    if (!cpu_here) {
        FILE *fp = cfg_open (fileName, "at", __func__);
        for (int f = 0 ; f < 12 ; ++f) {
            // Saves face, plane unrolling order, file names for result and each
            // of all fragments; records are separated with '---'.
            fprintf (fp, "%s\n", reg_printRanges (&face[f].face));
            fprintf (fp, "ep = %d, eq = %d\n", face[f].ep, face[f].eq);
            fprintf (fp, "%s\n", TFSF_filename (face + f,
                                                mc_TFSF_joinedStreamsCPU));
            for (int cpu = 0 ; cpu < cpu_total ; ++cpu)
                if (flags[cpu*12+f])
                    fprintf (fp, "%s\n", TFSF_filename (face + f, cpu));
            fprintf (fp, "---\n");
        }
        fclose (fp);
    }

    free (flags);
}

// ---------------------------------------------------------------------------
/// Prepares interface for writing data.
// ---------------------------------------------------------------------------
static void
TFSF_turnRec (void)
{
    const reg_t domain = {
        {cpu_min[0], cpu_min[1], cpu_min[2]},
        {cpu_max[0], cpu_max[1], cpu_max[2]}
    };

    // Offsets for E and H field domains.
    static const reg_t ghostE = {{0, 0, 0}, {-1, -1, -1}},
                       ghostH = {{1, 1, 1}, { 0,  0,  0}};

    // Fits E/H faces to node subdomain.
    TFSF_shapeFaces (tfsfE, &domain, &ghostE);
    TFSF_shapeFaces (tfsfH, &domain, &ghostH);

    // Pointer on coupled E and H faces.
    tfsfFace_t *faceE = tfsfE,
               *faceH = tfsfH,
               *end   = faceE + 12;
    for ( ; faceE < end ; ++faceE, ++faceH) {
        if (faceE->component >= 0) {
            faceE->fp = cfg_open (TFSF_filename (faceE, cpu_here),
                                  "wb", __func__);
            fwrite (&faceE->toCache, sizeof (reg_t), 1, faceE->fp);
        }

        if (faceH->component >= 0) {
            faceH->fp = cfg_open (TFSF_filename (faceH, cpu_here),
                                  "wb", __func__);
            fwrite (&faceH->toCache, sizeof (reg_t), 1, faceH->fp);
        }
    }

    // Writes consistency check parameters.
    if (!cpu_here) {
        FILE *fp = cfg_open (".EM_sources/TFSF/consistencyCheck.dat", "wt",
                             __func__);
        fprintf (fp, "@ %d	Interface size in X direction (nodes)\n",
                 (tfsfInterface.max[0] - tfsfInterface.min[0])*mc_have_x);
        fprintf (fp, "@ %d	Interface size in Y direction (nodes)\n",
                 (tfsfInterface.max[1] - tfsfInterface.min[1])*mc_have_y);
        fprintf (fp, "@ %d	Interface size in Z direction (nodes)\n",
                 (tfsfInterface.max[2] - tfsfInterface.min[2])*mc_have_z);
        fprintf (fp, "@ %.12e	Time step\n", tau);
        fprintf (fp, "@ %.12e	Spatial step along X\n", h1);
        fprintf (fp, "@ %.12e	Spatial step along Y\n", h2);
        fprintf (fp, "@ %.12e	Spatial step along Z\n", h3);
        fclose (fp);
    }

    // Resets frame counter.
    tfsfLastFrame = 0;

    // Exports information for defragmentator.
    TFSF_defragInit (tfsfE, ".EM_sources/TFSF/TFSF_defrag.map");
    TFSF_defragInit (tfsfH, ".EM_sources/TFSF/TFSF_defrag.map");

    say ("  - interface is initialized in 'record' mode;");
}

// ---------------------------------------------------------------------------
/// \brief Turns faces off if requested (using tag [TFSF: open faces]).
// ---------------------------------------------------------------------------
static void
TFSF_openFaces (void)
{
    int faces[2][3];

    // Gets opening marks.
    FILE *fp = fopen (".EM_sources/TFSF/openFaces.cfg", "rt");
    if (!fp)
        return;

    for (int axis = 0 ; axis < 3 ; ++axis)
        for (int top = 0 ; top < 2 ; ++top)
            faces[top][axis] = cfg_readInt (fp);
    fclose (fp);

    // Removes faces across degenerated axises.
    say ("  Removing faces:");
    for (int axis = 0 ; axis < 3 ; ++axis) {
        if (!ACTIVATOR[axis])
            continue;

        for (int top = 0  ; top  < 2 ; ++top)
        for (int comp = 0 ; comp < 2 ; ++comp)
            if (faces[top][axis]) {
                tfsfE[mf_face(axis, top, comp)].component = -1;
                tfsfH[mf_face(axis, top, comp)].component = -1;
                say ("  - '%s' and '%s'",
                            tfsfE[mf_face(axis, top, comp)].name,
                            tfsfH[mf_face(axis, top, comp)].name);
            }
    }
}


// ---------------------------------------------------------------------------
/// \brief Prepares interface for using written data. <b>File streams for E and
/// H are exchanged because of magnetic field is a source for electric and vice
/// versa.</b>
// ---------------------------------------------------------------------------
static void
TFSF_turnPlay (void)
{
    // Removes opened faces from futher processing.
    TFSF_openFaces ();

    // Shifts TFSF from boundary to avoid a mess with absorbing condition.
    const reg_t global = {
        { dmn_min[0] + 5*mc_have_x,
          dmn_min[1] + 5*mc_have_y,
          dmn_min[2] + 5*mc_have_z },
        { dmn_max[0] - 5*mc_have_x,
          dmn_max[1] - 5*mc_have_y,
          dmn_max[2] - 5*mc_have_z }
    };
    for (int f = 0 ; f < 12 ; ++f) {
        ENSURE (tfsfH[f].component < 0 ||
                   reg_isInside (&tfsfH[f].face, &global),
                   "TFSF interface is too close to a boundary (less than 5 cell)");
        ENSURE (tfsfE[f].component < 0 ||
                   reg_isInside (&tfsfE[f].face, &global),
                   "TFSF interface is too close to a boundary (less than 5 cell)");
    }

    // Exctacts local piece of TFSF face to handle.
    const reg_t domain = {
        { cpu_min[0] - mc_have_x,
          cpu_min[1] - mc_have_y,
          cpu_min[2] - mc_have_z },
        { cpu_max[0] + mc_have_x,
          cpu_max[1] + mc_have_y,
          cpu_max[2] + mc_have_z
        }
    };
    TFSF_shapeFaces (tfsfH, &domain, NULL);
    TFSF_shapeFaces (tfsfE, &domain, NULL);

    // Opens recorded source.
    for (int f = 0 ; f < 12 ; ++f) {
        //
        // IMPORTANT:
        //     File streams for E and H are exchanged because of magnetic field
        //     is a source for electric field and vice versa.
        //
        if (tfsfE[f].component >= 0) {
            const char *file = TFSF_filename (tfsfH + f, mc_TFSF_joinedStreamsCPU);
            tfsfE[f].fp = cfg_open (file, "rb", __func__);
        }

        if (tfsfH[f].component >= 0) {
            const char *file = TFSF_filename (tfsfE + f, mc_TFSF_joinedStreamsCPU);
            tfsfH[f].fp = cfg_open (file, "rb", __func__);
        }
    }

    FILE *fp = fopen (".EM_sources/TFSF/consistencyCheck.dat", "rt");
    ENSURE (fp, "there are no recorded data");

    for (int axis = 0 ; axis < 3 ; ++axis) {
        if (cfg_readInt (fp) != tfsfInterface.max[axis] - tfsfInterface.min[axis]
           && ACTIVATOR[axis])
            DIE ("inconsistent size along axis '%c'", 'X' + axis);
    }

    double oldTau = cfg_readDouble(fp);
    ENSURE (fabs (1 - oldTau/tau) < 1e-11,
            "inconsistent time step (τ = %e vs %e, δτ = %e)",
            tau, oldTau, tau - oldTau);

    ENSURE (mc_have_x*fabs (1 - cfg_readDouble (fp)/h1) < 1e-11,
            "inconsistent spatial step along X");

    ENSURE (mc_have_y*fabs (1 - cfg_readDouble (fp)/h2) < 1e-11,
            "inconsistent spatial step along Y");

    ENSURE (mc_have_z*fabs (1 - cfg_readDouble (fp)/h3) < 1e-11,
            "inconsistent spatial step along Z");
    fclose (fp);

    if (!cpu_here) {
        // Gets reference time T to inverse t -> T - t.
        fp = cfg_open (".EM_sources/TFSF/lastFrameTime.bin", "rb", __func__);
        fread (&tfsfLastFrame, sizeof (int), 1, fp);
        fclose (fp);
    }

    MPI_Bcast (&tfsfLastFrame, sizeof (int), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Resets counter of frames.
    tfsfFrame = (int) (Time/tau + 0.1);

    say ("  - interface is initialized in 'play%s' mode with "
              "T = %f (%d frames);",
              (tfsfInterface.barcode == mc_TFSF_playForward) ? "Forward" :
                                                               "Backward",
              tfsfLastFrame*tau, tfsfLastFrame);
}


// ---------------------------------------------------------------------------
/// Initial setup. If mode is \b inactive function still reads and reports all
/// parameters from file '.EM_sources/TFSF.cfg' (for usage in multi-pass setup
/// sessions and for debug reason - to help catch errors before actual run).
// ---------------------------------------------------------------------------
void
TFSF_init (void)
{
    ENSURE (tfsfInterface.barcode == mc_TFSF_uninitialized ||
               tfsfInterface.barcode == mc_TFSF_inactive,
               "double initialization");

    ENSURE (cpu_total < mc_TFSF_joinedStreamsCPU,
               "to small 'unrealistic' cpu number (appears to be realistic)");

    ENSURE (access (".EM_sources/TFSF/TFSF_defrag.map", F_OK),
               "input data is not compiled (use 'TFSFdefrag')");

    // Master reads and broadcasts parameters.
    TFSF_readParameters ();
    if (tfsfInterface.barcode == mc_TFSF_inactive)
        return;

    // Deactivation of all faces.
    for (int f = 0 ; f < 12 ; ++f)
        tfsfE[f].component = tfsfH[f].component = -1;

    const int* sides[2] = {			  // Top/bottom sides of interface.
        tfsfInterface.min,
        tfsfInterface.max
    };
    const double c[3] = {tau/h1, tau/h2, tau/h3}; // Transport coefficients.

    // Sets faces geometry.
    for (int axis = 0 ; axis < 3 ; ++axis) {
        // Skips degenerated axises.
        if (!ACTIVATOR[axis])
            continue;

        // Tangential frame.
        int ep = (axis + 1) % 3,
            eq = (axis + 2) % 3;
        // To keep 2D inner cycles active if possible.
        if (!ACTIVATOR[eq]) {
            int tmp = ep;
            ep = eq;
            eq = tmp;
        }

        for (int top = 0 ; top < 2 ; ++top)
            for (int comp = 0 ; comp < 2 ; ++comp) {
                tfsfFace_t *faceE = tfsfE + mf_face(axis, top, comp);			// Pointer on electric field face.
                faceE->face           = tfsfInterface;					// Sets bounding box sizes to start from.
                faceE->face.min[axis] = faceE->face.max[axis] = sides[top][axis];	// Sets position along axis.
                faceE->sourceFactor = (1 - 2*top)*(1 - 2*comp)*c[axis];			// Deducted from old transport coefficients.
                faceE->component = (axis + 1 + comp) % 3;				// Sets component to cache/update.
                faceE->face.min[faceE->component] += ACTIVATOR[faceE->component];	// Accounts for shift of the E-nodes.
                faceE->ep = ep;
                faceE->eq = eq;

                tfsfFace_t *faceH = tfsfH + mf_face(axis, top, comp);			// Pointer on magnetic field face.
                faceH->face = faceE->face;						// Magnetic field supports the same set of nodes.
                faceH->face.min[axis] += top;						// Accounts shift of magnetic field nodes.
                faceH->face.max[axis] += top;
                faceH->sourceFactor = faceE->sourceFactor;				// Deducted from old transport coefficients.
                faceH->component = 3 - faceE->component - axis;				// Sets component to cache/update.
                faceH->ep = ep;
                faceH->eq = eq;

                SAY_DEBUG ("    o adding '%s': comp %d, reg %s ...", faceE->name, faceE->component, reg_printRanges (& faceE->face));
                SAY_DEBUG ("    o adding '%s': comp %d, reg %s ...", faceH->name, faceH->component, reg_printRanges (& faceH->face));
            }
    }

    // Sets files and shapes 'toCache'.
    switch (tfsfInterface.barcode)
    {
        case  mc_TFSF_rec:
            TFSF_turnRec ();
        break;

        case mc_TFSF_playForward:
        case mc_TFSF_playBackward:
            TFSF_turnPlay ();
        break;
    }

    // Submits memory cleaner/file finalizator.
    ENSURE (!atexit (TFSF_shutDown), "cannot submit TFSF_shutDown");

    say("  All done.");
}

// ---------------------------------------------------------------------------
/// Writes field from mesh to the storage array (should be called at the
/// beginning of the main loop to cache \f$ E^n \f$ and \f$ H^{n-1/2} \f$.
// ---------------------------------------------------------------------------
static void
TFSF_writeMesh (tfsfFace_t *face, meshVecI_p F)
{
    /// \todo Strides in unrolled array -> mesh/reg/wrap.
    const int stride[3] = {
        mc_have_x*F->width_yz,
        mc_have_y*F->width_z,
        mc_have_z
    };

    for (const tfsfFace_t * const end = face + 12 ; face < end ; ++face) {
        if (face->fp) {
            const int c  = face->component,
                      ep = face->ep,
                      eq = face->eq;
            const int dp = face->toCache.max[ep] - face->toCache.min[ep],
                      dq = face->toCache.max[eq] - face->toCache.min[eq];
            double *buffer = face->buffer;
            const double coefficient = face->sourceFactor;

            int pos0 = mf_offset(F, face->toCache.min[0],
                                    face->toCache.min[1],
                                    face->toCache.min[2]);
            for (int p = 0 ; p <= dp ; ++p, pos0 += stride[ep]) {
                int pos = pos0;
                for (int q = 0 ; q <= dq ; ++q, pos += stride[eq], ++buffer)
                    *buffer = mv_unrl(F, pos).r[c]*coefficient;
            }

            fwrite (face->buffer, sizeof (double),
                    reg_volume (&face->toCache), face->fp);
        }
    }
}

static void TFSF_cacheFrame (void);

// ---------------------------------------------------------------------------
/// \brief Writes \f$ E^n \f$ and \f$ H^{n-1/2} \f$ to the data file or caches data for play time-step. <b>Should be called before any EM timestep activity</b>
// ---------------------------------------------------------------------------
void
TFSF_preHStep (meshVecI_p E, meshVecI_p H)
{
    TFSF_cacheFrame ();

    if (tfsfInterface.barcode != mc_TFSF_rec)
        return;

    // Writes caches on the surface of TFSF data and counts the frame.
    TFSF_writeMesh (tfsfE, E);
    TFSF_writeMesh (tfsfH, H);
    ++tfsfLastFrame;
}

// ---------------------------------------------------------------------------
/// \brief Gets data from file and writes it to the unrolled array.
// ---------------------------------------------------------------------------
static void
TFSF_cacheSource (tfsfFace_t *face, int frame)
{
    for (const tfsfFace_t * const end = face + 12 ; face < end ; ++face)
    {
        if (face->fp) {
            const reg_t * const f     = &face->face;
            const reg_t * const cache = &face->toCache;
            const int ep = face->ep,
                      eq = face->eq;
            const int dp = cache->max[ep] - cache->min[ep];
            const int strideFile  = (    f->max[eq] -     f->min[eq] + 1);								/// Reads subset of all data. \todo Span_t ?
            const int strideCache = (cache->max[eq] - cache->min[eq] + 1);
            long long int filePos =
                    sizeof (reg_t) + frame*sizeof (double)*reg_volume (f) +
                    sizeof (double)*           (cache->min[eq] - f->min[eq]) +
                    sizeof (double)*strideFile*(cache->min[ep] - f->min[ep]);
        /*      SAY_DEBUG ("Caching %s / %s (%s) [%d, %d]...", face->name, reg_printRanges (cache), reg_printRanges (f), dp, strideCache);*/
            if (face->component < 0)
                continue;

            double *buffer = face->buffer;
            for (int p = 0 ; p <= dp ; ++p) {
                fseek64 (face->fp, filePos, SEEK_SET);										// Gets data cached for given face line.
                fread (buffer, sizeof (double), strideCache, face->fp);
                buffer  += strideCache;
                filePos += strideFile*sizeof (double);
            }
        }
    }
}

// ---------------------------------------------------------------------------
/// Gets data from file and adds written source contribution (from given time
/// frame) to the mesh provided in the region asked.
// ---------------------------------------------------------------------------
static void
TFSF_addSource (tfsfFace_t *face, meshVecI_p F, const reg_t *toApply)
{
    /// \todo Strides in unrolled array -> mesh/reg/wrap.
    const int stride[3] = {
        mc_have_x*F->width_yz,
        mc_have_y*F->width_z,
        mc_have_z
    };

    /*  SAY_DEBUG ("Adding sources to mesh %s / subreg %s", F->name, reg_printRanges (toApply));
    */
    for (const tfsfFace_t * const end = face + 12 ; face < end ; ++face) {
        if (face->fp) {
            const reg_t *cache = &face->toCache;
            const int c  = face->component,
                      ep = face->ep,
                      eq = face->eq,
                      ea = 3 - ep - eq;

            // Skips if plane of cache is not in toApply.
            if (cache->min[ea] > toApply->max[ea] ||
                cache->max[ea] < toApply->min[ea])
                continue;

            // Gets subset of cached region.
            const int p1 = (cache->min[ep] < toApply->min[ep]) ? toApply->min[ep] : cache->min[ep];
            const int p2 = (cache->max[ep] > toApply->max[ep]) ? toApply->max[ep] : cache->max[ep];
            const int q1 = (cache->min[eq] < toApply->min[eq]) ? toApply->min[eq] : cache->min[eq];
            const int q2 = (cache->max[eq] > toApply->max[eq]) ? toApply->max[eq] : cache->max[eq];

            // Skips if there are no overlap.
            if (p1 > p2 || q1 > q2)
                continue;

/*      SAY_DEBUG ("  o %s (%s) -> %c = %d, %d <= %c <= %d, %d <= %c <= %d ...", face->name, reg_printRanges (cache),
                'X' + ea, cache->min[ea], p1, 'X' + ep, p2, q1, 'X' + eq, q2);
*/
            // Gets to the point (axis,p1,q1).
            int pos0 = cache->min[ea]*stride[ea] + p1*stride[ep] +
                                                   q1*stride[eq];
            const int dQ = cache->max[eq] - cache->min[eq] + 1;
            double *buffer0 = face->buffer + (q1 - cache->min[eq]) +
                                             (p1 - cache->min[ep])*dQ;
            for (int p = p1 ; p <= p2 ; ++p, pos0 += stride[ep], buffer0 += dQ){
                int pos = pos0;
                double *buffer = buffer0;
                for (int q = q1 ; q <= q2 ; ++q, pos += stride[eq], ++buffer)
                    mv_unrl(F, pos).r[c] += *buffer;
            }
        }
    }
}

// ---------------------------------------------------------------------------
/// Gets data from file and adds inverted written source contribution (from
/// given time frame) to the mesh provided in the region asked.
// ---------------------------------------------------------------------------
static void
TFSF_subSource (tfsfFace_t *face, meshVecI_p F, const reg_t *toApply)
{
    /// \todo Strides in unrolled array -> mesh/reg/wrap.
    const int stride[3] = {
        mc_have_x*F->width_yz,
        mc_have_y*F->width_z,
        mc_have_z
    };

    for (const tfsfFace_t * const end = face + 12 ; face < end ; ++face) {
        if (face->fp) {
            // Alias for shortening sources below.
            const reg_t *cache = &face->toCache;
            const int c  = face->component,
                      ep = face->ep,
                      eq = face->eq,
                      ea = 3 - ep - eq;

            // Skips if plane of cache is not in toApply.
            if (cache->min[ea] > toApply->max[ea] ||
                cache->max[ea] < toApply->min[ea])
                continue;

            // Gets subset of cached region.
            const int p1 = (cache->min[ep] < toApply->min[ep]) ? toApply->min[ep] : cache->min[ep];
            const int p2 = (cache->max[ep] > toApply->max[ep]) ? toApply->max[ep] : cache->max[ep];
            const int q1 = (cache->min[eq] < toApply->min[eq]) ? toApply->min[eq] : cache->min[eq];
            const int q2 = (cache->max[eq] > toApply->max[eq]) ? toApply->max[eq] : cache->max[eq];

            // Skips if there are no overlap.
            if (p1 > p2 || q1 > q2)
                continue;

            // Gets to the point (axis,p1,q1).
            int pos0 = cache->min[ea]*stride[ea] +
                                   p1*stride[ep] +
                                   q1*stride[eq];
            const int dQ = cache->max[eq] - cache->min[eq] + 1;
            double *buffer0 = face->buffer + (q1 - cache->min[eq]) +
                                             (p1 - cache->min[ep])*dQ;
            for (int p = p1 ; p <= p2 ; ++p, pos0 += stride[ep], buffer0 += dQ){
                int pos = pos0;
                double *buffer = buffer0;
                for (int q = q1 ; q <= q2 ; ++q, pos += stride[eq], ++buffer)
                    mv_unrl(F, pos).r[c] -= *buffer;
            }
        }
    }
}

// ---------------------------------------------------------------------------
/// \brief Cached data from disk into tfsfFace_t::buffer using toCache region
/// for unrolling.
/// Many updates may be requested so caching should be scheduled by client.
// ---------------------------------------------------------------------------
static void
TFSF_cacheFrame (void)
{
    if (tfsfInterface.barcode != mc_TFSF_playForward &&
        tfsfInterface.barcode != mc_TFSF_playBackward)
        return;

    // Checks if frame exists.
    if (tfsfFrame < tfsfLastFrame - 1) {
        if (tfsfInterface.barcode == mc_TFSF_playForward)
            TFSF_cacheSource (tfsfH, tfsfFrame);
        else
            TFSF_cacheSource (tfsfH, tfsfLastFrame - 1 - tfsfFrame);
    }

    // Checks if frame exists.
    if (tfsfFrame + 1 < tfsfLastFrame - 1) {
        if (tfsfInterface.barcode == mc_TFSF_playForward)
            TFSF_cacheSource (tfsfE, tfsfFrame + 1);
        else
            TFSF_cacheSource (tfsfE, tfsfLastFrame - 1 - tfsfFrame);
    }

    /*  SAY_DEBUG ("Frame %d is cached.\n", tfsfFrame);
    */
    tfsfCachedFrame = tfsfFrame;
}

// ---------------------------------------------------------------------------
/// Fixes <b>external H</b> by removing contribution of the <b>source E</b> into
/// external region nodes.
// ---------------------------------------------------------------------------
void
TFSF_postHStep (meshVecI_p H, const reg_t *toApply)
{
    if (tfsfInterface.barcode != mc_TFSF_playForward &&
        tfsfInterface.barcode != mc_TFSF_playBackward)
        return;

    // Checks if frame exists.
    if (tfsfFrame >= tfsfLastFrame - 1)
        return;

    ENSURE (tfsfFrame == tfsfCachedFrame, "frame is not cached");

    TFSF_addSource (tfsfH, H, toApply);
}

// ---------------------------------------------------------------------------
/// Adds to <b>E on the boundary</b> of interface solution incoming from outside
/// proportional to <b>source H</b>.
// ---------------------------------------------------------------------------
void
TFSF_postEStep (meshVecI_p E, const reg_t *toApply)
{
    if (tfsfInterface.barcode != mc_TFSF_playForward &&
        tfsfInterface.barcode != mc_TFSF_playBackward)
        return;

    // Checks if frame exists.
    if (tfsfFrame + 1 >= tfsfLastFrame - 1) {
        tfsfInterface.barcode = mc_TFSF_inactive;
        return;
    }

    ENSURE (tfsfFrame == tfsfCachedFrame, "frame is not cached");

    if (tfsfInterface.barcode == mc_TFSF_playForward)
        TFSF_addSource (tfsfE, E, toApply);
    else
        TFSF_subSource (tfsfE, E, toApply);
}

// ---------------------------------------------------------------------------
/// Increases internal counter (used to invalidate cache) and checks if frames
/// are out.
// ---------------------------------------------------------------------------
void
TFSF_completeFrame (void)
{
    ++tfsfFrame >= tfsfLastFrame - 1;
}

// ---------------------------------------------------------------------------
/// Closes all files and deactivates player.
// ---------------------------------------------------------------------------
static void
TFSF_shutDown (void)
{
    for (int f = 0 ; f < 12 ; ++f) {
        if (tfsfE[f].fp) {
            free (tfsfE[f].buffer);
            tfsfE[f].buffer = NULL;
            fclose (tfsfE[f].fp);
            tfsfE[f].fp = 0;
        }

        if (tfsfH[f].fp) {
            free (tfsfH[f].buffer);
            tfsfH[f].buffer = NULL;
            fclose (tfsfH[f].fp);
            tfsfH[f].fp = NULL;
        }
    }

    if (!cpu_here && tfsfInterface.barcode == mc_TFSF_rec) {
        FILE *fp = cfg_open (".EM_sources/TFSF/lastFrameTime.bin", "wb",
                             __func__);
        fwrite (&tfsfLastFrame, sizeof (int), 1, fp);
        fclose (fp);
        say ("TFSF_shutDown: %d frames written.", tfsfLastFrame);
    }
}
