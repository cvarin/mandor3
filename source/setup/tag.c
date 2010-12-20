/** \file tag.c
  * \brief Scans config file looking for the first line with '[' in the beginning of line. Than line is analized and tag ID is returned.
  */

#include <stdlib.h>
#include <string.h>

#include "setup/tag.h"

#include "log.h"
#include "misc_cfgReader.h"

// ---------------------------------------------------------------------------
/// Array of tag names (used to map index <-> name).
// ---------------------------------------------------------------------------
static const char* tags[TAG_TOTAL_NUM] =
{
    [TAG_UNITS]            = "units",
    [TAG_MESH]             = "mesh",
    [TAG_BOUNDARY]         = "boundaryConditions",

    [TAG_TFSF]             = "TFSF",
    [TAG_TFSF_FACES]       = "TFSF: open faces",
    [TAG_SRC_MIRROR]       = "focused laser",

    [TAG_EM_POINT]         = "point",
    [TAG_EM_WAVE]          = "EMWave",
    [TAG_EM_GAUSS_SPOT]    = "gaussSpot",
    [TAG_EM_PLANE_PULSE]   = "planePulse",
    [TAG_EM_RESONATOR]     = "EMResonator",

    [TAG_PLASMA]           = "plasma",
    [TAG_DF_UNIFORM]       = "DF:uniform",
    [TAG_TWO_STREAM]       = "DF:two-stream",
    [TAG_FOIL]             = "DF:foil",
    [TAG_PHOTOELECTRONS]   = "DF:photoelectrons",
    [TAG_MAXWELL]          = "DF:Maxwell",

    [TAG_TRIANGLE_PRIZM]   = "trianglePrizm",		// XXX: remove
    [TAG_CLUSTER]          = "cluster",			// XXX: remove
    [TAG_RING_DF]          = "DF:ringDF",		// XXX: remove

    [TAG_VELOCITY_SHIFT]   = "velocity shift",
    [TAG_CURRENT_BALANCER] = "current balancer",
    [TAG_MEAN_V_NOISE]     = "<V> perturbation",
    [TAG_PLASMA_WAVE]      = "plasmaWave",
    [TAG_SEED_PITS]        = "seed:PiTS",
    [TAG_SEED_WEIBEL]      = "seed:Weibel",

    [TAG_SCALES]           = "show scales",
    [TAG_SCISSORS]         = "scissors",

    [TAG_GRADIENT]         = "gradient",
};

// Incoming buffer.
static char line[1000] = "[NO TAG NAME SCANNED YET!";

// ---------------------------------------------------------------------------
/// Returns id of the next tag in the config file (>= 0), '-1' to signal the end
/// of file or -2 to signal unknown tag.
// ---------------------------------------------------------------------------
int
getTag (FILE *fp)
{
    while (fp && !feof (fp)) {
        line[0] = line[999] = 0;
        fscanf (fp, " %1000[^\n]%*[\n] ", line);
        ENSURE (line[999] == 0, "too long line");

        // Tag line detected.
        if (line[0] == '[') {
            for (int i = 0 ; i < TAG_TOTAL_NUM ; ++i)
                if (strncmp (line + 1, tags[i], strlen (tags[i])) == 0
                    && line[1+strlen (tags[i])] == ']') {
                    return i;
                }
            return TAG_UNKNOWN_TAG;
        }
    }
    return TAG_EOF;
}

// ---------------------------------------------------------------------------
/// After 'get_tag' we have the name of a last tag in buffer ('line').
// ---------------------------------------------------------------------------
const char *
getLastTagName (void)
{
    return line;
}
