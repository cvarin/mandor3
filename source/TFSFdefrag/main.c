/** \file main.c
  * \brief 'Total Field / Scattered Field' defragmentation module.
  *
  * Source data on TFSF interface are written in parallel. Real pass will use
  * other set of cpus so to reduce code complexity the data is compiled into
  * final file which makes replay code simplier.
  *
  * Here 'frame' means the one side of the interface for one component (field
  * components are allocated on shifted meshes).
  */

#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "type_reg.h"

#include "log.h"
#include "misc_cfgReader.h"

#define MC_FILES_MAX 	(3000)		///< Max number of opened file streams.

static int  frames_to_skip = 0;

// Data outside of the filter is neglected (that works like open face but
// allows me to merge many sources written in separate passes with different
// filters).
static reg_t  filter = {
    .min = {INT_MIN, INT_MIN, INT_MIN},
    .max = {INT_MAX, INT_MAX, INT_MAX},
};

// That is final result.
static FILE*   out_file;		// Defragmented result file.
static reg_t   face;			// Face region.
static int     ep, eq;			// Face orientation (unrolling order).
static double *full = NULL;		// Assembled frame.

// That is pieces of result written by separate parallel run.
static FILE*   files[MC_FILES_MAX];	// Inpute files (fragments).
static char*   names[MC_FILES_MAX];	// Names of files.
static int     piecesN = 0;		// Number of streams to merge.
static reg_t  *pieces  = NULL;		// Pieces assigned to src cpus.

static int     piece_size = 0;		// Input buffer.
static double *piece      = NULL;

/**
  * Reads single word from file (including special chars).
  */
static char*
read_word (FILE *fp)
{
    static char name[1000];
    fscanf (fp, "%1000[^ \t\n\r] ", name);
    name[sizeof (name) - 1] = 0;
    ENSURE (strlen (name) < sizeof (name) - 1, "too long name");
    return name;
}

/**
  * Prepares information for defragmentation pass: regions, face, files.
  */
static int
files_in_queue (FILE *fp)
{
    // Gets face frame orientation.
    fscanf (fp, "ep = %d, eq = %d\n", &ep, &eq);

    ENSURE (!pieces && !piecesN && !out_file, "unfinished previous job");

    // Opens file for updating.
    const char *root_name = read_word (fp);
    out_file = cfg_open (root_name,
                         access (root_name, F_OK) ? "w+b" : "r+b",
                         __func__);
    fwrite (&face, sizeof (reg_t), 1, out_file);

    say ("Processing face '%s' ...", reg_printRanges (&face));

    // Opens input file streams and get decomposition geometry.
    while (1) {
        char *filename = read_word (fp);

        if (filename[0] == '-')
            break;

        files[piecesN] = cfg_open (filename, "rb", __func__);
        names[piecesN] = strdup (filename);
        ++piecesN;

        SAY_DEBUG ("  - adding file '%s' to list...", filename);
    }

    if (piecesN) {
        SAY_DEBUG ("Allocating geometry (%d subregions) ...", piecesN);
        SAY_DEBUG ("Unrolling order: ep = %d, eq = %d\n", ep, eq);

        piece_size = 0;
        pieces     = (reg_t*) malloc (piecesN*sizeof (reg_t));
        for (int f = 0 ; f < piecesN ; ++f) {
            fread (pieces + f, sizeof (reg_t), 1, files[f]);
            piece_size = fmax (piece_size, reg_volume (pieces + f)) + 0.1;
            SAY_DEBUG ("    + %s", reg_printRanges (pieces + f));
        }

        full  = (double*) realloc (full,  reg_volume (&face)*sizeof (double));
        piece = (double*) realloc (piece, piece_size*sizeof (double));

        SAY_DEBUG ("buffers: %d ≤ %d", piece_size, reg_volume (&face));
    } else {
        SAY_DEBUG ("  no files listed - looks like degenerated axis at work.");
        fclose (out_file);
        out_file = NULL;
    }

    return piecesN;
}


/**
  * Inserts 'frames_skip' empty frames in the start of file if it is empty.
  */
static int
pad_with_empty_frames (int frames_to_skip)
{
    ENSURE (frames_to_skip >= 0, "unsupported shift (%d)", frames_to_skip);

    fseek (out_file, 0L, SEEK_END);
    long int tail    = ftell (out_file);
    size_t   to_skip = frames_to_skip*reg_volume (&face)*sizeof (double);
    if (tail < sizeof (reg_t) + to_skip) {
        fseek (out_file, 0L, SEEK_SET);
        fwrite (&face, sizeof (face), 1, out_file);
        memset (full, 0, reg_volume (&face)*sizeof (double));
        while (frames_to_skip-- >= 0)
            fwrite (full, sizeof (double), reg_volume (&face), out_file);
    }
    return (tail - sizeof (reg_t))/(reg_volume (&face)*sizeof (double));
}


/**
  * Joins all files to assemble a complete face and writes the frame.
  */
static int
defrag_join_files (int frameN)
{
    long int nodes = reg_volume (&face);
    fseek (out_file,
           sizeof (reg_t) + frameN*sizeof (double)*nodes,
           SEEK_SET);
    size_t got = fread (full, sizeof (double), nodes, out_file);
    if (got == 0) {
        memset (full, 0, nodes*sizeof (double));
    } else {
        ENSURE (got == nodes, "incomplete destination file (%d < %d)",
                                 got, nodes);
    }
    fseek (out_file,
           sizeof (reg_t) + frameN*sizeof (double)*nodes,
           SEEK_SET);

    memset (piece, 0, piece_size*sizeof (double));

    int got_data = 0;
    for (int f = 0 ; f < piecesN ; ++f) {
        // Buffers the input fragment.
        reg_t *reg = pieces + f;
        size_t got = fread (piece, sizeof (double), reg_volume (reg), files[f]);
        if (got == 0)
            continue;
        ENSURE (got == reg_volume (reg), "incomplete fragment file");
        got_data = 1;

        // Puts the piece into its place.
        int out_span = face.max[eq] - face.min[eq] + 1;
        int in_span  = reg->max[eq] - reg->min[eq] + 1,
           p0 = (reg->min[ep] < filter.min[ep]) ? filter.min[ep] : reg->min[ep],
           p1 = (reg->max[ep] > filter.max[ep]) ? filter.max[ep] : reg->max[ep],
           q0 = (reg->min[eq] < filter.min[eq]) ? filter.min[eq] : reg->min[eq],
           q1 = (reg->max[eq] > filter.max[eq]) ? filter.max[eq] : reg->max[eq],
           er = 3 - ep - eq,
           r  = reg->min[er];
        if (r < filter.min[er] || r > filter.max[er])
            continue;
        for (int p = p0 ; p <= p1 ; ++p)
            for (int q = q0 ; q <= q1 ; ++q)
                full[q - face.min[eq] + out_span*(p - face.min[ep])] +=
                      piece[q - reg->min[eq] + in_span*(p - reg->min[ep])];
    }

    if (got_data) {
        fwrite (full, sizeof (double), reg_volume (&face), out_file);
    }

    return got_data;
}

/**
  * Main function.
  */
int
main (int argc, char *argv[])
{
    log_open (0, 0, "output/defrag.log");

    FILE *fp = cfg_open (".EM_sources/TFSF/filter.cfg", "rt", __func__);
    fscanf (fp, "@ %d       I0.\n", &(filter.min[0]));
    fscanf (fp, "@ %d       I1.\n", &(filter.max[0]));
    fscanf (fp, "@ %d       ", &frames_to_skip);
    fclose (fp);

    int frames_written;
    fp = cfg_open (".EM_sources/TFSF/lastFrameTime.bin", "rb", __func__);
    fread (&frames_written, sizeof (int), 1, fp);
    fclose (fp);
    say ("Incoming frames: %d", frames_written);

    say ("Filter: [%d %d]", filter.min[0], filter.max[0]);
    say ("Shift:  %d",      frames_to_skip);

    // Opens list of jobs.
    fp = fopen (".EM_sources/TFSF/TFSF_defrag.map", "rt");
    if (!fp) {
        SAY_DEBUG ("File '.EM_sources/TFSF/TFSF_defrag.map' is removed, "
                   "no defragmentation necessary.");
        return 0;
    }

    while (1) {
        while (1) {
            // Gets face region.
            int read = fscanf (fp, "[%d %d] [%d %d] [%d %d]\n",
                                   &face.min[0], &face.max[0],
                                   &face.min[1], &face.max[1],
                                   &face.min[2], &face.max[2]);
            if (read != 6)
                goto exit;

            if (files_in_queue (fp))
                break;
        }

        int frames_in_file = pad_with_empty_frames (frames_to_skip);
        frames_written = fmax (frames_written, frames_in_file) + 0.1;

        // Defrags file streams.
        int frames = 0;
        while (defrag_join_files (frames_to_skip + frames))
            frames++;
        SAY_DEBUG ("    %d frames compiled", frames);
        frames_written = fmax (frames + frames_to_skip, frames_written) + 0.1;

        // Frees all used resources.
        for (int f = 0 ; f < piecesN ; ++f) {
            fclose (files[f]);
            remove (names[f]);
            free   (names[f]);
        }

        fclose (out_file);
        out_file = NULL;

        free (pieces);
        pieces  = NULL;
        piecesN = 0;
    }
exit:
    fclose (fp);

    // Removes list of jobs to signal the end of defragmentation.
    remove (".EM_sources/TFSF/TFSF_defrag.map");

    fp = cfg_open (".EM_sources/TFSF/lastFrameTime.bin", "wt", __func__);
    fwrite (&frames_written, sizeof (int), 1, fp);
    fclose (fp);

    say ("Done (%d frames total)", frames_written);

    return EXIT_SUCCESS;
}
