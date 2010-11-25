/** \file tag.h
  * \brief Predefined tag ID for switch in the main.c of the setup.out project.
  */

#ifndef MC_TAG_HEADER
#define MC_TAG_HEADER

#include <stdio.h>

// Bunch of prototypes.
#include "setup/tag_units.h"
#include "setup/tag_mesh.h"
#include "setup/tag_boundary.h"

#include "setup/tag_TFSF.h"
#include "setup/tag_TFSFFaces.h"
#include "setup/tag_laser.h"

#include "setup/tag_plasma.h"
#include "setup/tag_DF_uniform.h"

#include "setup/tag_point.h"
#include "setup/tag_EMWave.h"
#include "setup/tag_ringDF.h"
#include "setup/tag_vShift.h"
#include "setup/tag_photoelectrons.h"
#include "setup/tag_jBalancer.h"
#include "setup/tag_twoStream.h"
#include "setup/tag_gaussSpot.h"
#include "setup/tag_maxwell.h"
#include "setup/tag_planePulse.h"
#include "setup/tag_plasmaWave.h"
#include "setup/tag_EMResonator.h"
#include "setup/tag_scissors.h"
#include "setup/tag_meanVNoise.h"
#include "setup/tag_seedWeibel.h"
#include "setup/tag_seed_PITS.h"
#include "setup/tag_cluster.h"
#include "setup/tag_foil.h"
#include "setup/tag_scales.h"
#include "setup/tag_trianglePrizm.h"

// Bunch of unique IDs for switch statement in 'setup/main.c'.

enum {
    TAG_UNITS = 0,		///< Sets units: tag_units.h.
    TAG_MESH,			///< Sets size of the mesh and domain: tag_mesh.h.
    TAG_BOUNDARY,		///< Sets boundary conditions: tag_boundary.h.

    TAG_TFSF,			///< 'Total Field/Scattered Field' interface configuration: tag_TFSF.h.
    TAG_TFSF_FACES,		///< Sets opening marks for TFSF players: tag_TFSFFaces.h.
    TAG_SRC_MIRROR,		///< Soft source of EM waves::beam focused by parabolic mirror: tag_laser.h.

    TAG_EM_POINT,		///< Adds point perturbation of field to tests solvers: tag_point.h.
    TAG_EM_WAVE,		///< Adds plane EM wave: tag_EMWave.h.
    TAG_EM_GAUSS_SPOT,		///< Sets parameters of the gaussian hard source: tag_gaussSpot.h.
    TAG_EM_PLANE_PULSE,		///< Adds plane EM pulse to probe boundary conditions: tag_planePulse.h.
    TAG_EM_RESONATOR,		///< Excites standing wave in domain-resonator: tag_EMResonator.h.

    /// All this tag are going to be replaced by analytic setup in future.
    TAG_PLASMA,			///< Allocates arbitrary shaped piece of plasma with no velocity distribution: tag_plasma.h.
    TAG_DF_UNIFORM,		///< Sets uniform velocity distribution (to test BC implementations): tag_DF_uniform.h.
    TAG_TWO_STREAM,		///< Creates two cold beams to study two-stream instability: tag_twoStream.h.
    TAG_FOIL,			///< Foil with model (linear or exponential) charge density profile: tag_foil.h.
    TAG_PHOTOELECTRONS,		///< Allocates photo-ionization DF (for paper with V.U.Bychenkov): tag_photoelectrons.h.
    TAG_MAXWELL,		///< Uniform Maxwellian DF allocator: tag_maxwell.h.
    TAG_CLUSTER,		///< Cold cluster for Andrey and Valery: tag_cluster.h.
    TAG_TRIANGLE_PRIZM,		///< Cold plasma triangle prizm target: tag_trianglePrizm.h.
    TAG_RING_DF,		///< Creates model DF (on-going work): tag_ringDF.h.

    TAG_VELOCITY_SHIFT,		///< Shifts component of plasma or all particles in velocity space: tag_velocityShift.h.
    TAG_CURRENT_BALANCER,	///< Removes mean current to avoid zero mode plasma waves: tag_currentBalancer.h.
    TAG_MEAN_V_NOISE,		///< Adds mean velocity perturbation with white-noise spectrum in given subregion of k-space: tag_meanVNoise.h.
    TAG_PLASMA_WAVE,		///< Small kick to the mean velocity of the particles (dirty but fast plasma waves): tag_plasmaWave.h.
    TAG_SEED_PITS,		///< Precisely excites PiTS unstable mode: tag_seed_PITS.h.
    TAG_SEED_WEIBEL,		///< Sloppy excitation of PiW unstable mode: tag_seedWeibel.h.

    TAG_SCALES,			///< Prints scales of the simulation: tag_scales.h.
    TAG_SCISSORS,		///< Sets subregion where plasma must be created: tag_scissors.h.

    TAG_TOTAL_NUM,		///< Total number of tags (used as array size in tag.c).
    TAG_EOF         = -1,	///< Return code to signal the end of file.
    TAG_UNKNOWN_TAG = -2,	///< Return code to signal an unknown tag.
};

int 	    getTag (FILE *fp);
const char* getLastTagName (void);

#endif
