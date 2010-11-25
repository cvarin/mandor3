/** \file em_mirror.h
  * Soft source: laser beam focused by parabolic mirror.
  * <pre>
  * CREDITS
  * =======
  * Mirror code was taken from Kostya Popov's SCPIC code under his permission.
  * Thanks to him we had the first production version, which is always a major
  * breakthrough. You can find his code there:
  *     http://www.phys.ualberta.ca/~kpopov/scpic.html.
  *
  * Code was improved, documented and optimized (see corresponding section
  * below). Basically, it became 2-8 times faster with minor imporvements in
  * accuracy.
  *
  * \warning The mirror source works correctly only for 3D for now.
  *
  * CONFIG FILE
  * ===========
  * Example of the config file entry is
  *
  * [focused laser]
  * @ 1.0e+22   Poynting vector at the best focus: positive => I [W/cm²] / negative => e*E/(m*omega*c).
  * @ 0.9       Mirror f-number (defined as F_mirror/D_mirror; the smaller the number, the tighter the focusing). The value has to be greater then 0.25.
  * @ 10.0      Focus X coordinate [micron].
  * @ 10.0      Focus Y coordinate [micron].
  * @ 10.0      Focus Z coordinate [micron].
  * @ 20.0      Pulse duration, [fs] (defined by half of the maximum intensity), the envelope is Gaussian.
  * @ 40.0      Time of pulse peak arrival to focus [fs].
  * @ 0.0       Phase shift.
  * @ 10.0      Ratio between the radius of the pulse incident onto mirror and radius of the mirror (-1 means incident wave is a plane wave).
  * @ +1        Direction of propagation (positive or negative).
  *
  * INDEXING OF MIRROR REFLECTION
  * =============================
  * Lets assume that we have (optionally staggered) mesh, node with index 'i',
  * and mirror with coordinate 'xm = j*h/2' (j is integer). What would be the
  * index of reflection of this node?
  *
  * Answer is derived as follow ('shift' is an integer number):
  *
  *   x  = i*h + shift*h/2  =>   i  = int(x/h  - shift/2)
  *   x' = 2*xm - x         =>   i' = int(x'/h - shift/2)
  *                                 = int(2*xm/h) - i - shift.
  *
  * SOFT SOURCE IMPLEMENTATION: HACKER'S POINT OF VIEW
  * ==================================================
  * MESH LAYOUT
  *          Without source           With source           Without source
  *       ___________________^^^^^^^^^^^^^^^^^^^^^^^^^^^^^___________________
  *      ------x------|------x------|------x------|------x------|------x------
  *  Node:             I-1            I            I+1           I+2
  * Field:                   Hy[I]  Ey[I]
  *                          Hz[I]  Ez[I]
  *
  * Legend: region '_' : background solution only
  *         region '^' : background solution plus injected wave
  *
  * The result: we keep known solution of Maxwell equation (incoming focused
  * laser beam) only in part of the domain (marked '^'); nevertheless, any
  * other solution ('_') lives in the entire domain like nothing happend.

  * How do we compensate the dicontinuity in node I?
  *
  * Lets take a look on 'em.c::em_HSmallStep' and 'em.c::em_EStep_start':
  *    Hy[i,n+1/2] = Hy[i,n-1/2] + c1*(Ez[i,  n]     - Ez[i-1,n])
  *    Ez[i,n+1]   = Ez[i,n]     + c1*(Hy[i+1,n+1/2] - Hy[i,  n+1/2])
  *
  * Lets mark two EM waves as Hy_, Ez_, Hy^, Ez^, and drop index for clarity.
  * Both waves are solutions of the system, so
  *    H_[i,n+1/2] = H_[i,n-1/2] + c1*(E_[i,  n]     - E_[i-1,n])          (1a)
  *    E_[i,n+1]   = E_[i,n]     + c1*(H_[i+1,n+1/2] - H_[i,  n+1/2])
  *
  *    H^[i,n+1/2] = H^[i,n-1/2] + c1*E^[i,  n]     - E^[i-1,n])           (1b)
  *    E^[i,n+1]   = E^[i,n]     + c1*H^[i+1,n+1/2] - H^[i,  n+1/2])
  *
  * We have discontinyity in node I, so we add sources to compentsate it:
  *    H*[I,  n+½] = H*[I,n-½] + c1*(E*[I,n]   - E*[I-1,n])   =
  *                = H*[I,n-½] + c1*(E*[I,n]   - E_[I-1,n])   - c1*E^[I-1,n]
  *
  *    E_[I-1,n+1] = E_[I-1,n] + c1*(H_[I,n+½] - H_[I-1,n+½]) =
  *                = E_[I-1,n] + c1*(H*[I,n+½] - H_[I-1,n+½]) - c1*H^[I,n+½]
  *    ~~~~~~~~~~~   ~~~~~~~~~       ~~~~~~~~~   ~~~~~~~~~~~
  * where H* = H_ + H^, Ey* = Ey_ + Ey^.
  *
  * Please note that '~'-undescored values are already stored on the mesh, and
  * we add sources only to post-correct one of the terms in 'c1*(..)' to match
  * the other one to look like (1a) or (1b).
  *
  * Back to indexed form, full timestep now looks like (see site for the rest)
  * +--------------------------------------------+--------------------------+
  * | Hy[I]   +=  c1*(Ez*[I] - Ez[I-1])          | Yee timestep, part 1     |
  * | Hy[I]   += -c1*Ez^(Δx·(I-1), τ·n)          | Source of magnetic field |
  * | Ez[I-1] +=  c1*(Hy[I] - Hy[I-1])           | Yee timestep, part 2     |
  * | Ez[I-1] += -c1*Hy^(Δx·(I-1/2), τ·(n+1/2))  | Source of electric field |
  * +--------------------------------------------+--------------------------+
  * where E^(i, n) and H^(i, n) are KNOWN FUNCTIONS of time and space.
  *
  * \warning It is important to update magnetic field just after the first Yee
  *          step because of it is used immediately in second Yee step.
  *
  * \warning All signs of the source are included into the precomputed field!
  *
  *
  * MAIN IMPROVEMENTS COMPARED TO ORIGINAL (BORROWED) CODE
  * ======================================================
  * 1) Twice faster Stratton-Chu integration due to simultaneous computing of
  *    real and imaginary parts of the integral (reuses ALL heavy functions).
  *
  * 2) Account for (anti)symmetry of field distribution allows to compute only
  *    independent values of E and H (another 2-4 times integration speedup).
  *
  * 3) Analitical calibration of Ey amplitude in focal spot. Proof is on site,
  *    but key point is that 'Ey' depends on phase shift as 'exp(i*τ)', so we
  *    get exact amplitude and phase from just two calls to 'focused_ey'.
  *    Saves about 5-10 seconds.
  *
  * 4) Exact symmetry between waves ejected from left and right boundaries.
  *
  *    Now on my box the integration takes 16.2 seconds instead of 156.4.
  * </pre>
  */

#ifndef EM_MIRROR_HEADER
#define EM_MIRROR_HEADER

#include "type_mesh.h"

/// Complete set of tangent fields in X plane (enough to define incoming wave).
typedef struct {
   double y_cos,
          y_sin,
          z_cos,
          z_sin;
} wave_x_t;

/// Parameters of the parabolic mirror.
typedef struct {
   double    A;		// Amplitude of field.
   double    f0;	// Focal length.
   double    rm;	// Mirror radius.
   int       NN;	// Number of integration points.
   double    rratio;	// Ratio between the radius of the incident pulse and
                        // the mirror.
   double    xf;	// Coordinates of focus: X
   double    yf;	//                       Y
   double    zf;	//                       Z.
   double    dt;	// Pulse duration.
   double    T0;	// Time of arrival of the pulse pike to the focus.
   double    delay;	// Time to travel from the source to the focus.
   double    phi0;	// Phase shift.
   int       X_forward;	// Direction of propagation (±1).
   wave_x_t *E,		// Wave info (phase and amplitude).
            *H;
} mirror_t;

mirror_t  read_mirror_parameters  (FILE *fp);
void      print_mirror_parameters (mirror_t src);
void      save_mirror_source      (int mirror_number, mirror_t mirror);
int       mirror_source_is_loaded (int mirror_number, mirror_t *dest);

void add_mirror_H (mirror_t *src, meshVec_p H, double time, reg_t *to_update);
void add_mirror_E (mirror_t *src, meshVec_p E, double time, reg_t *to_update);

#endif
