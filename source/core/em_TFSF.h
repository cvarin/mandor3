/** \file em_TFSF.h
  * \brief Theory of "Total Field / Scattered Field" interface (part of Yee Maxwell solver).
  * This file is used by \b setup.out and \b core.out for consistency of all predefined constants. Data structures are in em_TFSF.c.
  *
  * <b> TF/SF interface introduction.</b>
  *
  * For nice description of the TF/SF interface, see <b>Taflove, Allen. "Computational electrodynamics:
  * the finite-difference time-domain method" // Boston: Artech House. 1995. 17 chapters, 599 p</b>.
  * Nevertheless, I arrived to the same method using another way which may be easier for coders.
  *
  * <b>1. Intro and problem description.</b>
  *
  * Problem description: a known solution of Maxwell equations should be added inside of the subregion
  * of computational domain without altering either this solution or any already existing fields.
  * Solution to introduce is known completely (it will be called \b source). \b Source should
  * land <b>without any single assignment to the mesh fields</b> - it is absolutely necessary in order
  * to keep all waves in the system uncorrupted. Any assignment will turn assignment point into mirror
  * which reflects all incident waves (\b hard sources in terms of Taflove). Only addition of
  * something is permitted (specific sources).
  *
  * <b>2. Solution of the issue, or how Yee solver really works.</b>
  *
  * Instead of thinking about solution of the Maxwell equation in terms of waves and potentials lets treat
  * all nodes as small generator of signal which they deliver to the neighbours. If all this inputs may
  * cancel each other to represent propagating wave packet or field of the point charge - who cares? Point
  * is that at each time-step each node which, for example, holds field component \f$ Ex_{\ i-1/2,\, j,\, k}\f$
  * does only four basic linear operations:
  * - Add \f$ +\tau/dz\cdot Ex \f$ to the \f$ Hy_{\ i-1/2,\, j,\, k+1/2} \f$.
  * - Add \f$ -\tau/dz\cdot Ex \f$ to the \f$ Hy_{\ i-1/2,\, j,\, k-1/2} \f$.
  * - Add \f$ +\tau/dy\cdot Ex \f$ to the \f$ Hz_{\ i-1/2,\, j-1/2,\, k} \f$.
  * - Add \f$ -\tau/dy\cdot Ex \f$ to the \f$ Hz_{\ i-1/2,\, j+1/2,\, k} \f$.
  *
  * That means that to safely remove \b source solution from external region two additional steps must
  * be performed:
  * - Removing of the contribution to the external region made by \b source solution in the boundary nodes
  *   (transmission coefficients for Ex are shown as example).
  * - Addition of the contribution to the boundary nodes which should normally be done by external nodes
  *   with \b source term presented.
  *
  * Once \b source is known exactly, any of this contributions can be founded and used (technically it's an
  * additional source in RHS).
  *
  * <b>3. Formal theoretical derivation.</b>
  *
  * <b>Problem formulation in 1D.</b> System of Maxwell equations is solved using Yee scheme and solver is used to study scattering of
  * waves on the object. Solution which represents incident wave is known completely, target is surrounded by convex region. What is the method to
  * support solution \b f<sub>Incident</sub> + \b f<sub>Scatter</sub> inside of the region and keep only \b f<sub>Scatter</sub>
  * solution outside of it without changing solutions at all? Term <b>"Without changing"</b> means that for any shape of region the results inside
  * or outside of the region are the same.
  *
  * <b>Solution.</b> Lets maintain solution \b f<sub>Incident</sub> discontinuity on the left boundary with index \f$ i_0 \f$. That means that on the left
  * only solution \b f<sub>Scatter</sub> remains while on the right both solutions stay unmodified (like in case of infinite domain and no discontonuity -
  * just wave coming from infinity (\b f<sub>Incident</sub>) which scatters on the target to transform part of energy into \b f<sub>Scatter</sub>).
  *
  * Before any Yee time step solution in mesh nodes is:
  * \f[
  * \left\{
  * \begin{array}{ll}
  * E{\,}_{i_0}^{n} = E_S{\,}_{i_0}^{n} + E_I{\,}_{i_0}^{n} & E{\,}_{i_0-1}^{n} = E_S{\,}_{i_0}^{n}\\[4mm]
  * H{\,}_{i_0+\frac 12}^{n} = H_S{\,}_{i_0+\frac 12}^{n} + H_I{\,}_{i_0+\frac 12}^{n} \qquad & H{\,}_{i_0-\frac 12}^{n} = H_S{\,}_{i_0/2}^{n}
  * \end{array}
  * \right.
  * \f]
  *
  * Time step is done using Yee scheme:
  * \f[
  * \left\{
  * \begin{array}{l}
  * \displaystyle \frac{E{\,}_{i}^{n+1} - E{\,}_{i}^{n}}\tau =
  * \frac{H{\,}_{i+\frac 12}^{n+\frac 12} - H{\,}_{i-\frac 12}^{n+\frac 12}}{h} \\[4mm]
  *   \displaystyle \frac{H{\,}_{i-\frac 12}^{n+\frac 12} - H{\,}_{i-\frac 12}^{n-\frac 12}}\tau =
  *   \frac{E{\,}_{i}^{n} - E{\,}_{i-1}^{n}}h
  * \end{array}
  * \right.
  * \f]
  *
  * Lets put full solution \f$ E = E_S + E_I \f$, \f$ H = H_S + H_I \f$ into system and regroup terms:
  * \f[
  * \frac{(E_S + E_I){\,}_{i_0}^{n+1} - (E_S + E_I){\,}_{i_0}^{n}}\tau = \frac{(H_S + H_I){\,}_{i_0+\frac 12}^{n+\frac 12} - (H_S + H_I) {\,}_{i_0-\frac 12}^{n+\frac 12}}{h} =
  * \frac{E{\,}_{i_0}^{n+1} - E{\,}_{i_0}^{n}}\tau = \frac{H{\,}_{i_0+\frac 12}^{n+\frac 12} - H{\,}_{i_0-\frac 12}^{n+\frac 12}}{h} - \frac 1h H_I {\,}_{i_0-\frac 12}^{n+\frac 12}
  * \f]
  * \f[
  * \frac{H_S{\,}_{i_0-\frac 12}^{n+\frac 12} - H_S{\,}_{i_0-\frac 12}^{n-\frac 12}}\tau = \frac{E_S{\,}_{i_0}^{n} - E_S{\,}_{i_0-1}^{n}}h =
  * \frac{H{\,}_{i_0-\frac 12}^{n+\frac 12} - H{\,}_{i_0-\frac 12}^{n-\frac 12}}\tau = \frac{E{\,}_{i_0}^{n} - E{\,}_{i_0-1}^{n}}h - \frac1h E_I{\,}_{i_0}^{n}
  * \f]
  *
  * Final expression may be interpreted as usual Yee time step (with values stored on the mesh) plus additional sources in the right hand side. This sources
  * maintain discontinuity of \b f<sub>Incident</sub> in requested form without corrupting it: at the end of time step the same pattern
  * \f[
  * \left\{
  * \begin{array}{ll}
  * E{\,}_{i_0}^{n} = E_S{\,}_{i_0}^{n} + E_I{\,}_{i_0}^{n} & E{\,}_{i_0-1}^{n} = E_S{\,}_{i_0}^{n}\\[4mm]
  * H{\,}_{i_0+\frac 12}^{n} = H_S{\,}_{i_0+\frac 12}^{n} + H_I{\,}_{i_0+\frac 12}^{n} \qquad & H{\,}_{i_0-\frac 12}^{n} = H_S{\,}_{i_0/2}^{n}
  * \end{array}
  * \right.
  * \f]
  * is reproduced for the next time index \b n.
  */

#ifndef EM_TFSF_HEADER
#define EM_TFSF_HEADER

#include "type_mesh.h"

#ifndef SGI
  // My Linux doesn't support files bigger 2 Gb.
  #define ftell64 ftell						///< Alias to map ftell64 -> ftell on Linux boxes.
  #define fseek64 fseek						///< Alias to map fseek64 -> fseek on Linux boxes.
#endif

void TFSF_init          (void);
void TFSF_preHStep      (meshVecI_p E, meshVecI_p H);
void TFSF_postHStep     (meshVecI_p H, const reg_t *toApply);
void TFSF_postEStep     (meshVecI_p E, const reg_t *toApply);
void TFSF_completeFrame (void);

/*
 * Possible regimes of the interface. Used by setup for consistency.
 */
#define mc_TFSF_playForward	0	///< <b>Playing forward</b> regime.
#define mc_TFSF_rec		1	///< <b>Recording source</b> regime.
#define mc_TFSF_playBackward	2	///< <b>Playing backward</b> regime.
#define mc_TFSF_inactive	3	///< Internal state "deactivated explicitly".
#define mc_TFSF_uninitialized	4	///< Internal state "uninitialized".

#endif
