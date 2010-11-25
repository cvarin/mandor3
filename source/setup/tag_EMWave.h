/** \file tag_EMWave.h
  * \brief Theory and declarations for plane EM wave test.
  *
  * Uniform plane EM-wave is excited in the domain to test periodic boundary
  * conditions and dispersion properties of the scheme.
  * Wave parameters:
  * - wave vector (expressed in wave numbers \f$ m_x\f$, \f$ m_y\f$, \f$ m_z\f$);
  * - amplitude of the E field;
  * - polarisation vector and polarisation type;
  *
  * Example of the config file entry is
  * <pre>
  * [EMWave]
  * @ physical		Type of the initialization (physical or Yee_2nd_order).
  * @ +0		mx = k_x*L_x/2pi.
  * @ +1		my = k_y*L_y/2pi.
  * @ +0		mz = k_z*L_z/2pi.
  * @ 5.0		E0 (electric field amplitude).
  * @ -1.0e0		nx:
  * @ -1.0e0		ny:  components of vector n: (E,[n,k]) = 0
  * @ +1.0e0		nz:
  * @ plane		Polarization of the wave (plane, left or right).
  * </pre>
  *
  * <h3>Physical initialization.</h3>
  *
  * EM wave is a solution of the Maxwell equations in the form
  * \f[ \vec E = \vec E_0\cdot e^{i(\vec k \vec r - \omega t)} \f]
  * \f[ \vec H = \vec H_0\cdot e^{i(\vec k \vec r - \omega t)} \f]
  * Vectors \f$ \vec E_0 \f$ and \f$ \vec H_0 \f$ are complex and include all
  * initial phase shifts. Charge and current density are constantly zero.
  * Resulting divergency-free condition requires vectors \f$ \vec E_0 \f$ and
  * \f$ \vec H_0 \f$ to be perpendicular to the \f$ \vec k \f$ (they are also
  * perpendicular to each other and have equal magnitudes).
  *
  * Wave vector \f$ \vec k = \left(\frac{2\pi m_x}{L_x}, \frac{2\pi m_y}{L_y},
  * \frac{2\pi m_z}{L_z} \right) \f$ is given and frequency \f$ \omega \f$ is
  * derived from dispersion equation \f[ \omega^2 = k^2c^2. \f]
  * Using \f$ \vec k \f$ and polarization vector guess \f$ \vec n \f$, it is
  * possible to find vector \f$ \vec E_0 \f$ which rests in the plane
  * \f$ \vec k, \vec n \f$ and have given magnitude. Vector \f$ \vec H \f$ is
  * defined as \f$ \vec H = [\vec k, \vec E ] /\omega \f$.
  *
  * <h3>Numerical initialization.</h3>
  *
  * Basic steps are the same like for obtaining dispersion equation for EM
  * waves. Solution in form \f[ \vec E = \vec E_0\cdot e^{i(\vec k\vec r - \omega t)} \f]
  * \f[ \vec H = \vec H_0\cdot e^{i(\vec k\vec r - \omega t)} \f] is substituted
  * into finite difference system of equations (see core/em.c).
  *
  * Finite differencing introduces two more parameters:
  * \f[ \vec D_r = \left(\frac 2{h_1}\cdot \sin \left(\frac{k_xh_1}{2}\right),
  *                      \frac 2{h_2}\cdot \sin \left(\frac{k_yh_2}{2}\right),
  *                      \frac 2{h_3}\cdot \sin \left(\frac{k_zh_3}{2}\right) \right) \f]
  * \f[ D_t = \frac 2{\tau}\cdot \sin \left(\frac{\omega \tau}{2}\right) \f]
  *
  * Vector \f$ \vec D_r \f$ replaces \f$ \vec k \f$ in all differential
  * operators. For example, condition of the divergency-free field means
  * \f[ (\vec D_r, \vec E_0) = 0 \f]
  * And magnetic field is found as
  * \f[ \vec H = [\vec D_r, \vec E ] / D_t \f]
  *
  * Now \f$ \omega \f$ and \f$ \vec k \f$ obey numerical dispersion equation:
  * \f[ \frac 4{\tau^2} \sin^2 (0.5\cdot \omega\tau) = \frac 4{h_1^2} \sin^2 (0.5\cdot k_xh_1) +
  *     \frac 4{h_2^2} \sin^2 (0.5\cdot k_yh_2) + \frac 4{h_3^2} \sin^2 (0.5\cdot k_zh_3)\f]
  *
  * <h3>Conclusion.</h3>
  *
  * This tag helps to estimate dispersion introduced by Yee solver for given
  * EM wave.
  *
  * Main points to note:
  * - Solver conserves energy of each field component with numerical accuracy.
  * - Energy of magnetic field is smaller than energy of electric field.
  * - Sin/cos wave doesn't change shape (linear solver has no nonlinearity due
  *   to implementation).
  * - Mesh acts like crystall with phase velocity dependent on angle with
  *   respect to axises.
  *
  * \sa tag_EMWave.c.
  */

#ifndef TAG_EMWAVE_HEADER
#define TAG_EMWAVE_HEADER

#include "type_mesh.h"

void tag_EMWave (FILE *fp, meshVec_p E, meshVec_p H);

#endif
