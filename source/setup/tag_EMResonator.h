/** \file tag_EMResonator.h
  * Theory and headers for standing EM wave test.
  *
  * <h3>Test of the mirror boundary conditions and dispertion properties of the Yee scheme.</h3>
  *
  * Physically solution of the Maxwell equations for the rectangular domain with mirror walls is standing wave.
  * Polarization of basic solution is one non-zero component of electric field and two (perpendicular) components
  * of magnetic field changing in the plane perpendicular to the electric field as shown below for X mode:
  * \f[ \vec E = E_{x\,0} \vec e_x \sin \left(\frac{\pi y m_y}{L_y}\right) \sin \left(\frac{\pi z m_z}{L_z}\right) \sin (\omega t) \f]
  *
  * Example of the config file entry is
    <pre>
    [EMResonator]
    @ physical		Type of the dispersion equation solved ('physical' or 'Yee_2nd_order').
    @ +0		mx = k_x*L_x/pi.
    @ +1		my = k_y*L_y/pi.
    @ +0		mz = k_z*L_z/pi.
    @ Ex		Polarization of the wave (non-zero component of the electric field).
    @ 5.0		E0 (electric field amplitude).
    </pre>
  *
  * Dispersion properties are described in details in the tag_EMWave.h file.
  *
  * Notes:
  * - Component of electric field is constant along shift direction so it is not necessary to account for '-0.5*h1' shift for X-mode,
  *   '-0.5*h2' shift for Y-mode, etc.
  * - Wave is \b standing so it is possible to choose phase of wave to cancel electric or magnetic field on given time layer. In this
  *   implementation magnetic field contribution is zero which means initial time shift is equal to \f$ 0.5\tau \f$.
  */

#ifndef tag_EMResonator_header
#define tag_EMResonator_header		///< \internal Guard.

#include "type_mesh.h"

void tag_EMResonator (FILE *fp, meshVec_p E);

#endif
