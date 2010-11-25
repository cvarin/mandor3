/** \file mainPage.c
  * Doxygen front page for the \b spectr.out diagnostic (no sources here).
  */

/** \mainpage Spectral energy diagnostic of the "Mandor".
  *
  * Deals with spatial spectum analysis of the EM-field.
  *
  * This is the spectral energy diagnostic - part of the standart diagnostic
  * package. During the simulation full distribution of the EM-field in
  * space is dumped into <b>./tmp</b> folder - mainly to avoid parallel FFT at
  * the run-time in the main simulation module. Than this binary files are
  * processed and saved in the compact (single accuracy) format to save disk
  * space. It is possible to remove high frequency modes by reducing the
  * size of the spectral domain - main interest is always in the physically
  * adequate modes with \f$ \lambda \geq 5-7\cdot h \f$.
  *
  * Please note that only spatial Fourier transformation is performed thus
  * there are no way to tell about sign (or directrion) of the wave vector
  * it there are any propagating waves.
  *
  * This diagnostic is supposed to be used interactively thus time limit for
  * interactive programs (1 hour) is installed to quit safely. For the sake
  * of data safety dump files are removed only if special options is provided
  * in the command line as shown below:
  *
  * &nbsp;&nbsp;&nbsp;&nbsp;<b>spectr.out -r</b>\n
  * &nbsp;&nbsp;&nbsp;&nbsp;<b>spectr.out --remove-dumps</b>
  *
  * For details on dumping of the raw data see spectr_dump.h. For fourier
  * transformation stuff see spectr_process.h.
  *
  * Example of the config file \b diag_spectr.cfg:
  <pre>
  @ 0           Start record.
  @ -1          Finish record (if negative then use the last one).
  @ 1           Step size (>0).
  @ 1           System of units (1 - micron/fs, 0 - lambda/t0).
  @ +6          max X - wave number (negative - whole range).
  @ +6          max Y - wave number (negative - whole range).
  @ +6          max Z - wave number (negative - whole range).
  </pre>
  */
