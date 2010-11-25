/** \file mainPage.c
  * Doxygen front page for the \b tecplot.out diagnostic (no sources here).
  */

/** \mainpage Tecplot energy diagnostic of the "Mandor".
  *
  * Deals with spatial EM-field snapshots. Data are saved just before particles
  * time-step (unlike system check-points) so both electric and magnetic fields
  * are defined on the integer time layer.
  *
  * Examples of the config file \b diag_tecplot.cfg:
  * - outputs whole domain
  <pre>
  @ 0             Start record.
  @ -1            Finish record (if negative then use the last one).
  @ 1             Step size (>0).
  @ 1             System of units (1 - micron/fs, 0 - lambda/t0).
  @ domain        Visualization region (may be 'domain', 'slice' or 'subdomain').
  </pre>
  *
  * - outputs only subdomain
  <pre>
  @ 0             Start record.
  @ -1            Finish record (if negative then use the last one).
  @ 1             Step size (>0).
  @ 1             System of units (1 - micron/fs, 0 - lambda/t0).
  @ subdomain     Visualization region (may be 'domain', 'slice' or 'subdomain').
  > 0             i-min
  > 0             j-min
  > 0             k-min
  > 10            i-max
  > 10            j-max
  > 10            k-max
  </pre>
  *
  * - outputs slice of the domain
  <pre>
  @ 0             Start record.
  @ -1            Finish record (if negative then use the last one).
  @ 1             Step size (>0).
  @ 1             System of units (1 - micron/fs, 0 - lambda/t0).
  @ slice         Visualization region (may be 'domain', 'slice' or 'subdomain').
  > 0             Axis to slice across.
  > 2             Position.
  </pre>
  */
