/** \file misc_markerPacker.h
  * \brief Type definitions for the markers pack/unpack routines.
  * \sa misc_markerPacker.c
  *
  * Main structure is just a wrapper around pointer to piece of memory to support
  * continuous packing process. All stuff is \b static to help compiler to optimize
  * (hopefully inline) all this routines. Header is required only to pass types to
  * pointers.
  *
  * \note You should exclude this file from the <b><tt>megaMake</tt></b> heiristic
  * target generation algorithm - bodies of the functions included themselves and
  * <b>there are no need to compile</b> misc_markerPacker.c.
  *
  * \note Dirty (\b finalized) flag allows multiple usage of the buffer_finalize
  * calls without messing up the mixed pack/read/read/pack-more access sessions for
  * cache-type buffers.
  */

#ifndef MC_BUFFERLIBRARY_HEADER
#define MC_BUFFERLIBRARY_HEADER						///< Guard against multiple include.

#endif
