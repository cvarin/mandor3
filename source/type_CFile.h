/** \file type_CFile.h
  * Fundamental type to power-up complex module IO (see type_CFile.c for methods).
  *
  * <h2>Motivation.</h2>
  * Very often diagnostic produces many data arrays which are weakly related to
  * each other (for example, different fields plus header with mesh parameters).
  * IO of such files requires strict layout – fixed positioning and ordering. In
  * case of development code that means that different versions of the code will
  * be inconsistent even if part of the information we probably would like to
  * get is still in old form.
  *
  * It also forces developer to review all input and output routines to check
  * consistency.
  *
  * This library designed to introduce some flexibility to the file IO. Data is
  * stored in \b chunks. Each chunk is sequential piece of data storage with
  * name. Order of chunks can be arbitrary, all chunks are referenced only by
  * name to simplify development – all information thus is in human-readable
  * form.
  *
  * Overhead is insignificant but flexibility provided helps to ignore ordering
  * in the file as long as structure of the chunk left unmodified. It also
  * permits adding of the new chunks to the file as a part of post-processing.
  *
  * Wrapping of the fread and fwrite routines help to ease data consistency
  * check implementations.
  *
  * Chunks are stored in form of
  * - file type identifier (simple piece of text)
  * - data chunks:
  *   - chunk length: long long integer (64 bits)
  *   - chunk deleted flag: integer, used to invalidate entire chunk
  *   - chunk name length: integer (\b not including terminating NULL)
  *   - line break
  *   - chunk name:   \b not NULL terminated string
  *   - line break
  *   - chunk body:   data written by client
  *
  * \attention Deleting chunks is \b not an option for very big data-sets! It is
  * done only to ease processing, i.e. support packed descriptions instead of
  * file trees - helps to use perl script for txt file scan instead of searching
  * entire file tree looking only for small files with different descriptions).
  * Goal was to give an option to enjoy fast searching/probing/updating of the
  * small content file rather than doing a replacement of entire filesystem.
  *
  * \todo MD4 for data consistency check (to ensure that check-point is written
  *       normally, sometimes job may be terminated with no confidence in the
  *       data quality - IBM cluster under NFS did this).
  *
  * \todo Deleting of the chunk.
  */

#ifndef MC_TYPE_CHUNKFILE_HEADER
#define MC_TYPE_CHUNKFILE_HEADER

#include <stdio.h>

/**
  * Target of the code is big clusters so we need to put big files (> 4Gb) from
  * the beginning.
  */
typedef long long int INT64;

/**
  * Basic element of the sectioned file - \b chunk.
  */
typedef struct {
    INT64  pos;			///< Position of the beginning of chunk.
    INT64  capacity;		///< Chunk's capacity.
    int    headerSize;		///< Size of the header (includes variable parameter – length of the chunk name).
    int    name;		///< Chunk name offset in reallocatable pool.
    int    deleted;		///< Flag to mark that chunk is deleted.
} chunk_t;

/**
  * Wrapper for the file pointer to form chunk file.
  */
typedef struct {
    FILE    *f;			///< File pointer.
    char    *name;		///< File name (used for debug outputs).
    int      mode;		///< Mode of the file access (read, write or append).

    int      chunksN;		///< Number of chunks.
    chunk_t *chunks;		///< Array of chunk descriptors.
    char    *namesPool;		///< One solid array with all names (to remove fragmentation of the memory and simplify deallocation).
    int      namesPoolN;	///< Size of the \b cNames array (for realloc).
    INT64    bytesProcessed;	///< Flag to indicate chunk opening - accumulates total size of data written or data read.
    int      currentChunk;	///< Current chunk's number (used by seek, might be used in error reports).
} CFile_t;

CFile_t *CF_openUpdate  (const char *nameFormat, ...);
CFile_t *CF_openWrite   (const char *nameFormat, ...);
CFile_t *CF_openRead    (const char *nameFormat, ...);
void     CF_close       (CFile_t *file);
void     CF_seek        (CFile_t *file, INT64 offset, int whence);

void     CF_openChunk   (CFile_t *file, const char *nameFormat, ...);
void     CF_closeChunk  (CFile_t *file);
int      CF_findChunk   (CFile_t *file, const char *nameFormat, ...);
int      CF_probeChunk  (CFile_t *file, const char *nameFormat, ...);
void     CF_deleteChunk (CFile_t *file, const char *nameFormat, ...);
CFile_t* CF_defrag      (CFile_t *file);

void     CF_write       (const void *data, size_t size, size_t nmemb, CFile_t *file);
void     CF_read        (void *data, size_t size, size_t nmemb, CFile_t *file);
void     CF_print       (CFile_t *file, const char *format, ...);
void     CF_scan        (CFile_t *file, const char *format, ...);

#endif
