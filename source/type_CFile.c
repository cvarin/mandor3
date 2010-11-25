/** \file type_CFile.c
  * Extension of the file IO with flexible layout and simple human-readable interface for developers (please see type_CFile.h).
  */

#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "log.h"
#include "type_CFile.h"

#define MC_TEMPORARY_NAME_BUFFER_SIZE		(500)			///< Size of the temporary buffer (for the file name).

// ---------------------------------------------------------------------------
/// Enumerator for the file mode.
// ---------------------------------------------------------------------------
enum {mode_write = 1, mode_read = 2, mode_update = 3};

static const char *magic = "============\n Chunk file\n============\n";								// File type indicator (magic).

// ---------------------------------------------------------------------------
/// Adds new chunk description to the cache tables and return number of chunk.
// ---------------------------------------------------------------------------
static int
CF_addChunkInfo (CFile_t *file, const char *name, INT64 pos, int headerSize, INT64 capacity, int deleted)
{
  assert (file && name);

  // Checks if name is unique.
  ENSURE (!CF_probeChunk (file, name),
          "chunk '%s' already exists in file '%s'", name, file->name);

  ++file->chunksN;														// Counts new chunk.

  file->chunks = (chunk_t*) realloc (file->chunks, sizeof (chunk_t)*file->chunksN);						// Adds new chunk.
  assert (file->chunks);
  file->chunks[file->chunksN-1].pos = pos;											// Writes field as told.
  file->chunks[file->chunksN-1].capacity = capacity;										// Total size of chunk.
  file->chunks[file->chunksN-1].headerSize = headerSize;
  file->chunks[file->chunksN-1].name = file->namesPoolN;									// Name added to the end of pool.
  file->chunks[file->chunksN-1].deleted = deleted;

  const int l = strlen (name) + 1;												// Adds NULL terminated name.
  file->namesPoolN += l;
  file->namesPool = (char*) realloc (file->namesPool, sizeof (char)*file->namesPoolN);
  assert (file->namesPool);
  memcpy (file->namesPool + file->namesPoolN - l, name, l);

  return file->chunksN - 1;
}

// ---------------------------------------------------------------------------
/// Imports chunk structure to the cache (used on opening in 'read' or 'append' mode).
// ---------------------------------------------------------------------------
static void
CF_importChunks (CFile_t *file)
{
  fseek (file->f, 0, SEEK_END);													// Gets file size.
  const INT64 tail = ftell (file->f);
  assert (tail >= strlen (magic));												// Checks if magic can fit the file.
  fseek (file->f, 0, SEEK_SET);													// Rewinds cursor back to start.

  int size = strlen (magic);
  char *buffer = (char*) malloc (size + 1);
  assert (buffer);
  fread (buffer, sizeof (char), size, file->f);
  buffer[size] = 0;														// Terminates string in buffer.
  assert (!strcmp (magic, buffer));												// Checks that magic is correct.

  INT64 pos = ftell (file->f);													// Remembers position of start.
  while (pos < tail - 1)													// Extracts chunk infos.
  {
    int deleted;
    INT64 capacity;
    fread (&capacity, sizeof (INT64), 1, file->f);										// Gets size of chunk.
    fread (&deleted, sizeof (int), 1, file->f);											// Gets chunk status.
    fread (&size, sizeof (int), 1, file->f);											// Gets length of chunk name.
    buffer = (char*) realloc (buffer, sizeof (char)*(size + 3));								// +2 characters are line breaks.
    assert (buffer);
    fread (buffer, sizeof (char), size + 2, file->f);										// Gets chunk name.
    buffer[size+1] = 0;														// Terminates chunk name.

    CF_addChunkInfo (file, buffer + 1, pos, ftell (file->f) - pos, capacity, deleted);						// Adds chunk info.
    fseek (file->f, pos + capacity, SEEK_SET);											// Jumps to the next chunk.
    pos += capacity;
  }
  file->currentChunk = -1;													// No chunks opened.
  free (buffer);
}

// ---------------------------------------------------------------------------
/// Alias for big common operation - expansion of the 'nameFormat + args' into the 'name' string.
// ---------------------------------------------------------------------------
#define MF_NAMEFORMAT_TO_NAME_EXPANDER										\
  va_list  argptr;												\
  char name[MC_TEMPORARY_NAME_BUFFER_SIZE+1];  									\
  name[MC_TEMPORARY_NAME_BUFFER_SIZE] = 0;		/* Marks end. */					\
  														\
  va_start (argptr, nameFormat);			/* Prints name. */					\
  vsnprintf (name, MC_TEMPORARY_NAME_BUFFER_SIZE + 1, nameFormat, argptr);					\
  va_end (argptr);												\
  														\
  if (name[MC_TEMPORARY_NAME_BUFFER_SIZE])		/* Checks end mark. */					\
    DIE("name '%s' is longer that %d chars", name, MC_TEMPORARY_NAME_BUFFER_SIZE);


// ---------------------------------------------------------------------------
/// Opens file for writing and returns pointer to the descriptor.
// ---------------------------------------------------------------------------
CFile_t*
CF_openWrite (const char *nameFormat, ...)
{
  CFile_t* file = (CFile_t*) calloc (1, sizeof (CFile_t));									// Allocs empty file struct.
  assert (file);

  MF_NAMEFORMAT_TO_NAME_EXPANDER;												// Expands arguments.

  file->name = (char*) malloc (sizeof (char)*(strlen (name) + 1));								// Makes copy of the file name.
  assert (file->name);
  strcpy (file->name, name);

  file->f = fopen (name, "wb");
  assert (file->f);
  int l = strlen (magic);
  fwrite (magic, sizeof (char), l, file->f);											// Writes header string.
  file->mode = mode_write;
  file->currentChunk = -1;													// No chunks opened.
  return file;
}

// ---------------------------------------------------------------------------
/// Opens file for reading and returns pointer to the descriptor or NULL if file cannot be opened.
// ---------------------------------------------------------------------------
CFile_t*
CF_openRead (const char *nameFormat, ...)
{
  CFile_t* file = (CFile_t*) calloc (1, sizeof (CFile_t));									// Allocates empty file struct.
  assert (file);

  MF_NAMEFORMAT_TO_NAME_EXPANDER;												// Expands arguments.

  file->name = (char*) malloc (sizeof (char)*(strlen (name) + 1));								// Makes copy of the file name.
  assert (file->name);
  strcpy (file->name, name);

  file->f = fopen (name, "rb");
  if (!file->f)
  {
    CF_close (file);
    return NULL;
  }

  CF_importChunks (file);													// Imports content of the file.
  file->mode = mode_read;

  return file;
}

// ---------------------------------------------------------------------------
/// Opens file for writing but imports all old chunks and keeps them.
// ---------------------------------------------------------------------------
CFile_t*
CF_openUpdate (const char *nameFormat, ...)
{
  CFile_t* file = (CFile_t*) calloc (1, sizeof (CFile_t));									// Allocates empty file struct.
  assert (file);

  MF_NAMEFORMAT_TO_NAME_EXPANDER;												// Expands arguments.

  file->name = (char*) malloc (sizeof (char)*(strlen (name) + 1));								// Makes copy of the file name.
  assert (file->name);
  strcpy (file->name, name);

  file->f = fopen (name, "rb");													// Checks if file exists.
  if (!file->f)															// Empty file starts as write only.
  {
    CF_close (file);
    file = CF_openWrite (name);													// Creates empty initialized file.
    fclose (file->f);
    file->f = fopen (name, "rb");												// Now file exists.
    assert (file->f);
  }

  CF_importChunks (file);													// Imports content of the file.
  fclose (file->f);

  file->f = fopen (name, "r+b");												// Reopens file for writing.
  assert (file->f);
  assert (file->bytesProcessed == 0);
  file->mode = mode_update;
  file->currentChunk = -1;													// No chunks opened.

  return file;
}

// ---------------------------------------------------------------------------
/// Closes file.
// ---------------------------------------------------------------------------
void
CF_close (CFile_t *file)
{
  assert (file);

  if (file->bytesProcessed > 0)
  {
    SAY_WARNING ("unfinalized session in file '%s' on closing", file->name);
    CF_closeChunk (file);
  }

  #define RELEASE(FUNC, POINTER) if ((POINTER)) { FUNC (POINTER); }								///< Alias for safe deallocator.
  RELEASE (free,   file->name);
  RELEASE (free,   file->chunks);
  RELEASE (free,   file->namesPool);
  RELEASE (fclose, file->f);
  RELEASE (free,   file);
  #undef RELEASE
}

// ---------------------------------------------------------------------------
/// Removes deleted chunks from file (import during opening already imports full map of chunks).
// ---------------------------------------------------------------------------
CFile_t *
CF_defrag (CFile_t *file)
{
  assert (file);
  assert (file->mode == mode_update && file->bytesProcessed <= 0);

  char buffer[4096];
  char name[500], oldName[500];
  name[499] = 0;
  snprintf (name, 500, "%s___defrag___", file->name);
  assert (name[499] == 0);

  CFile_t *newFile = CF_openWrite (name);
  for (int n = 0 ; n < file->chunksN ; ++n)
  {
    chunk_t *chunk = file->chunks + n;
    if (chunk->deleted)														// Skips deleted files.
      continue;

    CF_findChunk (file, file->namesPool + chunk->name);
    CF_openChunk (newFile, file->namesPool + chunk->name);
    while (file->bytesProcessed <= - 4096)											// Copies data using 4Kb pages.
    {
      CF_read (buffer, 1, 4096, file);												// Reads page.
      CF_write (buffer, 1, 4096, newFile);											// Writes page.
    }
    int remainer = - file->bytesProcessed;											// Computes tail size.
    if (remainer)
    {
      CF_read (buffer, 1, remainer, file);											// Reads the tail of the file.
      CF_write (buffer, 1, remainer, newFile);											// Writes the tail of the file.
    }
    CF_closeChunk (newFile);
  }
  CF_close (newFile);

  strcpy (oldName, file->name);
  CF_close (file);
  rename (name, oldName);

  return CF_openUpdate (oldName);
}

// ---------------------------------------------------------------------------
/// Fast checks if chunk exists without any \b vsnprintf or perturbations in reading/writing process.
// ---------------------------------------------------------------------------
static int
CF_identifyChunk (CFile_t *file, const char *name)
{
  for (int n = 0 ; n < file->chunksN ; ++n)											// Checks if name is unique.
    if (!strcmp (name, file->namesPool + file->chunks[n].name) && !file->chunks[n].deleted)
      return n;

  return -1;
}

// ---------------------------------------------------------------------------
/// Checks if chunk exists (does no perturbation in reading/writing process).
// ---------------------------------------------------------------------------
int
CF_probeChunk (CFile_t *file, const char *nameFormat, ...)
{
  assert (file && nameFormat);
  assert (file->f);														// Checks if file is opened.

  MF_NAMEFORMAT_TO_NAME_EXPANDER;												// Expands arguments.

  return (CF_identifyChunk (file, name) >= 0);
}

// ---------------------------------------------------------------------------
/// Searches for the chunk with given name and positions file cursor on the chunk's body (returns 0 on success).
// ---------------------------------------------------------------------------
int
CF_findChunk (CFile_t *file, const char *nameFormat, ...)
{
  assert (file && nameFormat);
  assert (file->mode & mode_read);
  assert (file->bytesProcessed <= 0);												// No opened chunks for writing.

  MF_NAMEFORMAT_TO_NAME_EXPANDER;												// Expands arguments.

  int n = CF_identifyChunk (file, name);

  if (n < 0)
  {
    file->currentChunk = -1;													// No chunks opened.
    return 1;
  }

  file->currentChunk = n;													// Remembers chunk's number.
  fseek (file->f, file->chunks[n].pos + file->chunks[n].headerSize, SEEK_SET);							// Positions cursor to the data.
  file->bytesProcessed = - file->chunks[n].capacity + file->chunks[n].headerSize;						// Saves total amount of data to read.
  return 0;
}

// ---------------------------------------------------------------------------
/// Starts writing of the new chunk.
// ---------------------------------------------------------------------------
void
CF_openChunk (CFile_t *file, const char *nameFormat, ...)
{
  assert (file && nameFormat);
  assert (file->mode & mode_write);
  assert (file->bytesProcessed <= 0);												// No opened chunks.

  MF_NAMEFORMAT_TO_NAME_EXPANDER;												// Expands arguments.

  int chunkPresented = (CF_identifyChunk (file, name) >= 0);
  assert (!chunkPresented);

  fseek (file->f, 0, SEEK_END);													// Rewinds file-pos to the end.
  const INT64 pos = ftell (file->f);

  INT64 size = -1;
  fwrite (&size, sizeof (INT64), 1, file->f);											// Writes dummy size (to be finalized).

  char lineBreak = '\n';
  int l = 0;
  fwrite (&l, sizeof (int), 1, file->f);											// Writes 'deleted' flag.
  l = strlen (name);
  fwrite (&l, sizeof (int), 1, file->f);											// Writes chunk name length.
  fwrite (&lineBreak, sizeof (char), 1, file->f);										// Breaks line.
  fwrite (name, sizeof (char), l, file->f);											// Writes chunk name.
  fwrite (&lineBreak, sizeof (char), 1, file->f);										// Breaks line.

  file->bytesProcessed = ftell (file->f) - pos;											// Gets header size (marks session).
  file->currentChunk = CF_addChunkInfo (file, name, pos, file->bytesProcessed, file->bytesProcessed, 0);			// Adds chunk info to cache table.
}

// ---------------------------------------------------------------------------
/// Finishes writing of the new chunk (file size is double checked to make sure client didn't write anything using direct write to the file->f).
// ---------------------------------------------------------------------------
void
CF_closeChunk (CFile_t *file)
{
  assert (file);
  assert (file->bytesProcessed > 0 && (file->mode & mode_write));								// Checks if chunk was opened.

  INT64 capacity = ftell (file->f) - file->chunks[file->currentChunk].pos;							// Calculates size of chunk.
  assert (capacity == file->bytesProcessed);											// Double checks the chunk size.

  fseek (file->f, file->chunks[file->currentChunk].pos, SEEK_SET);								// Rewinds to the chunk's origin.
  fwrite (&capacity, sizeof (INT64), 1, file->f);										// Writes size.
  fseek (file->f, 0, SEEK_END);													// Rewinds file-pos to the end.
  file->chunks[file->currentChunk].capacity = capacity;										// Updates size of the chunk.
  file->bytesProcessed = 0;
  file->currentChunk = -1;													// No chunks opened.
}

// ---------------------------------------------------------------------------
/// Marks chunk as deleted (use CF_defrag to remove empty chunks from file to save space).
// ---------------------------------------------------------------------------
void
CF_deleteChunk (CFile_t *file, const char *nameFormat, ...)
{
  assert (file && nameFormat);
  assert (file->mode == mode_update);
  assert (file->bytesProcessed <= 0);												// No opened chunks.

  MF_NAMEFORMAT_TO_NAME_EXPANDER;												// Expands arguments.

  int n = CF_identifyChunk (file, name);
  if (n >= file->chunksN)
    return;

  file->chunks[n].deleted = 1;
  fseek (file->f, file->chunks[n].pos + sizeof (INT64), SEEK_SET);								// Rewinds to the chunk's flag.
  fwrite (&file->chunks[n].deleted, sizeof (int), 1, file->f);									// Writes flag.
  fseek (file->f, 0, SEEK_END);													// Rewinds file-pos to the end.
  file->bytesProcessed = 0;
  file->currentChunk = -1;													// No chunks opened.
}

// ---------------------------------------------------------------------------
/// Writes data to the new chunk (later it is possible to send it through MD4).
// ---------------------------------------------------------------------------
void
CF_write (const void *data, size_t size, size_t nmemb, CFile_t *file)
{
  assert (file);
  assert (file->bytesProcessed > 0 && (file->mode & mode_write));

  fwrite (data, size, nmemb, file->f);												// Writes size.
  file->bytesProcessed += size*nmemb;
}

// ---------------------------------------------------------------------------
/// Reads data to the new chunk (later it is possible to send it through MD4).
// ---------------------------------------------------------------------------
void
CF_read (void *data, size_t size, size_t nmemb, CFile_t *file)
{
  assert (data);
  assert (file && (file->mode & mode_read));

  ENSURE (file->bytesProcessed + size*nmemb <= 0,
          "violation of the chunk's boundary in file '%s'", file->name);

  fread (data, size, nmemb, file->f);												// Reads data.
  file->bytesProcessed += size*nmemb;												// Counts readed data size.
}

// ---------------------------------------------------------------------------
/// Prints string and writes as raw data into file.
// ---------------------------------------------------------------------------
void
CF_print (CFile_t *file, const char *format, ...)
{
  assert (format && file);
  assert (file->bytesProcessed > 0 && (file->mode & mode_write));

  INT64 oldPos = ftell (file->f);												// Remembers old position.
  va_list  argptr;
  va_start (argptr, format);													// Prints line to file.
  vfprintf (file->f, format, argptr);
  va_end (argptr);

  file->bytesProcessed += ftell (file->f) - oldPos;										// Counts written bytes.
}

// ---------------------------------------------------------------------------
/// Reads text data from file.
// ---------------------------------------------------------------------------
void
CF_scan (CFile_t *file, const char *format, ...)
{
  assert (format && file);
  assert ((file->mode & mode_read));

  ENSURE (file->bytesProcessed < 0,
          "[file '%s'] violation of the chunk's boundary or mixed read/write",
          file->name);

  INT64 oldPos = ftell (file->f);
  va_list  argptr;
  va_start (argptr, format);													// Prints line to buffer.
  vfscanf (file->f, format, argptr);												// Reads data.
  va_end (argptr);

  file->bytesProcessed += ftell (file->f) - oldPos;										// Counts volume of readed data.
}

// ---------------------------------------------------------------------------
/// Reads text data from file.
// ---------------------------------------------------------------------------
void
CF_seek (CFile_t *file, INT64 offset, int whence)
{
  assert (file && whence == SEEK_SET && offset >= 0);
  assert ((file->mode & mode_read) && file->currentChunk >= 0);

  const chunk_t *chunk = file->chunks + file->currentChunk;
  fseek (file->f, chunk->pos + chunk->headerSize + offset, SEEK_SET);								// Positions cursor to the data segment.
  file->bytesProcessed = - chunk->capacity + chunk->headerSize + offset;							// Saves total amount of data to read.
  ENSURE (file->bytesProcessed < 0, "seek points outside of the chunk");
}
