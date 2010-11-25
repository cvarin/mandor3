// ---------------------------------------------------------------------------
/// Offsets of the fields in the packed partition header (eIntro_size \b must be the last one to hold array size properly).
// ---------------------------------------------------------------------------
enum {
   eIntro_imin = 0, eIntro_jmin, eIntro_kmin, eIntro_imax, eIntro_jmax, eIntro_kmax,
   eIntro_spltX, eIntro_spltY, eIntro_spltZ, eIntro_cpuN, eIntro_prtnRootN,
   eIntro_size	// <== Must be the last one: it is used as array size with the rest working as indices.
};

// ---------------------------------------------------------------------------
/// Offsets of the fields in the packed partition record (eRec_size \b must be the last one to hold array size properly).
// ---------------------------------------------------------------------------
enum {eRec_depth = 0, eRec_cpuMin, eRec_cpuMax, eRec_splitDir, eRec_splitM, eRec_childrenPos, eRec_nextRecPos, eRec_size};

// ---------------------------------------------------------------------------
/// Reallocates buffer if more room is asked for.
// ---------------------------------------------------------------------------
static void
partitionPack_growBuffer (int needSize, int *size, int **buffer)
{
  if (needSize <= *size)
    return;

  *size = needSize;
  *buffer = (int*) realloc (*buffer, needSize*sizeof (int));
  assert (*buffer);
}

// ---------------------------------------------------------------------------
/// \brief Packs entire partitioning: we cannot save or send pointers, so we need to marshall partitioning
/// into array of integers (keeping all informations necessary to unmarshall data structures).
/// Than we can send it to other CPUs (or save it).
// ---------------------------------------------------------------------------
static void
partition_marshal (int *packedSize, int **packedDataArray)
{
  int bufferSize = 0, *buffer = NULL;
  // Initial memory allocation guess.
  partitionPack_growBuffer (eIntro_size + prtnRootN*eRec_size + 1024, &bufferSize, &buffer);

  buffer[eIntro_cpuN] = prtnRootN - prtnNodes;						// Packs global information.
  buffer[eIntro_imin] = dmn_min[0];
  buffer[eIntro_jmin] = dmn_min[1];
  buffer[eIntro_kmin] = dmn_min[2];
  buffer[eIntro_imax] = dmn_max[0];
  buffer[eIntro_jmax] = dmn_max[1];
  buffer[eIntro_kmax] = dmn_max[2];
  buffer[eIntro_spltX] = prtnSplitted[0];
  buffer[eIntro_spltY] = prtnSplitted[1];
  buffer[eIntro_spltZ] = prtnSplitted[2];
  buffer[eIntro_prtnRootN] = prtnRootN;

  int bufferPos = eIntro_size;
  for (int pos = 0 ; pos < prtnRootN ; ++pos)						// Packs data.
  {
    partitionPack_growBuffer (bufferPos + eRec_size + prtnRoot[pos].splitSet.M,		// Ensures capacity of buffer.
                              &bufferSize, &buffer);

    buffer[bufferPos + eRec_depth] = prtnRoot[pos].depth;				// Packs partition node parameters.
    buffer[bufferPos + eRec_cpuMin] = prtnRoot[pos].cpuMin;
    buffer[bufferPos + eRec_cpuMax] = prtnRoot[pos].cpuMax;
    buffer[bufferPos + eRec_splitDir] = prtnRoot[pos].splitSet.axis;
    buffer[bufferPos + eRec_splitM] = prtnRoot[pos].splitSet.M;
    if (prtnRoot[pos].offsetChildren)
      buffer[bufferPos + eRec_childrenPos] = prtnRoot[pos].offsetChildren[0];
    buffer[bufferPos + eRec_nextRecPos] = bufferPos + eRec_size + prtnRoot[pos].splitSet.M - (prtnRoot[pos].splitSet.M != 0);

    for (int i = 0 ; i < prtnRoot[pos].splitSet.M - 1 ; ++i)				// Packs splitter planes.
      buffer[bufferPos + eRec_size + i] = prtnRoot[pos].splitTicks[i];

    bufferPos = buffer[bufferPos + eRec_nextRecPos];					// Places cursor to the next record.
  }

  *packedSize = bufferPos;								// Returns packed array with sizes.
  *packedDataArray = buffer;
}

// ---------------------------------------------------------------------------
/// Unpacks data structures and reinitializes all local arrays: \b completely resets the partitioning!
// ---------------------------------------------------------------------------
static void
partition_unmarshal (const int *buffer)
{
  for (int pos = 0 ; pos < prtnNodes ; ++pos)						// Releases dynamic memory.
  {
    free (prtnRoot[pos].splitTicks);
    free (prtnRoot[pos].offsetChildren);
  }

  const int cpuN = buffer[eIntro_cpuN];							// Extracts global parameters.
  _dmn_min[0] = buffer[eIntro_imin];
  _dmn_min[1] = buffer[eIntro_jmin];
  _dmn_min[2] = buffer[eIntro_kmin];
  _dmn_max[0] = buffer[eIntro_imax];
  _dmn_max[1] = buffer[eIntro_jmax];
  _dmn_max[2] = buffer[eIntro_kmax];
  prtnSplitted[0] = buffer[eIntro_spltX];
  prtnSplitted[1] = buffer[eIntro_spltY];
  prtnSplitted[2] = buffer[eIntro_spltZ];
  prtnRootN = buffer[eIntro_prtnRootN];
  prtnNodes = prtnRootN - cpuN;

  prtnRoot = (volume_t *) realloc (prtnRoot, prtnRootN*sizeof (volume_t));		// Reallocates partition tree.
  assert (prtnRoot);

  int bufferPos = eIntro_size;								// Positions pointer to the first record.
  for (int pos = 0 ; pos < prtnRootN ; ++pos)						// Extracts tree and sew pointers.
  {
    prtnRoot[pos].depth = buffer[bufferPos + eRec_depth];				// Extracts partition node parameters.
    prtnRoot[pos].cpuMin = buffer[bufferPos + eRec_cpuMin];
    prtnRoot[pos].cpuMax = buffer[bufferPos + eRec_cpuMax];
    const int ea = buffer[bufferPos + eRec_splitDir];
    const int M = buffer[bufferPos + eRec_splitM];
    const int offsetBase = buffer[bufferPos + eRec_childrenPos];
    const int ep = (ea + 1)%3, eq = (ea + 2)%3;						// Finds complementary basis vectors.

    if (pos == 0)
    {
      // Root node points to static global domain_<pos> variables.
      prtnRoot->pos1[mc_x] = dmn_min + 0;
      prtnRoot->pos1[mc_y] = dmn_min + 1;
      prtnRoot->pos1[mc_z] = dmn_min + 2;
      prtnRoot->pos2[mc_x] = dmn_max + 0;
      prtnRoot->pos2[mc_y] = dmn_max + 1;
      prtnRoot->pos2[mc_z] = dmn_max + 2;
      prtnRoot->offsetParent = -1;
    }

    int recPos = bufferPos + eRec_size;							// The offset to access splitting ticks.
    bufferPos = buffer[bufferPos + eRec_nextRecPos];					// Comes to the next record.

    prtnRoot[pos].splitSet.M = M;							// Prepares room for ticks and offsets.
    prtnRoot[pos].splitSet.axis = ea;
    if (M > 0)
    {
      prtnRoot[pos].splitTicks = (int*) malloc (M*sizeof (int));
      prtnRoot[pos].offsetChildren = (int*) malloc (M*sizeof (int));
      assert (prtnRoot[pos].splitTicks && prtnRoot[pos].offsetChildren);
    }
    else
    {
      prtnRoot[pos].splitTicks = prtnRoot[pos].offsetChildren = NULL;			// No child: do not allocate arrays.
      continue;
    }

    int i = 0;
    for ( ; i < M - 1 ; ++i)								// Gets positions of the separators.
    {
      prtnRoot[pos].splitTicks[i] = buffer[recPos + i];
      prtnRoot[pos].offsetChildren[i] = offsetBase + i;
    }
    prtnRoot[pos].offsetChildren[i] = offsetBase + i;					// Inits the last record to merge cycles.

    const int *axisMin = prtnRoot[pos].pos1[ea];					// Low boundary is shared with parent.
    for (i = 0 ; i < M ; ++i)								// Updates child.
    {
      const int *axisMax = prtnRoot[pos].splitTicks + i;
      prtnRoot[offsetBase + i].pos1[ea] = axisMin;					// Inits children's pointers.
      prtnRoot[offsetBase + i].pos1[ep] = prtnRoot[pos].pos1[ep];
      prtnRoot[offsetBase + i].pos1[eq] = prtnRoot[pos].pos1[eq];
      prtnRoot[offsetBase + i].pos2[ea] = axisMax;
      prtnRoot[offsetBase + i].pos2[ep] = prtnRoot[pos].pos2[ep];
      prtnRoot[offsetBase + i].pos2[eq] = prtnRoot[pos].pos2[eq];
      prtnRoot[offsetBase + i].offsetParent = pos;
      axisMin = axisMax;								// The separator step.
    }
    prtnRoot[offsetBase+M-1].pos2[ea] = prtnRoot[pos].pos2[ea];				// Top boundary is shared with parent.
  }
}
