/** \file misc_partition.c
  * \brief Partitioning and balancing manager. Takes domain and constructs hierarcy of separations
  * which theoretically allows perfect balancing for any loading density. Tree used to compute
  * loading efficiency and to define which partitioning needs to be adjusted.
  */

#define MC_MISC_PARALLEL_SOURCE		// Removes read-only attribute from global cpu parameters.

#define MC_SKEW_DEFAULT		3	///< Default amplitude of the boundary position jitter (see partition_skew()).

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "type_mesh.h"

#include "log.h"
#include "commandLine.h"
#include "misc_partition.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

/// Global sizes of the domain, local copy.
static int  _dmn_min[6],
           *_dmn_max = _dmn_min + 3;

/// Global parameters of the cpu.
static int cpu_mesh[6] = { 0,  0,  0, -1, -1, -1};
static int cpu_bc[6]   = {-1, -1, -1, -1, -1, -1};

/// Pointers to access a cpu boundaries.
int *cpu_min = cpu_mesh, *cpu_max = cpu_mesh + 3;
int *cpu_bc_min   = cpu_bc,   *cpu_bc_max   = cpu_bc   + 3;

// ---------------------------------------------------------------------------
/// Axis to slice across and number of slices - the simplest group of uniform subdivisioning.
// ---------------------------------------------------------------------------
typedef struct
{
  int axis;					///< Axis to slice across.
  int M;					///< Number of slices.
} prtnCutSet_t;

// ---------------------------------------------------------------------------
/// Tree structure to store order of decompositon: root holds domain and children, children are subdivided and
/// contain their children, etc. Children share the same external boundaries their parents do, so to keep it consistent
/// I store not coordinates but POINTERS to variable(s) which hold real position of wall.
// ---------------------------------------------------------------------------
typedef struct volume_s
{
  /* Tree-structure fields. */
  prtnCutSet_t splitSet;			///< Axis along which region is separated and a number of the subdivisions.
  int *splitTicks;				///< Array of the positions of the separators/boundaries.
  int *offsetChildren;				///< Offset of the beginning of the array of the sub-volumes.
  int offsetParent;				///< Parent leaf.

  /* Data fields. */
  int cpuMin, cpuMax, depth;			///< Range of the CPUs in the group.
  int const * pos1[3];				///< Low-left-rear corner of the domain.
  int const * pos2[3];				///< Up-right-forward corner of the domain.
} volume_t;

static int prtnSplitted[3] = {0, 0, 0};		///< Boolean flags: set if axis was splitted (to turn periodic Bc into splitter BC).
static int prtnMayCut[3] = {1, 1, 1};		///< Flag: gives permission to slice across axis (altered using mandor.rc).
static int prtnSkewed[3] = {0, 0, 0};		///< Boolean flags: set if user wants to perturb axis splitters' positions.

static int prtnSequenceN = 0;			///< Lenght of the sequence of split groups.
static prtnCutSet_t *prtnSequence = NULL;	///< Sequence of split groups readed from command line.

static int prtnRootN = 0, prtnNodes = 0;	///< Size of the tree and offset of the subarray of nodes.
static volume_t *prtnRoot = NULL;		///< Root of the domain partition tree.

#include "misc_partitionPack.c"			// Packing and unpacking of the partitions.

// ---------------------------------------------------------------------------
/// Returns list of region (each corresponds to one cpu, no wrapped copies).
// ---------------------------------------------------------------------------
regList_t
partition_getMapping (void)
{
  assert (prtnRoot);

  regList_t map;
  map.N = prtnRootN - prtnNodes;							// Allocates non-initialized list.
  map.list = (reg_t *) malloc (sizeof (reg_t)*(prtnRootN - prtnNodes));
  assert (map.list);

  volume_t *node = prtnRoot + prtnNodes;
  for (reg_t *reg = map.list, *end = reg + map.N ; reg < end ; ++reg, ++node)		// Packs parameters of the cpus.
  {
    reg->cpu = node->cpuMin;
    reg->min[mc_x] = *node->pos1[mc_x];
    reg->min[mc_y] = *node->pos1[mc_y];
    reg->min[mc_z] = *node->pos1[mc_z];
    reg->max[mc_x] = *node->pos2[mc_x];
    reg->max[mc_y] = *node->pos2[mc_y];
    reg->max[mc_z] = *node->pos2[mc_z];
  }

  return map;
}

// ---------------------------------------------------------------------------
/// Gets parameters of the domain and split it among CPUs (including changes of BC to BC_SPLITTER).
// ---------------------------------------------------------------------------
static void
node_init (const char *msg)
{
  static const char names[][20] = {[BC_PERIODIC] = "periodic", [BC_MIRROR] = "mirror",
                                   [BC_OPEN] = "open (Mur)",  [BC_SPLITTER] = "splitter" };

  // Gets cpu_here number and checks that partitioned array is correctly generated.
  assert (cpu_here + prtnNodes < prtnRootN);
  assert (prtnRoot[cpu_here+prtnNodes].cpuMin == cpu_here && prtnRoot[cpu_here+prtnNodes].cpuMax == cpu_here);

  cpu_min[0] = *prtnRoot[cpu_here+prtnNodes].pos1[0]*mc_have_x;
  cpu_min[1] = *prtnRoot[cpu_here+prtnNodes].pos1[1]*mc_have_y;
  cpu_min[2] = *prtnRoot[cpu_here+prtnNodes].pos1[2]*mc_have_z;
  cpu_max[0] = *prtnRoot[cpu_here+prtnNodes].pos2[0]*mc_have_x;
  cpu_max[1] = *prtnRoot[cpu_here+prtnNodes].pos2[1]*mc_have_y;
  cpu_max[2] = *prtnRoot[cpu_here+prtnNodes].pos2[2]*mc_have_z;

  if (cpu_max[0] - cpu_min[0] < 5*mc_have_x || cpu_max[1] - cpu_min[1] < 5*mc_have_y || cpu_max[2] - cpu_min[2] < 5*mc_have_z)
  {
    partition_show ();
    DIE ("cpu size is < 5 nodes, ghost cells span over few nodes");
  }

  // Checks if we should set a decomposition boundary condition.
  for (int b = 0 ; b < 6 ; ++b)
    cpu_bc_min[b] = (dmn_min[b] != cpu_min[b] || (dmn_bc_min[b] == BC_PERIODIC && prtnSplitted[b % 3])) ? BC_SPLITTER : dmn_bc_min[b];

  SAY_DEBUG ("node_init: (%s)", msg);
  SAY_DEBUG ("  - cpu #%d of %d", cpu_here, cpu_total);
  SAY_DEBUG ("  - mesh work corners are (%d, %d, %d) - (%d, %d, %d) ", cpu_min[0], cpu_min[1], cpu_min[2], cpu_max[0], cpu_max[1], cpu_max[2]);
  SAY_DEBUG ("  - boundary conditions are");
  SAY_DEBUG ("    X[%s, %s]", names[cpu_bc_min[0]], names[cpu_bc_max[0]]);
  SAY_DEBUG ("    Y[%s, %s]", names[cpu_bc_min[1]], names[cpu_bc_max[1]]);
  SAY_DEBUG ("    Z[%s, %s]", names[cpu_bc_min[2]], names[cpu_bc_max[2]]);
  SAY_DEBUG ("  - random jitter amplitudes: [%d, %d, %d]", prtnSkewed[0], prtnSkewed[1], prtnSkewed[2]);
  SAY_DEBUG ("  - decomposition: [%s, %s, %s]", (prtnMayCut[0]) ? "slice" : "keep",
             (prtnMayCut[1]) ? "slice" : "keep", (prtnMayCut[2]) ? "slice" : "keep");
}

// ---------------------------------------------------------------------------
/// Checks if direct domain decomposition is given in command line or in resource file.
// ---------------------------------------------------------------------------
static void
partition_importSequence (void)
{
  assert (prtnSequence == NULL && prtnSequenceN == 0);									// Should be nothing.

  char *sequence;
  if (!cl_findSequence ("decomposition->direct", &sequence) || cpu_total == 1)						// Checks if option was given.
    return;

  int good, M, cpuN = 1;
  char *cursor = sequence, ea;
  good = sscanf (cursor, "%d", &prtnSequenceN);										// Gets number of elements in the sequence.
  assert (good && prtnSequenceN > 0);											// Checks that integer is scanned.

  #define MF_NEXT	{while (*cursor != '\n' && *cursor) ++cursor;}							///< Jumps to next element in the sequence.
  MF_NEXT;

  prtnSequence = (prtnCutSet_t *) malloc (prtnSequenceN*sizeof (prtnCutSet_t));						// Allocates containers for all groups.
  assert (prtnSequence);
  for (int N = 0 ; N < prtnSequenceN ; ++N)										// Reads all members of the sequence.
  {
    assert (*cursor == '\n');												// Checks that element is ready.
    good = sscanf (++cursor, "%c%d", &ea, &M);
    MF_NEXT;

    assert (good && M > 0 && (ea == 'x' || ea == 'y' || ea == 'z'));							// Checks that M is scanned and item is valid.
    prtnSequence[N].axis = ea - 'x';											// Adds group of slices.
    prtnSequence[N].M = M;
    cpuN *= M;
  }
  assert (*cursor == 0);												// Checks that end of sequence is reached.
  free (sequence);
  #undef MF_NEXT

  if (cpuN == cpu_total)												// Checks applicability of sequence.
    return;

  SAY_WARNING ("bad sequence for key 'decomposition->direct'");
  free (prtnSequence);
  prtnSequence = NULL;
  prtnSequenceN = 0;
}

// ---------------------------------------------------------------------------
/// Scans command line and resource file to get domain decomposition hints.
// ---------------------------------------------------------------------------
static void
partition_getHints (void)
{
  assert (prtnSequence == NULL && prtnSequenceN == 0);				// Should be nothing.

  prtnMayCut[0] = prtnMayCut[1] = prtnMayCut[2] = 1;				// Sets default: no blocking;
  memset (prtnSkewed, 0, sizeof (int)*3);					//               no random shifts;
  memset (prtnSplitted, 0, sizeof (int)*3);					//               nothing splitted.

  if (cl_findVec3i ("decomposition->rand shifts", prtnSkewed))			// Gets amplitudes of the noise.
    SAY_DEBUG ("%s: key 'decomposition->rand shifts' is used.", __func__);

  if (cl_findVec3i ("decomposition->cut across", prtnMayCut))			// Gets slicing restrictions.
    SAY_DEBUG ("%s: key 'decomposition->cut across' is used.", __func__);

  // Sets default if all axises are blocked.
  if (! (ACTIVATOR[0]*prtnMayCut[0] || ACTIVATOR[1]*prtnMayCut[1] || ACTIVATOR[2]*prtnMayCut[2]))
  {
    prtnMayCut[0] = prtnMayCut[1] = prtnMayCut[2] = 1;
    SAY_WARNING ("all axises are blocked => slicing hints are ignored.");
  }
}

// ---------------------------------------------------------------------------
/// Takes number of CPUs avaliable and creates partion sequence. For test only, no balancing mensioned.
// ---------------------------------------------------------------------------
static void
partition_makeSequence (void)
{
  assert (prtnSequence == NULL && prtnSequenceN == 0);									// Should be nothing here.
  assert (ACTIVATOR[0]*prtnMayCut[0] || ACTIVATOR[1]*prtnMayCut[1] || ACTIVATOR[2]*prtnMayCut[2]);			// Checks that slicing is possible.

  int ea = 2, cpu_left = cpu_total;
  while (cpu_left > 1)
  {
    prtnSequence = (prtnCutSet_t *) realloc (prtnSequence, (++prtnSequenceN)*sizeof (prtnCutSet_t));			// Adds storage for new group.
    assert (prtnSequence);

    while (! (ACTIVATOR[ea]*prtnMayCut[ea]))										// Looks for next avaliable axis.
      ea = (ea + 1)%3;

    int M = 2;														// Searches for prime divisor of 'cpu_left'.
    while (cpu_left % M) ++M;

    prtnSequence[prtnSequenceN-1].M = M;										// Adds a cut group.
    prtnSequence[prtnSequenceN-1].axis = ea;

    ea = (ea + 1)%3;													// Next iteration.
    cpu_left /= M;
  }
}

// ---------------------------------------------------------------------------
/// \brief Perturbs splitting to introduce good randomness into tests; high bits are used as having more randomness
/// (see man rand ())). Random number initialization uses oscillation of startup time.
// ---------------------------------------------------------------------------
static int
partition_skew (void)
{
  assert (prtnRoot);													// Checks that partition is ready.
  srand (MPI_Wtime ()*1e6);												// Updates chaos.
  int perturbed = 0;
  for (int pos = 0 ; pos < prtnNodes ; ++pos)
  {
    const int axis = prtnRoot[pos].splitSet.axis, M = prtnRoot[pos].splitSet.M;
    for (int l = 0 ; l < M && prtnSkewed[axis] ; ++l)
    {
      prtnRoot[pos].splitTicks[l] += rand ()*(2*prtnSkewed[axis] + 1.0)/RAND_MAX - prtnSkewed[axis];
      perturbed = 1;
    }
  }
  return perturbed;
}

// ---------------------------------------------------------------------------
/// Creates the partition of the domain.
/// NOTE: all cpus at the end of array are real domains of the nodes and they are sorted already.
/// NOTE: this partitioning scheme is the simplest (used just to check || start-up etc).
// ---------------------------------------------------------------------------
static void
partition_buildTree (void)
{
  assert (prtnRoot == NULL);
  assert ((prtnSequence && prtnSequenceN) || cpu_total == 1);

  // Copies domain size to static variables to reference them from root of the tree by using pointer.
  memcpy (_dmn_min, dmn_min, 6*sizeof (int));

  prtnCutSet_t noCut = {0, 0};							// Empty group alias.

  prtnRootN = 1;								// Initial setup of the tree root.
  prtnRoot = (volume_t*) malloc (sizeof (volume_t));
  assert (prtnRoot);

  // Sets pointers to the local static variables with the positions of domain boundaries.
  prtnRoot->pos1[mc_x] = dmn_min + 0;
  prtnRoot->pos1[mc_y] = dmn_min + 1;
  prtnRoot->pos1[mc_z] = dmn_min + 2;
  prtnRoot->pos2[mc_x] = dmn_max + 0;
  prtnRoot->pos2[mc_y] = dmn_max + 1;
  prtnRoot->pos2[mc_z] = dmn_max + 2;
  prtnRoot->depth      = 0;
  prtnRoot->cpuMin     = 0;
  prtnRoot->cpuMax     = cpu_total - 1;
  prtnRoot->splitSet   = noCut;												// Default (no children).
  prtnRoot->splitTicks = prtnRoot->offsetChildren = NULL;
  prtnRoot->offsetParent = -1;

  for (int pos = 0 ; prtnRoot[pos].cpuMin != prtnRoot[pos].cpuMax ; ++pos)
  {
    assert (prtnRoot[pos].depth < prtnSequenceN);

    const int ea = prtnSequence[prtnRoot[pos].depth].axis,
              M = prtnSequence[prtnRoot[pos].depth].M,
              ep = (ea + 1)%3, eq = (ea + 2)%3;
    prtnSplitted[ea] += 1;												// Marks splitted axis.

    prtnRoot[pos].splitSet = prtnSequence[prtnRoot[pos].depth];								// Allocates array of the subregions.
    prtnRoot[pos].splitTicks = (int*) malloc (M*sizeof (int));
    prtnRoot[pos].offsetChildren = (int*) malloc (M*sizeof (int));
    assert (prtnRoot[pos].splitTicks && prtnRoot[pos].offsetChildren);

    const int offsetBase = prtnRootN;
    prtnRootN += M;
    prtnRoot = (volume_t*) realloc (prtnRoot, sizeof (volume_t)*prtnRootN);

    int cpu = prtnRoot[pos].cpuMin;											// Finds the cpu-group sizes.
    const int dCpu = (prtnRoot[pos].cpuMax - prtnRoot[pos].cpuMin + 1)/M;
    assert ((prtnRoot[pos].cpuMax - prtnRoot[pos].cpuMin + 1) % M == 0);

    const int *axisMin = prtnRoot[pos].pos1[ea];									// Low boundary is shared with parent.
    int axisStep = (*(prtnRoot[pos].pos2[ea]) - *(prtnRoot[pos].pos1[ea]))/M;
    for (int l = 0 ; l < M ; ++l)											// Places splitters.
    {
      prtnRoot[pos].splitTicks[l] = (*axisMin) + axisStep*(l + 1);
      prtnRoot[pos].offsetChildren[l] = offsetBase + l;
    }

    for (int l = 0 ; l < M ; ++l)											// Initializes children.
    {
      const int *axisMax = prtnRoot[pos].splitTicks + l;

      prtnRoot[offsetBase + l].pos1[ea] = axisMin;									// Inits boundary pointers.
      prtnRoot[offsetBase + l].pos1[ep] = prtnRoot[pos].pos1[ep];
      prtnRoot[offsetBase + l].pos1[eq] = prtnRoot[pos].pos1[eq];
      prtnRoot[offsetBase + l].pos2[ea] = axisMax;
      prtnRoot[offsetBase + l].pos2[ep] = prtnRoot[pos].pos2[ep];
      prtnRoot[offsetBase + l].pos2[eq] = prtnRoot[pos].pos2[eq];

      prtnRoot[offsetBase + l].depth = prtnRoot[pos].depth + 1;								// Inits range of CPUs inside.
      prtnRoot[offsetBase + l].cpuMin = cpu;
      prtnRoot[offsetBase + l].cpuMax = cpu + dCpu - 1;

      prtnRoot[offsetBase + l].splitSet = noCut;									// Default is no children.
      prtnRoot[offsetBase + l].splitTicks = prtnRoot[offsetBase + l].offsetChildren = NULL;
      prtnRoot[offsetBase + l].offsetParent = pos;

      axisMin = axisMax;												// Steps to the next chunk.
      cpu += dCpu;
    }

    prtnRoot[offsetBase+M-1].pos2[ea] = prtnRoot[pos].pos2[ea];								// Top boundary is shared with parent.
  }

  free (prtnSequence);													// Releases split sequence.
  prtnSequence = NULL;
  prtnSequenceN = 0;

  /* Checks structrure: last and only last 'cpu_total' elements must be 'single' cpu group. */
  for (int pos = 0 ; pos < prtnRootN ; ++pos)
    assert ((pos >= prtnRootN - cpu_total) == (prtnRoot[pos].cpuMin == prtnRoot[pos].cpuMax));
  prtnNodes = prtnRootN - cpu_total;
}

// ---------------------------------------------------------------------------
/// \brief Partition is created on the master node, packed and distributed to all cpus.
/// Than it is unpacked and local nodes are initialized. This ensures complete consistency of the decomposition over cpus.
// ---------------------------------------------------------------------------
void
partition_init (void)
{
    if (prtnRoot)	free (prtnRoot);
    if (prtnSequence)	free (prtnSequence);
    prtnRoot     = NULL;
    prtnSequence = NULL;
    prtnRootN    = prtnSequenceN = 0;

    // Gets blocked axises and noise options.
    partition_getHints ();
    if (!cpu_here) {
        // Checks if user requested his choise.
        partition_importSequence ();
        if (!prtnSequence)
            partition_makeSequence ();
        partition_buildTree ();

        // Adds noise if asked for it.
        if (partition_skew ())
            SAY_WARNING ("decomposition is perturbed for debug purpose");
    }

    int size, *buffer;
    if (!cpu_here)									// cpu#0 broadcasts marchalled partition.
        partition_marshal (&size, &buffer);
    MPI_Bcast (&size, 1, MPI_INT, 0, MPI_COMM_WORLD);					// Broadcasts size of the packed buffer.

    if (cpu_here)
    {											// Receivers prepare room for data.
        buffer = (int*) malloc (size*sizeof (int));
        assert (buffer);
    }
    MPI_Bcast (buffer, size, MPI_INT, 0, MPI_COMM_WORLD);					// Broadcasts packed partition.

    if (cpu_here)										// Receivers unpack the partition.
        partition_unmarshal (buffer);

    free (buffer);

    char msg[100];
    sprintf (msg, "%d cpus parition_init decomposition.", cpu_total);			// Final initialization of local node.
    node_init (msg);
}

// ---------------------------------------------------------------------------
/// Saves partition to file.
// ---------------------------------------------------------------------------
void
partition_save (const char *name)
{
  assert (prtnRoot);													// Checks if domain is partitioned.
  if (cpu_here)														// Only master can save partition.
    return;

  int size, *buffer;
  partition_marshal (&size, &buffer);
  FILE *fp = cfg_open (name, "wb", __func__);
  fwrite (&size, sizeof (int), 1, fp);
  fwrite (buffer, sizeof (int), size, fp);
  fclose (fp);
  free (buffer);
}

// ---------------------------------------------------------------------------
/// Loads partitioning from file.
// ---------------------------------------------------------------------------
void
partition_load (const char *name)
{
  int size, *buffer = NULL;

  if (!cpu_here)
  {
    FILE *fp = cfg_open (name, "rb", __func__);										// Opens file.
    fread (&size, sizeof (int), 1, fp);											// Gets buffer size.
    buffer = (int *) malloc (sizeof (int)*size);									// Allocates memory.
    assert (buffer);
    fread (buffer, sizeof (int), size, fp);										// Gets packed partition.
    fclose (fp);
  }

  MPI_Bcast (&size, 1, MPI_INT, 0, MPI_COMM_WORLD);				// Broadcasts size of the packed buffer.

  if (cpu_here)									// Receivers prepare room for data.
  {
    buffer = (int*) malloc (size*sizeof(int));
    assert (buffer);
  }

  MPI_Bcast (buffer, size, MPI_INT, 0, MPI_COMM_WORLD);				// Broadcasts packed partition.

  partition_unmarshal (buffer);							// Unpacks partition.
  free (buffer);

  char msg[100];
  sprintf (msg, "partitioning is loaded from '%s' file.", name);		// Initializes node parameters.
  node_init (msg);
}

// ---------------------------------------------------------------------------
/// Prints partitioning structure (all level).
// ---------------------------------------------------------------------------
void
partition_show (void)
{
  assert (prtnRoot);								// Makes sure partition is here.
  if (!cpu_here)
  {
    FILE *fp = cfg_open ("output/nodes.dat", "wt", __func__);			// Draws all nodes.
    for (volume_t *leaf = prtnRoot ; leaf < prtnRoot + prtnRootN ; ++leaf)
    {
      reg_t reg = { .min = {*(leaf->pos1)[0], *(leaf->pos1)[1], *(leaf->pos1)[2]},
                    .max = {*(leaf->pos2)[0], *(leaf->pos2)[1], *(leaf->pos2)[2]} };
      fprintf (fp, "'%d-%d: %s'\n", leaf->cpuMin, leaf->cpuMax, reg_printCorners (&reg));
    }
    fclose (fp);
  }
}
