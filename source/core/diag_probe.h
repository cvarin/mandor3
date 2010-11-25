/** \file diag_probe.h
  * \brief Probe diagnostic.
  *
  * Stort facts:
  * - \b probe is a point where all fields are recorded with high resolution in time;
  * - this diagnostic is useful to study mainly the high frequency part of the information;
  * - probes are located at the mesh nodes (coordinates are rounded), all fields from staggered meshes are interpolated linearly;
  * - config file is analyzed at start-up, new probes may be added and old probes may be removed safely;
  * - allocation primitives in file \b diag_probes.cfg (in the root folder of the variant) are:
  *   - \b Point: single point specified by coordinates;
  *   - <b>Line of probes:</b> line is specified by ends, line ends are included and by default and additional equally spaced probes are allocated in between;
  *   - <b>2D subset of probes:</b> specified by origin vector and 2 vertices located on the edges with origin vertex on the other side.
  *     All 4 vertices are included by default plus additional probes in between;
  *   - <b>3D subset of probes:</b> specified by origin vector and 3 vertices located on the edges with origin vertex on the other side.
  *   - <b>terminator:</b> used to terminate the input; very handy - if you want to keep something just place it behind the terminator.
  * - Results may be viewed using \b vProbes.sh script (see pictures below).
  * \image html probe1.png "Ez component of the recorded field for gaussian laser pulse."
  * \image latex probe1.eps "Ez component of the recorded field for gaussian laser pulse." width=10cm
  * \image html probe2.png "Map of the probe locations."
  * \image latex probe2.eps "Map of the probe locations." width=10cm
  *
  * Exmaple of the config file is:
  <pre>
  @ point         Type: single probe
  > 10            X coordinate [micron]
  > 10            Y coordinate [micron]
  > 10            Z coordinate [micron]
  @ volume        Type: probes spaced in 3D subdomain.
  > 2.0           XYZ coordinates of the origin vertex [micron]
  > 2.0
  > 2.0
  > 19.0          XYZ coordinates of the vertex shifted along edge 1 [micron]
  > 2.0
  > 2.0
  > 2.0           XYZ coordinates of the vertex shifted along edge 2 [micron]
  > 19.0
  > 2.0
  > 2.0           XYZ coordinates of the vertex shifted along edge 3 [micron]
  > 2.0
  > 19.0
  > 5             Number of additional points along the 1/2/3 edges.
  > 5
  > 5
  @ end           Terminator record (end of file).

  @ line          Type: probes along line with equal spacing and N additional points in between.
  > 2.2           X coordinate of the start [micron]
  > 2.2           Y coordinate of the start [micron]
  > 2.2           Z coordinate of the start [micron]
  > 18.2          X coordinate of the start [micron]
  > 18.2          Y coordinate of the start [micron]
  > 18.2          Z coordinate of the start [micron]
  > 11            Number of additional points.

  @ rectangle     Type: probes spaced on 2D rectangle.
  > 2.0           XYZ coordinates of the origin vertex [micron]
  > 2.0
  > 2.0
  > 19.0          XYZ coordinates of the vertex shifted along one edge [micron]
  > 2.0
  > 2.0
  > 2.0           XYZ coordinates of the vertex shifted along another edge [micron]
  > 2.0
  > 19.0
  > 5             Number of additional points along the first edge.
  > 5             Number of additional points along the second edge.
  </pre>
  *
  * \attention The most sutable time to call to this diagnostic is after the first half-step for magnetic field - at this moment
  *   all field components are defined on the same time layer.
  */

#ifndef mc_diag_probe_header
#define mc_diag_probe_header				///< \internal Guard

#include "probe.h"
#include "type_mesh.h"

void probe_touch (void);				// MPI_Barrier inside!
void probe_allocate (int cont);				// MPI_Barrier inside!
void probe_postData (double time, meshVec_RO_p E, meshVec_RO_p H);

int probe_loadData (probeHeader_t *header, probeSample_t **data);

#endif
