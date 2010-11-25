/** \file em_caps.h
  * Electromagnetic module boundary cells management (see em.h for Maxwell
  * solver structure).
  *
  * <h2>Motivation.</h2>
  * Parallel exchange means that new values of magnetic and electric field in
  * the boundary cells to be evaluated and sent/received. This procedure
  * includes optional usage of the boundary conditions (BCs) in a transversal
  * direction. As a result, in order to keep all boundary condition code in
  * one place the time step is performed for the entire boundary layer at once.
  * Result is cached for subsequent use and stored in <b>boundary caps</b>
  * (see cap_t).
  *
  * This cached data is used to send data to the neighbour cpus (data sent
  * <b>has to be written back to the main meshes</b> to ensure consistency of
  * all values over the cluster). "What we send is what we have" rule is
  * important, because of the disregarding of it may lead to numerical
  * instability (see <a href="#instability">example</a> below).
  *
  * <h2>Data logic.</h2>
  * The shortest way to introduce design ideas is to list key topics:
  * -# Each boundary has associated cap to store electric and magnetic field on
  *    new time layer. All caps are stored in static arrays (em_caps.c::capsE
  *    and em_caps.c::capsH) with 6 elements in each, so array element is
  *    associated with boundary at compile time (storage order is <b>{xMin,
  *    xMax, yMin, yMax, zMin, zMax}</b>).
  * -# Cap covers entire boundary and ghost cells, so caps overlap with each
  *    other.
  * -# Cap consists of mesh to store new field and few geometrical parameters
  *    stored as reg_t. \n This regions are:
  *   - region to \b update - that is region where new field is evaluated
  *     during the caching. Field in this region doesn't depend on boundary
  *     conditions (update uses only fields which are known on the old time
  *     layer). This region is subregion of cap_t::toFlush. Sometimes timestep
  *     region may be extended for particular boundary condition.
  *   - region to \b flush - that is region where new field is evaluated after
  *     time step and after boundary condition usage. This domain is bigger
  *     than cap_t::toUpdate and all values from this region will be used in
  *     parallel exchange and also they will be written back to main mesh
  *     (flushed).
  *   - region to \b claim - that is region which should be retrieved from
  *     neighbour cpus because of local cpu doesn't possesse all information
  *     to evaluate new field in this cells.
  * -# Size of all regions are defined by boundary condition combination. All
  *    offsets for different BCs are documented below. New value of magnetic
  *    field is used to find new electric field so domains for magnetic field
  *    caps depends on region of cached new electric field.
  * -# Caps overlap with each other and to ensure binary equivalence of data
  *    an explicit syncronization is performed (list of regions cap_t::toCopy
  *    holds regions to copy from reference owner).
  *
  * <h2>Boundary conditions offsets.</h2>
  *
  * <h3>Field domains.</h3>
  *
  * There are two ways to treat boundary condition problem.
  *
  * <h4>Partitioning in physical space.</h4>
  *
  * That means to use direct geometrical view and take into account shifts of
  * all mesh nodes (see em.h and type_reg.h for notation of regions). This way
  * is straightforward but it have big disadvantage: all components of EM-field
  * will have different domains of definition. For example, components of
  * electric field will be defined in regions (<i>defined</i> means that node
  * is located inside or on the boundary of the domain's bounding box):
  * - Ex: <b>[1, I] x [0, J] x [0, K]</b>
  * - Ey: <b>[0, I] x [1, J] x [0, K]</b>
  * - Ez: <b>[0, I] x [0, J] x [1, K]</b>
  *
  * Each boundary has only 2 component of field defined on it:
  * - Ex (boundary values): <b>[1, I] x [0, J] x {0, K}</b> and <b>[1, I] x {0, J} x [0, K]</b>
  * - Ey (boundary values): <b>[0, I] x [1, J] x {0, K}</b> and <b>{0, I} x [1, J] x [0, K]</b>
  * - Ez (boundary values): <b>[0, I] x {0, J} x [1, K]</b> and <b>{0, I} x [0, J] x [1, K]</b>
  *
  * That leads to separate treatment of each component = code complication.
  *
  * <h4>Partitioning in index space.</h4>
  *
  * Another way is partitioning of the mesh in the <i>index</i> space. In this
  * case each cpu is responsible for all mesh nodes in the region
  * <b>[0, I] x [0, J] x [0, K]</b> regardless to their positions with respect
  * to the cpu's bounding box. It makes code  for parallel decomposition and
  * boundary layer management really straighforward and simple.
  *
  * \b Note: for proper interpolation of the Lorentz force all fields must be
  * defined in extended domain <b>[0, I+1] x [0, J+1] x [0, K+1]</b> before PIC
  * step.
  *
  * <h3>1D single cpu time step.</h3>
  *
  * Lets consider case of the single cpu simulation. All boudnary conditions are local and field in region <b>[0, I]</b> is
  * supposed to be known and valid before time step. Timestep looks as follow:
  * - E and H is defined in <b>[0, I]</b>
  * - magnetic field is advanced and defined properly in <b>[1, I]</b>
  * - BC for H-field is used to update H in <b>[0, 0]</b>
  * - electric field is advanced and defined properly <b>[0, I-1]</b>
  * - BC for E-field is used to update E in <b>[I, I]</b>
  *
  * Domain of definition is reduced by one in accordance with direction of derivative in Yee solver (upward/downward). So after timestep
  * the boundary conditions are used to restore field in the rest of the domain nodes (and to correct boundary values if necessary).
  *
  * <h3>1D parallel time step.</h3>
  *
  * Prior to timestep field is supposed to be valid and known in the extended domain:
  * - E is defined in <b>[-1, I+1]</b>, H is defined in <b>[0, I+1]</b>
  * - magnetic field is advanced and defined properly in <b>[0, I+1]</b>
  * - electric field is advanced and defined properly <b>[0, I-1]</b>
  * - E and H field are received to extend domain of definition back to original one.
  *
  * Purpose of domain extension is that immidiately after timestep the magnetic field is valid in the region <b>[0, I+1]</b> and
  * can be used in Lorentz force interpolation withour parallel exchange. Neverless, field in ghost cells is treated as temporary substitute and
  * actual value for the new field is requested from the neighbour. In this case boundary condition (parallel exchange) is used to update field in bigger region.
  *
  * <h3>2D/3D time step.</h3>
  *
  * In 2D/3D case decomposition boundary is limited by another 2/4 boundaries and it makes situation more complex. Lets say domain size is <b>[0, I] x [0, J]</b>,
  * decomposition is done along X axis and along Y axis local (periodic, absorbing or mirror) boundary condition is set. Main question is <i>who is responsible for
  * updating field in the ghost cells <b>{-1, I, I+1} x {-1, J, J+1}</b>?</i>
  *
  * It is easy to show that there are two ways to update ghost cells. For example, for E-field:
  * - receive field in <b>{-1, I-1, I, I+1} x [0, J-1]</b> and use local BC to update E-field in <b>[-1, I+1] x {J}</b>.
  * - update field in <b>[0, I-1] x {J}</b> and send/receive field in <b>{-1, I-1, I, I+1} x [0, J]</b>.
  *
  * In present algorithm second way is used. During parallel exchange new electric field should be found in advance to send it using non-blocking communications.
  * That requires cached magnetic field to be defined properly in whole region which means BC usage. To keep uniformity of the code BCs are also applied to E-field
  * and it removes any additional passes to apply local BC after parallel exchange.
  *
  * <h3>Sizes of cached regions.</h3>
  *
  * Along axises with local (periodic, mirror or absorbing) BC only cells required by PIC interpolation are requested from neighbours cpus. If BC is decomposition
  * BC then region is extended to keep data for advancing H-field in bigger region.
  *
  * For periodic BC field is updated using new values and translation symmetry (see emCap_periodic.c).
  *
  * <table border=0>
  * <caption align=top>Periodic boundary condition.</caption>
  * <tr> <td class="indexkey"> Region     </td> <td class="indexkey">      E                </td> <td class="indexkey">          H            </td> </tr>
  * <tr> <td class="indexkey"> To request:</td> <td class="indexvalue"> [0, I+1]            </td> <td class="indexvalue"> [0, I+1]            </td> </tr>
  * <tr> <td class="indexkey"> To update: </td> <td class="indexvalue"> [0, 1] + {I-1}      </td> <td class="indexvalue"> [0, 2] + [I-1, I]   </td> </tr>
  * <tr> <td class="indexkey"> Cached:    </td> <td class="indexvalue"> [0, 1] + [I-1, I+1] </td> <td class="indexvalue"> [0, 2] + [I-1, I+1] </td> </tr>
  * </table>
  *
  * For mirror BC field is updated using new values and mirror symmetry (see emCap_mirror.c). Time step on boundary is performed to update values attributed
  * to boundary but physically shifted inward domain.
  *
  * <table border=0>
  * <caption align=top>Mirror boundary condition.</caption>
  * <tr> <td class="indexkey"> Region     </td> <td class="indexkey">      E                </td> <td class="indexkey">          H            </td> </tr>
  * <tr> <td class="indexkey"> To request:</td> <td class="indexvalue"> [0, I+1]            </td> <td class="indexvalue"> [0, I+1]            </td> </tr>
  * <tr> <td class="indexkey"> To update: </td> <td class="indexvalue"> {1} + [I-1, I]      </td> <td class="indexvalue"> [1, 2] + [I-1, I]   </td> </tr>
  * <tr> <td class="indexkey"> Cached:    </td> <td class="indexvalue"> [0, 1] + [I-1, I+1] </td> <td class="indexvalue"> [0, 2] + [I-1, I+1] </td> </tr>
  * </table>
  *
  * For Mur ABC field is updated using new values and Mur formula (see emCap_Mur.c). Time step on boundary is performed to update values attributed
  * to boundary but physically shifted inward domain. External values are updated using direct copy from the closest inner mesh node.
  *
  * <table border=0>
  * <caption align=top>Mur absorbing boundary condition.</caption>
  * <tr> <td class="indexkey"> Region     </td> <td class="indexkey">      E                </td> <td class="indexkey">          H            </td> </tr>
  * <tr> <td class="indexkey"> To request:</td> <td class="indexvalue"> [0, I+1]            </td> <td class="indexvalue"> [0, I+1]            </td> </tr>
  * <tr> <td class="indexkey"> To update: </td> <td class="indexvalue"> {1} + [I-1, I]      </td> <td class="indexvalue"> [0, 2] + [I-1, I]   </td> </tr>
  * <tr> <td class="indexkey"> Cached:    </td> <td class="indexvalue"> [0, 1] + [I-1, I+1] </td> <td class="indexvalue"> [0, 2] + [I-1, I+1] </td> </tr>
  * </table>
  *
  * For decomposition BC no boundary conditions are used at cache stage but field, requested from other nodes, will replace "speculative" values used for
  * PIC interpolation. Field sended is field cached and flushed so at the end of the time step there will be no discrepancy between values of field in the
  * same mesh node stored on different processor units (see emCap_decomposition.c).
  *
  * <table border=0>
  * <caption align=top>Decomposition boundary condition.</caption>
  * <tr> <td class="indexkey"> Region     </td> <td class="indexkey">      E                </td> <td class="indexkey">          H            </td> </tr>
  * <tr> <td class="indexkey"> To request:</td> <td class="indexvalue"> [-1, I+1]           </td> <td class="indexvalue"> [0, I+1]            </td> </tr>
  * <tr> <td class="indexkey"> To update: </td> <td class="indexvalue"> [0, 1] + [I-1, I-1] </td> <td class="indexvalue"> [0, 2] + [I-1, I+1] </td> </tr>
  * <tr> <td class="indexkey"> Cached:    </td> <td class="indexvalue"> [0, 1] + [I-1, I+1] </td> <td class="indexvalue"> [0, 2] + [I-1, I+1] </td> </tr>
  * </table>
  *
  * <h3>Inner time step.</h3>
  *
  * Cached values of the new field are used to update field in the entire boundary layer so region to advance field in Yee solver is reduced to
  * - <b>[3, I-2] x [3, J-2] x [3, K-2]</b> for magnetic field
  * - <b>[2, I-2] x [2, J-2] x [2, K-2]</b> for electric field
  *
  * <a name="instability"><h3>Example of the unstable decomposition boundary condition.</h3></a>
  *
  * Boundary condition is used to find out fields in the ghost cell region. In case of the decomposition boundary condition this values is stored on
  * the other cpu, but to untie first Yee half-step from parallel syncronization the electric and magnetic fields are defined in extended domains.
  * In scheme below magnetic field is defined in region <b>[0, I+1]</b>, electric field is defined in region <b>[-1, I+2]</b>, and only electric field
  * on the boundaries <b>{-1, I+2}</b> is exchanged and updated. This electric field is mapped to the inner point <b>{2, I-1}</b> of the neighbours and can be found
  * and send prior to main time step (to overlap computations and exchange). Field sended is forgotten and actual field on the sender cpu is found as
  * usual in order not to treat inner points <b>{2, I-1}</b> differently. Theoretically enough volume of information is transmitted and scheme must
  * work. It works but small round-off errors start numerical instability (looks like radiation source is the point where 4 cpus contact each other):
  * \image html splitterOld.png
  * \image latex splitterOld.eps "Numerical instability for impicit syncronisation protocol." width=10cm
  * Instability source is the difference between different values of field in the same mesh nodes stored on the separate cpus.
  *
  * <h3>TF/SF interface.</h3>
  *
  * Full theoretical review of the TF/SF interface is given in em_TFSF.h. Here I want to mention few implementation details:
  * - TF/SF correction is added to cached fields. Reason is that cached field is send to the other cpus <i> and cached magnetic field is used
  *   to find cached electric field used in exchange</i>. That means that TF/SF noise in magnetic field must be canceled before electric field computations.
  * - Once TF/SF contribution is accounted new received field doesn't need any correction and TF/SF should be done for inner regions only:
  *   - <b>[3, I-2] x [3, J-2] x [3, K-2]</b> for magnetic field
  *   - <b>[2, I-2] x [2, J-2] x [2, K-2]</b> for electric field
  * - \f$ H^n \f$ in inner region is found using half step for magnetic step with no TF/SF correction (if particles cross Tf/SF interface solution is unphysical
  *   anyway and no correction may fix it). In cap region solution is found as \f$ H^n = 0.5\cdot (H^{n-1/2} + H^{n+1/2}) \f$ using TF/SF corrected fields.
  *   It leads to visual artifact in magnetic field visualization by tecplot (example is shown below):
  *   \image html TFSF_tecplot.png
  *   \image latex TFSF_tecplot.eps "Visual artifact for Tecplot H-field visualization on TF/SF interface with parallel decomposition." width=10cm
  *   That is only the result of different treatment of TF/SF correction in inner and cached regions. Final result \f$ H^{n+1/2} \f$ is always correct which
  *   was verified by direct dump of data in capped and ghost regions on all nodes using megaDump() function (see /335 snapshot of the code for exact function body).
  * - TF/SF correction is added to the field before using list cap_t::toCopy because of at the time step execution the region cap_t::toFlush is truncated already.
  */

#ifndef EM_CAPS_HEADER
#define EM_CAPS_HEADER

#include "type_mesh.h"

void cap_reportNaNs (const char *msg, meshVec_p H, const reg_t *reg);
void cap_init       (meshVec_p E, meshVec_p H);

void cap_cacheEH           (meshVec_RO_p E, meshVec_p H);
void cap_catchParallelMsgs (meshVec_p    E, meshVec_p H);

void cap_flushHalfH (meshVec_p H);
void cap_flushH     (meshVec_p H);
void cap_flushE     (meshVec_p H);

void cap_resetGhosts (meshVec_p E, meshVec_p H);

#endif
