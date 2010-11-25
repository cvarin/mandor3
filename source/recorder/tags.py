# -*- coding: utf8 -*-

from __future__ import with_statement

r"""
Tools to convert 'recorder' config file into set of Mandor '.input' and '.cfg'
files. Script to run all recording passes is created as well.
"""

import re, os, stat

from math import *

__all__ = ["make_all"]

interface = {}
sources   = []

# I can import units in many writings (unicode, text, human text) and this dict
# reduces all this zoo to a canonical internal name of the unit.
dealias = {
    "μm"           : "μm",
    "um"           : "μm",
    "micron"       : "μm",
    "microns"      : "μm",
    "fs"           : "fs",
    "femtosec"     : "fs",
    "femtoseconds" : "fs",
    "r0"           : "r0",
    "t0"           : "t0",
}

# Constants in CGS.
CGS = {
    "c"     : 2.9979e+10,	# [cm/sec], speed of light in vacuum.
    "e"     : 4.8032e-10,	# [statcoulomb (statcoul)],
                                #     absolute value of the charge of electron.
    "m"     : 9.1094e-28,	# [gr], mass of the electron.
    "eV"    : 1.6022e-12,	# [erg], energy assiciated with 1 electron-volt.
    "W/cm²" : 1.0e7,		# [erg/(sec⋅cm²)].
}


def set_units (Lambda):
    """Finds all derived scales used to get dimensionless parameters.

    See 'source/setup/tag_units.c' ('[units]' tag)."""
    if Lambda < 0:
        # Computes λ[μm] from n_cr: ω = 2πc/λ = ω_pe(n_e = n_cr).
        n_cr   = -Lambda
        omega  = sqrt (4*pi*n_cr*CGS['e']**2/CGS['m'])
        Lambda = 2*pi*CGS['c']/omega/1e-4

    # Converts Lambda from μm to cm.
    Lambda *= 1e-4

    global UNIT, TIME, LENGTH
    UNIT = {}
    UNIT['r'] = Lambda
    UNIT['λ'] = Lambda
    UNIT['v'] = CGS['c']				# v = c
    UNIT['t'] = Lambda/UNIT['v']			# t = r0/v
    UNIT['ω'] = 2*pi/UNIT['t']				# ω = 2π/t0
    UNIT['A'] = CGS['m']*CGS['c']*UNIT['ω']/CGS['e']	# A0 = mcω/e
    UNIT['E'] = UNIT['A']/(2*pi)			# E0 = mcω/(2πe)
    UNIT['n'] = CGS['m']*UNIT['ω']**2/(16*pi**3*CGS['e']**2) # n = mω²/(16π³e²)
    UNIT['ρ'] = CGS['e']*UNIT['n']			# ρ = en

    # Dimensionless values of main temporal/spatial values.
    UNIT['c'] = CGS['c']/UNIT['v']

    TIME = {
        "t0" : 1.0,
        "fs" : 1.0e-15/UNIT['t'],
        "μm" : 1.0e-4/(CGS['c']*UNIT['t']),
        "r0" : UNIT['r']/(CGS['c']*UNIT['t']),
    }

    LENGTH = {
        "t0" : CGS['c']*UNIT['t']/UNIT['r'],
        "fs" : CGS['c']*1.0e-15/UNIT['r'],
        "μm" : 1.0e-4/UNIT['r'],
        "r0" : 1.0,
    }


def TFSF (I, J, K):
    """Sets Total Field / Scatered Field interface (called from config file)."""
    interface.update ({
        "I0" : I[0], "I1" : I[1],
        "J0" : J[0], "J1" : J[1],
        "K0" : K[0], "K1" : K[1],
    })

def source (I, polarization, focus_X, focus_Y, focus_Z, emitter,
            spot_FWHM, length_FWHM,
            t_front, t_back, t_plato = "0 fs", t_arrival = "0 fs"):
    """Adds new beam to the set (called from config file)."""
    sources.append ({ "I"            : I,
                      "polarization" : polarization,
                      "focus_X"      : focus_X,
                      "focus_Y"      : focus_Y,
                      "focus_Z"      : focus_Z,
                      "emitter"      : emitter,
                      "spot_FWHM"    : spot_FWHM,
                      "length_FWHM"  : length_FWHM,
                      "t_arrival"    : t_arrival,
                      "t_back"       : t_back,
                      "t_front"      : t_front,
                      "t_plato"      : t_plato  })


def cfg_read_int (line, field, tag):
    """Mimics tag scanner in Mandor 'setup'."""
    tokens = re.match ("@\s+([+\-]?\d+)\s+", line)
    if tokens:
        return int (tokens.group (1))
    raise ValueError ("cannot get integer parameter [%s].%s from %s" \
                      % (tag, field, line))

def cfg_read_double (line, field, tag):
    """Mimics tag scanner in Mandor 'setup'."""
    tokens = re.match ("@\s+([+\-]?(\d*)(\.)?\d+([eE][+\-]?\d+)?)\s+", line)
    if tokens:
        return float (tokens.group (1))
    raise ValueError ("cannot get float parameter [%s].%s from %s" \
                      % (tag, field, line))

def parse_config (config):
    """Extracts parameters of the mesh and scales."""

    # Imports all 'config' definitions into global namespace.
    exec config in globals ()

    # Gets [units] and [mesh] to get scales and convertion factors.
    cfg = [x for x in re.split ("[\r\n]", CONFIG_FILE)]

    # Sets units to get all characteristic scales to convert between fs, μm, r0.
    if '[units]' not in cfg:
        raise ValueError ("tag [units] is missing")
    i      = cfg.index ('[units]')
    Lambda = cfg_read_double (cfg[i+1], 'unit', 'units')
    set_units (Lambda)

    # Extracts parameters of the mesh simular to 'source/setup/tag_mesh.c'.
    if '[mesh]' not in cfg:
        raise ValueError ("tag [mesh] is missing")

    global h1, h2, h3, tau, imax, jmax, kmax, Lx, Ly, Lz
    mesh = cfg.index ('[mesh]')
    imax = cfg_read_int    (cfg[mesh+1], 'imax', 'mesh')
    jmax = cfg_read_int    (cfg[mesh+2], 'jmax', 'mesh')
    kmax = cfg_read_int    (cfg[mesh+3], 'kmax', 'mesh')
    Lx   = cfg_read_double (cfg[mesh+4], 'Lx',   'mesh')
    Ly   = cfg_read_double (cfg[mesh+5], 'Ly',   'mesh')
    Lz   = cfg_read_double (cfg[mesh+6], 'Lz',   'mesh')
    tau  = cfg_read_double (cfg[mesh+7], 'tau',  'mesh')

    before_mesh_size = "\n".join (cfg[:mesh+1])
    after_mesh_size  = "\n".join (cfg[mesh+7:])

    h1 = Lx if imax <= 1 else Lx/float (imax)
    h2 = Ly if jmax <= 1 else Ly/float (jmax)
    h3 = Lz if kmax <= 1 else Lz/float (kmax)
    if tau < 0:
        tau /= -sqrt ((imax > 1)/h1**2 + (jmax > 1)/h2**2 + (kmax > 1)/h3**2)

    return (before_mesh_size, after_mesh_size)


def value (s, type = None):
    """Extracts value from string like '1.0e14   W/cm²'."""
    return s.split ()[0]

def is_unit (s, units):
    """Tests if unit in string like '1.0e14   W/cm²' belongs to the group."""
    return s.split (None, 1)[1] in units

def L(x, unit):
    """Converts length/time to length in given units (with speed of light)."""
    unit_x = dealias[x.split (None, 1)[1]]
    unit   = dealias[unit]
    return float (value (x))*LENGTH[unit_x]/LENGTH[unit]

def T(t, unit):
    """Converts length/time to time in given units (with speed of light)."""
    unit_t = dealias[t.split (None, 1)[1]]
    unit   = dealias[unit]
    return float (value (t))*TIME[unit_t]/TIME[unit]


def make_gauss_tag (s, shift = 0):
    """Makes gauss tag with all parameters in correct units."""
    dX = shift/LENGTH["μm"]
    params = {
        "fX"     : repr (L (s["focus_X"],     "μm") + dX),
        "fY"     : repr (L (s["focus_Y"],     "μm")),
        "fZ"     : repr (L (s["focus_Z"],     "μm")),
        "spot"   : repr (L (s["spot_FWHM"],   "μm")),
        "L"      : repr (T (s["length_FWHM"], "fs")),
        "t_back" : repr (T (s["t_back"],      "fs")),
        "Ey"     : int ('Ey' in s["polarization"]),
        "Ez"     : int ('Ez' in s["polarization"]),
        "I"      : value (s["I"])
    }

    # Sign of I choses the units.
    assert is_unit (s["I"], ["W/cm²", "W/cm^2", "eE/mcω", "eE/mc\omega"]), \
            "unit of 'I = %s' is unknown" % s["I"]
    if is_unit (s["I"], ["eE/mcω", "eE/mc\omega"]):
        S["I"] = "-" + S["I"]

    tag = """
[gaussSpot]
@ %(I)s   I: positive => [W/cm²], negative => eE/mcω.
@ 1.0      Cyclic frequency [ω_0]. Must be 1 or you understand all that.
@ %(fX)5s    Focus X coordinate [micron].
@ %(fY)5s    Focus Y coordinate [micron].
@ %(fZ)5s    Focus Z coordinate [micron].
@ %(spot)5s    Gauss-spot width (FWHM for laser intensity) [micron].
@ %(L)5s    Pulse duration (FMHW for laser intensity) [fs].
@ %(t_back)5s    t_front (t_back if for reverse mode) [fs].
@ %(Ey)5d    Ey != 0
@ %(Ez)5d    Ez != 0

""" % params

    return tag


def make_TFSF_tag (F0, F1, add_mode = False, play_mode = False, shift = 0):
    """Prepares TFSF tag with filtering of the region."""
    params = dict (interface)
    params["mode"]   = "playBackward" if play_mode else "record"
    params["filter"] = "add" if add_mode else "new"
    params["F0"]     = F0 + shift
    params["F1"]     = F1 + shift
    params["I0"]    += shift
    params["I1"]    += shift

    tag = """
[TFSF]
@ %(I0)3d		TF/SF interface position (Imin).
@ %(I1)3d		TF/SF interface position (Imax).
@ %(J0)3d		TF/SF interface position (Jmin).
@ %(J1)3d		TF/SF interface position (Jmax).
@ %(K0)3d		TF/SF interface position (Kmin).
@ %(K1)3d		TF/SF interface position (Kmax).
@ %(mode)s	TF/SF interface regime (record, playForward, playBackward).
> %(filter)s		Mode of information accumulation ('new' or 'add').
> %(F0)d		I0    Only nodes with i ∊ [I0, I1] will be actually accumulated.
> %(F1)d		I1
""" % params

    return tag


def make_all (config):
    """Creates all files to record and blend few laser beams."""
    (before_mesh_size, after_mesh_size) = parse_config (config)
    print before_mesh_size
    print "... new mesh/domain ..."
    print after_mesh_size

    print "Key elements to verify:"
    print "    ∘ h1:", h1
    print "    ∘ h2:", h2
    print "    ∘ h3:", h3
    print "    ∘ τ: ", tau

    tags   = []
    times  = []
    meshes = []

    for s_num, s in enumerate (sources):
        # Prints the source input parameters.
        l = max (map (len, s.keys ()))
        f = "%%%ds : %%s" % l
        print "--------- source %d ---------" % s_num
        for k, v in sorted (s.items ()):
            print f % (k, v)

        # Makes emitter/receiver tags.
        (X_l, X_r) = [interface[k]*h1 for k in ("I0", "I1")]
        focus      = L (s["focus_X"], "r0")
        length     = L (s["t_back"],  "r0") + L (s["t_front"], "r0")

        # Computes new imin and imax so the wall is 'length/2' far away from
        # the emitter to cut-off the reflected signal.
        if s["emitter"] == "right":
            assert focus < X_r, "focal plane on the wrong side of interface"
            I    = interface["I1"]
            imin = int (min (focus, X_l)/h1)    - 10
            imax = int (X_r/h1 + 0.5*length/h1) + 10
        else:
            assert focus > X_l, "focal plane on the wrong side of interface"
            I    = interface["I0"]
            imin = int (X_l/h1 - 0.5*length/h1) - 10
            imax = int (max (focus, X_r)/h1)    + 10

        tag   = make_gauss_tag (s,
                                shift = -imin*h1)
        tag  += "\n\n"
        tag  += make_TFSF_tag (F0    = I - 5,
                               F1    = I + 5,
                               shift = -imin,
                               add_mode = bool (tags))

        meshes.append ([imin, imax])

        tags.append (tag)

        # Computes times.
        times.append ({
            "rec" : abs (I*h1 - L (s["focus_X"], "r0"))/UNIT["c"] +
                    T(s["t_back"], "t0") + T(s["t_front"], "t0"),
            "pad" : - T(s["t_back"], "t0") - T(s["t_arrival"], "t0"),
        })

    # Computes global time offset and writes config files.
    start_time = min ([t["pad"] for t in times])
    script     = open ("laser.rec_all.sh", "wt")
    script.write ("""#!/bin/sh

NCPU=4

bad_exit ()
{
    echo "Cannot record source $num"
    exit -1;
}

""")

    for num, (time, tag, mesh) in enumerate (zip (times, tags, meshes)):
        with open ("laser_%d.rec.input" % num, "wt") as fp:
            fp.write ("""
%(before)s
@ %(imax)d\t\timax
@ %(jmax)d\t\tjmax
@ %(kmax)d\t\tkmax
@ %(Lx)s\t\tLx
@ %(Ly)s\t\tLy
@ %(Lz)s\t\tLz
%(after)s
""" %       { "before" : before_mesh_size,
              "imax"   : mesh[1] - mesh[0],
              "jmax"   : jmax,
              "kmax"   : kmax,
              "Lx"     : repr ((mesh[1] - mesh[0])*h1),
              "Ly"     : repr (Ly),
              "Lz"     : repr (Lz),
              "after"  : after_mesh_size,
            })
            fp.write (tag)
            fp.write ("> %d		Frames to pad.\n\n" % \
                      ((time["pad"] - start_time)/tau))

        with open ("laser_%d.run_mandor.cfg" % num, "wt") as fp:
            params = {
                "all" : time["rec"]/tau,
                "chk" : time["rec"]/tau*2,
                "tec" : time["rec"]/tau/20,
            }
            fp.write ("""@ %(all)d\t\tTotal number of steps.
@ %(chk)d\t	Period between checkpoints.
@ %(tec)d\t	Period between TecPlot snapshots (negative means 'None').
@ 0		Continue run from last checkpoint (1 - yes, 0 - no)
@ 20		Check Stop Flag every $< steps
@ -25		Dump period (time steps) for spectral density visualization.
""" % params)

        script.write (r"""
# Records source #%(num)d.
num="%(num)d"
cp -f laser_${num}.run_mandor.cfg run_mandor.cfg
mpirun -np $NCPU ./setup.out laser_${num}.rec.input create 	&& \
    mpirun -np $NCPU ./core.out 				&& \
    ./TFSFdefrag.out 						|| \
    bad_exit;
""" % { "num" : num})

    script.close ()
    os.chmod ("laser.rec_all.sh",
              stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)

    with open ("laser.play.input", "wt") as fp:
        fp.write (CONFIG_FILE + "\n\n")
        fp.write (make_TFSF_tag (F0 = 0,
                                 F1 = 1e6,
                                 play_mode = True))
