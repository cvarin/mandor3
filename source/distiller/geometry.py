# -*- coding: utf8 -*-

r"""
Filter instance represents filter of particles, or a subdomain in
(x, y, z, px, py, pz, q/M) parameter space.

Filter is convinient to imagine as a simple body in phase space (box, cylinder,
sphere, etc.). All filters can be concatenated, intersected, inversed and
subtracted (scissored) with each other.

Special filters work with other parameters like charge, momentum, energy --
very useful to select given component, or given energy range, etc.

You can use Python constants 'None' and 'True' as 'empty' or 'infinite' domain.

Supported operations:
    '+' - concatenation
    '*' - intersection
    '^' - scissors: first object without any intersection with the second
    '~' - inversion

Function 'list_filters()' returns list of installed filters, together with list
of expected arguments.

Example of filter declaration file. First four lines are special and select
records to scan and optional filter block. Parameters are in fixed format:

    @ <parameter> <comment>,

one parameter per line.

------------------------    Example of the file    ------------------------

    @ 0		Start record
    @ 10000	End record
    @ 1		Step
    @ 0		Use filter

    # Lets think that beam propagates, reflects from the plane
    # and we want to see only the beam vicinity:
    #
    #                \     /
    #                 \   /
    #                θ \ /     θ - angle between incident and
    # ----->----->------/          reflected beams in degrees.
    #                  /       δ - width of layer with n > n_cr
    #                 /            we want to grab

    (R, X0, Y0, Z0) = (0.2, 1.0, 1.0, 1.0)
    (th, d) = (40.0, 0.05)

    th *= acos (-1)/180.0		# Degrees to radians conversion.

    incident  = cylinder (r0  = ( 0, Y0, Z0),
                          dir = ( 1,  0,  0),
                          r   = R)

    reflected = cylinder (r0  = (       X0,       Y0, Z0),
                          dir = (-cos (th), sin (th),  0),
                          r   = R)

    plane_of_target = plane (r0 = (X0 + d*cos (th/2), Y0 - d*sin (th/2), Z0),
                             n  = (     - cos (th/2),        sin (th/2),  0))

    return (incident + reflected)^plane_of_target
"""

import types, re

__all__ = ["box", "cylinder", "sphere", "plane", "plasma", "list_filters"]

namestack = []

C_NAMES = {}
C_STACK = []
C_CONST = []

# Width of source file to honor while splitting long expressions.
WIDTH = 60

def ADD_TO_STACK(prefix, expr, counter=0):
    """Creates temporary variable which can be used as prefix."""
    while [x for x in C_NAMES if x.startswith(prefix + str(counter))]:
        counter = 0 if isinstance(counter, str) else counter + 1
    C_NAMES[prefix + str(counter)] = expr
    C_STACK.append(prefix + str(counter))
    return prefix + str(counter)


def ADD_TO_CONST(name, value):
    """Declares constant to use in main loop."""
    assert name not in C_NAMES, "name collision; '%s': '%s' vs '%s'" % (name,
                                                         value, C_NAMES[name])
    if isinstance(value, (int, float)):
        res = "double %s = %s;" % (name, repr(value))
    elif isinstance(value, (tuple, list)):
        assert len(value) == 3, "only 3D vectors, sorry"
        res = "double %s[3] = {%s};" % (name, ", ".join(repr(x) for x in value))
    else:
        raise TypeError("unknown type of constant '%s'" % value)
    C_NAMES[name] = res
    C_CONST.append(name)
    return name


# Recursive conversion to C-source, functions.
def my_repr(x):
    func = eval("do_" + x[0])
    return func(*(x[1:]))

def do_or(*argc):
    res = "(%s)" % " || ".join(my_repr(x) for x in argc)
    if len(res) > WIDTH:
        return ADD_TO_STACK("union", res)
    return res

def do_and(*argc):
    res = "(%s)" % " && ".join(my_repr(x) for x in argc)
    if len(res) > WIDTH:
        return ADD_TO_STACK("cross", res)
    return res

def do_strip(x, y):
    res = "(%s && !%s)" % (my_repr(x), my_repr(y))
    if len(res) > WIDTH:
        return ADD_TO_STACK("leftover", res)
    return res

def do_inverse(x):
    return "(!%s)" % my_repr(x)

def do_plane(r0, n):
    name = ADD_TO_STACK("p", None)
    C_NAMES[name] = "is_under_plane (p, %(n)s_r0, %(n)s_n)" % {'n' : name}
    for l in ["r0", "n"]:
        ADD_TO_CONST(name + "_" + l, eval(l))
    return name

def do_box(r0, e1, e2, e3):
    name = ADD_TO_STACK("b", None)
    C_NAMES[name] = "is_in_box (p, %(n)s_r0, %(n)s_e1, %(n)s_e2, %(n)s_e3)" % {'n' : name}
    for l in ["r0", "e1", "e2", "e3"]:
        ADD_TO_CONST(name + "_" + l, eval(l))
    return name

def do_cylinder(r0, axis, R):
    name = ADD_TO_STACK("c", None)
    C_NAMES[name] = "is_in_cylinder (p, %(n)s_r0, %(n)s_axis, %(n)s_R)" % {'n' : name}
    for l in ["r0", "axis", "R"]:
        ADD_TO_CONST(name + "_" + l, eval(l))
    return name

def do_sphere(r0, R):
    name = ADD_TO_STACK("s", None)
    C_NAMES[name] = "is_in_sphere (p, %(n)s_r0, %(n)s_R)" % {'n' : name}
    for l in ["r0", "R"]:
        ADD_TO_CONST(name + "_" + l, eval(l))
    return name

def do_plasma(qDivM):
    name = ADD_TO_STACK("plasma", None)
    C_NAMES[name] = "is_plasma (p, %s)" % repr(qDivM)
    return name

#def do_energy (emin, emax, mc2):
    #name = ADD_TO_STACK ("energy", None)
    #C_NAMES[name] = "energy_in (u, v, w, %(n)s_emin, %(n)s_emax, %(n)s_mc2)" % {'n' : name}
    #for l in ["emin", "emax", "mc2"]:
        #ADD_TO_CONST (name + "_" + l, eval (l))
    #return name


def do_repr(arg):
    return repr(arg)


class Filter:
    def __init__(self, shape, name=None):
        assert name is None or name not in namestack, \
               "name '%s' is used; see all names: %s" % (name, namestack)
        assert isinstance(shape, tuple), "tuple ('func', arg1, ..) expected"
        self.name = name
        self._    = shape
        namestack.append(name)

    def __repr__(self):
        return my_repr(self._)

    def C(self, varname="res"):
        """Returns C source: constants, tmp variables, final expression."""
        C_STACK[:] = []
        C_CONST[:] = []
        C_NAMES.clear()
        res   = ADD_TO_STACK(varname, my_repr(self._), counter="")
        const = [C_NAMES[k] for k in C_CONST]
        stack = ["int %s = %s;" % (k, C_NAMES[k]) for k in C_STACK]
        return (const, stack, res)

    def __add__(self, othr):
        """Logical OR: returns concatenation of both objects."""
        if isinstance (othr, Filter):
            stack = ['or',]
            stack.extend(self._[1:] if self._[0] == 'or' else (self._,))
            stack.extend(othr._[1:] if othr._[0] == 'or' else (othr._,))
            return Filter(tuple(stack))
        if othr is None:
            return self
        if othr is True:
            return Filter(("repr", 1))
        raise NotImplemented("cannot join spatial object with %s" % othr)

    def __mul__(self, othr):
        """Logical AND: returns intersection of both objects."""
        if isinstance(othr, Filter):
            stack = ['and',]
            stack.extend(self._[1:] if self._[0] == 'and' else (self._,))
            stack.extend(othr._[1:] if othr._[0] == 'and' else (othr._,))
            return Filter(tuple(stack))
        if othr is None:
            return Filter(("repr", 0))
        if othr == True:
            return self
        raise NotImplemented("cannot overlap spatial object with %s" % othr)

    def __xor__(self, other):
        """Scissors: returns part of the object which is not clipped by the other."""
        if isinstance(other, Filter):
            return Filter(('strip', self._, other._))
        if other is None:
            return self
        if other == True:
            return Filter(("repr", 0))
        raise NotImplemented("cannot overlap spatial object with %s" % other)

    def __invert__(self):
        """Inversion: space not occupied by the object."""
        return Filter(('inverse', self._))

    memo = "only '+' (concatenation), '*' (intersection), '^' (scissor), and '~' (inversion) operators defined."
    def __sub__(self, other): raise NotImplemented(Filter.memo)
    def __neg__(self, other): raise NotImplemented(Filter.memo)
    def __pos__(self, other): raise NotImplemented(Filter.memo)


def box(center, side1, side2, side3, name=None):
    return Filter(("box", center, side1, side2, side3), name)

def cylinder(r0, dir, r, name=None):
    return Filter(("cylinder",  r0, dir, r), name)

def sphere(r0, r, name=None):
    return Filter(("sphere", r0, r), name)

def plane(r0, n, name=None):
    return Filter(("plane", r0, n), name)

def plasma(qDivM, name=None):
    return Filter(("plasma", qDivM), name)


def list_filters():
    """Shows prototypes of all filters."""
    import inspect
    print "Installed filters:"
    for l in __all__:
        if l != 'list_filters':
            print "    %s %s" % (l, \
                inspect.formatargspec(*inspect.getargspec(eval(l))))
