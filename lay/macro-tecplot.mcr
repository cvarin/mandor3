#!MC 1000

# ================================================================================
# This script is used to manage set of the preplotted files by reloading them
# automatically using this macro-menu to browse data. Main job of file preparation
# is done by vTec_preplot.pl script. See the tail of this script to get details.
# ================================================================================

#
# Parameters to set range of zones (actually, symlinks to prepared data-sets).
#
$!VarSet |zone_max| = 1					# Reads parameters from the shell
$!VarSet |zone_now| = 1
$!VarSet |zone_loaded| = -1

#
# Different flag to set state of plot when new zone is activated.
#
$!VarSet |flag_resetAxis| = 0				# Useful for 1D tecplot meshes.
$!VarSet |flag_commonLevels| = 0			# Useful for moderate variations case.
$!VarSet |flag_commonLevelsAreSet| = 0			# Used to avoid unnecessary resettings.
$!VarSet |flag_rescaleContours| = 1			# Useful for big variations of magnitude.
$!VarSet |flag_killMesh| = 1				# Useful to deactivate mesh (tec default).
$!VarSet |flag_nameIsSet| = 0				# Useful to remind about proper name for export.

#
# Different flag to set state of plot when new zone is activated.
#
$!VarSet |param_skipI| = 5				# Skip steps to make mesh coarse.
$!VarSet |param_skipJ| = 5
$!VarSet |param_skipK| = 5

$!VarSet |param_globContourVar| = 4			# Tecplot cannot tell which variable is in charge - have to ask user.
$!VarSet |param_levelsN| = 15				# Number of contour levels (smaller for isosurfaces).
$!VarSet |param_prefix| = "tec3D"			# Export file prefix.

# Maximum and minimum for common contour variable in the dataset over all files.
$!VarSet |param_minCommonVar| = -10
$!VarSet |param_maxCommonVar| = +10

$!VarSet |slices_N| = -1
$!VarSet |slices_axis| = 1

# ------------------------------------------------------------------------------------------
# Sometimes we have mixed together slices and isosurfaces. Tecplot disables isosurfaces so
# this macros is necessary to reanimate them back. I also make resolution very small to save
# adjusting time.
# ------------------------------------------------------------------------------------------
$!Macrofunction Name = "=== Normalize View ==="
  ShowInMacroPanel = True
  $!FIELD  [1]
    POINTS {
      IJKSKIP {
        I = 10
        J = 10
        K = 10
      }
    }
    MESH {  
      COLOR = RED  
    }
    CONTOUR {
      CONTOURTYPE = FLOOD
      COLOR = RED
      USELIGHTINGEFFECT = YES
    }
    VECTOR {
      COLOR = RED
    }
    SCATTER {
      COLOR = RED
    }
    SHADE {
      COLOR = WHITE
    }
    BOUNDARY {
      SHOW = NO
      COLOR = RED
    }
    POINTS {
      POINTSTOPLOT = SURFACENODES
    }
    SURFACES {
      SURFACESTOPLOT = EXPOSEDCELLFACES
    }
    VOLUMEMODE {
      VOLUMEOBJECTSTOPLOT {
        SHOWISOSURFACES = YES
        SHOWSLICES = YES
        SHOWSTREAMTRACES = YES
      }
    }
    
  $!BLANKING 						# Disables blanking.
    IJK {
      INCLUDE = NO
    }
    
  $!VarSet |param_skipI| = 10				# Makes grid very coarse.
  $!VarSet |param_skipJ| = 10
  $!VarSet |param_skipK| = 10
  $!VarSet |zone_loaded| = -1				# Marks dataset as dirty (will reload).
  
############################################################  $!RESET3DSCALEFACTORS 				# ___3D_view_reset___ (perl marker).
############################################################  $!RESET3DORIGIN 					# ___3D_view_reset___ (perl marker).
############################################################    ORIGINRESETLOCATION = DATACENTER			# ___3D_view_reset___ (perl marker).
  $!VIEW FIT
  
  $!RunMacroFunction "====== Redraw ======"
$!Endmacrofunction

#
# ================================================================================
#
$!Macrofunction Name = "X-resolution"
  ShowInMacroPanel = True
  $!PromptForTextString |param_skipI|
    Instructions = "Enter the number of nodes to skip (now it is |param_skipI|)."
    
  $!VarSet |zone_loaded| = -1				# Marks dataset as dirty (will reload).
  
  $!If |param_skipI| < 1				# Checks range.
    $!VarSet |param_skipI| = 1
  $!EndIf
  $!If |param_skipI| > 30
    $!VarSet |param_skipI| = 30
  $!EndIf
  
  $!RunMacroFunction "====== Redraw ======"
$!Endmacrofunction

#
# ================================================================================
#
$!Macrofunction Name = "Y-resolution"
  ShowInMacroPanel = True
  $!PromptForTextString |param_skipJ|
    Instructions = "Enter the number of nodes to skip (now it is |param_skipJ|)."
    
  $!VarSet |zone_loaded| = -1				# Marks dataset as dirty (will reload).
  
  $!If |param_skipJ| < 1				# Checks range.
    $!VarSet |param_skipJ| = 1
  $!EndIf
  $!If |param_skipJ| > 30
    $!VarSet |param_skipJ| = 30
  $!EndIf
  
  $!RunMacroFunction "====== Redraw ======"
$!Endmacrofunction

#
# ================================================================================
#
$!Macrofunction Name = "Z-resolution"
  ShowInMacroPanel = True
  $!PromptForTextString |param_skipK|
    Instructions = "Enter the number of nodes to skip (now it is |param_skipK|)."
    
  $!VarSet |zone_loaded| = -1				# Marks dataset as dirty (will reload).
  
  $!If |param_skipK| < 1				# Checks range.
    $!VarSet |param_skipK| = 1
  $!EndIf
  $!If |param_skipK| > 30
    $!VarSet |param_skipK| = 30
  $!EndIf
  
  $!RunMacroFunction "====== Redraw ======"
$!Endmacrofunction

# ---------------------------------------------------------------------------
#   Zone navigation functions
# ---------------------------------------------------------------------------

$!Macrofunction Name = "====== Redraw ======"
  ShowInMacroPanel = True
  $!If |zone_now| < 1					# Checks ranges.
    $!VarSet |zone_now| = 1
  $!EndIf
  
  $!If |zone_now| > |zone_max|
    $!VarSet |zone_now| = (|zone_max|)
  $!EndIf
  
  $!If |zone_now| != |zone_loaded|
    $!DRAWGRAPHICS FALSE
    $!READDATASET  "../output/tecplot/link.tecplot_|zone_now|.plt"
      READDATAOPTION = NEW
      IJKSKIP
      {
        I = |param_skipI|
        J = |param_skipJ|
        K = |param_skipK|
      }
      RESETSTYLE = NO
      INCLUDETEXT = YES
      INCLUDEGEOM = NO
      INCLUDECUSTOMLABELS = NO
      VARLOADMODE = BYNAME
      VARNAMELIST = '"x [<greek>l</greek><sub>0</sub>]" "y [<greek>l</greek><sub>0</sub>]" "z [<greek>l</greek><sub>0</sub>]" "E<sub>x</sub> [a<sub>0</sub>]" "E<sub>y</sub> [a<sub>0</sub>]" "E<sub>z</sub> [a<sub>0</sub>]" "E<sup>2</sup> [a<sub>0</sub><sup>2</sup>]" "H<sub>x</sub> [a<sub>0</sub>]" "H<sub>y</sub> [a<sub>0</sub>]" "H<sub>z</sub> [a<sub>0</sub>]" "H<sup>2</sup> [a<sub>0</sub><sup>2</sup>]" "j<sub>x</sub> [a.u.]" "j<sub>y</sub> [a.u.]" "j<sub>z</sub> [a.u.]" "<greek>r</greek> [|e|n<sub>e crit</sub>]"' 
  
    $!If |flag_killMesh| != 0
      $!FIELDLAYERS SHOWMESH = NO			# Disables mesh (sometimes tecplot set it as default).
    $!EndIf
    
    $!VarSet |zone_loaded| = (|zone_now|)
  $!EndIf
  
  $!RunMacroFunction "ReSet Levels"
  $!RunMacroFunction "ReSet Title"
  $!If |flag_resetAxis| == 1
    $!RunMacroFunction "Update axis range"  
  $!EndIf
  $!DRAWGRAPHICS TRUE
  $!ReDraw
$!Endmacrofunction

#
# ================================================================================
#
$!Macrofunction Name = "Activate"
  ShowInMacroPanel = True
  $!PromptForTextString |zone_now|
    Instructions = "Enter the number of zone (1 to |zone_max|)."
  $!RunMacroFunction "====== Redraw ======"
$!Endmacrofunction

#
# ================================================================================
#
$!Macrofunction Name = "Next data set"
  ShowInMacroPanel = True
  $!VarSet |zone_now| = (|zone_now| + 1)
  $!RunMacroFunction "====== Redraw ======"
$!Endmacrofunction

#
# ================================================================================
#
$!Macrofunction Name = "Prev data set"
  ShowInMacroPanel = True
  $!VarSet |zone_now| = (|zone_now| - 1)
  $!RunMacroFunction "====== Redraw ======"
$!Endmacrofunction

#---------------------------------------------------------------------------
#     Animation of the dataset
#---------------------------------------------------------------------------

$!Macrofunction Name = "Animate frames"
  ShowInMacroPanel = True
  $!Loop |zone_max|
    $!VarSet |zone_now| = (|LOOP|)
    $!RunMacroFunction  "====== Redraw ======"
  $!EndLoop
$!Endmacrofunction


#---------------------------------------------------------------------------
#     Color levels management
#---------------------------------------------------------------------------

$!Macrofunction Name = "ReSet Levels"
    ShowInMacroPanel = False
  $!If |flag_rescaleContours| == 1		# Checks if level are not locked.
    $!If |flag_commonLevels| == 1
      $!CONTOURLEVELS DELETERANGE
        CONTOURGROUP = 1
        RANGEMIN = -1.0e10
        RANGEMAX = +1.0e10
      $!VarSet |delta| = ((|param_maxCommonVar| - |param_minCommonVar|)/(|param_levelsN|))
      $!VarSet |var| = (|param_minCommonVar| + 0.5*|delta|)
      $!Loop |param_levelsN|
        $!CONTOURLEVELS ADD
          CONTOURGROUP = 1
          RAWDATA
          1
          |var|
        $!VarSet |var| = (|var| + |delta|)
      $!EndLoop
      $!RemoveVar |var|
      $!RemoveVar |delta|
    $!Else
      $!ContourLevels Reset
        NumValues = |param_levelsN|
        ContourGroup = 1
      $!VarSet |flag_commonLevelsAreSet| = 0
    $!EndIf
  $!EndIf
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "=== Contour lines N ==="
  ShowInMacroPanel = True
  $!PromptForTextString |param_levelsN|
    Instructions = "Enter the number of contour lines to generate (2 - 30)."
  $!If |param_levelsN| < 2
    $!VarSet |param_levelsN| = 2
  $!EndIf
  $!If |param_levelsN| > 30
    $!VarSet |param_levelsN| = 30
  $!EndIf
  $!RunMacroFunction "====== Redraw ======"
$!Endmacrofunction


#---------------------------------------------------------------------------
#     Flag setting interface functions.
#---------------------------------------------------------------------------

$!Macrofunction Name = "Use separate levels"
    ShowInMacroPanel = True
  $!VarSet |flag_commonLevels| = 0
  $!RunMacroFunction "====== Redraw ======"
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Use common levels"
    ShowInMacroPanel = True
    
  $!If |flag_commonLevels| == 0
    $!VarSet |flag_commonLevels| = 1
    $!RunMacroFunction "====== Redraw ======"
  $!EndIf
$!Endmacrofunction

# -------------------------------------------------------------------------------
# Trick to get common level amplitude - Perl uses registry and updates this func.
# -------------------------------------------------------------------------------
$!Macrofunction Name = "Common levels variable"
    ShowInMacroPanel = True
    
  $!VarSet |OldVar| = (|param_globContourVar|)
  $!PromptForTextString |param_globContourVar|
    Instructions = "Confirm contour variable code (|param_globContourVar| now) (1/Ex, 2/Ey, 3/Ez, 4/E2, 5/Hx, 6/Hy, 7/Hz, 8/H2, 9/jx, 10/jy, 11/jz, 12/rho)."
    
  $!If |param_globContourVar| <= 1
    $!VarSet |param_minCommonVar| = -5.7456e-01		# ___Ex_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +5.7456e-01		# ___Ex_max___ (that is Perl search mark)
  $!EndIf
  $!If |param_globContourVar| == 2
    $!VarSet |param_minCommonVar| = -1.3085e+00		# ___Ey_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +1.3085e+00		# ___Ey_max___ (that is Perl search mark)
  $!EndIf
  $!If |param_globContourVar| == 3
    $!VarSet |param_minCommonVar| = -9.6164e-01		# ___Ez_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +9.6164e-01		# ___Ez_max___ (that is Perl search mark)
  $!EndIf
  $!If |param_globContourVar| == 4
    $!VarSet |param_minCommonVar| = +3.3884e-01		# ___E2_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +2.1115e+00		# ___E2_max___ (that is Perl search mark)
  $!EndIf

  $!If |param_globContourVar| == 5
    $!VarSet |param_minCommonVar| = -5.1706e-01		# ___Hx_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +5.1706e-01		# ___Hx_max___ (that is Perl search mark)
  $!EndIf
  $!If |param_globContourVar| == 6
    $!VarSet |param_minCommonVar| = -1.3038e+00		# ___Hy_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +1.3038e+00		# ___Hy_max___ (that is Perl search mark)
  $!EndIf
  $!If |param_globContourVar| == 7
    $!VarSet |param_minCommonVar| = -9.6783e-01		# ___Hz_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +9.6783e-01		# ___Hz_max___ (that is Perl search mark)
  $!EndIf
  $!If |param_globContourVar| == 8
    $!VarSet |param_minCommonVar| = +2.9004e-01		# ___H2_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +2.0552e+00		# ___H2_max___ (that is Perl search mark)
  $!EndIf

  $!If |param_globContourVar| == 9
    $!VarSet |param_minCommonVar| = +0.0000e+00		# ___jx_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +2.2251e-308		# ___jx_max___ (that is Perl search mark)
  $!EndIf
  $!If |param_globContourVar| == 10
    $!VarSet |param_minCommonVar| = +0.0000e+00		# ___jy_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +2.2251e-308		# ___jy_max___ (that is Perl search mark)
  $!EndIf
  $!If |param_globContourVar| == 11
    $!VarSet |param_minCommonVar| = +0.0000e+00		# ___jz_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +2.2251e-308		# ___jz_max___ (that is Perl search mark)
  $!EndIf
  $!If |param_globContourVar| >= 12
    $!VarSet |param_minCommonVar| = +0.0000e+00		# ___rho_min___ (that is Perl search mark)
    $!VarSet |param_maxCommonVar| = +2.2251e-308		# ___rho_max___ (that is Perl search mark)
  $!EndIf
  
  $!If |param_globContourVar| != |oldVar|
    $!RunMacroFunction "====== Redraw ======"
  $!EndIf
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Reset levels on step"
    ShowInMacroPanel = True
  
  $!If |flag_rescaleContours| == 0
    $!VarSet |flag_rescaleContours| = 1
    $!RunMacroFunction "====== Redraw ======"
  $!EndIf
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Keep levels on step"
    ShowInMacroPanel = True
  $!VarSet |flag_rescaleContours| = 0
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Reset axis on step"
    ShowInMacroPanel = True
  $!If |flag_resetAxis| == 0
    $!VarSet |flag_resetAxis| = 1
    $!RunMacroFunction "====== Redraw ======"
  $!EndIf
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Keep axis on step"
    ShowInMacroPanel = True
  $!VarSet |flag_resetAxis| = 0
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Update axis range"
    ShowInMacroPanel = False
  $!View AxisFit
    Axis = 'Z'
$!Endmacrofunction


#---------------------------------------------------------------------------
#     Zone name management
#---------------------------------------------------------------------------

$!Macrofunction Name = "ReSet Title"
    ShowInMacroPanel = False
  $!RunMacroFunction "Remove title"
  $!RunMacroFunction "Add title"
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Add title"
    ShowInMacroPanel = True
  $!ATTACHTEXT
    XYPOS
      {
      X = 64
      Y = 95
      }
    TEXTSHAPE
      {
      FONT = TIMES
      HEIGHT = 25
      }
    TEXT = '&(ZONENAME:1)'
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Remove title"
    ShowInMacroPanel = True
  $!PICK ADD
    X = 8.1838028169
    Y = 0.593306933011
  $!PICK CUT
$!Endmacrofunction



#---------------------------------------------------------------------------
#     Export of the figures
#---------------------------------------------------------------------------

$!Macrofunction Name = "Change Output Name"
  ShowInMacroPanel = True
  $!PromptForTextString |param_prefix|
    Instructions = "Enter the new file prefix (old is |param_prefix|)."
  $!VarSet |flag_nameIsSet| = 1
$!Endmacrofunction

#---------------------------------------------------------------------------
#     Export of the PNG-figures
#---------------------------------------------------------------------------

$!Macrofunction Name = "Export (PNG)"
  ShowInMacroPanel = True
  $!VarSet |Fix| = ""
  $!If |zone_now| < 100
    $!VarSet |Fix| = "0"
  $!EndIf
  $!If |zone_now| < 10
    $!VarSet |Fix| = "00"
  $!EndIf

  $!ExportSetup ExportFormat = PNG
  $!ExportSetup ImageWidth = 1000
  $!ExportSetup ExportFName = "../tmp/|param_prefix|_|Fix||zone_now|.png"
  $!Export
  
  $!RemoveVar |FIX|
  
  $!RunMacroFunction "Next data set"
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Export All (PNG)"
  ShowInMacroPanel = True
  
  $!If |flag_nameIsSet| != 1
    $!RunMacroFunction "Change Output Name"
  $!Endif
  
  $!VarSet |zone_now| = 1
  $!RunMacroFunction "====== Redraw ======"
  $!Loop |zone_max|
    $!RunMacroFunction "Export (PNG)"
  $!EndLoop
  $!VarSet |flag_nameIsSet| = 0			# Protects exported figures.
$!Endmacrofunction

#---------------------------------------------------------------------------
#     Export of the PNG slice series.
#---------------------------------------------------------------------------

$!Macrofunction Name = "Configure slice export"
  ShowInMacroPanel = True

  $!PromptForTextString |slices_N|
    Instructions = "Enter number of slices in sequence:"
    
  $!PromptForTextString |slices_axis|
    Instructions = "Enter slice axis (1 => x, 2 => y, 3 => z):"
    
  $!If |slices_axis| <= 1			# With clamp against bad axis number.
    $!VarSet |slices_axis| = 1
  $!EndIf
  
  $!If |slices_axis| == 2
    $!GLOBALSLICE SLICESURFACE = YPLANES
  $!EndIf
  
  $!If |slices_axis| >= 3			# With clamp against bad axis number.
    $!VarSet |slices_axis| = 3
  $!EndIf
  
  $!If |slices_N| <= 5				# Clamp against bad number of slices.
    $!VarSet |slices_N| = 5
  $!EndIf
  
  $!If |slices_N| >= 50				# Clamp against bad number of slices.
    $!VarSet |slices_N| = 50
  $!EndIf
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Export slices for current zone"
  ShowInMacroPanel = True
  $!VarSet |Fix| = ""				# Prefix to keep width of the zone number the same.
  $!If |zone_now| < 100				# It is necessary for proper sorting in the slideshow.
    $!VarSet |Fix| = "0"
  $!EndIf
  $!If |zone_now| < 10
    $!VarSet |Fix| = "00"
  $!EndIf

  $!If |slices_N| < 5
    $!RunMacroFunction "Configure slice export"  
  $!Endif
  
  $!If |slices_axis| <= 1			# With clamp against bad axis number.
    $!GLOBALSLICE SLICESURFACE = XPLANES
    $!VarSet |MinOfPos| = |AXISMINX|
    $!VarSet |MaxOfPos| = |AXISMAXX|
  $!EndIf
  
  $!If |slices_axis| == 2
    $!GLOBALSLICE SLICESURFACE = YPLANES
    $!VarSet |MinOfPos| = |AXISMINY|
    $!VarSet |MaxOfPos| = |AXISMAXY|
  $!EndIf
  
  $!If |slices_axis| >= 3			# With clamp against bad axis number.
    $!GLOBALSLICE SLICESURFACE = ZPLANES
    $!VarSet |MinOfPos| = |AXISMINZ|
    $!VarSet |MaxOfPos| = |AXISMAXZ|
  $!EndIf
  
  $!GLOBALSLICE SHOW = YES
  $!GLOBALSLICE SHOWINTERMEDIATESLICES = NO
  $!GLOBALSLICE SHOWPOSITION2 = NO

  $!VarSet |dPos| = ((|MaxOfPos| - |MinOfPos|)/|slices_N|)
  $!VarSet |pos| = (|MinOfPos| + 0.5*|dPos|)
  
  $!RunMacroFunction "====== Redraw ======"
    
  $!Loop |slices_N|
    $!If |slices_axis| <= 1			# With clamp against bad axis number.
      $!GLOBALSLICE POSITION1{X = |pos|}
    $!EndIf
    
    $!If |slices_axis| == 2
      $!GLOBALSLICE POSITION1{Y = |pos|}
    $!EndIf
    
    $!If |slices_axis| >= 3			# With clamp against bad axis number.
      $!GLOBALSLICE POSITION1{Z = |pos|}
    $!EndIf
    
    $!VarSet |Fix2| = ""			# Prefix to keep width of the zone number the same.
    $!If |LOOP| < 100				# It is necessary for proper sorting in the slideshow.
      $!VarSet |Fix2| = "0"
    $!EndIf
    $!If |LOOP| < 10
      $!VarSet |Fix2| = "00"
    $!EndIf
    
    $!ExportSetup ExportFormat = PNG
    $!ExportSetup ImageWidth = 1000
    $!ExportSetup ExportFName = "../tmp/|param_prefix|_|Fix||zone_now|_|Fix2||LOOP|.png"
    $!Export
    
    $!VarSet |pos| = (|pos| + |dPos|)
  $!EndLoop
  
  $!RemoveVar |MinOfPos|
  $!RemoveVar |MaxOfPos|
  $!RemoveVar |Fix|
  $!RemoveVar |Fix2|
  $!RemoveVar |pos|
  $!RemoveVar |dPos|
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Export slices for all zones"
  ShowInMacroPanel = True

  $!If |flag_nameIsSet| != 1
    $!RunMacroFunction "Change Output Name"
  $!Endif
  
  $!VarSet |zone_now| = 1
  $!RunMacroFunction "====== Redraw ======"

  $!Loop |zone_max|
    $!RunMacroFunction "Export slices for current zone"
    $!RunMacroFunction "Next data set"
  $!EndLoop
  
  $!VarSet |flag_nameIsSet| = 0			# Protects exported figures.
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Disable automesh"
    ShowInMacroPanel = True
  $!VarSet |flag_killMesh| = 1
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Enable automesh"
    ShowInMacroPanel = True
  $!VarSet |flag_killMesh| = 0
  $!RunMacroFunction "====== Redraw ======"
$!Endmacrofunction
