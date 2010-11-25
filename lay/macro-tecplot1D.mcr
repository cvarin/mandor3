#!MC 1000

# ================================================================================
# This script is used to manage set of the preplotted files by reloading them
# automatically using this macro-menu to browse data. Main job of file preparation
# is done by vTec_preplot.pl script. See the tail of this script to get details.
# ================================================================================

#
# Parameters to set range of zones (actually, symlinks to prepared data-sets).
#
$!VarSet |zone_now| = 1

#
# Different flag to set state of plot when new zone is activated.
#
$!VarSet |flag_resetAxis| = 0				# Useful for 1D tecplot meshes.
$!VarSet |flag_nameIsSet| = 0				# Useful to remind about proper name for export.
$!VarSet |param_prefix| = "tec3D"			# Export file prefix.

# ---------------------------------------------------------------------------
#   Zone navigation functions
# ---------------------------------------------------------------------------

$!Macrofunction Name = "====== Redraw ======"
  ShowInMacroPanel = True
  $!If |zone_now| < 1					# Checks ranges.
    $!VarSet |zone_now| = 1
  $!EndIf
  
  $!If |zone_now| > |NUMZONES|
    $!VarSet |zone_now| = (|NUMZONES|)
  $!EndIf
  
  $!ACTIVEFIELDZONES += [|zone_now|]			# Activates only chosen zone.
  $!Loop |NumZones|
    $!If |Loop| != |zone_now|
      $!ACTIVEFIELDZONES -= [|Loop|]
    $!EndIf
  $!Endloop

  $!RunMacroFunction "ReSet Title"
  $!If |flag_resetAxis| == 1
    $!RunMacroFunction "Update axise range"  
  $!EndIf
  $!ReDraw
$!Endmacrofunction

#
# ================================================================================
#
$!Macrofunction Name = "Activate"
  ShowInMacroPanel = True
  $!PromptForTextString |zone_now|
    Instructions = "Enter the number of zone (1 to |NUMZONES|)."
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
  $!Loop |NUMZONES|
    $!VarSet |zone_now| = (|LOOP|)
    $!RunMacroFunction  "====== Redraw ======"
  $!EndLoop
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Reset axis range on step"
    ShowInMacroPanel = True
  $!VarSet |flag_resetAxis| = 1
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Keep axis range on step"
    ShowInMacroPanel = True
  $!VarSet |flag_resetAxis| = 0
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Update axise range"
    ShowInMacroPanel = True
  $!View AxisFit
    Axis = 'Y'
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
    TEXT = '&(ZONENAME:|zone_now|)'
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
  $!Loop |NUMZONES|
    $!RunMacroFunction "Export (PNG)"
  $!EndLoop
  $!VarSet |flag_nameIsSet| = 0			# Protects exported figures.
$!Endmacrofunction
