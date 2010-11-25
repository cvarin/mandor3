#!MC 900

$!VarSet |CurrentZoneNum| = 1
$!VarSet |doContourLevels| = 0
$!VarSet |numberOfContourLevels| = 15
$!VarSet |commonLevels| = 1
$!VarSet |commonLevelsAreSet| = 0
$!VarSet |nextPrevUpdate| = 1
$!VarSet |nextPrevAxisScale| = 0
$!VarSet |IAmInMovie| = 0
$!VarSet |prefix| = 'tec'

# ---------------------------------------------------------------------------
#   Zone navigation functions
# ---------------------------------------------------------------------------

$!Macrofunction Name = "reActivate"
  ShowInMacroPanel = False
  $!If |CurrentZoneNum| < 1
    $!VarSet |CurrentZoneNum| = 1
  $!EndIf
  $!If |CurrentZoneNum| > |NumZones|
    $!VarSet |CurrentZoneNum| = (|NumZones|)
  $!EndIf
  $!RunMacroFunction "ReSet Levels"
  $!RunMacroFunction "ReSet Title"
  $!If |nextPrevAxisScale| == 1
    $!RunMacroFunction "Update axises range to var min/max"  
  $!EndIf
  $!ReDraw
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Activate"
  ShowInMacroPanel = True
  $!PromptForTextString |CurrentZoneNum|
    Instructions = "Enter the number of zone to activate."
  $!RunMacroFunction "reActivate"
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Next data set"
  ShowInMacroPanel = True
  $!VarSet |CurrentZoneNum| = (|CurrentZoneNum| + 1)
  $!RunMacroFunction "reActivate"
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Prev data set"
  ShowInMacroPanel = True
  $!VarSet |CurrentZoneNum| = (|CurrentZoneNum| - 1)
  $!RunMacroFunction "reActivate"
$!Endmacrofunction



#---------------------------------------------------------------------------
#     Animation of the dataset
#---------------------------------------------------------------------------

$!Macrofunction Name = "Animate frames"
  ShowInMacroPanel = True
  $!VarSet |IAmInMovie| = 1
  $!Loop |NumZones|
    $!VarSet |CurrentZoneNum| = (|LOOP|)
    $!RunMacroFunction  "reActivate"
  $!EndLoop
  $!VarSet |IAmInMovie| = 0
$!Endmacrofunction



#---------------------------------------------------------------------------
#     Color levels management
#---------------------------------------------------------------------------

$!Macrofunction Name = "ReSet Levels"
    ShowInMacroPanel = False
  $!If |CurrentZoneNum| > |NumZones|
    $!VarSet |CurrentZoneNum| = (|CurrentZoneNum| - 1)
  $!EndIf

  # Stupid "OR" logic
  $!VarSet |INeedToReset| = 0
  $!If |nextPrevUpdate| == 1
    $!VarSet |INeedToReset| = 1
  $!EndIf
  $!If |IAmInMovie| == 1
    $!VarSet |INeedToReset| = 1
  $!EndIf

  $!If |doContourLevels| == 0
    $!VarSet |INeedToReset| = 0
  $!EndIf
  
  $!If |INeedToReset| == 1
    $!If |commonLevels| == 1
      $!If |commonLevelsAreSet| == 0
        $!Loop |NumZones|
          $!ACTIVEFIELDZONES += [|Loop|]
        $!Endloop
        $!ContourLevels Reset
          NumValues = |numberOfContourLevels|
          ContourGroup = 1
        $!VarSet |commonLevelsAreSet| = 1
      $!EndIf
    $!EndIf
  $!EndIf

  $!ACTIVEFIELDZONES += [|currentZoneNum|]
  $!Loop |NumZones|
    $!If |Loop| != |currentZoneNum|
      $!ACTIVEFIELDZONES -= [|Loop|]
    $!EndIf
  $!Endloop

  $!If |INeedToReset| == 1
    $!If |commonLevels| == 0
      $!ContourLevels Reset
        NumValues = |numberOfContourLevels|
        ContourGroup = 1
      $!VarSet |commonLevelsAreSet| = 0
    $!EndIf
  $!EndIf
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Do contours"
    ShowInMacroPanel = True
  $!VarSet |doContourLevels| = 1
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Coutour lines N"
  ShowInMacroPanel = True
  $!PromptForTextString |numberOfContourLevels|
    Instructions = "Enter the number of contour lines to generate."
  $!If |numberOfContourLevels| < 2
    $!VarSet |numberOfContourLevels| = 2
  $!EndIf
  $!If |numberOfContourLevels| > 30
    $!VarSet |numberOfContourLevels| = 30
  $!EndIf
  $!RunMacroFunction "reActivate"
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Do NOT contours"
    ShowInMacroPanel = True
  $!VarSet |doContourLevels| = 0
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Use common levels"
    ShowInMacroPanel = True
  $!VarSet |commonLevels| = 1
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Use separate levels"
    ShowInMacroPanel = True
  $!VarSet |commonLevels| = 0
  $!RunMacroFunction "reActivate"
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Update levels on next/prev"
    ShowInMacroPanel = True
  $!VarSet |nextPrevUpdate| = 1
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Don't update levels on next/prev"
    ShowInMacroPanel = True
  $!VarSet |nextPrevUpdate| = 0
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Update axises range on next/prev"
    ShowInMacroPanel = True
  $!VarSet |nextPrevAxisScale| = 1
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Don't update axises range on next/prev"
    ShowInMacroPanel = True
  $!VarSet |nextPrevAxisScale| = 0
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Update axises range to var min/max"
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
      X = 74
      Y = 95
      }
    TEXTSHAPE
      {
      FONT = TIMES
      HEIGHT = 28
      }
    TEXT = '&(ZONENAME:|CurrentZoneNum|)'
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
#     Export of the EPS-figures
#---------------------------------------------------------------------------

$!Macrofunction Name = "Change Output Name"
  ShowInMacroPanel = True
  $!PromptForTextString |prefix|
    Instructions = "Enter the new file(s) prefix."
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Export (EPS)"
  ShowInMacroPanel = True
  $!VarSet |Fix| = ""
  $!If |currentZoneNum| < 100
    $!VarSet |Fix| = "0"
  $!EndIf
  $!If |currentZoneNum| < 10
    $!VarSet |Fix| = "00"
  $!EndIf
  $!ExportSetup EpsPreviewImage{IMAGETYPE = NONE}
  $!PrintSetup Precision = 1
  $!PrintSetup Palette = Color
  $!ExportSetup ExportFormat = EPS
  $!ExportSetup ExportFName = "../tmp/|prefix|_|Fix||currentZoneNum|.eps"
  $!Export
  $!RunMacroFunction "Next data set"
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Export All (EPS)"
  ShowInMacroPanel = True
  $!VarSet |CurrentZoneNum| = 1
  $!RunMacroFunction "reActivate"
  $!Loop |NumZones|
  $!RunMacroFunction "Export (EPS)"
  $!EndLoop
$!Endmacrofunction



#---------------------------------------------------------------------------
#     Export of the PNG-figures
#---------------------------------------------------------------------------

$!Macrofunction Name = "Export (PNG)"
  ShowInMacroPanel = True
  $!VarSet |Fix| = ""
  $!If |currentZoneNum| < 100
    $!VarSet |Fix| = "0"
  $!EndIf
  $!If |currentZoneNum| < 10
    $!VarSet |Fix| = "00"
  $!EndIf

  $!ExportSetup ExportFormat = PNG
  $!ExportSetup ImageWidth = 1000
  $!ExportSetup ExportFName = "../tmp/|prefix|_|Fix||currentZoneNum|.png"
  $!Export
  
  $!RemoveVar |FIX|
  
  $!RunMacroFunction "Next data set"
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Export All (PNG)"
  ShowInMacroPanel = True
  $!VarSet |CurrentZoneNum| = 1
  $!RunMacroFunction "reActivate"
  $!Loop |NumZones|
  $!RunMacroFunction "Export (PNG)"
  $!EndLoop
$!Endmacrofunction


#---------------------------------------------------------------------------
#     Export of the PNG-figures
#---------------------------------------------------------------------------

$!Macrofunction Name = "Export to AVI"
  ShowInMacroPanel = True
  $!EXPORTSETUP EXPORTFORMAT = AVI
  $!EXPORTSETUP BITDUMPREGION = ALLFRAMES
  $!EXPORTSETUP ANIMATIONSPEED = 4
  $!EXPORTSETUP IMAGEWIDTH = 800
  $!EXPORTSETUP EXPORTFNAME = '../tmp/|prefix|.avi'
  $!EXPORTSTART

  $!VarSet |CurrentZoneNum| = 1
  $!RunMacroFunction "reActivate"

  $!Loop |NumZones|
  $!ExportNextFrame
  $!RUNMACROFUNCTION  "Next data set"
  $!EndLoop

  $!EXPORTFINISH
$!Endmacrofunction


#---------------------------------------------------------------------------
#     Export of the PNG-slices
#---------------------------------------------------------------------------

$!VarSet |NumSlices| = 10
$!VarSet |SliceAxis| = 1

$!Macrofunction Name = "Configure slice export"
  ShowInMacroPanel = True

  $!PromptForTextString |NumSlices|
    Instructions = "Enter number of slices in sequence:"
    
  $!PromptForTextString |SliceAxis|
    Instructions = "Enter slice axis (1 => x, 2 => y, 3 => z):"
    
  $!If |SliceAxis| <= 1			# With clamp against bad axis number.
    $!VarSet |SliceAxis| = 1
  $!EndIf
  
  $!If |SliceAxis| == 2
    $!GLOBALSLICE SLICESURFACE = YPLANES
  $!EndIf
  
  $!If |SliceAxis| >= 3			# With clamp against bad axis number.
    $!VarSet |SliceAxis| = 3
  $!EndIf
  
  $!If |NumSlices| <= 5			# Clamp against bad number of slices.
    $!VarSet |NumSlices| = 5
  $!EndIf
  
  $!If |NumSlices| >= 50		# Clamp against bad number of slices.
    $!VarSet |NumSlices| = 50
  $!EndIf
$!Endmacrofunction

#---------------------------------------------------------------------------

$!Macrofunction Name = "Export slices for current zone"
  ShowInMacroPanel = True
  $!VarSet |Fix| = ""			# Prefix to keep width of the zone number the same.
  $!If |currentZoneNum| < 100		# It is necessary for proper sorting in the slideshow.
    $!VarSet |Fix| = "0"
  $!EndIf
  $!If |currentZoneNum| < 10
    $!VarSet |Fix| = "00"
  $!EndIf

  $!FIELDLAYERS SHOWCONTOUR = NO	# Disables contours.
  
  $!If |SliceAxis| <= 1			# With clamp against bad axis number.
    $!GLOBALSLICE SLICESURFACE = XPLANES
    $!VarSet |MinOfPos| = |AXISMINX|
    $!VarSet |MaxOfPos| = |AXISMAXX|
  $!EndIf
  
  $!If |SliceAxis| == 2
    $!GLOBALSLICE SLICESURFACE = YPLANES
    $!VarSet |MinOfPos| = |AXISMINY|
    $!VarSet |MaxOfPos| = |AXISMAXY|
  $!EndIf
  
  $!If |SliceAxis| >= 3			# With clamp against bad axis number.
    $!GLOBALSLICE SLICESURFACE = ZPLANES
    $!VarSet |MinOfPos| = |AXISMINZ|
    $!VarSet |MaxOfPos| = |AXISMAXZ|
  $!EndIf
  
  $!GLOBALSLICE SHOW = YES
  $!GLOBALSLICE SHOWINTERMEDIATESLICES = NO
  $!GLOBALSLICE SHOWPOSITION2 = NO

  $!VarSet |dPos| = ((|MaxOfPos| - |MinOfPos|)/|NumSlices|)
  $!VarSet |pos| = (|MinOfPos| + 0.5*|dPos|)
  
  $!RunMacroFunction "reActivate"
    
  $!Loop |NumSlices|
    $!If |SliceAxis| <= 1		# With clamp against bad axis number.
      $!GLOBALSLICE POSITION1{X = |pos|}
    $!EndIf
    
    $!If |SliceAxis| == 2
      $!GLOBALSLICE POSITION1{Y = |pos|}
    $!EndIf
    
    $!If |SliceAxis| >= 3		# With clamp against bad axis number.
      $!GLOBALSLICE POSITION1{Z = |pos|}
    $!EndIf
    
    $!VarSet |Fix2| = ""		# Prefix to keep width of the zone number the same.
    $!If |LOOP| < 100			# It is necessary for proper sorting in the slideshow.
      $!VarSet |Fix2| = "0"
    $!EndIf
    $!If |LOOP| < 10
      $!VarSet |Fix2| = "00"
    $!EndIf
    
    $!ExportSetup ExportFormat = PNG
    $!ExportSetup ImageWidth = 1000
    $!ExportSetup ExportFName = "../tmp/|prefix|_|Fix||currentZoneNum|_|Fix2||LOOP|.png"
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
  
  $!RunMacroFunction "Configure slice export"
  
  $!VarSet |CurrentZoneNum| = 1
  $!RunMacroFunction "reActivate"
  
  $!Loop |NumZones|
    $!RunMacroFunction "Export slices for current zone"
    $!RunMacroFunction "Next data set"
  $!EndLoop
$!Endmacrofunction
