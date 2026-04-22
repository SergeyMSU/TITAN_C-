#!MC 1410
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A5\Tecplot_Tecplot_print_2D__2d_(0, 0, 1, 0)_.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"xx" "yy" "X" "Y" "Z" "rho" "p" "Vx" "Vy" "Vz" "Bx" "By" "Bz" "Q" "BB/8pi" "zone_geo" "Mach" "Mach_Alf" "|B|/8pi" "rho_Th" "p_Th" "T_Th" "J" "Jxx" "Jyy" "Jx" "Jy" "Jz" "J_an" "Jxx_an" "Jyy_an" "Bxx_pot" "Byy_pot"'
$!TwoDAxis XDetail{VarNum = 3}
$!TwoDAxis YDetail{VarNum = 4}
$!RedrawAll 
$!Pick AddAtPosition
  X = 0.945643153527
  Y = 4.03423236515
  ConsiderStyle = Yes
$!FrameControl ActivateByNumber
  Frame = 1
$!Pick Shift
  X = -2.58921161826
  Y = 0
  PickSubposition = Left
$!Pick Shift
  X = 2.57593360996
  Y = 0
  PickSubposition = Right
$!RedrawAll 
$!GlobalRGB RedChannelVar = 6
$!GlobalRGB GreenChannelVar = 1
$!GlobalRGB BlueChannelVar = 11
$!SetContourVar 
  Var = 1
  ContourGroup = 1
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 2
  ContourGroup = 2
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 6
  ContourGroup = 3
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 7
  ContourGroup = 4
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 8
  ContourGroup = 5
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 9
  ContourGroup = 6
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 10
  ContourGroup = 7
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 11
  ContourGroup = 8
  LevelInitMode = ResetToNice
$!FieldLayers ShowContour = Yes
$!SetContourVar 
  Var = 6
  ContourGroup = 1
  LevelInitMode = ResetToNice
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0
$!ContourLevels New
  ContourGroup = 1
  RawData
5
0
0.75
1.5
2.25
3
$!RedrawAll 
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.00106712061805
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.000362486897137
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.000198839245055
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
8.97859236857E-05
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.000251796739412
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.000272303727288
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.00029203172981
$!TwoDAxis 
  XDetail
    {
    RangeMin = -373.75283845476684519
    RangeMax = 478.23494920757076443
    }
  YDetail
    {
    RangeMin = -271.83743465232402059
    RangeMax = 222.16806852185135313
    }
$!GlobalTwoDVector UVar = 8
$!GlobalTwoDVector VVar = 9
$!StreamtraceLayers Show = Yes
$!Streamtrace ResetDeltaTime
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -211.6623573498199
    Y = 122.9250596393097
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -211.0707800436531
    Y = 97.60555093537084
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -210.4792027374863
    Y = 72.28604223143199
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -209.8876254313195
    Y = 46.96653352749313
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -209.2960481251527
    Y = 21.64702482355428
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -208.7044708189859
    Y = -3.672483880384604
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -208.1128935128191
    Y = -28.99199258432344
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -207.5213162066524
    Y = -54.31150128826229
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -206.9297389004855
    Y = -79.63100999220114
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -206.3381615943188
    Y = -104.95051869614
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 189.7820026149678
    Y = 89.38262637965238
    }
  AltStartPos
    {
    X = -206.3381615943188
    Y = -104.95051869614
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 190.669368574218
    Y = 70.21552165984819
    }
  AltStartPos
    {
    X = -206.3381615943188
    Y = -104.95051869614
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 191.5567345334682
    Y = 51.04841694004401
    }
  AltStartPos
    {
    X = -206.3381615943188
    Y = -104.95051869614
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 192.4441004927184
    Y = 31.88131222023983
    }
  AltStartPos
    {
    X = -206.3381615943188
    Y = -104.95051869614
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 193.3314664519686
    Y = 12.71420750043565
    }
  AltStartPos
    {
    X = -206.3381615943188
    Y = -104.95051869614
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 194.2188324112188
    Y = -6.452897219368552
    }
  AltStartPos
    {
    X = -206.3381615943188
    Y = -104.95051869614
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 195.106198370469
    Y = -25.62000193917272
    }
  AltStartPos
    {
    X = -206.3381615943188
    Y = -104.95051869614
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 195.9935643297192
    Y = -44.78710665897691
    }
  AltStartPos
    {
    X = -206.3381615943188
    Y = -104.95051869614
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 196.8809302889694
    Y = -63.95421137878108
    }
  AltStartPos
    {
    X = -206.3381615943188
    Y = -104.95051869614
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 197.7682962482196
    Y = -83.12131609858528
    }
  AltStartPos
    {
    X = -206.3381615943188
    Y = -104.95051869614
    }
$!Linking BetweenFrames{LinkAxisPosition = Yes}
$!Linking BetweenFrames{LinkYAxisRange = Yes}
$!Linking BetweenFrames{LinkXAxisRange = Yes}
$!Linking BetweenFrames{LinkSolutionTime = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!Linking BetweenFrames{LinkContourLevels = Yes}
$!CreateNewFrame 
  XYPos
    {
    X = -1.0129
    Y = 0.28983
    }
  Width = 10.861
  Height = 7.5353
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A5\TS_srez__2d_(0, 0, 1, 0)_.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"xx" "yy" "X" "Y"'
$!FieldLayers ShowMesh = Yes
$!RedrawAll 
$!TwoDAxis XDetail{VarNum = 3}
$!TwoDAxis YDetail{VarNum = 4}
$!Linking BetweenFrames{LinkAxisPosition = Yes}
$!Linking BetweenFrames{LinkYAxisRange = Yes}
$!Linking BetweenFrames{LinkXAxisRange = Yes}
$!Linking BetweenFrames{LinkSolutionTime = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!FrameLayout IsTransparent = Yes
$!HOESettings NumSubdivisionLevels = 2
$!FieldMap [1]  Mesh{LineThickness = 0.400000000000000022}
$!RedrawAll 
$!CreateNewFrame 
  XYPos
    {
    X = -1.2452
    Y = 0.51556
    }
  Width = 12.7
  Height = 7.2896
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A5\HP_srez__2d_(0, 0, 1, 0)_.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"xx" "yy" "X" "Y"'
$!FieldLayers ShowMesh = Yes
$!TwoDAxis XDetail{VarNum = 3}
$!TwoDAxis YDetail{VarNum = 4}
$!RedrawAll 
$!Linking BetweenFrames{LinkAxisPosition = Yes}
$!Linking BetweenFrames{LinkYAxisRange = Yes}
$!Linking BetweenFrames{LinkXAxisRange = Yes}
$!Linking BetweenFrames{LinkSolutionTime = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!FieldMap [1]  Mesh{LineThickness = 0.400000000000000022}
$!FrameName = 'HP-1'
$!FrameLayout IsTransparent = Yes
$!RedrawAll 
$!CreateNewFrame 
  XYPos
    {
    X = -1.544
    Y = 0.37614
    }
  Width = 14.041
  Height = 7.8407
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A5\BS_srez__2d_(0, 0, 1, 0)_.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"xx" "yy" "X" "Y"'
$!FieldLayers ShowMesh = Yes
$!TwoDAxis XDetail{VarNum = 3}
$!TwoDAxis YDetail{VarNum = 4}
$!Linking BetweenFrames{LinkAxisPosition = Yes}
$!Linking BetweenFrames{LinkYAxisRange = Yes}
$!Linking BetweenFrames{LinkXAxisRange = Yes}
$!Linking BetweenFrames{LinkSolutionTime = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!FieldMap [1]  Mesh{LineThickness = 0.400000000000000022}
$!FrameLayout IsTransparent = Yes
$!FrameName = 'BS-1'
$!RedrawAll 
$!Pick SetMouseMode
  MouseMode = Select
$!TwoDAxis 
  XDetail
    {
    RangeMin = -113.614402480321246
    RangeMax = 218.09651323312522209
    }
  YDetail
    {
    RangeMin = -121.00215177684879109
    RangeMax = 71.332785646376180466
    }
$!FrameControl ActivateByNumber
  Frame = 2
$!FieldLayers ShowMesh = No
$!FieldLayers ShowMesh = Yes
$!FieldLayers ShowMesh = No
$!FieldLayers ShowMesh = Yes
$!FrameName = 'TS-1'
$!FrameControl ActivateByNumber
  Frame = 1
$!RedrawAll 
$!FrameName = 'rho 2D'
$!RedrawAll 
$!FrameControl MoveToTopByNumber
  Frame = 1
$!TwoDAxis 
  XDetail
    {
    RangeMin = -92.687027812395314186
    RangeMax = 197.16913856519929027
    }
  YDetail
    {
    RangeMin = -108.86789539119212122
    RangeMax = 59.198529260719482181
    }
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
3.81647571949
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
3.99906547675
$!ContourLevels DeleteNearest
  ContourGroup = 1
  RangeMin = 3.99906547675
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
3.99808609389
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
4.29644023227
$!FrameControl MoveToBottomByNumber
  Frame = 4
$!FrameControl ActivateAtPosition
  X = 3.10995850622
  Y = 5.94626556017
$!Pick SetMouseMode
  MouseMode = Select
$!Pick AddAtPosition
  X = 3.10995850622
  Y = 5.94626556017
  CollectingObjectsMode = HomogeneousAdd
  ConsiderStyle = Yes
$!Pick SetMouseMode
  MouseMode = Select
$!TwoDAxis 
  XDetail
    {
    RangeMin = -143.81658171874340724
    RangeMax = 248.29869247154732648
    }
  YDetail
    {
    RangeMin = -138.51419005684641661
    RangeMax = 88.844823926373862832
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -261.29504641477285531
    RangeMax = 365.77715716757677455
    }
  YDetail
    {
    RangeMin = -206.63137307836203149
    RangeMax = 156.96200694788942087
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -320.44550704861956092
    RangeMax = 424.92761780142348016
    }
  YDetail
    {
    RangeMin = -240.92840584286042827
    RangeMax = 191.25903971238781764
    }
$!FrameControl ActivateByNumber
  Frame = 1
$!FieldLayers ShowScatter = Yes
$!RedrawAll 
$!RedrawAll 
$!FieldLayers ShowScatter = No
$!RedrawAll 
