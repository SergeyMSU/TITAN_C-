#!MC 1410
$!Pick AddAtPosition
  X = 0.646887966805
  Y = 3.19771784232
  ConsiderStyle = Yes
$!Pick AddAtPosition
  X = 1.04522821577
  Y = 3.66244813278
  ConsiderStyle = Yes
$!FrameControl ActivateByNumber
  Frame = 1
$!Pick Shift
  X = -2.6622406639
  Y = 0
  PickSubposition = Left
$!Pick Shift
  X = 2.59585062241
  Y = 0
  PickSubposition = Right
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A6\Tecplot_Tecplot_print_2D__2d_(0, 0, 1, 0)_.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"xx" "yy" "X" "Y" "Z" "rho" "p" "Vx" "Vy" "Vz" "Bx" "By" "Bz" "Q" "rho_He" "rho_H1" "Vx_H1" "Vy_H1" "Vz_H1" "p_H1" "rho_H2" "Vx_H2" "Vy_H2" "Vz_H2" "p_H2" "rho_H3" "Vx_H3" "Vy_H3" "Vz_H3" "p_H3" "rho_H4" "Vx_H4" "Vy_H4" "Vz_H4" "p_H4" "MK_n_H" "MK_IVx_H" "MK_IVy_H" "MK_IVz_H" "MK_IT_H" "MK_n_H1" "MK_Vx_H1" "MK_Vy_H1" "MK_Vz_H1" "MK_T_H1" "MK_IVx_H1" "MK_IVy_H1" "MK_IVz_H1" "MK_IT_H1" "MK_n_H2" "MK_Vx_H2" "MK_Vy_H2" "MK_Vz_H2" "MK_T_H2" "MK_IVx_H2" "MK_IVy_H2" "MK_IVz_H2" "MK_IT_H2" "MK_n_H3" "MK_Vx_H3" "MK_Vy_H3" "MK_Vz_H3" "MK_T_H3" "MK_IVx_H3" "MK_IVy_H3" "MK_IVz_H3" "MK_IT_H3" "MK_n_H4" "MK_Vx_H4" "MK_Vy_H4" "MK_Vz_H4" "MK_T_H4" "MK_IVx_H4" "MK_IVy_H4" "MK_IVz_H4" "MK_IT_H4" "BB/8pi" "zone_geo" "Mach" "Mach_Alf" "|B|/8pi" "rho_Th" "p_Th" "T_Th" "J" "Jxx" "Jyy" "Jx" "Jy" "Jz" "J_an" "Jxx_an" "Jyy_an" "Bxx_pot" "Byy_pot"'
$!TwoDAxis XDetail{VarNum = 3}
$!TwoDAxis YDetail{VarNum = 4}
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
1.25
2.5
3.75
5
$!RedrawAll 
$!Pick AddAtPosition
  X = 0.427800829876
  Y = 2.28817427386
  CollectingObjectsMode = HomogeneousAdd
  ConsiderStyle = Yes
$!TwoDAxis 
  XDetail
    {
    RangeMin = -527.48293990321712954
    RangeMax = 792.0930035940682501
    }
  YDetail
    {
    RangeMin = -402.38491965166792852
    RangeMax = 357.75308654833213495
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -76.587376445986450335
    RangeMax = 341.197440136837713
    }
  YDetail
    {
    RangeMin = -142.64780008700680014
    RangeMax = 98.015966983671120261
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -144.02562555777848274
    RangeMax = 273.75919102504576585
    }
  YDetail
    {
    RangeMin = -114.1162331550947755
    RangeMax = 126.54753391558317333
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -69.453042234963760393
    RangeMax = 199.18660770223115719
    }
  YDetail
    {
    RangeMin = -71.158906535488682721
    RangeMax = 83.590207295977137392
    }
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.233498777003
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0860772965454
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0326262988274
$!ContourLevels DeleteNearest
  ContourGroup = 1
  RangeMin = 0.0326262988274
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0340132494337
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0222655255119
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0108183537673
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0918610461574
$!ContourLevels DeleteNearest
  ContourGroup = 1
  RangeMin = 0.0918610461574
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0931179094189
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.210007234955
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
1.76569043045
$!ContourLevels DeleteNearest
  ContourGroup = 1
  RangeMin = 1.76569043045
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
1.96986753205
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.61552326142
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.75775543623
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.81587295924
$!ContourLevels DeleteNearest
  ContourGroup = 1
  RangeMin = 2.81587295924
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.81562708382
$!ContourLevels DeleteNearest
  ContourGroup = 1
  RangeMin = 2.81562708382
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.82295874084
$!ContourLevels DeleteNearest
  ContourGroup = 1
  RangeMin = 2.82295874084
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.8255770208
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.82963103458
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.39844031818
$!ContourLevels DeleteNearest
  ContourGroup = 1
  RangeMin = 0.39844031818
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.457681323758
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
1.97119435471
$!GlobalTwoDVector UVar = 8
$!GlobalTwoDVector VVar = 9
$!StreamtraceLayers Show = Yes
$!Streamtrace ResetDeltaTime
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 18.54493251973827
    Y = 6.970341052322041
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 18.21136758180333
    Y = 5.265453591765697
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 17.87780264386839
    Y = 3.560566131209352
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 17.54423770593345
    Y = 1.855678670653006
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 17.21067276799852
    Y = 0.1507912100966617
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 16.87710783006358
    Y = -1.554096250459684
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 16.54354289212864
    Y = -3.258983711016029
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 16.2099779541937
    Y = -4.963871171572372
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 15.87641301625877
    Y = -6.668758632128718
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 15.54284807832383
    Y = -8.373646092685064
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 41.56091323724894
    Y = 9.472078086834074
    }
  AltStartPos
    {
    X = 15.54284807832383
    Y = -8.373646092685064
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 40.94937751770156
    Y = 6.692370270709596
    }
  AltStartPos
    {
    X = 15.54284807832383
    Y = -8.373646092685064
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 40.33784179815417
    Y = 3.912662454585119
    }
  AltStartPos
    {
    X = 15.54284807832383
    Y = -8.373646092685064
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 39.72630607860679
    Y = 1.132954638460642
    }
  AltStartPos
    {
    X = 15.54284807832383
    Y = -8.373646092685064
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 39.11477035905941
    Y = -1.646753177663836
    }
  AltStartPos
    {
    X = 15.54284807832383
    Y = -8.373646092685064
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 38.50323463951202
    Y = -4.426460993788314
    }
  AltStartPos
    {
    X = 15.54284807832383
    Y = -8.373646092685064
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 37.89169891996463
    Y = -7.206168809912789
    }
  AltStartPos
    {
    X = 15.54284807832383
    Y = -8.373646092685064
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 37.28016320041725
    Y = -9.985876626037268
    }
  AltStartPos
    {
    X = 15.54284807832383
    Y = -8.373646092685064
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 36.66862748086987
    Y = -12.76558444216175
    }
  AltStartPos
    {
    X = 15.54284807832383
    Y = -8.373646092685064
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 36.05709176132248
    Y = -15.54529225828622
    }
  AltStartPos
    {
    X = 15.54284807832383
    Y = -8.373646092685064
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -140.7291386219202991
    RangeMax = 270.46270408918780959
    }
  YDetail
    {
    RangeMin = -112.2173011389150048
    RangeMax = 124.64860189940360158
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -184.12753654886341792
    RangeMax = 227.06430616224469077
    }
  YDetail
    {
    RangeMin = -101.49534400402316692
    RangeMax = 135.37055903429543946
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -58.87944102583847
    Y = 65.57572378002337
    }
  AltStartPos
    {
    X = 36.05709176132248
    Y = -15.54529225828622
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -66.28269714278761
    Y = 46.42937175343083
    }
  AltStartPos
    {
    X = -58.87944102583847
    Y = 65.57572378002337
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -69.34611346704241
    Y = 26.77245033946247
    }
  AltStartPos
    {
    X = -66.28269714278761
    Y = 46.42937175343083
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -71.64367571023351
    Y = 9.668375862373125
    }
  AltStartPos
    {
    X = -69.34611346704241
    Y = 26.77245033946247
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -67.814405304915
    Y = -9.222691470531515
    }
  AltStartPos
    {
    X = -71.64367571023351
    Y = 9.668375862373125
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -57.09244817002318
    Y = -27.09262002868456
    }
  AltStartPos
    {
    X = -67.814405304915
    Y = -9.222691470531515
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -26.45828492747512
    Y = -32.19831390244258
    }
  AltStartPos
    {
    X = -57.09244817002318
    Y = -27.09262002868456
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 6.473440558264087
    Y = -16.11537820010483
    }
  AltStartPos
    {
    X = -26.45828492747512
    Y = -32.19831390244258
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 9.536856882518862
    Y = 18.34805544776177
    }
  AltStartPos
    {
    X = 6.473440558264087
    Y = -16.11537820010483
    }
$!Pick SetMouseMode
  MouseMode = Select
$!TwoDAxis 
  XDetail
    {
    RangeMin = -82.059858082336660345
    RangeMax = 124.99662769571818899
    }
  YDetail
    {
    RangeMin = -42.699545255689919543
    RangeMax = 76.57476028596235551
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -72.804321089289359747
    RangeMax = 134.25216468876547538
    }
  YDetail
    {
    RangeMin = -59.28238236823304419
    RangeMax = 59.991923173419223758
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -41.386790732792604786
    RangeMax = 102.83463433226873462
    }
  YDetail
    {
    RangeMin = -41.184402481132622142
    RangeMax = 41.893943286318830133
    }
$!Linking BetweenFrames{LinkAxisPosition = Yes}
$!Linking BetweenFrames{LinkYAxisRange = Yes}
$!Linking BetweenFrames{LinkXAxisRange = Yes}
$!Linking BetweenFrames{LinkSolutionTime = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!Linking BetweenFrames{LinkContourLevels = Yes}
$!FrameName = 'rho 2D - 1'
$!RedrawAll 
$!CreateNewFrame 
  XYPos
    {
    X = -0.32905
    Y = 0.76784
    }
  Width = 11.91
  Height = 6.6523
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A6\HP_srez__2d_(0, 0, 1, 0)_.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"xx" "yy" "X" "Y"'
$!FieldLayers ShowMesh = Yes
$!TwoDAxis XDetail{VarNum = 3}
$!TwoDAxis YDetail{VarNum = 4}
$!FieldMap [1]  Mesh{LineThickness = 0.400000000000000022}
$!FieldMap [1]  Mesh{Color = Red}
$!Linking BetweenFrames{LinkAxisPosition = Yes}
$!Linking BetweenFrames{LinkYAxisRange = Yes}
$!Linking BetweenFrames{LinkXAxisRange = Yes}
$!Linking BetweenFrames{LinkSolutionTime = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!FrameLayout IsTransparent = Yes
$!FrameControl ActivateByNumber
  Frame = 2
$!FrameName = 'HP - 1'
$!RedrawAll 
$!CreateNewFrame 
  XYPos
    {
    X = -1.5838
    Y = 0.29647
    }
  Width = 14.141
  Height = 7.8871
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A6\HP_srez_init_setka_2d_(0, 0, 1, 0)_init.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"X" "Y" "Z"'
$!FieldLayers ShowMesh = Yes
$!Linking BetweenFrames{LinkAxisPosition = Yes}
$!Linking BetweenFrames{LinkXAxisRange = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!Linking BetweenFrames{LinkSolutionTime = Yes}
$!Linking BetweenFrames{LinkYAxisRange = Yes}
$!FieldMap [1]  Mesh{LineThickness = 0.400000000000000022}
$!FrameName = 'HP - init'
$!RedrawAll 
$!FrameLayout IsTransparent = Yes
$!RedrawAll 
$!Pick SetMouseMode
  MouseMode = Select
$!TwoDAxis 
  XDetail
    {
    RangeMin = -98.198049419400646798
    RangeMax = 159.64589301887667716
    }
  YDetail
    {
    RangeMin = -73.910368061634670767
    RangeMax = 74.619908866820821913
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -99.168198032253670249
    RangeMax = 160.61604163172978588
    }
  YDetail
    {
    RangeMin = -74.4692194380376975
    RangeMax = 75.178760243223905491
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -46.838022896768279679
    RangeMax = 108.28586649624460847
    }
  YDetail
    {
    RangeMin = -44.324569418001516397
    RangeMax = 45.034110223187880706
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -29.80855446923221308
    RangeMax = 91.256398068708563187
    }
  YDetail
    {
    RangeMin = -34.514792086925339731
    RangeMax = 35.224332892111732463
    }
$!Page Name = '2D'
$!RedrawAll 
$!Page Name = 'Untitled'
$!PageControl Create
$!Page Name = 'HP'
$!RedrawAll 
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A6\Tecplot_print_all_gran_in_surface_HP.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"X" "Y" "Z"'
$!FieldLayers ShowMesh = Yes
$!PlotType = Cartesian3D
$!Pick AddAtPosition
  X = 1.12489626556
  Y = 4.87074688797
  ConsiderStyle = Yes
$!Pick AddAtPosition
  X = 1.01867219917
  Y = 4.57863070539
  ConsiderStyle = Yes
$!FrameControl ActivateByNumber
  Frame = 1
$!Pick Shift
  X = -2.50954356846
  Y = 0
  PickSubposition = Left
$!Pick Shift
  X = 2.52946058091
  Y = 0
  PickSubposition = Right
$!RedrawAll 
$!ThreeDView 
  PSIAngle = 93.5379
  ThetaAngle = -86.2919
  AlphaAngle = 4.33558E-213
  ViewerPosition
    {
    X = 1131.47188341708
    Y = -61.55990694342515
    Z = -70.88364124913414
    }
  ViewWidth = 203.065
$!ThreeDView 
  PSIAngle = 89.7074
  ThetaAngle = -86.0365
  AlphaAngle = 4.33544E-213
  ViewerPosition
    {
    X = 1133.30529946205
    Y = -66.8267405394313
    Z = 6.130884643434256
    }
  ViewWidth = 203.065
$!ThreeDView 
  ViewerPosition
    {
    X = 1133.30529946205
    Y = -66.82674053943131
    Z = 6.130884643434256
    }
  ViewWidth = 316.857
$!ThreeDView 
  PSIAngle = 87.324
  ThetaAngle = -85.2704
  ViewerPosition
    {
    X = 1130.899585667656
    Y = -82.09165414169473
    Z = 54.0596416212735
    }
  ViewWidth = 316.857
$!PageControl SetCurrentByNumber
  PageNum = 1
$!Pick SetMouseMode
  MouseMode = Select
$!PageControl SetCurrentByNumber
  PageNum = 1
$!Page Name = 'Untitled'
$!PageControl Create
$!Pick AddAtPosition
  X = 1.01867219917
  Y = 4.39937759336
  ConsiderStyle = Yes
$!FrameControl ActivateByNumber
  Frame = 1
$!Pick Shift
  X = -2.52946058091
  Y = 0
  PickSubposition = Left
$!Pick Shift
  X = 2.56265560166
  Y = 0
  PickSubposition = Right
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A6\Tecplot_Tecplot_print_1D_(1, 0, 0)_.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"R" "X" "Y" "Z" "rho" "p" "Vx" "Vy" "Vz" "Bx" "By" "Bz" "Q" "rho_He" "rho_H1" "Vx_H1" "Vy_H1" "Vz_H1" "p_H1" "rho_H2" "Vx_H2" "Vy_H2" "Vz_H2" "p_H2" "rho_H3" "Vx_H3" "Vy_H3" "Vz_H3" "p_H3" "rho_H4" "Vx_H4" "Vy_H4" "Vz_H4" "p_H4" "MK_n_H" "MK_IVx_H" "MK_IVy_H" "MK_IVz_H" "MK_IT_H" "MK_n_H1" "MK_Vx_H1" "MK_Vy_H1" "MK_Vz_H1" "MK_T_H1" "MK_IVx_H1" "MK_IVy_H1" "MK_IVz_H1" "MK_IT_H1" "MK_n_H2" "MK_Vx_H2" "MK_Vy_H2" "MK_Vz_H2" "MK_T_H2" "MK_IVx_H2" "MK_IVy_H2" "MK_IVz_H2" "MK_IT_H2" "MK_n_H3" "MK_Vx_H3" "MK_Vy_H3" "MK_Vz_H3" "MK_T_H3" "MK_IVx_H3" "MK_IVy_H3" "MK_IVz_H3" "MK_IT_H3" "MK_n_H4" "MK_Vx_H4" "MK_Vy_H4" "MK_Vz_H4" "MK_T_H4" "MK_IVx_H4" "MK_IVy_H4" "MK_IVz_H4" "MK_IT_H4" "BB/8pi" "zone_geo" "BB_8Pi" "rho_Th" "T_Th"'
$!Pick DeselectAll
$!Pick AddAllInRect
  SelectText = Yes
  SelectGeoms = Yes
  SelectMaps = Yes
  ConsiderStyle = Yes
  X1 = 3.23609958506
  X2 = 3.26929460581
  Y1 = 2.90560165975
  Y2 = 3.25746887967
$!ActiveLineMaps -= [1]
$!ActiveLineMaps += [4]
$!ActiveLineMaps += [5]
$!ActiveLineMaps += [12]
$!ActiveLineMaps += [77]
$!LineMap [2-77]  Lines{LineThickness = 0.400000000000000022}
$!LineMap [5]  Lines{Color = Blue}
$!LineMap [4]  Lines{Color = Black}
$!XYLineAxis YDetail 1 {CoordScale = Log}
$!RedrawAll 
$!GlobalLinePlot Legend{Show = Yes}
$!GlobalLinePlot Legend{Box{BoxType = None}}
$!Pick AddAtPosition
  X = 2.63195020747
  Y = 3.87489626556
  CollectingObjectsMode = HomogeneousAdd
  ConsiderStyle = Yes
$!Pick SetMouseMode
  MouseMode = Select
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = -82.798360936264458587
    RangeMax = 300.8400112704972571
    }
  YDetail 1
    {
    RangeMin = 0.040007248533265545032
    RangeMax = 153.52316199918956841
    }
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = -143.66255199011118293
    RangeMax = 361.70420232434400987
    }
  YDetail 1
    {
    RangeMin = 0.01080268814472651849
    RangeMax = 568.56582504537880141
    }
$!View Zoom
  X1 = -25.6128938952
  Y1 = 0.0274033241007
  X2 = 86.7635777079
  Y2 = 273.694544313
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = -46.146271850950100202
    RangeMax = 107.29695566369527171
    }
  YDetail 1
    {
    RangeMin = 0.0050935808032429603226
    RangeMax = 1472.4690924009769333
    }
$!View Zoom
  X1 = -2.97783018963
  Y1 = 0.0260586960918
  X2 = 80.5886834465
  Y2 = 291.942057739
$!Pick SetMouseMode
  MouseMode = Select
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = -3.8177111346419159155
    RangeMax = 79.74880250150677341
    }
  YDetail 1
    {
    RangeMin = 0.019866376851472713344
    RangeMax = 222.56796416058790555
    }
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = 13.52436875829671159
    RangeMax = 62.406722608568145461
    }
  YDetail 1
    {
    RangeMin = 0.13754849129306837008
    RangeMax = 32.145892764889509863
    }
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = 2.0711803799750589938
    RangeMax = 50.953534230246496861
    }
  YDetail 1
    {
    RangeMin = 0.048882426271856675715
    RangeMax = 11.42411100442155103
    }
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = 15.684556498203225061
    RangeMax = 37.340158112018329462
    }
  YDetail 1
    {
    RangeMin = 0.22325944338576764192
    RangeMax = 2.5012973938586036127
    }
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = 20.056394923951152975
    RangeMax = 32.9683196862704051
    }
  YDetail 1
    {
    RangeMin = 0.36362521795579055839
    RangeMax = 1.535752297474147765
    }
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = 13.637827442806518263
    RangeMax = 39.386887167415039812
    }
  YDetail 1
    {
    RangeMin = 0.17767749899738857766
    RangeMax = 3.1429880938573484528
    }
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = -0.93356659418444642995
    RangeMax = 53.958281204406020493
    }
  YDetail 1
    {
    RangeMin = 0.034958526312364790467
    RangeMax = 15.974307924348183363
    }
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = -18.407626801715359477
    RangeMax = 71.432341411936931763
    }
  YDetail 1
    {
    RangeMin = 0.0049753157261939020065
    RangeMax = 112.24177411598081733
    }
$!XYLineAxis 
  XDetail 1
    {
    RangeMin = -7.5724430727500697813
    RangeMax = 82.267525140902336034
    }
  YDetail 1
    {
    RangeMin = 0.019408584452943817428
    RangeMax = 437.85240414175001433
    }
$!Page Name = '1D'
$!Page Name = 'Untitled'
$!PageControl Create
$!PageControl SetCurrentByNumber
  PageNum = 3
$!RedrawAll 
$!PageControl SetCurrentByNumber
  PageNum = 2
$!TwoDAxis 
  XDetail
    {
    RangeMin = -11.807643158926119042
    RangeMax = 73.255486758402525993
    }
  YDetail
    {
    RangeMin = -24.14541781890521932
    RangeMax = 24.85495862409166179
    }
