#!MC 1410
$!Page Name = 'HP'
$!Pick AddAtPosition
  X = 1.03858921162
  Y = 4.15373443983
  ConsiderStyle = Yes
$!FrameControl ActivateByNumber
  Frame = 1
$!Pick Shift
  X = -2.54273858921
  Y = 0
  PickSubposition = Left
$!Pick Shift
  X = 2.53609958506
  Y = 0
  PickSubposition = Right
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A4\Tecplot_print_all_gran_in_surface_HP.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"X" "Y" "Z"'
$!PlotType = Cartesian3D
$!FieldLayers ShowMesh = Yes
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!Linking BetweenFrames{Link3DView = Yes}
$!RedrawAll 
$!FrameName = 'HP-1'
$!CreateNewFrame 
  XYPos
    {
    X = -1.0925
    Y = 0.39606
    }
  Width = 13.112
  Height = 7.449
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A4\Tecplot_print_all_gran_in_surface_HP-.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"X" "Y" "Z"'
$!PlotType = Cartesian3D
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!Linking BetweenFrames{Link3DView = Yes}
$!FieldLayers ShowMesh = Yes
$!FrameName = 'HP-init'
$!Page Name = 'Untitled'
$!PageControl Create
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A4\Tecplot_print_all_gran_in_surface_TS-.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"X" "Y" "Z"'
$!PlotType = Cartesian3D
$!FieldLayers ShowMesh = Yes
$!Linking BetweenFrames{LinkGroup = 2}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!Linking BetweenFrames{Link3DView = Yes}
$!Page Name = 'TS'
$!FrameName = 'TS-init'
$!Pick SetMouseMode
  MouseMode = Select
$!Pick AddAtPosition
  X = 0.925726141079
  Y = 3.90145228216
  ConsiderStyle = Yes
$!FrameControl ActivateByNumber
  Frame = 1
$!Pick Shift
  X = -2.56265560166
  Y = 0
  PickSubposition = Left
$!Pick Shift
  X = 2.60248962656
  Y = 0
  PickSubposition = Right
$!RedrawAll 
$!CreateNewFrame 
  XYPos
    {
    X = -0.52822
    Y = 0.54876
    }
  Width = 12.501
  Height = 6.8382
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A4\Tecplot_print_all_gran_in_surface_TS.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"X" "Y" "Z"'
$!PlotType = Cartesian3D
$!FieldLayers ShowMesh = Yes
$!Linking BetweenFrames{LinkGroup = 2}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!Linking BetweenFrames{Link3DView = Yes}
$!FrameName = 'TS-1'
$!Page Name = 'Untitled'
$!PageControl Create
$!Page Name = '2D XY'
$!Pick SetMouseMode
  MouseMode = Select
$!Pick AddAtPosition
  X = 1.01867219917
  Y = 4.33298755187
  ConsiderStyle = Yes
$!FrameControl ActivateByNumber
  Frame = 1
$!Pick Shift
  X = -2.58921161826
  Y = 0
  PickSubposition = Left
$!Pick Shift
  X = 2.59585062241
  Y = 0
  PickSubposition = Right
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A4\Tecplot_Tecplot_print_2D__2d_(0, 0, 1, 0)_.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"xx" "yy" "X" "Y" "Z" "rho" "p" "Vx" "Vy" "Vz" "Bx" "By" "Bz" "Q" "rho_He" "rho_H1" "rho_H2" "rho_H3" "rho_H4" "rho_H5" "rho_H6" "rho_H7" "rho_H8" "rho_H9" "S_k1__p_H1" "S_k1__p_H2" "S_k1__p_H3" "S_k1__p_H4" "S_k1__Pui_1_H3" "S_k1__Pui_1_H4" "S_k1__Pui_2_H3" "S_k1__Pui_2_H4" "S_k2__p_H1" "S_k2__p_H2" "S_k2__p_H3" "S_k2__p_H4" "S_k2__Pui_1_H3" "S_k2__Pui_1_H4" "S_k2__Pui_2_H3" "S_k2__Pui_2_H4" "Mach" "Mach_Alf" "|B|/8pi" "rho_Th" "p_Th" "T_Th" "J" "Jxx" "Jyy" "Jx" "Jy" "Jz" "J_an" "Jxx_an" "Jyy_an" "Bxx_pot" "Byy_pot"'
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
$!RedrawAll 
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
7
0
0.333333333333
0.666666666667
1
1.33333333333
1.66666666667
2
$!RedrawAll 
$!Pick AddAtPosition
  X = 0.520746887967
  Y = 3.58278008299
  CollectingObjectsMode = HomogeneousAdd
  ConsiderStyle = Yes
$!TwoDAxis 
  XDetail
    {
    RangeMin = -294.85423882459451761
    RangeMax = 1011.1161004319594667
    }
  YDetail
    {
    RangeMin = -372.12428619087131665
    RangeMax = 384.04932360912846434
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -23.278244863827694644
    RangeMax = 739.54010647119287114
    }
  YDetail
    {
    RangeMin = -214.87829559261629697
    RangeMax = 226.80333301087367204
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -321.74715808728859656
    RangeMax = 441.07119324773196922
    }
  YDetail
    {
    RangeMin = -227.73102391324377436
    RangeMax = 213.95060469024619465
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -145.24126225898905318
    RangeMax = 264.56529741943262479
    }
  YDetail
    {
    RangeMin = -125.53184190467868575
    RangeMax = 111.75142268168127657
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -168.76883712771254409
    RangeMax = 241.03772255070910546
    }
  YDetail
    {
    RangeMin = -123.23023131969486599
    RangeMax = 114.05303326666509633
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -100.65694355235152102
    RangeMax = 172.92582897534825292
    }
  YDetail
    {
    RangeMin = -83.792569615590636545
    RangeMax = 74.615371562561037422
    }
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0951876076992
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0338650343189
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0175622109172
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0111646260848
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0405412126604
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.058202111814
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0424863973573
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0617653227446
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.052060341697
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0368328055221
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0314267805835
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.97125737848
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.81937261198
$!ContourLevels DeleteNearest
  ContourGroup = 1
  RangeMin = 2.81937261198
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.81937261198
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.51438642494
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.45650769019
$!TwoDAxis 
  XDetail
    {
    RangeMin = -231.01891295777087976
    RangeMax = 303.28779838076769693
    }
  YDetail
    {
    RangeMin = -159.27382521749757416
    RangeMax = 150.09662716446803188
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -245.02285648317106848
    RangeMax = 289.28385485536750821
    }
  YDetail
    {
    RangeMin = -150.93814454761650268
    RangeMax = 158.43230783434910336
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -167.23408527358969877
    RangeMax = 211.49508364578616693
    }
  YDetail
    {
    RangeMin = -105.89744633318325384
    RangeMax = 113.39160961991591137
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
    X = 120.8777838549822
    Y = 53.03010195119624
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 121.0616046823037
    Y = 42.70987550300742
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 121.2454255096251
    Y = 32.38964905481861
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 121.4292463369465
    Y = 22.06942260662979
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 121.613067164268
    Y = 11.74919615844097
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 121.7968879915894
    Y = 1.428969710252147
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 121.9807088189108
    Y = -8.891256737936672
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 122.1645296462323
    Y = -19.21148318612548
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 122.3483504735537
    Y = -29.53170963431431
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 122.5321713008751
    Y = -39.85193608250313
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -76.821939787849657932
    RangeMax = 121.08293816004621135
    }
  YDetail
    {
    RangeMin = -53.54765218653054859
    RangeMax = 61.041815473263248748
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 18.64988347155449
    Y = 0.2304193686316296
    }
  AltStartPos
    {
    X = 122.5321713008751
    Y = -39.85193608250313
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 17.90888351902041
    Y = -3.845080370305908
    }
  AltStartPos
    {
    X = 18.64988347155449
    Y = 0.2304193686316296
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 6.42338425474189
    Y = -8.908580045955574
    }
  AltStartPos
    {
    X = 17.90888351902041
    Y = -3.845080370305908
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -12.22511455069957
    Y = -8.291080085510501
    }
  AltStartPos
    {
    X = 6.42338425474189
    Y = -8.908580045955574
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -21.36411396528679
    Y = -0.5105805839024669
    }
  AltStartPos
    {
    X = -12.22511455069957
    Y = -8.291080085510501
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -21.11711398110876
    Y = 9.245918791129824
    }
  AltStartPos
    {
    X = -21.36411396528679
    Y = -0.5105805839024669
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -5.185615001625649
    Y = 17.02641829273784
    }
  AltStartPos
    {
    X = -21.11711398110876
    Y = 9.245918791129824
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 6.793884231008931
    Y = 11.71591863291015
    }
  AltStartPos
    {
    X = -5.185615001625649
    Y = 17.02641829273784
    }
$!Pick SetMouseMode
  MouseMode = Select
$!TwoDAxis 
  XDetail
    {
    RangeMin = -140.65132875056687567
    RangeMax = 184.9123271227635712
    }
  YDetail
    {
    RangeMin = -90.505688684313795989
    RangeMax = 97.999851971046652466
    }
$!FrameName = 'rho'
$!CreateNewFrame 
  XYPos
    {
    X = -1.5639
    Y = 0.28983
    }
  Width = 14.134
  Height = 7.9137
$!FrameControl ActivateByNumber
  Frame = 1
$!Linking BetweenFrames{LinkGroup = 2}
$!Linking BetweenFrames{LinkGroup = 3}
$!Linking BetweenFrames{LinkAxisPosition = Yes}
$!Linking BetweenFrames{LinkYAxisRange = Yes}
$!Linking BetweenFrames{LinkXAxisRange = Yes}
$!Linking BetweenFrames{LinkSolutionTime = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!FrameControl ActivateByNumber
  Frame = 2
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A4\HP_srez__2d_(0, 0, 1, 0)_.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"xx" "yy" "X" "Y"'
$!FieldLayers ShowMesh = Yes
$!Linking BetweenFrames{LinkGroup = 2}
$!Linking BetweenFrames{LinkGroup = 3}
$!Linking BetweenFrames{LinkAxisPosition = Yes}
$!Linking BetweenFrames{LinkYAxisRange = Yes}
$!Linking BetweenFrames{LinkXAxisRange = Yes}
$!Linking BetweenFrames{LinkSolutionTime = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!TwoDAxis XDetail{VarNum = 3}
$!TwoDAxis YDetail{VarNum = 4}
$!FrameLayout IsTransparent = Yes
$!FieldMap [1]  Mesh{LineThickness = 0.400000000000000022}
$!RedrawAll 
$!FrameName = 'HP'
$!Pick SetMouseMode
  MouseMode = Select
$!TwoDAxis 
  XDetail
    {
    RangeMin = -81.863648296528978676
    RangeMax = 176.42749089285803166
    }
  YDetail
    {
    RangeMin = -65.29423333922075301
    RangeMax = 84.25965561922087943
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -132.31400907035254022
    RangeMax = 125.97713011903442748
    }
  YDetail
    {
    RangeMin = -55.139687879952745675
    RangeMax = 94.414201078488829921
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -174.55655771531277765
    RangeMax = 168.21967876399477859
    }
  YDetail
    {
    RangeMin = -79.598666415556863285
    RangeMax = 118.8731796140930328
    }
$!Page Name = 'Untitled'
$!PageControl Create
$!Pick AddAtPosition
  X = 0.965560165975
  Y = 4.5387966805
  ConsiderStyle = Yes
$!FrameControl ActivateByNumber
  Frame = 1
$!Pick Shift
  X = -2.54273858921
  Y = 0
  PickSubposition = Left
$!Pick Shift
  X = 2.48962655602
  Y = 0
  PickSubposition = Right
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A4\Tecplot_Tecplot_print_2D__2d_(0, 1, 0, 0)_.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"xx" "yy" "X" "Y" "Z" "rho" "p" "Vx" "Vy" "Vz" "Bx" "By" "Bz" "Q" "rho_He" "rho_H1" "rho_H2" "rho_H3" "rho_H4" "rho_H5" "rho_H6" "rho_H7" "rho_H8" "rho_H9" "S_k1__p_H1" "S_k1__p_H2" "S_k1__p_H3" "S_k1__p_H4" "S_k1__Pui_1_H3" "S_k1__Pui_1_H4" "S_k1__Pui_2_H3" "S_k1__Pui_2_H4" "S_k2__p_H1" "S_k2__p_H2" "S_k2__p_H3" "S_k2__p_H4" "S_k2__Pui_1_H3" "S_k2__Pui_1_H4" "S_k2__Pui_2_H3" "S_k2__Pui_2_H4" "Mach" "Mach_Alf" "|B|/8pi" "rho_Th" "p_Th" "T_Th" "J" "Jxx" "Jyy" "Jx" "Jy" "Jz" "J_an" "Jxx_an" "Jyy_an" "Bxx_pot" "Byy_pot"'
$!TwoDAxis XDetail{VarNum = 3}
$!TwoDAxis YDetail{VarNum = 5}
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
0.5
1
1.5
2
$!RedrawAll 
$!GlobalTwoDVector UVar = 8
$!GlobalTwoDVector VVar = 10
$!Streamtrace ResetDeltaTime
$!StreamtraceLayers Show = Yes
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 151.3617510668977
    Y = 42.79486503521048
    }
  AltStartPos
    {
    X = 6.793884231008931
    Y = 11.71591863291015
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 151.0897345018638
    Y = 32.54890775226843
    }
  AltStartPos
    {
    X = 6.793884231008931
    Y = 11.71591863291015
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 150.81771793683
    Y = 22.30295046932639
    }
  AltStartPos
    {
    X = 6.793884231008931
    Y = 11.71591863291015
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 150.5457013717961
    Y = 12.05699318638434
    }
  AltStartPos
    {
    X = 6.793884231008931
    Y = 11.71591863291015
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 150.2736848067622
    Y = 1.811035903442296
    }
  AltStartPos
    {
    X = 6.793884231008931
    Y = 11.71591863291015
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 150.0016682417284
    Y = -8.434921379499762
    }
  AltStartPos
    {
    X = 6.793884231008931
    Y = 11.71591863291015
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 149.7296516766945
    Y = -18.6808786624418
    }
  AltStartPos
    {
    X = 6.793884231008931
    Y = 11.71591863291015
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 149.4576351116607
    Y = -28.92683594538384
    }
  AltStartPos
    {
    X = 6.793884231008931
    Y = 11.71591863291015
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 149.1856185466268
    Y = -39.17279322832589
    }
  AltStartPos
    {
    X = 6.793884231008931
    Y = 11.71591863291015
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 148.9136019815929
    Y = -49.41875051126794
    }
  AltStartPos
    {
    X = 6.793884231008931
    Y = 11.71591863291015
    }
$!Pick SetMouseMode
  MouseMode = Select
$!TwoDAxis 
  XDetail
    {
    RangeMin = -454.18078474839683167
    RangeMax = 839.4364928584094514
    }
  YDetail
    {
    RangeMin = -399.05524237753951411
    RangeMax = 358.11646722246069885
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -65.96684062253200409
    RangeMax = 451.22254873254485119
    }
  YDetail
    {
    RangeMin = -171.82836403208477805
    RangeMax = 130.889588877006247
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -206.9100487453183348
    RangeMax = 310.27934060975849206
    }
  YDetail
    {
    RangeMin = -152.90543331189584819
    RangeMax = 149.81251959719517686
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -67.75061933663150171
    RangeMax = 171.11991120107171582
    }
  YDetail
    {
    RangeMin = -71.453534905540209365
    RangeMax = 68.360621190839509609
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -77.695906077334285555
    RangeMax = 161.17462446036890356
    }
  YDetail
    {
    RangeMin = -69.193242464471410358
    RangeMax = 70.620913631908308616
    }
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.159600229154
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0668072533332
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0289220991093
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0156652476386
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.00931727460604
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0396415267558
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0524936019696
$!ContourLevels DeleteNearest
  ContourGroup = 1
  RangeMin = 0.0524936019696
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0521628442372
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0331901147555
$!ContourLevels DeleteNearest
  ContourGroup = 1
  RangeMin = 0.0331901147555
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0331850148794
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.036283106009
$!TwoDAxis 
  XDetail
    {
    RangeMin = -120.18940396942792859
    RangeMax = 118.68112656827527474
    }
  YDetail
    {
    RangeMin = -68.590497813519718306
    RangeMax = 71.223658282860000668
    }
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
0.0382201708093
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.97583085106
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.66708449413
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.36946408334
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.34066691256
$!TwoDAxis 
  XDetail
    {
    RangeMin = -220.39472572685045293
    RangeMax = 218.8864483256979554
    }
  YDetail
    {
    RangeMin = -127.24202952964691349
    RangeMax = 129.87518999898725269
    }
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
2.15241857721
$!ContourLevels Add
  ContourGroup = 1
  RawData
1
1.89632610027
$!TwoDAxis 
  XDetail
    {
    RangeMin = -85.285249325507678009
    RangeMax = 83.77697192435523732
    }
  YDetail
    {
    RangeMin = -48.160623459357942977
    RangeMax = 50.793783928698324814
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -73.873791018106686579
    RangeMax = 95.18843023175622875
    }
  YDetail
    {
    RangeMin = -41.974973161888243567
    RangeMax = 56.979434226168024225
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 20.94169643922669
    Y = -1.293657541718972
    }
  AltStartPos
    {
    X = 148.9136019815929
    Y = -49.41875051126794
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = 19.34195929332935
    Y = 3.505553895973037
    }
  AltStartPos
    {
    X = 20.94169643922669
    Y = -1.293657541718972
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -16.49215277477099
    Y = 24.94203165099736
    }
  AltStartPos
    {
    X = 19.34195929332935
    Y = 3.505553895973037
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -146.2313801283195005
    RangeMax = 167.54601934196907109
    }
  YDetail
    {
    RangeMin = -84.326849870330519821
    RangeMax = 99.331310934610300478
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -34.58533971971914
    Y = -35.94054022031402
    }
  AltStartPos
    {
    X = -16.49215277477099
    Y = 24.94203165099736
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -51.60813632393848
    Y = -26.63738393661275
    }
  AltStartPos
    {
    X = -34.58533971971914
    Y = -35.94054022031402
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -55.76486572729436
    Y = -1.103189030283758
    }
  AltStartPos
    {
    X = -51.60813632393848
    Y = -26.63738393661275
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -57.54632118587546
    Y = 12.75257564756919
    }
  AltStartPos
    {
    X = -55.76486572729436
    Y = -1.103189030283758
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -59.92159513065025
    Y = 26.21246133462634
    }
  AltStartPos
    {
    X = -57.54632118587546
    Y = 12.75257564756919
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -60.51541361684394
    Y = -12.97955875415771
    }
  AltStartPos
    {
    X = -59.92159513065025
    Y = 26.21246133462634
    }
$!Streamtrace Add
  StreamType = TwoDLine
  StreamDirection = Both
  StartPos
    {
    X = -82.09081861521497
    Y = 40.46410500327508
    }
  AltStartPos
    {
    X = -60.51541361684394
    Y = -12.97955875415771
    }
$!Pick SetMouseMode
  MouseMode = Select
$!TwoDAxis 
  XDetail
    {
    RangeMin = -81.788100041512052485
    RangeMax = 103.1027392551616515
    }
  YDetail
    {
    RangeMin = -46.607325401834110323
    RangeMax = 61.611786466113905192
    }
$!TwoDAxis 
  XDetail
    {
    RangeMin = -84.587322636576743662
    RangeMax = 100.30351666009696032
    }
  YDetail
    {
    RangeMin = -42.058588684854015582
    RangeMax = 66.160523183093999933
    }
$!FrameName = 'rho'
$!Linking BetweenFrames{LinkGroup = 2}
$!Linking BetweenFrames{LinkGroup = 3}
$!Linking BetweenFrames{LinkGroup = 4}
$!Linking BetweenFrames{LinkAxisPosition = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = No}
$!Linking BetweenFrames{LinkSolutionTime = Yes}
$!Linking BetweenFrames{LinkXAxisRange = Yes}
$!Linking BetweenFrames{LinkYAxisRange = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!Page Name = '2D XZ'
$!CreateNewFrame 
  XYPos
    {
    X = -1.5241
    Y = 0.30975
    }
  Width = 13.982
  Height = 7.9071
$!ReadDataSet  '"D:\Science_code\TITAN_C-\TITAN_A4\HP_srez__2d_(0, 1, 0, 0)_.txt" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"xx" "yy" "X" "Y"'
$!FieldLayers ShowMesh = Yes
$!TwoDAxis XDetail{VarNum = 3}
$!TwoDAxis YDetail{VarNum = 4}
$!Linking BetweenFrames{LinkGroup = 2}
$!Linking BetweenFrames{LinkGroup = 3}
$!Linking BetweenFrames{LinkGroup = 4}
$!Linking BetweenFrames{LinkAxisPosition = Yes}
$!Linking BetweenFrames{LinkYAxisRange = Yes}
$!Linking BetweenFrames{LinkXAxisRange = Yes}
$!Linking BetweenFrames{LinkSolutionTime = Yes}
$!Linking BetweenFrames{LinkFrameSizeAndPosition = Yes}
$!FrameLayout IsTransparent = Yes
$!RedrawAll 
$!Pick SetMouseMode
  MouseMode = Select
$!TwoDAxis 
  XDetail
    {
    RangeMin = -98.032331387713199433
    RangeMax = 113.74852541123341609
    }
  YDetail
    {
    RangeMin = -49.928134367266586935
    RangeMax = 74.030068865506592601
    }
$!FrameControl ActivateByNumber
  Frame = 2
$!FrameName = 'HP'
$!FrameControl MoveToBottomByNumber
  Frame = 2
$!FrameControl ActivateAtPosition
  X = 4.11244813278
  Y = 5.67406639004
$!TwoDAxis 
  XDetail
    {
    RangeMin = -121.74877573721421697
    RangeMax = 137.46496976073467522
    }
  YDetail
    {
    RangeMin = -63.809690392064119635
    RangeMax = 87.91162489030426741
    }
