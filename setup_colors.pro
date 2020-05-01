;
; 16 Color Definitions mimicking Versterm-Pro Tek 4105/407 default color Palette
;
; Store new variables in Common Block to save local variable space
;
common EDGEDB_Colors,Red_Table,Green_Table,Blue_Table,Color_Name,$
              white,black,red,green,blue,cyan,magenta,yellow,$
              orange,Lime,TurquoiseGreen,TurquoiseBlue,Purple,Pink,darkgray,lightgray

common DECW_DISPLAY,DECW_DISPLAY,_xsize,_ysize,_window,_wtitle,_mag,_colors,_xpos,_ypos
if type_of(_colors) eq 0 then _colors=16
if type_of(_mag) eq 0 then _mag=1.0
if type_of(_ysize) eq 0 then _ysize=480
if type_of(_xsize) eq 0 then _xsize=640 
if type_of(_wtitle) eq 0 then _wtitle='X'
if type_of(_xpos) eq 0 then _xpos=550+_xsize-640
if type_of(_ypos) eq 0 then _ypos=420+_ysize-480
if type_of(_window) eq 0 then _window=0
;
; If DISPLAY is not defined and device is X reset it to TEK
; If DISPLAY is defined and device is TEK reset it to X
;
reset=1
if !d.name eq "X" or !d.name eq "TEK" then begin & $
   reset=0 & $
   if strlen(getenv('DISPLAY')) lt 1 then begin & $
      if !d.name eq "X" or type_of(red) eq 0 then begin & $
         set_plot,'tek' & $
         _colors=16 & $
         device,/tek4100,colors=_colors & $
         reset=1 & $
      endif & $
   endif else begin & $
      if !d.name eq "TEK" or type_of(red) eq 0 then begin & $
         set_plot,'X' & $ 
         device,retain=2,decompose=0,true=24 & $
         window,_window,colors=_colors,retain=2,xsize=_xsize,ysize=_ysize,title=_wtitle,xpos=_xpos,ypos=_ypos & $
;         wdelete & $
         reset=1 & $
      endif & $
   endelse & $
endif
;
; Versterm-Pro Tektronix should be setup NOT to ignore host color mapping (see color palette menu)
;
; All devices should paint WHITE for color index=0 and BLACK for color index=1
;
; color Index: 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
red_table=  [  0,255,255,  0,  0,  0,255,255,255,128,  0,  0,128,255, 85,170]
green_table=[  0,255,  0,255,  0,255,  0,255,128,255,255,128,  0,  0, 85,170]
blue_table= [  0,255,  0,  0,255,255,255,  0,  0,  0,128,255,255,128, 85,170]
;
; Modify the color table for the X terminal graphics: In order to get
; a white background, define the "black" color (no color) to paint white and the "white" color
; (first nonzero color index) to paint black
;
if !D.name eq 'X' or !D.name eq 'PS' then begin & Red_table(1)=0 & Red_Table(0)=255 & endif
if !D.name eq 'X' or !D.name eq 'PS' then begin & Green_table(1)=0 & Green_Table(0)=255 & endif
if !D.name eq 'X' or !D.name eq 'PS' then begin & Blue_table(1)=0 & Blue_Table(0)=255 & endif
;
; Load Color Tables
;
if reset eq 1 and _colors eq 16 then tvlct,red_table,green_table,blue_table
;
; To get the following definitions do:
; IDL>@Setup_Colors
;
; color Index: 0       1      2      3       4     5        6        7        8        9
Color_Name=['White','Black','Red','Green','Blue','Cyan','Magenta','Yellow','Orange','Lime']
; color Index:                10              11           12      13      14          15
Color_Name=[Color_Name,'TurquoiseGreen','TurquoiseBlue','Purple','Pink','DarkGray','LightGray']
white=0
black=1
red=2
green=3
blue=4
cyan=5
magenta=6
yellow=7
orange=8
Lime=9
TurquoiseGreen=10
TurquoiseBlue=11
Purple=12
Pink=13
darkgray=14
lightgray=15
!p.color=black
