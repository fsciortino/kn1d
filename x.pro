;+
pro x,xsize=xsize,ysize=ysize,title=title,mag=mag,window=window,colors=colors,ypos=ypos,xpos=xpos
common DECW_DISPLAY,DECW_DISPLAY,_xsize,_ysize,_window,_wtitle,_mag,_colors,_xpos,_ypos
if not keyword_set(colors) then colors=16
if not keyword_set(mag) then mag=1.0
if not keyword_set(xsize) then xsize=640
if not keyword_set(ysize) then ysize=480
if not keyword_set(title) then title='X'
if not keyword_set(window) then window=0
if type_of(xpos) eq 0 then xpos=550+xsize-640
if type_of(ypos) eq 0 then ypos=420+ysize-480
;-
_colors=colors
_mag=mag
_xsize=xsize*mag
_ysize=ysize*mag
_wtitle=title
_window=window
_xpos=xpos
_ypos=ypos

if !d.name eq 'PS' then device,/close
if strlen(getenv('DISPLAY')) gt 0 then begin
   set_plot,'X'
   device,retain=2,decompose=0,true=24
   window,_window,colors=colors,xsize=_xsize,ysize=_ysize,title=_wtitle,retain=2,xpos=_xpos,ypos=_ypos
;   wdelete
   !p.font=-1
   !x.thick=0.0
   !y.thick=0.0
   !p.thick=0.0
endif
@setup_colors
end
