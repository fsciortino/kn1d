;+
;_______________________________________________________________________________
; DBOPlot.pro
;_______________________________________________________________________________
; This procedure calls OPLOT and XYOUTS to display database data. Sets of
; data points can be displayed with different plot symbols and with a 
; corresponding legend
;_______________________________________________________________________________
   pro dboplot,xdata_in,ydata_in,sigma_in,psym=psym,thick=thick,$
         symthick=symthick,trajectory=trajectory,shots=shots,$
         color=color,fillcolor=fillcolor,drawcolor=drawcolor,$
         xcutoff=xcutoff,ycutoff=ycutoff,set=set,symsize=symsize,$
         inquire=inquire,indices=nindices,xlegend=xlegend,ylegend=ylegend,$
         legendsize=legendsize,xnorm=xnorm,ynorm=ynorm,$
         nolegend=nolegend,legendcolor=legendcolor,nocopy=nocopy,iactive=iactive,$
         lastyleg=lastyleg,nodrawcolor=nodrawcolor
;_______________________________________________________________________________
; Input:
;	xdata_in	- fltarr(*), X value of all data points which may be plotted
;	ydata_in	- fltarr(*), Y value of all data points which may be plotted
;	Sigma_in	- optional: fltarr(*), error bar (1 sigma) associated with Y
;
; Keywords:
;	xnorm		- float, normalize x variable by this factor
;	ynorm		- float, normalize y variable by this factor
;	xcutoff		- float-> largest x (not normalized) allowed on plot
;			  -or- fltarr(2)-> [smallest,largest] x allowed on plot
;	ycutoff		- float-> largest y (not normalized) allowed on plot
;			  -or- fltarr(2)-> [smallest,largest] y allowed on plot
;	indices		- lonarr, indices of x,y that lie within 'selection
;			  box'
;	symsize		- same as PLOT keyword
;	Legendsize	- character size of legend text
;	psym		- same as PLOT keyword (ignored if SET is defined, see below)
;	color		- same as PLOT keyword (ignored if SET is defined, see below)
;	FillColor	- color to fill solid symbols with (ignored if SET is defined, see below)
;	DrawColor	- color to draw outline of symbols with (ignored if SET is defined, see below)
;	NoDrawColor	- if set, override DRAWCOLOR -> do not draw any outlines
;	SymThick	- thickness of lines making up a symbol
;	shots		- shot numbers (input) corresponding to xdata_in, ydata_in
;	trajectory 	- if set, then lines are plotted connecting data points that have the
;			  same shot number
;
;	set		- struct array, of the form:
;  [{Subset:lonarr(n_elements(xdata)), PSym:integer, Color:integer, FillColor:integer,
;	    DrawColor:integer, SymSize:float, legend:string},
;   { },{ }, ...]
;		when SET exists, then data is group into n_elements(sets)
;		data sets, and each is plotted with specified PSYM, COLOR, FILLCOLOR,
;		DRAWCOLOR, SYMSIZE, SYMTHICK, and LEGEND. Data to be included in each set is indicated by the
;		array Set(i).SUBSET. Elements of this array equal to 1
;		indicate that the corresponding data point belongs to
;		this set. (0=does not belong)
;
;	All tag names except SUBSET are optional
;
;	xlegend		- sets x coordinate of legend (in NORMAL plot units)
;	ylegend		- sets y coordinate of legend (in NORMAL plot units)
;
;
;  History:
;	Written by B. LaBombard 5/7/94
;	04/21/95 added FillColor,DrawColor keywords and allowed code to look
;	         for PSYM, COLOR, FILLCOLOR, DRAWCOLOR, SYMSIZE, and LEGEND
;		 tag names of SET - B. LaBombard
;	3/25/96 added SymThick tag
;	3/27/98 added shots and trajectory keywords
;_______________________________________________________________________________
;-
   Common DBPLOT,nleg,xleg,yleg,legsize,Yn,Xn,Xcut,Ycut,LegendDir,Thk
   if type_of(color) eq 0 then color=!p.color
   if keyword_set(thick) then thk=thick else thk=!p.thick*0.8
   if keyword_set(symthick) then symthick=symthick else symthick=!p.thick*0.8
   if keyword_set(legendcolor) then legendcolor=legendcolor else legendcolor=0
   if keyword_set(nolegend) then nolegend=1 else nolegend=0
   if keyword_set(inquire) then inquire=1 else inquire=0
   if keyword_set(ycutoff) then ycutoff=ycutoff else ycutoff=ycut
   if keyword_set(xcutoff) then xcutoff=xcutoff else xcutoff=xcut
   if keyword_set(xnorm) then xn=xnorm 
   if keyword_set(ynorm) then yn=ynorm
   if keyword_set(psym) then psym=psym else psym=!p.psym
   if keyword_set(legendsize) then begin
      legsize=legendsize
      nleg=0
   endif
   if keyword_set(symsize) then symsize=symsize else symsize=1
   if keyword_set(trajectory) then trajectory=trajectory else trajectory=0
   if n_elements(shots) lt 2 then trajectory=0
   if keyword_set(title) then title=title else title=''
   if keyword_set(ytitle) then ytitle=ytitle else ytitle=''
   if keyword_set(xtitle) then xtitle=xtitle else xtitle=''
   if keyword_set(xlegend) then begin
      xleg=xlegend
      nleg=0
   endif
   if keyword_set(ylegend) then begin
      yleg=ylegend
      nleg=0
   endif
   if type_of(Fillcolor) ne 0 then use_fill_color=1 else use_fill_color=0
   if type_of(noDrawcolor) eq 0 then nodrawcolor=0
   if type_of(Drawcolor) ne 0 then use_draw_color=1 else use_draw_color=0
   if nodrawcolor then use_draw_color=0
   if n_params() gt 2 then dosigma=1 else dosigma=0
   if n_tags(set) gt 0 then doset=1 else doset=0
   if n_elements(xcutoff) eq 1 then xcut=[-1e32,xcutoff]
   if n_elements(ycutoff) eq 1 then ycut=[-1e32,ycutoff]

   iactive=[-1]
   xdata=float(xdata_in) & xdata(*)=1e32
   ok=where(xdata_in lt 1e32,count)
   if count gt 0 then xdata(ok)=xdata_in(ok)/xn
   bad=where(xdata_in ge 1e32,count)
   if count gt 0 then xdata(bad)=1e32
   ydata=float(ydata_in) & ydata(*)=1e32
   if dosigma then begin
      Sigma=float(Sigma_in)
      Sigma(*)=1e32
   endif
   ok=where(ydata_in lt 1e32,count)
   if count gt 0 then begin 
      ydata(ok)=ydata_in(ok)/yn
      if dosigma then Sigma(ok)=Sigma_in(ok)/yn
   endif
   bad=where(ydata_in ge 1e32,count)
   if count gt 0 then ydata(bad)=1e32
   xc=xcut
   if abs(xn) gt 1.0 then begin
      xc=xcut/xn
   endif else begin
      if 1e32*abs(xn) gt -xcut(0) then xc(0)=xcut(0)/xn else xc(0)=-1e32
      if 1e32*abs(xn) gt  xcut(1) then xc(1)=xcut(1)/xn else xc(1)=1e32
   endelse
   yc=ycut
   if abs(yn) gt 1.0 then begin
      yc=ycut/yn
   endif else begin
      if 1e32*abs(yn) gt -ycut(0) then yc(0)=ycut(0)/yn else yc(0)=-1e32
      if 1e32*abs(yn) gt  ycut(1) then yc(1)=ycut(1)/yn else yc(1)=1e32
   endelse
   inrange=xdata gt xc(0) and xdata lt xc(1) and ydata gt yc(0) and ydata lt yc(1)
   if doset eq 0 then begin
      anyactive=inrange
      indices=where(inrange,count)
      if count lt 1 then begin
         print,'nothing to plot'
         goto,exit
      endif
      symbol=psymbol(psym,color=color,thick=symthick)
      if symbol le 0 then oplot,xdata(indices),ydata(indices),thick=thk,psym=0,color=color
      if use_fill_color or use_draw_color then begin
         for ii=0,n_elements(indices)-1 do begin
            if use_fill_Color then fcolor=fillcolor else fcolor=color
            plots,xdata(indices(ii)),ydata(indices(ii)),psym=abs(psymbol(psym,color=fcolor,thick=symthick)),symsize=symsize,noclip=0,color=fcolor
            if use_draw_color then begin
               plots,xdata(indices(ii)),ydata(indices(ii)),psym=abs(psymbol(psym,color=drawcolor,thick=symthick,/nofill)),symsize=symsize,$
			noclip=0,color=drawcolor
            endif
         endfor
      endif else begin
         plots,xdata(indices),ydata(indices),psym=abs(symbol),symsize=symsize,noclip=0,color=color
      endelse
      if dosigma then begin
         color_save=!p.color
         !p.color=color
         width=.01*symsize
         errplot,xdata(indices),ydata(indices)-sigma(indices),ydata(indices)+sigma(indices),width=width
         !p.color=color_save
      endif
      if trajectory then begin
         Ushots=shots(unique_shots(shots))
         for ii=0,n_elements(Ushots)-1 do begin
            kk=where(Ushots(ii) eq shots(indices),count)
            if count gt 1 then plots,xdata(indices(kk)),ydata(indices(kk)),noclip=0,color=1
         endfor
      endif
   endif else begin
      SetTags=strupcase(tag_names(SET))
      Taglist=SetTags(0)
      for ii=1,n_tags(set)-1 do TagList=TagList+','+SetTags(ii)
      if strpos(Taglist,'SUBSET') ne -1 then do_subset=1 else do_subset=0
      if do_subset eq 0 then begin
         message,'SET must contain a SUBSET tag',/cont
         stop
      endif      
      if strpos(Taglist,'PSYM') ne -1 then do_psym=1 else do_psym=0
      if strpos(Taglist,'COLOR') ne -1 then do_COLOR=1 else do_COLOR=0
      if strpos(Taglist,'FILLCOLOR') ne -1 then do_FILLCOLOR=1 else do_FILLCOLOR=0
      if strpos(Taglist,'DRAWCOLOR') ne -1 then do_DRAWCOLOR=1 else do_DRAWCOLOR=0
      if nodrawcolor then do_DRAWCOLOR=0
      if strpos(Taglist,'SYMSIZE') ne -1 then do_SYMSIZE=1 else do_SYMSIZE=0
      if strpos(Taglist,'SYMTHICK') ne -1 then do_SYMTHICK=1 else do_SYMTHICK=0
      if strpos(Taglist,'LEGEND') ne -1 then do_LEGEND=1 else do_LEGEND=0

      if n_elements(xdata) ne n_elements(set(0).subset) then begin
         message,'Data arrays and SUBSET arrays must be the same size',/cont
         stop
      endif
;
; Setup overall plot scale range
;
      active=set(0).subset
      for i=1,n_elements(set)-1 do active=active or set(i).subset
      anyactive=active and inrange
      indices=where(anyactive,count)
      if count lt 1 then begin
         print,'nothing to plot'
         goto,exit
      endif
;
; Plot each set with legend
;
      dyleg=.035*legsize & dxleg=.02
      for i=0,n_elements(set)-1 do begin
         act=set(i).subset and inrange
         indices=where(act,count)
         if count gt 0 then begin
            if do_psym then symbol=set(i).psym else symbol=psym
            if do_color then ccolor=set(i).color else ccolor=color
            if do_fillcolor then fcolor=set(i).fillcolor else fcolor=ccolor
            if do_drawcolor then dcolor=set(i).drawcolor else dcolor=ccolor
            if do_symsize then ssize=set(i).symsize else ssize=symsize
            if do_symthick then sthick=set(i).symthick else sthick=symthick
            print,format='("",$)'
            if symbol le 0 then oplot,xdata(indices),ydata(indices),thick=thk,psym=0,color=ccolor
            if do_fillcolor or do_drawcolor then begin
               for ii=0,n_elements(indices)-1 do begin
                  plots,xdata(indices(ii)),ydata(indices(ii)),$
                        psym=abs(psymbol(symbol,color=fcolor,thick=sthick)),symsize=ssize,$
			noclip=0,color=fcolor
                  if do_drawcolor then begin
                     plots,xdata(indices(ii)),ydata(indices(ii)),$
                        psym=abs(psymbol(symbol,color=dcolor,/nofill,thick=sthick)),symsize=ssize,$
            		noclip=0,color=dcolor
                  endif
               endfor
            endif else begin
               plots,xdata(indices),ydata(indices),$
                       psym=abs(psymbol(symbol,color=ccolor,thick=sthick)),symsize=ssize,noclip=0,color=ccolor
            endelse
            if dosigma then begin
               color_save=!p.color
               !p.color=ccolor
               width=.01*ssize
               errplot,xdata(indices),ydata(indices)-sigma(indices),ydata(indices)+sigma(indices),width=width
               !p.color=color_save
            endif
            if do_legend then begin 
               if set(i).legend ne '' and nolegend eq 0 then begin
                  if legendcolor eq 0 then legcolor=ccolor else legcolor=legendcolor
                  if do_fillcolor or do_drawcolor then begin
                     plots,[xleg],[yleg-nleg*dyleg],/normal,psym=psymbol(symbol,color=fcolor,thick=sthick),$
                           color=fcolor,symsize=ssize,noclip=0
                     if do_drawcolor then plots,[xleg],[yleg-nleg*dyleg],/normal,$
                               psym=psymbol(symbol,color=dcolor,/nofill,thick=sthick),$
                              color=dcolor,symsize=ssize,noclip=0
                  endif else begin
                     plots,[xleg],[yleg-nleg*dyleg],/normal,psym=psymbol(symbol,color=ccolor,thick=sthick),$
                           color=ccolor,symsize=ssize,noclip=0
                  endelse
                  xyouts,[xleg+dxleg],[yleg-nleg*dyleg-Legsize*.01],/normal,$
                         set(i).legend,color=legcolor,charsize=legsize,$
                         noclip=0
                  nleg=nleg+legenddir
                  lastyleg=yleg-nleg*dyleg
               endif
            endif
         endif
      endfor
;
; Plot shot trajectories
;
      if trajectory then begin
         act=set(0).subset
         for i=1,n_elements(set)-1 do begin
            act=set(i).subset or act
         endfor
         act=act and inrange
         indices=where(act,count)
         if count gt 0 then begin
            Ushots=shots(unique_shots(shots))
            for ii=0,n_elements(Ushots)-1 do begin
               kk=where(Ushots(ii) eq shots(indices),count)
               if count gt 1 then plots,xdata(indices(kk)),ydata(indices(kk)),noclip=0,color=1
            endfor
         endif
      endif
   endelse
;
   if inquire then begin
Redo:
      nindices=[-1]
More:
;      if !D.name eq 'X' then begin
;         ind=where(anyactive,count)
;         box_cursor, x0, y0, dx, dy
;         pll = convert_coord(x0, y0, /device, /to_data)
;         pur = convert_coord(x0+dx, y0+dy, /device, /to_data)
;         found=(xdata ge pll(0)) and (xdata le pur(0)) and (ydata ge pll(1)) and (ydata le pur(1)) and anyactive
;         nindices=where(found)
;      endif else begin
         ind=where(anyactive,count)
         print,'Enter lower left corner of selection box'
         cursor,xll,yll,/down,/data
         wait,0.2
         print,'Enter upper right corner of selection box'
         cursor,xur,yur,/down,/data
         oplot,[xll,xll,xur,xur,xll],[yll,yur,yur,yll,yll],linestyle=1
         found=(xdata ge xll) and (xdata le xur) and (ydata ge yll) and (ydata le yur) and anyactive
         nindices=[nindices,where(found)]
         ans=' '
         on_ioerror,retry
Retry:
         print,format='("(R)edo, (M)ore points, (E)xit",$)'
         read,ans
         ans=strupcase(ans)
         if ans eq 'R' then goto,redo
         if ans eq 'M' then goto,more
         if ans eq 'E' then goto,done
         goto,retry
Done:
         on_ioerror,null
         if n_elements(nindices) gt 1 then nindices=nindices(1:n_elements(nindices)-1)
;      endelse
   endif
   iactive=[where(anyactive)]
Exit:
   return
   end
