function normal2data,norm,x=x,y=y
   if keyword_set(y) then begin
      scale=!Y.S(1)
      offset=!Y.S(0)
      type=!y.type
   endif else begin
      scale=!x.S(1)
      offset=!x.S(0)
      type=!x.type
   endelse
   data=(norm-offset)/scale
   if type ne 0 then data=10.0^data
   return,data
   end
