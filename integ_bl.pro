;
; Integ_BL.pro
;+
Function Integ_BL,X,Y,value_only=value_only,rev=rev
_Y=Y
_X=X
if keyword_set(rev) then begin
   _Y=reverse(Y)
   _X=-reverse(X)
endif
ans=_Y & ans(*)=0
for ii=1,n_elements(_Y)-1 do ans(ii)=ans(ii-1)+(_X(ii)-_X(ii-1))*0.5*(_Y(ii)+_Y(ii-1))
if keyword_set(value_only) then begin
   return,ans(n_elements(ans)-1)
endif else begin
   if keyword_set(rev) then return,reverse(ans)
   return,ans
endelse
end
;-
