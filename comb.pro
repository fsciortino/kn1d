;+
; Comb.pro
;
Function Comb,data,npts=npts
if not keyword_set(npts) then npts=2000
if n_elements(data) gt npts then begin
   ii=indgen(npts+1)*(n_elements(data)-1)/npts
endif else begin
   ii=indgen(n_elements(data))
endelse
return,data(ii)
end
;-
