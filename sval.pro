;
; sval.pro
;
function sval,val,len=len
s=strtrim(string(val),2)
if keyword_set(len) then s=STRMID(s,0,len)
return,s
end
