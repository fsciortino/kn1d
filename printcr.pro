;+
; PrintCR.pro
;
;   Prints a string with carriage return
;
pro printCR,string,rev=rev
if keyword_set(rev) then revon
print,string
if keyword_set(rev) then revoff
return
end
