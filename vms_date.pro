;
; VMS_Date.pro
;
; returns current date/time string in MMM-DD-YYYY HH:MM:SS.SS format
;
function vms_date,dummy
time = systime()
return,strupcase(strmid(time,8,2)+'-'+strmid(time,4,3)+'-'+strmid(time,20,4)+' '+strmid(time, 11,8))
end
