;+
; SigmaV_rec_H1s.pro
;
; Returns maxwellian averaged <sigma V) for electron-ion radiative
; recombination to the atomic hydrogen in the 1s state.
; Coefficients are taken from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.32.
;
;________________________________________________________________________________
   Function SigmaV_rec_H1s,Te
;________________________________________________________________________________
;  Input:
;	Te	- fltarr(*) or float, electron temperature (eV)
;
;  Output:
;	returns <sigma V> for 0.1 < Te < 2e4.
;	units: m^3/s
;________________________________________________________________________________
;-
   t=type_of(Te,nDim=nDim)
   _Te=[Te]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
;
; Data for nl = 1s
;
   n=1
   Ry=13.58
   Eion_n=Ry/n
   Anl=3.92
   Xnl=0.35
;
   Bn=Eion_n/_Te
   result=Anl*1e-14*sqrt(Eion_n/Ry)*(Bn^1.5/(Bn+Xnl))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
