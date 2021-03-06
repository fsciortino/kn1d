;+
; SigmaV_cx_HH.pro
;
; Returns maxwellian averaged <sigma V) for charge exchange of molecular
; hydrogen. Coefficients are taken
; from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.292.
;
;________________________________________________________________________________
   Function SigmaV_CX_HH,T,E
;________________________________________________________________________________
;  Input:
;	T	- fltarr(*) or float, molecular ion [neutral] temperature (eV)
;	E	- fltarr(*) or float, molecular neutral [ion] mono-energy (eV)
;
;  Output:
;	returns <sigma V> for 0.1 < Te < 2e4 and 0.1 < E < 2e4 
;	units: m^3/s
;________________________________________________________________________________
;-
   if n_elements(T) ne n_elements(E) then message,'number of elements of T and E are different!'
   ty=type_of(T,nDim=nDim)
   _T=double([T])
   _E=double([E])
   alpha=dblarr(9,9)
;               n,m
;               E,T
;
; E Index      0			1			2
   alpha[0:2,*]=$
     [[-2.013143517466D+01,	1.875197914224D-01,	6.865479604288D-02],$
      [ 2.643458086299D-01,    -1.177247941077D-01,    -6.758032286178D-03],$
      [ 7.295645990688D-02,	6.053418575149D-03,    -1.068656224307D-02],$
      [-1.022454343675D-02,     7.350954380641D-03,     6.814213275702D-04],$
      [-4.801198168030D-03,    -1.111612877392D-03,     8.373319888351D-04],$
      [ 1.141613586234D-03,    -1.371389288760D-04,    -1.733761953296D-04],$
      [-3.388853048483D-05,     4.426148343648D-05,     9.992317920676D-06],$
      [-6.418225985394D-06,    -3.652063962019D-06,     1.351312819077D-07],$
      [ 3.555592819527D-07,     1.012701361110D-07,    -1.993091213299D-08]]
;
; E Index      3			4			5
   alpha[3:5,*]=$
     [[ 6.246595384100D-03,    -5.017891372102D-03,    -3.907644829287D-04],$
      [ 8.585003721992D-03,    -3.261863407467D-04,    -3.322528542186D-04],$
      [-9.371235639464D-04,     9.735708783528D-04,    -9.933049259228D-05],$
      [-8.156435157073D-04,     2.903991825737D-05,     3.223596225946D-05],$
      [ 1.392977576749D-04,    -9.316910697276D-05,     8.814981236658D-06],$
      [ 1.602610140599D-05,     1.464235749797D-05,    -2.944711701791D-06],$
      [-5.333970870280D-06,    -2.999105886511D-07,     2.275612517364D-07],$
      [ 4.285396408056D-07,    -7.184302986068D-08,    -3.265552364687D-10],$
      [-1.131561847140D-08,     3.678869095972D-09,    -3.639982258214D-10]]
;
; E Index      6			7			8
   alpha[6:8,*]=$
     [[ 2.786239030986D-04,    -2.942576591004D-05,     9.352275354690D-07],$
      [ 6.015471216449D-05,    -4.039435357369D-06,     9.730479674748D-08],$
      [-6.786246802840D-06,     1.438327767305D-06,    -5.530742535057D-08],$
      [-5.199055182831D-06,     2.852443990256D-07,    -4.825480212106D-09],$
      [ 6.675626166047D-07,    -1.325441927019D-07,     5.012529587757D-09],$
      [ 6.365231650682D-08,     1.872976659964D-08,    -1.014883015867D-09],$
      [-1.173422836715D-08,    -1.364602870139D-09,     9.566404348683D-11],$
      [ 1.585228996542D-10,     6.431866226702D-11,    -4.507074278992D-12],$
      [ 2.056662091085D-11,    -1.804254277469D-12,     9.042973335167D-14]]

   _E=_E > 0.1
   _E=_E < 2.01e4
   _T=_T > 0.1
   _T=_T < 2.01e4
   _alogE=alog(_E)
   _alogT=alog(_T)
   result=_E & result(*)=0.0
   for n=0,8 do for m=0,8 do result=result+alpha(n,m)*_alogE^n*_alogT^m
   result=EXP(result)*1D-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
