;
; Plot_sigma_v_ionization.pro
;
;  Plots H ionization rate, <sigma V> using Johnson-Hinnov coefficients.
;
common JH_Coef,DKnot,TKnot,order,LogR_BSCoef,LogS_BSCoef,LogAlpha_BSCoef,A_Lyman,A_Balmer
Tmin=0.35
Tmax=706.0
Te=exp(alog(Tmin)+(alog(Tmax)-alog(Tmin))*findgen(1001)/1000)
density1=replicate(1.0e18,1001)
density2=replicate(1.0e19,1001)
density3=replicate(1.0e20,1001)
;
SV1=JHS_Coef(density1,Te)
SV2=JHS_Coef(density2,Te)
SV3=JHS_Coef(density3,Te)
!x.ticklen=1
!y.ticklen=1
dbplot,/logx,/logy,te,SV3,$
   title='e!U-!N + H -> H!U+!n + e!U-!N + e!U-!N',$
   xtitle='Te (eV)',ytitle='<sigma v> cm!U3!N s!U-1!N',/nodata,$
   yrange=[1.0e-11,1.0e-7],xrange=[0.1,2000.0],xstyle=1
dboplot,te,SV1*1.0e6,color=2
dboplot,te,SV2*1.0e6,color=3
dboplot,te,SV3*1.0e6,color=4
xyouts,/normal,.15,.85,'n!De!N = 10!U14!N cm!U-3!N',color=2
xyouts,/normal,.15,.8,'n!De!N = 10!U13!N cm!U-3!N',color=3
xyouts,/normal,.15,.75,'n!De!N = 10!U12!N cm!U-3!N',color=4

   yleg=.25
   dy=.03
   cs=0.8
   xyouts,.42,yleg,/normal,'data from collisional-radiative model',charsize=cs
   xyouts,.42,yleg-dy,/normal,'of L.C.Johnson and E. Hinnov, J. Quant. ',charsize=cs
   xyouts,.42,yleg-2*dy,/normal,'Spectrosc. Radiat. Transfer. vol. 13 pp.333-358',charsize=cs
end
