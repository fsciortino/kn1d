;
; Test_Lyman_Alpha.pro
;
tek
; 
; density scan test
;
density=10^(16.0+(21.-16)*findgen(101)/100)
Te=replicate(20,101)
n0=replicate(3.3e16,101)
ly=lyman_alpha(density,te,n0)
dbplot,density,ly,/logx,/logy,title='Lyman-alpha',xtitle='Density (m^-3)',ytitle='watts m^-3'
press_return

Lyman_Alpha_N0_Source,Density,Te,Ly,N0,Source,Error=Error
dbplot,density,N0,/logx,/logy,title='N0 from Lyman_Alpha',xtitle='Density (m^-3)',ytitle='m^-3'
press_return
N0_new=N0_FROM_LYMAN_ALPHA(ly,Density,Te,source=source2)
dboplot,density,N0_new,color=2
press_return
dbplot,density,Source,/logx,/logy,title='Source from Lyman_Alpha',xtitle='Density (m^-3)',ytitle='m^-3 s^-1'
dboplot,density,Source2,color=2

press_return

;
; Te scan test
;
Te=10^(alog10(0.35)+(alog10(700)-alog10(0.35))*findgen(101)/100)
density=replicate(1.0e19,101)
n0=replicate(3.3e16,101)
ly=lyman_alpha(density,te,n0)
dbplot,Te,ly,/logx,/logy,title='Lyman-alpha',xtitle='Te (eV)',ytitle='watts m^-3'
press_return

Lyman_Alpha_N0_Source,Density,Te,Ly,N0,Source,Error=Error
dbplot,Te,N0,/logx,/logy,title='N0 from Lyman_Alpha',xtitle='Te (eV)',ytitle='m^-3'
press_return
N0_new=N0_FROM_LYMAN_ALPHA(ly,Density,Te,source=source2)
dboplot,Te,N0_new,color=2
press_return
dbplot,Te,Source,/logx,/logy,title='Source from Lyman_Alpha',xtitle='Te (eV)',ytitle='m^-3 s^-1'
dboplot,Te,Source2,color=2

press_return

end
