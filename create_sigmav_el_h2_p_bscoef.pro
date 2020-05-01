;
; Create_SIGMAV_EL_H2_P_BSCOEF.PRO
;
;   Creates a save set storing Bi-cubic spline interpolation
; coefficients for parameters in Johnson-Hinov rate equations.
;
Pro Create_SIGMAV_EL_H2_P_BSCOEF
;
   Sigma_Function='sigma_el_P_HH'
;
; Particle is HH
;
; Target is P
;
   mE=5
   Emin=0.1 & Emax=1.0e3
   nT=5
   Tmin=0.1 & Tmax=1.0e3
   E_Particle=10^(alog10(Emin)+(alog10(Emax)-alog10(Emin))*findgen(mE)/(mE-1))
   Mu_Particle=2.0
   T_Target=10^(alog10(Tmin)+(alog10(Tmax)-alog10(Tmin))*findgen(nT)/(nT-1))
   Mu_Target=1.0
   SigmaV=fltarr(mE,nT)
   for iT=0,n_elements(T_Target)-1 do begin
      print,'Processing T='+sval(T_target(iT))
      SigmaV(*,iT)=Make_SigmaV(E_Particle,mu_particle,T_target(iT),mu_target,sigma_function)
   endfor

   print,'Computing B-Spline coefficients'
   order_EL_H2_P=4
   order=order_EL_H2_P
   LogE_Particle=alog(E_Particle)
   LogT_Target=alog(T_Target)
   bsnak,LogE_Particle,order,Eknot_EL_H2_P
   bsnak,LogT_Target,order,Tknot_EL_H2_P
   LogSigmaV_EL_H2_P=alog(SigmaV)
   bs2in,LogE_Particle,LogT_Target,LogSigmaV_EL_H2_P,order,order,EKnot_EL_H2_P,TKnot_EL_H2_P,LogSigmaV_EL_H2_P_BSCoef
   print,'Saving results in file: /home/labombard/edge/neutrals/sigmav_el_h2_p_bscoef.dat'
   save,file='/home/labombard/edge/neutrals/sigmav_el_h2_p_bscoef.dat',$
        EKnot_EL_H2_P,TKnot_EL_H2_P,order_EL_H2_P,LogSigmaV_EL_H2_P_BSCoef
   return
   end
