;
; Create_SIGMAV_EL_H2_P_DATA.PRO
;
;   Creates a 2D SigmaV data table in particle energy and target temperature
;   and saves it as a save set.
;
Pro Create_SIGMAV_EL_H2_P_DATA
;
   Sigma_Function='sigma_el_P_HH'
;
; Particle is HH
;
; Target is P
;
   mE=50
   Emin=0.1 & Emax=1.0e3
   nT=50
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
   Ln_E_Particle=alog(E_Particle)
   Ln_T_Target=alog(T_Target)

   print,'Saving results in file: /home/labombard/edge/neutrals/sigmav_el_h2_p_data.dat'
   save,file='/home/labombard/edge/neutrals/sigmav_el_h2_p_data.dat',$
        Ln_E_Particle,Ln_T_Target,SigmaV
   return
   end
