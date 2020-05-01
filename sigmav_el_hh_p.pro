;+
; SigmaV_EL_HH_P.pro
;
;  Evaluates the elastic scattering (momentum transfer) <sigma v> for
;  H2 molecules with energy E impacting on protons with temperature T.
;
;  Data generated from cross sections tablulated in:
;
; Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 -
; Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and
; Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
; Passovets, page 305.
;
;________________________________________________________________________________
Function SigmaV_EL_HH_P,T,E,use_Bspline=use_Bspline
;________________________________________________________________________________
;  Input:
;	T	- fltarr(*) or float, proton temperature (eV)
;	E	- fltarr(*) or float, molecule mono-energy (eV)
;
;  Output:
;	returns <sigma V> for 0.1 < T < 1e3 and 0.1 < E < 1e4 
;	units: m^3/s
;       if T and/or E is outside this range, the value on the boundary is returned
;________________________________________________________________________________
;-
   if n_elements(E) ne n_elements(T) then message,'Number of elements of E and T are different!'
   _E=E > .0001
   _T=T > .0001
   _E=_E < 1.0e5
   _T=_T < 1.0e5
   LEP=Alog(_E)
   LTH2=Alog(_T)
   if keyword_set(use_Bspline) then begin
      common JSigmaV_EL_HH_P,EKnot_EL_H2_P,TKnot_EL_H2_P,order_EL_H2_P,LogSigmaV_EL_H2_P_BSCoef
      if type_of(LogSigmaV_EL_H2_P_BSCoef) eq 0 then begin
         restore,'sigmav_el_h2_p_bscoef.dat'
      endif
      LEP=LEP > min(Eknot_EL_H2_P)
      LEP=LEP < max(Eknot_EL_H2_P)
      LTH2=LTH2 > min(Tknot_EL_H2_P)
      LTH2=LTH2 < max(Tknot_EL_H2_P)
      Result=exp(BS2DR(0,0,LEP,LTH2,order_EL_H2_P,order_EL_H2_P,Eknot_EL_H2_P,Tknot_EL_H2_P,LogSigmaV_EL_H2_P_BSCoef))
   endif else begin
      common SIGMAV_EL_H2_P_DATA,Ln_E_Particle,Ln_T_Target,SigmaV,nEP,nT
      if type_of(Ln_E_Particle) eq 0 then begin
         restore,'sigmav_el_h2_p_data.dat'
         nEP=n_elements(Ln_E_Particle)-1
         nT=n_elements(Ln_T_Target)-1
      endif
      LEP=LEP > min(Ln_E_Particle)
      LEP=LEP < max(Ln_E_Particle)
      LTH2=LTH2 > min(Ln_T_Target)
      LTH2=LTH2 < max(Ln_T_Target)
      iE=(LEP-Ln_E_Particle(0))*nEP/(Ln_E_Particle(nEP)-Ln_E_Particle(0))
      iT=(LTH2-Ln_T_Target(0))*nT/(Ln_T_Target(nT)-Ln_T_Target(0))
      result=interpolate(SigmaV,iE,iT)
   endelse
   return,result
   end
