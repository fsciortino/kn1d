Function N0_from_Balmer_Alpha,B_Alpha,Density,Te,Source=Source,$
            Ionization=Ionization,Recombination=Recombination,create=create
;________________________________________________________________________________
; Input:
;  	B_alpha 	- fltarr, local Balmer-alpha emissivity (watts m^-3)
;  	Density		- fltarr, electron density (=hydrogen ion density) (m^-3)
;  	Te		- fltarr, electron temperature (eV)
;
; Keywords:
;	Source       - returns the local net ionization source rate m^-3 s^-1
;	Ionization   - returns the local ionization rate m^-3 s^-1
;	Recombination- returns the local recombination rate m^-3 s^-1
;	create	- if set, then create bi-cubic spline coefficients for
;		  interpolation of r0(p) r1(p) and save them in the
;		  default save set. 
;________________________________________________________________________________
; History:
;    Coding by B. LaBombard  6/29/99
;    Coefficients from J. Terry's idl code JH_RATES.PRO
;   
;-
   common JH_Coef,DKnot,TKnot,order,LogR_BSCoef,LogS_BSCoef,LogAlpha_BSCoef,A_Lyman,A_Balmer
   if keyword_set(create) then create_JH_BSCoef
   if type_of(LogR_BSCoef) eq 0 then begin
      restore,'/home/labombard/edge/jh/jh_bscoef.dat'
   endif
;
;  From Johnson and Hinnov, eq (11):
;
;   n(3) =  ( r0(3) + r1(3)*n(1)/NHsaha(1) )*NHsaha(3)
;
;  Inverting:
;
;  n(1) =  ( n(3)/NHsaha(3) - r0(3) )*NHsaha(1)/r1(3)
;
   if n_elements(Density) ne n_elements(Te) then message,'Number of elements of Density and Te are different!'
   if n_elements(Density) ne n_elements(B_alpha) then message,'Number of elements of Density and B_alpha are different!'
   N0=Density & N0(*)=1.0e32
   n3=N0
   Source=N0
   Ionization=N0
   Recombination=N0
   r03=JHR_Coef(Density,Te,0,3)
   r13=JHR_Coef(Density,Te,1,3)
   NHSaha1=NHSaha(Density,Te,1)
   NHSaha3=NHSaha(Density,Te,3)
   S=JHS_Coef(Density,Te)
   Alpha=JHAlpha_Coef(Density,Te)
   ok=where(B_alpha lt 1e32 and B_alpha gt 0 and r03 lt 1.0e32 and r13 lt 1.0e32 and $
             NHSaha1 lt 1.0e32 and NHSaha3 lt 1.0e32 and S gt 0 and S lt 1.0e32 and Alpha gt 0 and Alpha lt 1.0e32,count)
   if count gt 0 then begin
      n3(ok)=B_alpha(ok)/(A_Balmer(0)*13.6057*(0.25-1.0/9.0)*1.6e-19)
      N0(ok)=(n3(ok)/NHSaha3(ok) - r03(ok))*NHSaha1(ok)/r13(ok)
;
; Evaluate ionization, recombination and net source from Johnson-Hinnov equation (12)
;
;   Ionization = n(1)*Density*S
;   Recombination = Density^2*Alpha
;   Source = Ionization - Recombination
;
      Ionization(ok)=N0(ok)*(Density(ok)*S(ok))
      Recombination(ok)=Density(ok)*(Density(ok)*Alpha(ok))
      Source(ok)=Ionization(ok)-Recombination(ok)
   endif
   return,N0
   end
