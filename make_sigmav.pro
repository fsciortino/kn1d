;
; Make_SigmaV.pro
;
Function Make_SigmaV,E_Particle,mu_particle,T_target,mu_target,sigma_function
;
;________________________________________________________________________________
;  Input:
;	E_Particle	- fltarr(*), eV
;	Mu_particle 	- float
;	T_Target	- float, eV
;	mu_target	- float
;	Sigma_Function	- string, name of function to return Sigma (m^2) for inputted E_Particle (eV)
;
; Result SigmaV (m^2 s^-1)
;________________________________________________________________________________
; Test particle has velocity Vxa
; Target particles are a maxwellian in Vrb and Vb
;
   Trange=[E_Particle,T_target]
   Trange=Trange(sort(Trange))
   nvb=100
   Create_VrVxMesh,nvb,Trange,vxb,vrb,Tnorm
   Make_dVr_dVx,vrb,vxb,Vr2pidVrb,VrVr4pidVrb,dVxb
   mH=1.6726231D-27
   q=1.602177D-19				
   vth=sqrt(2*q*Tnorm/(mu_target*mH))
;
; Set normalized particle velocities
;
   Vxa=sqrt(2*q*E_Particle/(mu_particle*mH))/vth
   nvxa=n_elements(vxa)
   nvrb=n_elements(vrb)
   nvxb=n_elements(vxb)
;
; Set Normalized Target Distribution Function
;
   fi_hat=dblarr(nvrb,nvxb)
   for i=0,nvrb-1 do begin
      fi_hat(i,*)=exp(-(vrb(i)^2+vxb(*)^2) * Tnorm/T_target)
   endfor
   fi_hat=fi_hat/(total(Vr2pidVrb*(Fi_hat#dVxb)))

;
; Compute relative velocity at each mesh point
;
   Vrel=dblarr(nvrb,nvxb,nvxa)
   for k=0,nvxa-1 do begin
      for i=0,nvrb-1 do begin
         Vrel(i,*,k)=sqrt(Vrb(i)^2+(Vxb(*)-Vxa(k))^2)
      endfor
   endfor
;
; Get sigma for inputted E_Particle
;
   istat=execute('sig='+Sigma_Function+'(0.5*Vth*vth*Vrel*Vrel*mu_particle*mH/q)')
;
; Compute Sigmav by integrating sigma x Vrel x Fi_hat over velocity space
;
   sigv=dblarr(nvxa)
   for k=0,nvxa-1 do sigv(k)=vth*total(Vr2pidVrb*((sig(*,*,k)*Vrel(*,*,k)*Fi_hat)#dVxb))
   return,sigv
   end
