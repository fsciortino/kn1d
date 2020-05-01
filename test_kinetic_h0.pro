;
; Test_kinetic_H0.pro
;
   tek
   key_default,name,'test'
   file='Kinetic_H0_'+name+'.dat'
   key_default,plot,0
   key_default,debug,0
   key_default,debrief,1
   key_default,test,1
   key_default,simple_cx,1
   compute_errors=1
   nBC=0.0
;
; Test #1 - 3 eV neutrals
;
   if test eq 1 then begin
      nx=100
      xa=0.0
      xb=0.05
      x=xa+(xb-xa)*findgen(nx)/(nx-1)

      Ti=10.0*exp(x/.025)

      nv=20
      Create_VrVxMesh,nv,Ti,vx,vr,Tnorm

      iz=where(vx eq 0)
      iz=iz(0)
      ip=where(vx gt 0)
      nvr=n_elements(vr)
      nvx=n_elements(vx)
      nx=n_elements(x)

      GammanBC=1.0e22
      n=1.0e19*exp(x/.030)
      Te=Ti
      fnBC=fltarr(nvr,nvx)
      Tneut=3.0
      for i=iz+1,nvx-1 do begin
         fnBC(*,i)=exp(-(vr(*)^2+vx(i)^2)/(Tneut/Tnorm))
      endfor
      fSH0=fltarr(nvr,nvx,nx)
  endif
  if test eq 2 then begin
;
; Test #2 no ionization, Ti=T0
;
      max_gen=100
      nx=70
      xa=0.0
      xb=0.05
      x=xa+(xb-xa)*findgen(nx)/(nx-1)

      Ti=replicate(10.0,nx)

      nv=20
      Create_VrVxMesh,nv,Ti,vx,vr,Tnorm

      iz=where(vx eq 0)
      iz=iz(0)
      ip=where(vx gt 0)
      nvr=n_elements(vr)
      nvx=n_elements(vx)
      nx=n_elements(x)

      GammanBC=1.0e22
;      n=1.0e19*exp(x/.030)
      n=replicate(5.0e19,nx)
      Te=replicate(0.1,nx)
      fnBC=fltarr(nvr,nvx)
      Tneut=Ti(0)
      for i=iz+1,nvx-1 do begin
         fnBC(*,i)=exp(-(vr(*)^2+vx(i)^2)/(Tneut/Tnorm))
      endfor
      fSH0=fltarr(nvr,nvx,nx)
   endif
  if test eq 3 then begin
;
; Test #3 with large ionization/charge exchange fraction, Ti=T0 
;
      max_gen=100
      nx=70
      xa=0.0
      xb=0.05
      x=xa+(xb-xa)*findgen(nx)/(nx-1)

      Ti=replicate(10.0,nx)

      nv=40
      Create_VrVxMesh,nv,Ti,vx,vr,Tnorm

      iz=where(vx eq 0)
      iz=iz(0)
      ip=where(vx gt 0)
      nvr=n_elements(vr)
      nvx=n_elements(vx)
      nx=n_elements(x)

      GammanBC=1.0e22
      n=replicate(5.0e19,nx)
      Te=replicate(30.0,nx)
      fnBC=fltarr(nvr,nvx)
      Tneut=Ti(0)
      for i=iz+1,nvx-1 do begin
         fnBC(*,i)=exp(-(vr(*)^2+vx(i)^2)/(Tneut/Tnorm))
      endfor
      fSH0=fltarr(nvr,nvx,nx)
   endif
;
   print,'Test ='+sval(test)
   multiplot,x,n,x,te,x,ti,title='Inputted profiles',xtitle='x',$
             ytitle=['n','Te','Ti'],color=[2,3,4],$
             varlabel=['Density (m^-3)','Te (ev)','Ti (eV)']
   press_return
   plot,vx,fnBC(0,*),/nodata,yrange=[0,max(fnBC)],title='Inputted fnBC'
   for i=0,nvr-1 do oplot,vx,fnBC(i,*),color=(i mod 8)+2
   press_return

   truncate=1.e-4
   print,'truncate=',truncate

   mu=2

   Simple_CX=1
   Kinetic_H0,vx,vr,x,Tnorm,mu,Ti,Te,n,fnBC,nBC,GammanBC,fSH0,$
        fn,n0,gammax0,vx0,p0,pi0_xx,pi0_yy,pi0_zz,$
 	T0,qx0,qx0_total,NetH0Source,Sion,Qin,Rxin,Qin_total,Albedo,truncate=truncate,$
        No_Johnson_Hinnov=No_Johnson_Hinnov,No_Recomb=No_Recomb,Simple_CX=Simple_CX,$
        max_gen=max_gen,Max_dx=Max_dx,$
        error=error,compute_errors=compute_errors,$
        vbar_error=vbar_error,mesh_error=mesh_error,max_mesh_error=max_mesh_error,$
        ave_mesh_error=ave_mesh_error,moment_error=moment_error,$
	max_moment_error=max_moment_error,qx0_total_error=qx0_total_error,$
        Qin_total_error=Qin_total_error,NetH0Source_Error=NetH0Source_Error,$
	plot=plot,debug=debug,debrief=debrief

   save,file=file,vx,vr,x,Tnorm,mu,Ti,Te,n,fnBC,GammanBC,fSH0,$
        fn,n0,gammax0,vx0,p0,pi0_xx,pi0_yy,pi0_zz,$
 	T0,qx0,qx0_total,NetH0Source,Qin,Rxin,Qin_total,Albedo,truncate,$
        No_Johnson_Hinnov,Simple_CX,$
        max_gen,Max_dx,$
        error,compute_errors,$
        vbar_error,mesh_error,max_mesh_error,$
        ave_mesh_error,moment_error,$
	max_moment_error,qx0_total_error,$
        Qin_total_error,NetH0Source_Error,$
	plot,debug,debrief

   end
