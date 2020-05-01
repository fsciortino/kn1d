;
; Test_INTERP_FVRVXX.PRO
;

      nx_a=100
      xa_a=0.0
      xb_a=0.02
      x_a=xa_a+(xb_a-xa_a)*findgen(nx_a)/(nx_a-1)

      Ti_a=1.0*exp(x_a/.025)

      nv_a=15
;      E0=[.01,0.1]
      Create_VrVxMesh,E0=E0,nv_a,Ti_a,vx_a,vr_a,Tnorm_a

      nvr_a=n_elements(vr_a)
      nvx_a=n_elements(vx_a)
      nx_a=n_elements(x_a)

      ip_a=where(vx_a gt 0)

      fa=dblarr(nvr_a,nvx_a,nx_a)
      Tneut_a=5.0+1.0*findgen(nx_a)/(nx_a-1)
      mH=1.6726231D-27
      q=1.602177D-19
      mu=1
      Vtha=sqrt(2*q*Tnorm_a/(mu*mH))
      Uxa=0.5*vtha

      for j=0,nvx_a-1 do begin
         for k=0,nx_a-1 do begin
            fa(*,j,k)=exp(-x_a(k)/.01)*exp(-(vr_a(*)^2+(vx_a(j)-Uxa/vtha)^2)/(Tneut_a(k)/Tnorm_a))
         endfor
      endfor

      nx_b=100
      xa_b=0.0
      xb_b=0.05
      x_b=xa_b+(xb_b-xa_b)*findgen(nx_b)/(nx_b-1)

      Ti_b=10.0*exp(x_b/.025)

      nv_b=10
      Create_VrVxMesh,nv_b,Ti_b,vx_b,vr_b,Tnorm_b

      nvr_b=n_elements(vr_b)
      nvx_b=n_elements(vx_b)
      nx_b=n_elements(x_b)

      ip_b=where(vx_b gt 0)

      fb_test=dblarr(nvr_b,nvx_b,nx_b)
      
      Interp_scalarX,Tneut_a,X_a,_Tneut_a,X_b,warn=warn,debug=debug

      debug=1
      correct=1
      Interp_fVrVxX,fa,Vr_a,Vx_a,X_a,Tnorm_a,fb,Vr_b,Vx_b,X_b,Tnorm_b,warn=warn,debug=debug,correct=correct

      levels=10^(-findgen(10))
      levels=reverse(levels)
      contour,fa(*,*,0),sqrt(Tnorm_a)*Vr_a,sqrt(Tnorm_a)*Vx_a,levels=levels
      contour,fb(*,*,0),sqrt(Tnorm_b)*Vr_b,sqrt(Tnorm_b)*Vx_b,color=2,/overplot,levels=levels


end
