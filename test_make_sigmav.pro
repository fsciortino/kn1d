;
; TEST_MAKE_SIGMAV.PRO
;
   E_Particle=[0.1,1,10,50,100,1000,2e4]
   mu_particle=01.0
   T_Target=10.0
   Mu_Target=1.0
   Sigma_Function='sigma_cx_h0'
   SigmaV=Make_SigmaV(E_Particle,mu_particle,T_target,mu_target,sigma_function)
   
   _T_Target=replicate(T_Target,n_elements(E_Particle))
   sigv=sigmav_cx_h0(_T_Target,E_Particle)
   
   print,sigmav,sigv
   end
