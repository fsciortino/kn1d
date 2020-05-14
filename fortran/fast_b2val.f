c=============================================================================
c  Wrapper Routine

        subroutine fast_b2val(ARGC, ARGV)
        INTEGER*4               ARGC    !Argument count
        INTEGER*8               ARGV(*) !Vector of pointers to argments

C       Local variables

        INTEGER                 ARG_CNT

C       The argument count is passed in by value. Get the location of
C       this value in memory (a pointer) and convert it into an
C       Fortran integer.

        ARG_CNT = LOC(ARGC)

C	Insure that we got the correct number of arguments

	IF(ARG_CNT .ne. 14)THEN
	   WRITE(*,*)'fast_b2val => incorrect number of arguments'
	   RETURN
	ENDIF

C       To convert the pointers to the IDL variables contained in ARGV
C       we must use the Fortran function %VAL. This funcion is used
C       in the argument list of a Fortran sub-program. Call the Fortran
C       subroutine that will actually perform the desired operations.
C       Set the return value to the value of this function.

        call fast_b2val_code( %val(ARGV(1)), %val(ARGV(2)),
     &                        %val(ARGV(3)), %val(ARGV(4)),
     &                        %val(ARGV(5)), %val(ARGV(6)),
     &                        %val(ARGV(7)), %val(ARGV(8)),
     &                        %val(ARGV(9)), %val(ARGV(10)),
     &                        %val(ARGV(11)), %val(ARGV(12)),
     &                        %val(ARGV(13)), %val(ARGV(14)) )


        RETURN

        END

c=============================================================================
c
c fast_b2val.for
c
c   Calls public domain routine b2val to evaluate the derivative of a two-dimensional 
c tensor-product spline given its the tensor-product B-spline representation. 
c Returns a vector of derivative evaluations at the positions specified by 
c vectors X and Y.
c
      SUBROUTINE fast_b2val_code(npts,IXD,IYD,Xd,Yd,Kx,Ky,Xk,Yk,
	1 NXCOEF,NYCOEF,BSC,RESULT,WORK)
c
c______________________________________________________________________________
c Input:
c   npts - integer*4, number of data points to evaluate 
c   IXD  - integer*4, Order of the derivative in X-direction
c   IYD  - integer*4, Order of the derivative in Y-direction 
c   Xd   - real*4, size npts, containing x coordinates at which the deriv is to be evaluated 
c   Yd   - real*4, size npts, containing y coordinates at which the deriv is to be evaluated 
c   KX	 - integer*4, Order of spline in X direction. 2=linear spline, 3=quadratic 
c   KY	 - integer*4, Order of spline in Y direction. 2=linear spline, 3=quadratic 
c   Xk   - real*4, size NXCOEF, X knot sequence -> from B2ink 
c   Yk   - real*4, Y knot sequence -> from B2ink
c   NXCOEF - integer*4, size of Xk
c   NYCOEF - integer*4, size of Yk
c   BSC   - Real*4 size NXCOEF*NYCOEF, tensor product B-spline coefficients -> from B2ink
c   WORK    Real 1D array (size 3*max(KX,KY) + KY), A working storage array.
c
c Output:
c   Result - real*4, size npts
c______________________________________________________________________________
c
	Integer*4 npts,IXD,IYD,KX,KY,NXCOEF,NYCOEF
	Real*4 Xd(npts),Yd(npts),Xk(NXCOEF),Yk(NYCOEF)
	Real*4 BSC(NXCOEF*NYCOEF),RESULT(npts),WORK(*)

        do i=1,npts
           XdP=Xd(i)
           YdP=Yd(i)
           RESULT(i)=B2VAL(XdP,YdP,IXD,IYD,Xk,Yk,NXCOEF,NYCOEF,
	1 Kx,Ky,BSC,WORK)
        end do
	return
	end
