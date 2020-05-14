c=============================================================================
c  Wrapper Routine

        subroutine call_b2ink(ARGC, ARGV)
        INTEGER*4               ARGC    !Argument count
        INTEGER*8               ARGV(*) !Vector of pointers to argments

C       Local variables

        INTEGER                 ARG_CNT

C       The argument count is passed in by value. Get the location of
C       this value in memory (a pointer) and convert it into an
C       Fortran integer.

        ARG_CNT = LOC(ARGC)

C	Insure that we got the correct number of arguments

	IF(ARG_CNT .ne. 13)THEN
	   WRITE(*,*)'call_b2ink => incorrect number of arguments'
	   RETURN
	ENDIF

C       To convert the pointers to the IDL variables contained in ARGV
C       we must use the Fortran function %VAL. This funcion is used
C       in the argument list of a Fortran sub-program. Call the Fortran
C       subroutine that will actually perform the desired operations.
C       Set the return value to the value of this function.

        call b2ink( %val(ARGV(1)), %val(ARGV(2)),
     &                        %val(ARGV(3)), %val(ARGV(4)),
     &                        %val(ARGV(5)), %val(ARGV(6)),
     &                        %val(ARGV(7)), %val(ARGV(8)),
     &                        %val(ARGV(9)), %val(ARGV(10)),
     &                        %val(ARGV(11)), %val(ARGV(12)),
     &                        %val(ARGV(13)) )


        RETURN

        END
