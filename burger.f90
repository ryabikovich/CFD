Program Pr
      Implicit none

      INTEGER, parameter:: IO = 12 ! input-output unit
      INTEGER NX,NT,I,J,ID,m, ISCHEME
      REAL(8),ALLOCATABLE :: U(:),UN(:),X(:)
      REAL(8) L,h,CFL,dt,t,Time
      REAL(8) C, C0, C1, pi, k,G,FE

      WRITE(*,*) 'Read input file'
      OPEN(IO,FILE='Input.txt')
      READ(IO,*) L
      READ(IO,*) m
      READ(IO,*) C
      READ(IO,*) C0, C1
      READ(IO,*) NX
      READ(IO,*) NT
      READ(IO,*) CFL
      READ(IO,*) ISCHEME
      CLOSE(IO)

      ALLOCATE(U(0:NX+1),UN(0:NX+1),X(0:NX+1))

      pi=3.14159265359d0
      k = m*pi/L
      h = 1.0d-2 !x  step(0<=x<=1)
      dt = 1.0d-2
      Time = dt*NT

      WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX
      WRITE(*,*) 'CFL=', CFL, 'dt=', dt, 'Time=', Time, 'NT=', NT

      X(0)=-h
      DO I=1, NX+1
       X(I)=X(I-1)+h
       write(*,*) X(I)
      END DO


      U(:)=0.0
      UN(:)=0.0

      t=0.0d0
      CALL InitValue(X,U,NX,k)
!      CALL BoundValue()
      OPEN(IO,FILE='Res.dat')
      call Output(U,X,NX,t)
!-------------------------  Solve equation ------------------
    do
        t = t + dt
        if (t .gt. Time) then
            exit
        end if
        do I = 1,NX
            if (U(I) .ge. 0.0d0) then
              UN(I) =  -((U(I))**2.0d0-(U(I-1))**2.0d0)*dt/(2.0d0*h) + U(I)
               write(IO,*) UN(I), X(I), t
            else
               UN(I) =  -((U(I+1))**2.0d0-(U(I))**2.0d0)*dt/(2.0d0*h) + U(I)
               write(IO,*) UN(I), X(I), t
            end if


        end do
        call BoundValue(UN,NX)
!        call Output(UN,X,NX,t)
        U = UN
!        call Output(U,X,NX,t)
    end do

!-------------------------  Output results ------------------

!      OPEN(IO,FILE='Res.dat')
!      call Output()
      CLOSE(IO)

      End program

!----------------------- Set Initial Value -----------------
      subroutine InitValue(X,U,NX,k)
      IMPLICIT NONE
      integer, intent(in):: NX
      real*8,intent(in)::k
      integer i
       REAL*8,intent(in out):: U(0:NX+1)
       REAL*8,intent(in out):: X(0:NX+1)
      U(0) = sin(k*X(0))
      do I = 1, NX+1
        U(I) = sin(k*X(I))
      end do

      end subroutine

!----------------------- Set Boundary Condition ------------
      subroutine BoundValue(UN,NX)
      IMPLICIT NONE
      integer, intent(in):: NX
      REAL*8,intent(in out):: UN(0:NX+1)
      UN(0) = UN(NX-1)
      UN(NX+1)=UN(2)
      end subroutine

!----------------------- Output results ------------
      subroutine Output(U,X,NX,t)
      IMPLICIT NONE
     INTEGER, parameter:: IO=12
      integer,intent(in)::NX
      real*8,intent(in out)::t
      integer i
      REAL*8,intent(in out):: U(0:NX+1)
      REAL*8,intent(in out):: X(0:NX+1)
      write(IO,*)U(0),X(0),t
      do i=1,NX+1
        write(IO,*) U(i),X(i),t
      end do


      end subroutine

