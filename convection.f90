      Program Pr
      Implicit none

      INTEGER, parameter:: IO = 12 ! input-output unit
      INTEGER NX,NT,I,J,ID,m, ISCHEME
      REAL(8),ALLOCATABLE :: U(:),UN(:),X(:)
      REAL(8):: X_theory(0:1000),theory(0:1000)
      REAL(8) L,h,CFL,dt,t,Time,p
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
      h = L/(NX-1)
      dt = CFL*h/c

      Time = dt*NT

      WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX, 'k=',k, 'beta=',k*h
      WRITE(*,*) 'CFL=', CFL, 'dt=', dt, 'Time=', Time, 'NT=', NT

      X(0)=-h
      DO I=1, NX+1
       X(I)=X(I-1)+h

      END DO
      write(*,*) X


      U(:)=0.0
      UN(:)=0.0

      t=0.0d0
      CALL InitValue(X,U,NX,k)
!    OPEN(1,FILE='Res1.dat')
!        p = (X(NX+1)-X(0))/1000
!        write(*,*)'p=',p
!        X_theory(0) = X(0)
!         DO I=1, 1000
!            X_theory(I)= X_theory(I-1)+p
!         END DO
!
!        theory(0) = sin(k*(x(0)-c*Time))
!
!        WRITE(1,*)THEORY(0),X(0)
!        do I = 1, 1000
!            THEORY(I) = sin(k*(X_theory(I)-c*Time))
!            WRITE(1,*)THEORY(i),X_theory(I)
!        end do
!    CLOSE(1)
!      CALL BoundValue()
      OPEN(1,FILE='Res1.dat')
        theory(0) = sin(k*(x(0)-c*Time))
        WRITE(1,*)THEORY(0),X(0)
        do I = 1, NX+1
            THEORY(I) = sin(k*(x(I)-c*Time))
            WRITE(1,*)THEORY(i),X(i)
        end do
      CLOSE(1)


      OPEN(IO,FILE='Res.dat')

      call Output(U,X,NX,t)
!-------------------------  Solve equation ------------------
    do
        t = t + dt
        if (t .gt. Time) then
            exit
        end if
        do I = 1,NX
            UN(I) = U(I)-CFL*(U(I)-U(I-1))
        end do
        call BoundValue(UN,NX)
        call Output(UN,X,NX,t)
        U = UN
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
