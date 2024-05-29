
Program Pr
      Implicit none

      INTEGER, parameter:: IO = 12 ! input-output unit
      INTEGER NX,NT,I,J,ID,m, ISCHEME
      REAL(8),ALLOCATABLE :: U(:),UN(:),X(:), U_add(:)
      REAL(8):: X_theory(101),theory(101)
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

      ALLOCATE(U(NX),UN(NX),X(NX),U_add(NX))

      pi=3.14159265359d0
!      k = m*pi/L
      h = L/(NX-1)
      dt = CFL*h/c
      NT = 0.5/(c*dt)

      Time = dt*NT

      WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX
      WRITE(*,*) 'CFL=', CFL, 'dt=', dt, 'Time=', Time, 'NT=', NT

      DO I=1, NX
        X(I)=(I-1)*h
      END DO


      U(:)=0.0
      UN(:)=0.0

      t=0.0d0
      CALL InitValue(X,U,NX)


    OPEN(1,FILE='Res1.dat')
        p = L/100
        write(*,*)'p=',p
         DO I=1, 101
            X_theory(I)= (I-1)*p
         END DO

        do I = 1, 101
                if ( (X_theory(I)-c*Time) .lt. 0.5) then
                    THEORY(i) = 1
                else
                    THEORY(i) = 3
                end if
            WRITE(1,*)THEORY(i),X_theory(I)
        end do
    CLOSE(1)
!      CALL BoundValue()
!      OPEN(1,FILE='Res1.dat')
!        theory(0) = sin(k*(x(0)-c*Time))
!        WRITE(1,*)THEORY(0),X(0)
!        do I = 1, NX+1
!            THEORY(I) = sin(k*(x(I)-c*Time))
!            WRITE(1,*)THEORY(i),X(i)
!        end do
!      CLOSE(1)


      OPEN(IO,FILE='Res.dat')

      call Output(U,X,NX,t)
!-------------------------  Solve equation ------------------
    do
        t = t + dt
        if (t .gt. (Time+0.0000001)) then
            exit
        end if
        do I = 2,NX-1
            UN(I) = U(I)-CFL*(U(I)-U(I-1))
        end do
        call BoundValue(UN,NX)
        call Output(UN,X,NX,t)
        U = UN
!        do I = 2,NX-1
!            U_add(I) = U(I)-CFL*(U(I+1)-U(I))
!        end do
!        call BoundValue(U_add,NX)
!        do I = 2,NX-1
!            UN(I) = 0.5*(U(I)+U_add(I)-CFL*(U_add(I)-U_add(I-1)))
!        end do
!        call BoundValue(UN,NX)
!
!        call Output(UN,X,NX,t) ! разобраться  с этим почему не Un
!        U = UN
    end do

!-------------------------  Output results ------------------

!      OPEN(IO,FILE='Res.dat')
!      call Output()
      CLOSE(IO)


      End program

!----------------------- Set Initial Value -----------------
      subroutine InitValue(X,U,NX)
      IMPLICIT NONE
      integer, intent(in):: NX
      integer i
       REAL*8,intent(in out):: U(NX)
       REAL*8,intent(in out):: X(NX)

      do I = 1, NX
        if (X(I) .lt. 0.5) then
            U(I) = 1
        else
            U(I) = 3
        end if
      end do
      end subroutine

!----------------------- Set Boundary Condition ------------
      subroutine BoundValue(UN,NX)
      IMPLICIT NONE
      integer, intent(in):: NX
      REAL*8,intent(in out):: UN(NX)
      UN(1) = 1
      UN(NX)=3
      end subroutine

!----------------------- Output results ------------
      subroutine Output(U,X,NX,t)
      IMPLICIT NONE
     INTEGER, parameter:: IO=12
      integer,intent(in)::NX
      real*8,intent(in out)::t
      integer i
      REAL*8,intent(in out):: U(NX)
      REAL*8,intent(in out):: X(NX)
      do i=1,NX
        write(IO,*) U(i),X(i),t
      end do


      end subroutine

