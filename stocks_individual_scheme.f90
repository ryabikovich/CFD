!program hello
!    implicit none
!
!      INTEGER, parameter:: IO = 12 ! input-output unit
!      INTEGER NX,NT,I,J,ID,m,ISCHEME
!      REAL,ALLOCATABLE :: U(:),UN(:),X(:), U_interval(:)
!      REAL L,h,VNM,dt,t,Time,time1,V1,in_time
!      REAL C0, C1, a, in_val, delta_1
!
!      WRITE(*,*) 'Read input file'
!      OPEN(IO,FILE='Input.txt')
!      READ(IO,*) L
!      READ(IO,*) NX
!      READ(IO,*) Time
!      READ(IO,*) VNM
!      READ(IO,*) a
!      READ(IO,*) C0, C1
!      READ(IO,*) ISCHEME
!      CLOSE(IO)
!
!      ALLOCATE(U(1:NX),UN(1:NX),X(1:NX), U_interval(1:NX))
!
!      h= L/(Nx-1)
!      dt=(h**2)/(vnm*2)
!      NT=Time/dt + 1
!
!
!      WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX, 'a=', a
!      WRITE(*,*) 'VNM=', VNM, 'dt=', dt, 'Time=', Time, 'NT=', NT
!
!      X(1)=0.0
!      DO I=2, NX-1
!       X(I)=X(I-1)+h
!      END DO
!      X(NX)=L
!
!      U(:)=0.0
!      UN(:)=0.0
!      U_interval(:)=0.0
!
!      in_val = 0.0
!      time1 = 0.0
!
!
!      CALL InitValue(U,NX, in_val)
!
!
!      OPEN(IO,FILE='Res.dat')
!      call Output(U,X,NX,time1)
!!-------------------------  Solve equation ------------------
!
!    do
!        time1 = time1 + dt
!        if(time1 .gt. time) then
!            exit
!        end if
!
!        do i=2,NX-1
!                U_interval(i)= U(i) + VNM*dt/(2*(h**2))*(U(i+1)-2*U(i)+U(i-1))
!        end do
!        call BoundValue(U_interval(1),U_interval(NX),C0,C1)
!        call BoundValue(UN(1),UN(NX),C0,C1)
!        UN(2) = VNM*dt/(12*h**2)*(-U_interval(4)+16*U_interval(3)-29*U_interval(2)+16*U_interval(1)-2)+U(2)
!        UN(NX-1) = VNM*dt/(12*h**2)*(16*U_interval(NX)-31*U_interval(NX-1)+16*U_interval(NX-2)-U_interval(NX-3))+U(NX-1)
!        do i=3,NX-2
!        UN(i)= U(i) + VNM*dt/(12*h**2)*(16*U_interval(i+1)-U_interval(i+2)-30*U_interval(i)+16*U_interval(i-1)-U_interval(i-2))
!        end do
!        call Output(UN,X,NX,time1)
!        U = UN
!
!
!    end do
!
!
!!---------------------------Output results
!
!      close(IO)
!
!      End program

!!----------------------- Set Initial Value -----------------
!      SUBROUTINE InitValue(U, NX, in_val)
!      IMPLICIT NONE
!      integer, intent(in):: NX
!      real,intent(in)::in_val
!       REAL,intent(in out):: U(NX)
!       U = in_val
!
!      END SUBROUTINE
!
!!----------------------- Set Boundary Condition ------------
!      SUBROUTINE BoundValue(a,b,C0,C1)
!      IMPLICIT NONE
!      real, intent(in)::C0,C1
!      real, intent(in out):: a, b
!      a = C1
!      b = C0
!
!      END SUBROUTINE
!
!!-----------------------
!      SUBROUTINE Output(U,X,NX,time)
!      IMPLICIT NONE
!      INTEGER, parameter:: IO=12
!      integer,intent(in)::NX
!      real,intent(in out)::time
!      integer i
!      REAL,intent(in out):: U(NX)
!      REAL,intent(in out):: X(NX)
!      do i=1,NX
!        write(IO,*) U(i),X(i),time
!      end do
!
!
!
!      END SUBROUTINE
      program hello
    implicit none

      INTEGER, parameter:: IO = 12 ! input-output unit
      INTEGER NX,NT,I,J,ID,m,ISCHEME
      REAL,ALLOCATABLE :: U(:),UN(:),X(:), U_interval(:)
      REAL L,h,VNM,dt,t,Time,time1,V1,in_time
      REAL C0, C1, a, in_val, delta_1

      WRITE(*,*) 'Read input file'
      OPEN(IO,FILE='Input.txt')
      READ(IO,*) L
      READ(IO,*) NX
      READ(IO,*) Time
      READ(IO,*) VNM
      READ(IO,*) a
      READ(IO,*) C0, C1
      READ(IO,*) ISCHEME
      CLOSE(IO)

      ALLOCATE(U(1:NX),UN(1:NX),X(1:NX), U_interval(1:NX))

      h= L/(Nx-1)
      !dt=(h**2)/(a*2)
      dt=(h**2)/(a*6)
      NT=Time/dt + 1
      VNM = a*dt/(h**2) !изменила число неймана


      WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX, 'a=', a
      WRITE(*,*) 'VNM=', VNM, 'dt=', dt, 'Time=', Time, 'NT=', NT

      X(1)=0.0
      DO I=2, NX-1
       X(I)=X(I-1)+h
      END DO
      X(NX)=L

      U(:)=0.0
      UN(:)=0.0
      U_interval(:)=0.0

      in_val = 0.0
      time1 = 0.0


      CALL InitValue(U,NX, in_val)


      OPEN(IO,FILE='Res1.dat')
      call Output(U,X,NX,time1)
!-------------------------  Solve equation ------------------

    do
        time1 = time1 + dt
        if(time1 .gt. time) then
            exit
        end if

        do i=2,NX-1
                U_interval(i)= U(i) + (VNM/2)*(U(i+1)-2*U(i)+U(i-1))
        end do
        call BoundValue(U_interval(1),U_interval(NX),C0,C1)
        call BoundValue(UN(1),UN(NX),C0,C1)
        UN(2) = (VNM/12)*(-U_interval(4)+16*U_interval(3)-29*U_interval(2)+16*U_interval(1)-2)+U(2)
        UN(NX-1) = (VNM/12)*(16*U_interval(NX)-31*U_interval(NX-1)+16*U_interval(NX-2)-U_interval(NX-3))+U(NX-1)
        do i=3,NX-2
        UN(i)= U(i) + (VNM/12)*(16*U_interval(i+1)-U_interval(i+2)-30*U_interval(i)+16*U_interval(i-1)-U_interval(i-2))
        end do
        call Output(UN,X,NX,time1)
        U = UN


    end do


!---------------------------Output results

      close(IO)

      End program

!----------------------- Set Initial Value -----------------
      SUBROUTINE InitValue(U, NX, in_val)
      IMPLICIT NONE
      integer, intent(in):: NX
      real,intent(in)::in_val
       REAL,intent(in out):: U(NX)
       U = in_val

      END SUBROUTINE

!----------------------- Set Boundary Condition ------------
      SUBROUTINE BoundValue(a,b,C0,C1)
      IMPLICIT NONE
      real, intent(in)::C0,C1
      real, intent(in out):: a, b
      a = C1
      b = C0

      END SUBROUTINE

!-----------------------
      SUBROUTINE Output(U,X,NX,time)
      IMPLICIT NONE
      INTEGER, parameter:: IO=12
      integer,intent(in)::NX
      real,intent(in out)::time
      integer i
      REAL,intent(in out):: U(NX)
      REAL,intent(in out):: X(NX)
      do i=1,NX
        write(IO,*) U(i),X(i),time
      end do



      END SUBROUTINE


