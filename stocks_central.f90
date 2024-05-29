
program hello
    implicit none

      INTEGER, parameter:: IO = 12 ! input-output unit
      INTEGER NX,NT,I,J,ID,m,ISCHEME
      REAL,ALLOCATABLE :: U(:),UN(:),X(:)
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

      ALLOCATE(U(1:NX),UN(1:NX),X(1:NX))

      h= L/(Nx-1)
!      dt=(h**2)/(a*2)
      dt=(h**2)/(a*6) !изменила шаг
      NT=Time/dt + 1
      VNM = a*dt/(h**2) !изменила число фон неймана


      WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX, 'a=', a
      WRITE(*,*) 'VNM=', VNM, 'dt=', dt, 'Time=', Time, 'NT=', NT

      X(1)=0.0
      DO I=2, NX-1
       X(I)=X(I-1)+h
      END DO
      X(NX)=L

      U(:)=0.0
      UN(:)=0.0

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
                UN(i)= U(i) + VNM*(U(i+1)-2*U(i)+U(i-1))
        end do
        call BoundValue(UN(1),UN(NX),C0,C1)
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



