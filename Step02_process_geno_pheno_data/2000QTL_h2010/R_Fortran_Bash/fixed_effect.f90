program aaaaaa
use kinds
implicit none
real(r8)::y,x1,yy,u
real(r8), allocatable::b1(:),b2(:)
integer::ia,j,k,nfix1,nfix2,i
!print*,'# of classes fix1 and fix2'
!read(*,*)nfix1,nfix2
nfix1=100
nfix2=4
allocate(b1(nfix1),b2(nfix2))

open(15,file='true_fixed')

x1=253723536.
      do i=1,nfix1
       call normal(x1,u)
       b1(i)=10+u*sqrt(5.0)
       write(15,15)i,b1(i)
       enddo
         do i=1,nfix2
       call normal(x1,u)
       b2(i)=2+u*sqrt(2.0)
       write(15,15)i,b2(i)
             enddo

open(1,file='pheno')

1  read(1,*,end=5)ia,y

call unif(x1,u)
j=1+u*nfix1
call unif(x1,u)
k=1+u*nfix2
yy=y+b1(j)+b2(k)
write(11,100)j,k,ia,yy
go to 1
5 continue
100 format(2i4,i8,f12.3)
15  format(i3,f12.3)
end
! ----------------------------------------------------------------------
subroutine normal(x1,z)
! generacion de un numero normal z -> n(0,1)
use kinds
implicit none

real(r8)  :: x1   ! semilla
real(r8)  :: z    ! valor N(0,1)
real(r8)  :: u1,u2

!   real*8 x1,z,u1,u2
    call unif(x1,u1)
    call unif(x1,u2)
    z=((-2.*log(u1))**0.5)*cos(2.*3.1416*u2)
    return
    end
! ----------------------------------------------------------------------
subroutine unif(x1,u)
!  generacion de un numero uniforme u[0,1]

use kinds
implicit none
real(r8)       :: x1,u ! semilla y U(0,1)
real(r8)       :: divis,trans,lsol,divid

      divis=2.**31.-1.
      trans=7.**5.
      divid=trans*x1
      lsol=int(divid/divis)
      x1=divid-lsol*divis
      u=x1/divis
      return
      end
