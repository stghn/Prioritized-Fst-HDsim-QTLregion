program MAF_S1S2pop
implicit none
integer::i,k,kk,id1,id2,id3,j,aa,idd
integer,parameter::n=1500,m=600000
real::bv,freq0(m),freq1(m),freq2(m),p(m),q(m)
integer::ff(m)
integer::sum0(m),sum1(m),sum2(m)

open(11,file="S1S2.gen9")

open(13,file="allele.freq.S1S2pop")

sum0=0
sum1=0
sum2=0

 do k=1,n
   read(11,*)idd,(ff(i),i=1,m)

     do i=1,m
   if (ff(i).eq. 0) sum0(i)=sum0(i)+1
   if (ff(i).eq. 1) sum1(i)=sum1(i)+1
   if (ff(i).eq. 2) sum2(i)=sum2(i)+1
     enddo
enddo


  do j=1,m
     kk=j
  aa=0
  aa=(sum0(kk)+sum1(kk)+sum2(kk))
  freq0(j)=float(sum0(kk))/aa
  freq1(j)=float(sum1(kk))/aa
  freq2(j)=float(sum2(kk))/aa

  p(j)=freq0(j)+(0.5*freq1(j))
   !print*,p(j)

  q(j)=1-p(j)
   !print*,q(j)

    ! if (min(p(j),q(j)).ge.0.05)then

  ! write(13,*)j,min(p(j),q(j)) !,min(p(j),q(j))
   write(13,*)j,max(p(j),q(j)),min(p(j),q(j)) !Allele Frequency

enddo

          end program MAF_S1S2pop

