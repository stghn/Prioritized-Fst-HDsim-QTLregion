program ReplaceSNPcode
       implicit none
       integer::igen(600000),igenn(600000),i,k,m,j

       open(1,file='new.geno')
       open(11,file='geno')


      do m=1,30000

     read(11,100)(igen(k),k=1,600000)

       do i=1,600000



           if(igen(i).eq.0)igenn(i)=0
           if(igen(i).eq.2)igenn(i)=2
           if(igen(i).eq.3)igenn(i)=1
           if(igen(i).eq.4)igenn(i)=1
       enddo
       write(1,200)(igenn(k),k=1,600000)
   rewind(12)
    enddo

100     format(600000(1x,i1))
200     format(600000(1x,i2))

         end

