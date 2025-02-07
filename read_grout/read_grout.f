      program read_grout

c this program will pick out information at certain depths from the ground.out
c file for comparison purposes

      implicit none

      integer mnodes
	parameter (mnodes = 1000)   !maximum number of nodes

      integer i,k,year(mnodes),day(mnodes),hour(mnodes),minute(mnodes)
	integer node(mnodes),tot_nodes,freq_id,lflag,io
	real*8 tdepth,depth(mnodes),temp(mnodes),moist(mnodes),ice(mnodes)
	real*8 num1,num2,ntemp,nmoist,nice,water(mnodes),nwater
	real*8 vapor(mnodes),nvapor,denom
	character*120 fname1,header,fname2
	character in_data*101,stype(mnodes)*2,nstype*2
	character*1 state(mnodes),nstate,answer


      write(*,'('' Enter the name of the input file ''\)')
	read(*,'(a)') fname1
	open(unit=5,file=fname1,status='old')

	write(*,*)' '
	write(*,'('' You can enter as many depths as you like. You will 
     &be asked for them '',/, ''one at a time. Each depth has its own 
     &output file. ''\)')
	write(*,*)' '

      answer = 'y'
      do while(answer.eq.'Y'.or.answer.eq.'y')
 	  write(*,*) ' '
	  write(*,'('' Enter the target depth in meters ''\)')
        read(*,*) tdepth

c open the output file
        write(*,*)' '
        write(*,'('' Enter the name of the output file ''\)')
	  read(*,'(a)') fname2
	  open(unit=6,file=fname2,status='unknown')

c read in the data
        read(5,'(a)') in_data
	  read(in_data(25:28),'(i4)') tot_nodes
	  read(in_data(40:49),'(i10)') freq_id
	  read(5,'(a)') header
	  read(5,'(a)') header

	  write(6,'(''freq_id: '',i10)') freq_id
        write(6,'('' Year   JD   Hr    M  USCS    Depth   Grtemp    Wate
     &r         Ice       W + I       Vapor     F/T'')')
        write(6,'(''                                 m      K'')')

        io = 0
        do while(io.ne.-1)
          do i=1,tot_nodes
            read(5,'(a)',iostat=io) in_data
            read(in_data(2:5),'(i4)') year(i)
	      read(in_data(8:10),'(i3)') day(i)
	      read(in_data(14:15),'(i2)') hour(i)
	      read(in_data(19:20),'(i2)') minute(i)
	      read(in_data(23:25),'(i3)') node(i)
	      read(in_data(29:30),'(a2)') stype(i)
            read(in_data(35:40),*) depth(i)
	      read(in_data(43:49),*) temp(i)
	      read(in_data(52:61),*) water(i)
	      read(in_data(64:73),*) ice(i)
	      read(in_data(76:85),*) moist(i)
	      read(in_data(88:97),*) vapor(i)
	      read(in_data(101:101),'(a)') state(i)
	    end do

          lflag = 0
          do while(lflag.eq.0)
            do k=2,tot_nodes
	        if((depth(k-1).le.tdepth).and.(depth(k).ge.tdepth)) then
                lflag = 1
	          if(depth(k).eq.tdepth) then
	            nstype = stype(k)
	            nstate = state(k)
	            ntemp = temp(k)
	            nwater = water(k)
	            nice = ice(k)
	            nmoist = moist(k)
	            nvapor = vapor(k)
	          else
	            nstype = stype(k-1)
	            nstate = state(k-1)
	            num1 = depth(k) - tdepth  !depth(k-1)
	            num2 = tdepth - depth(k-1)
	            denom = num1 + num2
                  ntemp = (temp(k-1)*num1 + temp(k)*num2)/denom
	            nwater = (water(k-1)*num1 + water(k)*num2)/denom
	            nice = (ice(k-1)*num1 + ice(k)*num2)/denom
	            nmoist = (moist(k-1)*num1 + moist(k)*num2)/denom
	            nvapor = (vapor(k-1)*num1 + vapor(k)*num2)/denom
	          end if
	          write(6,1001) year(k),day(k),hour(k),minute(k),nstype,
     &                        tdepth,ntemp,nwater,nice,nmoist,nvapor,
     &                        nstate
                if(lflag.eq.1) exit
	        end if
              if(k.eq.tot_nodes.and.lflag.eq.0) lflag = 1
	      end do
          end do
        end do

c Are there any more target depths
        write(*,*)' '
        write(*,'('' Do you want to pick any more target depths? ''\)')
	  read(*,'(a)') answer
	  if(answer.eq.'y'.or.answer.eq.'Y') then
	    rewind(5)
	    close(6)
	  end if
	end do     

      close(5)

 1001 format(4(i5),3x,a2,1x,2(f9.3),4(e12.4),3x,a1)

 	stop
	end