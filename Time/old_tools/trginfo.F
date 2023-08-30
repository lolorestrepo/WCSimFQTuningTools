
      subroutine trginfo(trgofst)
*
*     read trigger information from DSTRGOFF bank

      implicit none

      integer dsevtrig_ibuf(4) ! Store primary trigger info
      real dsevtrig_rbuf(4)
      real trgofst
      equivalence (dsevtrig_ibuf(1), dsevtrig_rbuf(1))
      integer lbank,ndata

      call kzbloc('DSTRGOFF',lbank)
      
      if (lbank.le.0) then
        write(*,*) 'No DSTRGOFF bank!'
        stop
      endif
      
      call kzget1('DSTRGOFF', 1, ndata, dsevtrig_ibuf)

*      write(*,*) '   primary trigger type: ',dsevtrig_ibuf(1)
*      write(*,*) '   primary trigger it0sk: ',dsevtrig_ibuf(2)
*      write(*,*) '   primary trigger dgeantt0: ',dsevtrig_rbuf(3)
*      write(*,*) '   primary trigger number: ',dsevtrig_ibuf(4)
      
      trgofst=dsevtrig_rbuf(3)
      
      return
      end
