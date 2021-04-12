      subroutine elmt47(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)
c
c     lineares 4-Knoten Element / 4 GAUSS-Punkte
c     HOOKE-Material
c     
c
c             eta
c              ^
c        4     |     3
c        o-----|-----o
c        | x   |   x |
c        |     |     |
c        |     +-----|--->xi
c        |           |
c        | x       x |
c        o-----------o
c        1           2  
c
c
      implicit double precision (a-h,o-z)
c
c
      include   'cdata.h'     ! numnp,numel,nummat,nen,neq,ipr
      include   'eltran.h'
      include   'iofile.h'
      include   'prstrs.h'
      include   'comblk.h'
      include   'tdata.h'
      include   'ptdat6.h'
      include   'eldata.h'   ! element number n for testing
      include   'counts.h'   ! iteration number

      
      real*8 d(*),ul(ndf,nen,*),xl(ndm,*),tl(*),s(nst,*),p(*)
      integer ix(*)
      integer ndf,ndm,nst,isw   

      

      real*8      lambda,mue,rho,damp_ska,sig(3),eps(3),xx,yy,c
      real*8      cc(3,3),sigm(2),upp(2),up(2)
      real*8      shp(3,nen),sg(3,100),ss(2),xsj
      real*8      Bi1(3,2),Bj1(3,2),hilf(3,2)
      real*8      kmat(2*nen,2*nen)
      real*8      mmat(2*nen,2*nen)
      real*8      dmat(2*nen,2*nen)
      real*8      valp(9),pel
      integer     i,j,i1,j1,k1,l,lint,ngp

c.... go to correct array processor
c                 4       8
      go to(1,2,3,3,2,3,2,3,2,2,2,2,3,3,3), isw
      return
c.... input material properties

1     continue


      call dinput(d(1),4)
      lambda   = d(1)
      mue      = d(2)
      rho      = d(3)
      damp_ska = d(4)
 

      write(iow,2000)  lambda,mue,rho

      return

2     continue
      return

3     continue

      lambda   = d(1)
      mue      = d(2)
      rho      = d(3)
      damp_ska = d(4)
      

      do i=1,2*nen
        do j=1,2*nen
          kmat(i,j)=0.d0
        enddo
      enddo

      do i = 1,2*nen
         do j = 1,2*nen
          mmat(i,j)=0.d0
         enddo
      enddo 


      do i = 1,2*nen
         do j = 1,2*nen
          dmat(i,j)=0.d0
         enddo
      enddo 

      ngp = 2
      call int2d(ngp,lint,sg)

c...  GP-loop:
      do 320 l = 1,lint

      ss(1)=sg(1,l)
      ss(2)=sg(2,l)

      call shp2d(ss,xl,shp,xsj,ndm,nen,ix,.false.)

      xsj = xsj*sg(3,l)
      xx = 0.d0
      yy = 0.d0
      do i=1,3
        eps(i) =0.d0
        sig(i) =0.d0
      enddo
      
      do i =1,2
        up(i)  =0.d0
        upp(i) =0.d0
      enddo 

      do j=1,nen
c       GP position
        xx    = xx    + shp(3,j)*xl(1,j)
        yy    = yy    + shp(3,j)*xl(2,j)
c       GP strain
        eps(1) = eps(1) + shp(1,j)*ul(1,j,1)
        eps(2) = eps(2) + shp(2,j)*ul(2,j,1)
        eps(3) = eps(3) + shp(1,j)*ul(2,j,1) + shp(2,j)*ul(1,j,1)
c       GP  velocities,acceleration
        up(1)  = up(1) + shp(3,j)*ul(1,j,4)
        up(2)  = up(2) + shp(3,j)*ul(2,j,4)

        upp(1) = upp(1) + shp(3,j)*ul(1,j,5)
        upp(2) = upp(2) + shp(3,j)*ul(2,j,5)
      enddo



c     C-Matrix
      cc(1,1) = lambda+2.d0*mue
      cc(1,2) = lambda
      cc(1,3) = 0.d0
      cc(2,1) = lambda
      cc(2,2) = lambda+2.d0*mue
      cc(2,3) = 0.d0
      cc(3,1) = 0.d0
      cc(3,2) = 0.d0
      cc(3,3) = mue
      do k = 1,3
         do j = 1,3
           sig(k) = sig(k)+cc(k,j)*eps(j)
         enddo
      enddo
     

c... Berechnung der Hauptspannungen (nur zum plotten)
      sigm(1)=(sig(1)+sig(2))/2.d0
     x       +dsqrt(((sig(1)-sig(2))/2.d0)**2+sig(3)**2)
      sigm(2)=(sig(1)+sig(2))/2.d0        
     x       -dsqrt(((sig(1)-sig(2))/2.d0)**2+sig(3)**2)

c... Berechnung von psi_elastisch
      pel = 0.d0
      do j=1,3
        pel = pel + eps(j)*sig(j)
      enddo
      pel = pel/2.d0


      if(mod(isw,3).eq.0) then

c.... I-loop
      do 300 i=1,nen
c       B_I
        Bi1(1,1)=shp(1,i)
        Bi1(1,2)=0.d0
        Bi1(2,1)=0.d0
        Bi1(2,2)=shp(2,i)
        Bi1(3,1)=shp(2,i)
        Bi1(3,2)=shp(1,i)


        if(isw.eq.3) then 

c....   J-loop
        do 310 j=1,nen
          Bj1(1,1)=shp(1,j)
          Bj1(1,2)=0.d0
          Bj1(2,1)=0.d0
          Bj1(2,2)=shp(2,j)
          Bj1(3,1)=shp(2,j)
          Bj1(3,2)=shp(1,j)
  
c         Kuu - terme einsortieren
          
          do i1=1,3
           do j1=1,2
            hilf(i1,j1)=0.d0
            do k1=1,3
              hilf(i1,j1)=hilf(i1,j1)+cc(i1,k1)*Bj1(k1,j1)
            enddo
           enddo
          enddo

          do i1=1,2
           do j1=1,2
            do k1=1,3
              kmat(2*(i-1)+i1,2*(j-1)+j1)=kmat(2*(i-1)+i1,2*(j-1)+j1)
     x                   +Bi1(k1,i1)*hilf(k1,j1)*xsj
            enddo
           enddo
          enddo
         
          do i1=1,2
           do j1=1,2
             mmat(2*(i-1)+1,2*(j-1)+1)=mmat(2*(i-1)+1,2*(j-1)+1)
     x                   +shp(3,i)*rho*shp(3,j)*xsj
           
             mmat(2*(i-1)+2,2*(j-1)+2)=mmat(2*(i-1)+2,2*(j-1)+2)
     x                   +shp(3,i)*rho*shp(3,j)*xsj

           enddo
          enddo  

          do i1=1,2
           do j1=1,2
             dmat(2*(i-1)+1,2*(j-1)+1)=dmat(2*(i-1)+1,2*(j-1)+1)
     x                   +shp(3,i)*damp_ska*shp(3,j)*xsj
           
             dmat(2*(i-1)+2,2*(j-1)+2)=dmat(2*(i-1)+2,2*(j-1)+2)
     x                   +shp(3,i)*damp_ska*shp(3,j)*xsj

           enddo
          enddo  

             

310     continue
        
        endif ! isw=3

c     compute residual -(Bt sigma+ NI rho upp+NI damp_ska up )
          
        
        do i1=1,2
          do k1=1,3
            p(2*(i-1)+i1)=p(2*(i-1)+i1)-Bi1(k1,i1)*sig(k1)*xsj
          enddo
            p(2*(i-1)+i1)=p(2*(i-1)+i1)-
     x      (shp(3,i) *rho*upp(i1)*xsj+shp(3,i)*up(i1)*damp_ska*xsj) 

        enddo


300   continue

      endif ! mod(isw,3)
       


      if(isw.eq.4) then
        write(iow,2500) xx,yy,c,eps,sig
        write(*,2500) xx,yy,c,eps,sig
      endif
      
c     im plot-modus wird durch stre,i geplottet, was in valp(i) steht

      if (isw.eq.8) then        
          valp(1)=sig(1)         ! sig_11 
          valp(2)=sig(2)         ! sig_22
          valp(3)=sig(3)         ! sig_12 
          valp(4)=sigm(1)        ! sig_1
          valp(5)=sigm(2)        ! sig_2
          valp(6)=0.d0           ! 
          valp(7)=0.d0           ! 
          valp(8)=0.d0           ! 
          valp(9)=0.d0           ! 

          call elmt47plot(ix,valp,shp,xsj,
     x                    hr(nph),hr(nph+numnp),hr(ner),
     x                    nen,numnp)  

      endif

c     for energy output with tplo,ener
      if(isw.eq.13) then
c     write energy data to epl
        epl(1) = epl(1) + pel*xsj
      endif

320   continue ! l=1,lint (GP-loop)

c     form overall element matrix

c        IF (n.EQ.5) THEN
        
c        WRITE(iow,'(A)') 'mmatr'
c        do i = 1,8 
c         WRITE(iow,100) mmat(i,1),mmat(i,2),mmat(i,3),mmat(i,4),
c     x          mmat(i,5), mmat(i,6),mmat(i,7),mmat(i,8)
c        enddo

c        WRITE(iow,'(A)') 'kmatr'
c        do i = 1,8 
c         WRITE(iow,100) kmat(i,1),kmat(i,2),kmat(i,3),kmat(i,4),
c     x          kmat(i,5), kmat(i,6),kmat(i,7),kmat(i,8)
c        enddo

c       PRINT*,'ctan'
c       PRINT*,ctan(1),ctan(2),ctan(3)      
c       ENDIF 
 
      do i=1,2*nen
        do j=1,2*nen
          s(i,j)=ctan(1)*kmat(i,j)+ctan(2)*dmat(i,j)+ctan(3)*mmat(i,j)
        enddo
      enddo


      return

9     continue
        write(*,3000) isw
        write(iow,3000) isw
      return

c     Formatangaben
2000  format(' 4-node-element for isotropic material ',/,
     x       ' for cracks',//,
     x       ' lambda    = ',1e12.4,/,
     x       ' mue       = ',1e12.4,/,
     x       ' rho       = ',1e12.4)

100   format(T2,E9.2,TR2,E9.2,TR2,E9.2,TR2,E9.2,TR2,E9.2,TR2,E9.2
     x       TR2,E9.2,TR2,E9.2)

2500  format(' Element stresses (x_1,x_2,c): ',e12.4,e12.4,e12.4,/,
     x       ' eps_11, eps_22 ,2eps_12     : ',3e12.4,/,
     x       ' sig_11, sig_22 , sig_12     : ',3e12.4,//)
      

3000  format(' Warning: Element called with isw=',i4,' ! ',/,
     x       '          No action implemented ! ',/)

      end

c     7-----------------------------------------------------------------
      subroutine elmt47plot(ix,value,shp,xsj,dt,st,ser,nel,numnp)
 
      implicit none
 
c      integer nel,nen,numnp
c      integer i,l,ll
      integer nel,numnp
      integer i,ll
      real*8  xg

      integer ix(*)
      real*8  dt(numnp),st(numnp,*),ser(*)
      real*8  xsj,shp(3,4),value(9)
 
      integer k

      save
 
c     Lumped and consistent projection routine
 
c       Compute lumped projection and assemble stress integrals
 
        do i = 1,nel
          ll = ix(i)
          if(ll.gt.0) then
 
            xg     = shp(3,i)*xsj
            dt(ll) = dt(ll) + xg
 
c           Projections of GP-Data
 
            do k=1,9
              st(ll,k) = st(ll,k) + value(k)*xg
            enddo
            
          endif
        end do
 
      end
c     7-----------------------------------------------------------------
