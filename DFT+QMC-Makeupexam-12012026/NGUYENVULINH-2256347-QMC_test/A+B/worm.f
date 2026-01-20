      program main
      implicit none
      include 'worm.h'
      integer iblk,istp,imoves
      call input
      call restart(1)
      call misure(1)
      do iblk=1,nblk
       call medie(1,iblk)
       do istp=1,nstp
        do imoves=1,nmoves
         call move(1)
        enddo
        call misure(2)
        call medie(2,iblk)
       enddo
       call medie(3,iblk)
       call misure(3)
       call restart(-1)
       call reset(iblk)
      enddo
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine add_misura(j,w,a)
      implicit none
      include 'worm.h'
      integer j
      real*8 w,a
      prop(j)=prop(j)+w*a
      wt(j)=wt(j)+w
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine add_tablename(word,nta,mta,tname)
      implicit none
      integer i,nta,mta
      character*48 word,tname(mta)
      do i=1,nta
       if(word.eq.tname(i))return
      enddo
      nta=nta+1
      tname(nta)=word
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine add_v(ndim,xp,x,el,eli,drt,drti,mgrid,ngrid,t,v,dv,ind)
      implicit none
      integer ndim,mgrid,ngrid,i,ind
      real*8 xp(ndim),x(ndim),el(ndim),eli(ndim),drt,drti,t(0:mgrid,4)
      real*8 r,rvec,rem,v,dv
      r=0.d0
      do i=1,ndim
       rvec=xp(i)-x(i)
       rvec=rvec-el(i)*nint(rvec*eli(i))
       r=r+rvec**2
      enddo
      r=sqrt(r)
      ind=int(r*drti)
      rem=(r-ind*drt)*drti
      ind=min(ngrid,ind)
      v=v+(t(ind,1)+rem*(t(ind,2)+rem*(t(ind,3)+rem* t(ind,4))))
      dv=dv+(t(ind,2)+rem*(2*t(ind,3)+rem*(3*t(ind,4))))*r*drti
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine advance
      implicit none
      include'worm.h'
      logical cemasha
      integer emme,it,itau,itaup,ip,ipp,j,irandom_number,next_emme_tau
      real*8 rnd,xnew(mdim,0:m),v,dv,p,delta_azione,quantov(0:m)
      if(zsector)return
      emme=irandom_number(emme_bar(iadvance))
      call get_ira(ip,itau,it)
      if(c_n.ne.0.and.np(next_emme_tau(emme,itau)).gt.c_n)return
      if(boltzmann)then
       j=ntau+1-emme
       if(.not.cemasha(j,ip,ipp,itau,itaup))return
      endif
      call diffuse(mdim,ndim,emme,x(1,ip,itau),xnew,sigma(it))
      call calcola_quantov(emme,quantov)
      delta_azione=0
      j=0
      call potential(xnew(1,j),v,dv,it,itau,ip_ira,quantov(j),.false.)
      delta_azione=delta_azione-tau*v
      itau=next_tau(itau)
      do j=1,emme
       call potential(xnew(1,j),v,dv,it,itau,0,quantov(j),.false.)
       call masha_correzione(xnew(1,j),v,it,itau)
       delta_azione=delta_azione-tau*v
       itau=next_tau(itau)
      enddo
      p=exp(delta_azione+emme*mu(it)*tau)
      att(iadvance)=att(iadvance)+1
      call random_number(rnd)
      if(p.lt.rnd)return
      acc(iadvance)=acc(iadvance)+1
      call get_ira(ip,itau,it)
      call metti(emme,ip,itau,it,xnew)
      call set_ira(ip,itau,it)
      if(verbose)write(*,*)'advance np(itau_ira) = ',np(itau)
      if(np(itau).gt.mnp)stop'np.gt.mnp in advance'
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wiggle
      implicit none
      include 'worm.h'
      integer itau0,itau1,itau,ip0,ip1,ip,it,i,j,emme,irandom_number
      real*8 rnd,xnew(mdim,0:m),v,dv,delta_azione,p
      if(gsector)return
      itau0=irandom_number(ntau)
      if(np(itau0).eq.0)return
      j=irandom_number(np(itau0))
      ip0=jp(j,itau0)
      emme=emme_bar(iwiggle) ! irandom_number(emme_bar(iwiggle))
      call next_emme_p_tau(emme,ip0,ip1,itau0,itau1)
      it=itype(ip0,itau0)
      call bridge(x(1,ip0,itau0),x(1,ip1,itau1),emme,sigma(it)
     &           ,el,eli,ndim,mdim,xnew)
      ip=ip0
      itau=itau0
      delta_azione=0
      do j=1,emme-1
       ip=next(ip,itau)
       itau=next_tau(itau)
       call potential(x(1,ip,itau),v,dv,it,itau,ip,1.d0,.false.)
       delta_azione=delta_azione+tau*v
       call potential(xnew(1,j),v,dv,it,itau,ip,1.d0,.false.)
       delta_azione=delta_azione-tau*v
      enddo
      p=exp(delta_azione)
      att(iwiggle)=att(iwiggle)+1
      call random_number(rnd)
      if(p.lt.rnd)return
      acc(iwiggle)=acc(iwiggle)+1
      ip=ip0
      itau=itau0
      do j=1,emme-1
       ip=next(ip,itau)
       itau=next_tau(itau)
       do i=1,ndim
        x(i,ip,itau)=xnew(i,j)
       enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bridge(x0,x1,m,sigma,el,eli,ndim,mdim,xnew)
      implicit none
      integer ndim,mdim,m,l1,l2,l3,i,j
      real*8 x0(ndim),x1(ndim),el(ndim),eli(ndim),xnew(mdim,0:m),sigma
      real*8 pi,d,s,xi
      data pi/3.14159265358979d0/
      l3=m
      do i=1,ndim
       xnew(i,0)=x0(i)
       xnew(i,l3)=x1(i)
      enddo
      do j=1,m-1
       l1=j-1
       l2=j
       s=sigma*(dble(l3-l2)/dble(l3-l1))**0.5d0
       do i=1,ndim
        d=xnew(i,l3)-xnew(i,l1)
        d=d-el(i)*nint(d*eli(i))
        xnew(i,j)=xnew(i,l1)+d/dble(l3-l1)+xi(s)
       enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine calcola_quantov(m,q)
      implicit none
      integer i,m
      real*8 q(0:m)
      q(0)=0.5d0
      do i=1,m-1
       q(i)=1.d0
      enddo
      q(m)=0.5d0
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function ceira(j,ip,itau)
      implicit none
      include 'worm.h'
      logical ceira
      integer j,ip,itau,emme
      ceira=.true.
      emme=j
      do j=0,emme-1
       if(ip.eq.ip_ira.and.itau.eq.itau_ira)return
       ip=next(ip,itau)
       itau=next_tau(itau)
      enddo
      if(ip.eq.ip_ira.and.itau.eq.itau_ira)return
      ceira=.false.
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function cemasha(j,ip,ipp,itau,itaup)
      implicit none
      include'worm.h'
      logical cemasha
      integer j,ip,ipp,itau,itaup,emme
      cemasha=.true.
      emme=j
      itaup=itau
      ipp=ip
      do j=0,emme-1
       if(ipp.eq.ip_masha.and.itaup.eq.itau_masha)return
       ipp=prev(ipp,itaup)
       itaup=prev_tau(itaup)
      enddo
      if(ipp.eq.ip_masha.and.itaup.eq.itau_masha)return
      cemasha=.false.
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine close
      implicit none
      include 'worm.h'
      integer emme,ip,itau,it,j
      real*8 v,dv,delta_azione,p,rnd,twosigma2,r2,rho
      real*8 xnew(mdim,0:m)
      if(zsector)return
      if(c_n.ne.0.and.np(itau_ira).ne.c_n)return
      emme=itau_masha-itau_ira
      if(emme.lt.0)emme=emme+ntau
      if(emme.eq.0)return
      if(emme.gt.emme_bar(iclose))return
      twosigma2=2*emme*sigma(it_ira)**2
      call gaussian(ndim,x(1,ip_ira,itau_ira),x(1,ip_masha,itau_masha)
     &             ,el,eli,twosigma2,r2,rho)
      if(r2.gt.4)return ! ekko
      call bridge(x(1,ip_ira,itau_ira),x(1,ip_masha,itau_masha)
     &           ,emme,sigma(it_ira),el,eli,ndim,mdim,xnew)
      delta_azione=0
      j=0
      call potential(xnew(1,j),v,dv,it_ira,itau_ira,ip_ira,.5d0,.false.)
      delta_azione=delta_azione-tau*v
      itau=next_tau(itau_ira)
      do j=1,emme-1
       call potential(xnew(1,j),v,dv,it_ira,itau,0,1.0d0,.false.)
       delta_azione=delta_azione-tau*v
       itau=next_tau(itau)
      enddo
      j=emme
      call potential(xnew(1,j),v,dv,it_ira,itau,ip_masha,0.5d0,.false.)
      delta_azione=delta_azione-tau*v
      p=rho*exp(delta_azione+emme*mu(it_ira)*tau)
     &        /(c*emme_bar(iclose)*ntau*np(itau_ira))
      att(iclose)=att(iclose)+1
      call random_number(rnd)
      if(p.lt.rnd)return
      acc(iclose)=acc(iclose)+1
      call get_ira(ip,itau,it)
      call metti(emme-1,ip,itau,it,xnew)
      next(ip,itau)=ip_masha
      prev(ip_masha,itau_masha)=ip
      gsector=.false.
      zsector=.true.
      ip_ira=0
      ip_masha=0
      itau_ira=0
      itau_masha=0
      it_ira=0
      it_masha=0
      if(verbose)write(*,*)'close np(1) = ',np(1)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine diffuse(m,n,emme,x0,x,sigma)
      implicit none
      integer m,n,i,emme,j
      real*8 x0(n),x(m,0:emme),sigma,xi
      do i=1,n
       x(i,0)=x0(i)
      enddo
      do j=1,emme
       do i=1,n
        x(i,j)=x(i,j-1)+xi(sigma)
       enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gaussian(n,x1,x2,el,eli,twosigma2,r2,rho)
      implicit none
      integer i,n
      real*8 x1(n),x2(n),el(n),eli(n),twosigma2,rho,r2,d,pi
      data pi/3.14159265358979d0/
      r2=0
      do i=1,n
       d=x1(i)-x2(i)
       d=d-el(i)*nint(d*eli(i))
       r2=r2+d*d
      enddo
      r2=r2/twosigma2
      rho=((twosigma2*pi)**(-0.5*n))*exp(-r2)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_ira(ip,itau,it)
      implicit none
      include'worm.h'
      integer ip,itau,it
      itau=itau_ira
      ip=ip_ira
      it=it_ira
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_masha(ip,itau,it)
      implicit none
      include'worm.h'
      integer ip,itau,it
      itau=itau_masha
      ip=ip_masha
      it=it_masha
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dmix(r0,drt,i0,nt,mt,f,vshift,vtail)
c 1/r3 potential
      implicit none
      integer i,i0,nt,mt,j
      real*8 r0,drt,f(0:mt,4),r,rcut,vshift,vtail
      call r_set(4*(mt+1),f(0,1),0.d0)
      vshift=1.d0/(nt*drt)**3
      rcut=0.27d0
      do i=i0+1,nt
       r=r0+drt*i
       f(i,1)=1.d0/r**3-vshift
       f(i,2)=-3.d0/r**4
      enddo
      j=int(rcut/drt)
      do i=i0,j
       f(i,1)=1.d0/rcut**3-vshift
       f(i,2)=0.d0
       f(i,3)=0.d0
       f(i,4)=0.d0
      enddo
      call spline(mt,nt,drt,f)
      vtail=0.5d0*2*acos(-1.d0)/(nt*drt)
      write(*,*)'rcut = ',nt*drt
      write(*,*)' 1/r3 vshift = ',vshift,' vtail (2d) = rho * ',vtail
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine hspotential(r0,drt,i0,nt,mt,f,vshift,vtail)   ! RIC
c hard-spheres potential discontinuous
      implicit none
      integer i,i0,nt,mt
      real*8 r0,drt,f(0:mt,4),r,sigma,epsr,vshift,vtail
      parameter(epsr=1.d20,sigma=1.d0)
      call r_set(4*(mt+1),f(0,1),0.d0)
      do i=i0,nt
       r=r0+drt*i
       if(r.lt.sigma)then
        f(i,1)=epsr
        f(i,2)=0.d0
        f(i,3)=0.d0
        f(i,4)=0.d0
       else
        f(i,1)=0.d0
        f(i,2)=0.d0
        f(i,3)=0.d0
        f(i,4)=0.d0
       endif
      enddo
      vshift=0.d0
      vtail=0.d0
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dipolar_smoothed(r0,drt,i0,nt,mt,f,vshift,vtail)
c 1/r3 potential
      implicit none
      integer i,i0,nt,mt
      real*8 r0,drt,f(0:mt,4),r,lmr,vshift,vtail
      call r_set(4*(mt+1),f(0,1),0.d0)
      do i=i0+1,nt
       r=r0+drt*i
       lmr=drt*nt*2-r
       f(i,1)=1.d0/r**3+1.d0/lmr**3
       f(i,2)=-3.d0/r**4+3.d0/lmr**4
      enddo
      f(i0,1)=f(i0+1,1)
      f(i0,2)=0.d0
      do i=i0,nt
       f(i,1)=f(i,1)-f(nt,1)
      enddo
      call spline(mt,nt,drt,f)
      f(i0,3)=0.d0
      f(i0,4)=0.d0
      vshift=0.d0
      vtail=0.d0
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dipolar(r0,drt,i0,nt,mt,f,vshift,vtail)
c 1/r3 potential
      implicit none
      integer i,i0,nt,mt
      real*8 r0,drt,f(0:mt,4),r,s,vshift,vtail
      call r_set(4*(mt+1),f(0,1),0.d0)
      vshift=1.d0/(nt*drt)**3
      do i=i0+1,nt-1
       r=r0+drt*i
       f(i,1)=1.d0/r**3-vshift
       f(i,2)=-3.d0/r**4
      enddo
      f(i0,1)=f(i0+1,1)
      f(i0,2)=0.d0
      call spline(mt,nt,drt,f)
      f(i0,3)=0.d0
      f(i0,4)=0.d0
      vtail=0.5d0*2*acos(-1.d0)/r
      write(*,*)'rcut = ',nt*drt
      write(*,*)' 1/r3 vshift = ',vshift,' vtail (2d) = rho * ',vtail
      write(*,*)'vtail out 2d = rho * ',s
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine zero(r0,drt,i0,nt,mt,f,vshift,vtail)
c non-interacting
      implicit none
      integer i,i0,nt,mt
      real*8 r0,drt,f(0:mt,4),r,v,dv,vshift,vtail
      call r_set(4*(mt+1),f(0,1),0.d0)
      vshift=0.d0
      vtail=0.d0
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine hfdhe2(r0,drt,i0,nt,mt,f,vshift,vtail)
c aziz potential
      implicit none
      integer i,i0,nt,mt
      real*8 r0,drt,f(0:mt,4),r,v,dv,vshift,vtail
      external v_hfdhe2
      call v_hfdhe2(nt*drt,vshift,dv)
      call r_set(4*(mt+1),f(0,1),0.d0)
      do i=i0,nt-1
       r=r0+drt*i
       r=max(r,drt*0.01d0)
       call v_hfdhe2(r,f(i,1),f(i,2))
       f(i,1)=f(i,1)-vshift
      enddo
      call spline(mt,nt,drt,f)
      call how_many_dim(i)
      call tail_correction(r,drt,i,v_hfdhe2,vtail)
      write(*,*)'rcut = ',nt*drt
      write(*,*)'hfdhe2 vshift = ',vshift,' vtail = rho * ',vtail
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lj(r0,drt,i0,nt,mt,f,vshift,vtail,p)
c aziz potential
      implicit none
      integer i,i0,nt,mt
      real*8 r0,drt,f(0:mt,4),r,v,dv,vshift,vtail,p(*),epsilon,sigma
      real*8 f0,df0
      common /c_lj/epsilon,sigma
      external v_lj
      epsilon=p(1)
      sigma=p(2)
      call v_lj(nt*drt,vshift,dv)
      call r_set(4*(mt+1),f(0,1),0.d0)
      call v_lj(0.5d0*sigma,f0,df0)
c f(x)=4*epsilon*((sigma/x)**12-(sigma/x)**6)
c df(x)=4*epsilon*(-12*(sigma/x)**12+6*(sigma/x)**6)/x
c f0 = f(sigma/2)
c df0 = df(sigma/2)
c g(x)=f0-df0*(sigma/2-x)
      do i=i0,nt-1
       r=r0+drt*i
c      r=max(r,drt*0.01d0)
       if(r.le.0.5d0*sigma)then
        f(i,1)=f0-df0*(0.5d0*sigma-r)
        f(i,2)=df0
       else
        call v_lj(r,f(i,1),f(i,2))
       endif
       f(i,1)=f(i,1)-vshift
      enddo
      call spline(mt,nt,drt,f)
      call how_many_dim(i)
      call tail_correction(r,drt,i,v_lj,vtail)
      write(*,*)'rcut = ',nt*drt
      write(*,*)'lj vshift = ',vshift,' vtail = rho * ',vtail
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine how_many_dim(i)
      implicit none
      include 'worm.h'
      integer i
      i=ndim
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tail_correction(r,drt,ndim,f,tail)
      implicit none
      integer i,ndim
      real*8 r,drt,tail,v,dv,s
c     external f
      r=r+0.5d0*drt
      s=0.d0
      do i=1,10000
       call f(r,v,dv)
       s=s+v*r**(ndim-1)
       r=r+drt
      enddo
      if(ndim.eq.1)tail=0.5d0                 *drt*s
      if(ndim.eq.2)tail=0.5d0*2*acos(-1.d0)   *drt*s
      if(ndim.eq.3)tail=0.5d0*4*acos(-1.d0)**2*drt*s
      write(*,*)'rout = ',r,' v = ',v
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine v_hfdhe2(r,f,df)
c aziz potential
      implicit none
      real*8 f,df,r,x,xi,a0,a1,a2,b0,b1,b2
     &      ,aux0,aux1,aux2,epsilon,a,alpha,c6,c8,c10,rmi,d
      data epsilon/10.8d0/
      data a/.5448504d6/
      data alpha/13.353384d0/
      data c6/1.3732412d0/
      data c8/0.4253785d0/
      data c10/0.1781d0/
      data rmi/0.337006706d0/
      data d/1.241324d0/
      x=r*rmi
      xi=1.d0/x
c parte attrattiva
      f=a*exp(-alpha*x)
      df=-alpha*rmi*f
c parte repulsiva
      aux1=xi**6
      aux2=xi*xi
      b0=-aux1*(c6+aux2*(c8+aux2*c10))
      aux1=-aux1*xi*rmi
      b1=-aux1*(c6*6.d0+aux2*(c8*8.d0+aux2*c10*10.d0))
      aux1=-aux1*xi*rmi
      b2=-aux1*(c6*42.d0+aux2*(c8*72.d0+aux2*c10*110))
      if(x.ge.d)then
       f=epsilon*(f+b0)
       df=epsilon*(df+b1)
      else
       aux0=rmi*xi*xi
       aux1=aux0*(d*d  *xi    *2.d0  -d           *2.d0 )
       aux0=aux0*rmi*xi
       aux2=aux0*(d           *4.d0  -d*d  *xi    *6.d0 )
       aux0=aux0*rmi*xi
       a0=exp(-(d*xi-1)*(d*xi-1))
       a1=aux1*a0
       a2=aux2*a0+aux1*a1
       f=epsilon*(f+a0*b0)
       df=epsilon*(df+a0*b1+a1*b0)
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine v_lj(r,f,df)
c LJ potential
      implicit none
      real*8 r,f,df,epsilon,sigma,rm1,f6,f12
      common /c_lj/epsilon,sigma
c parte attrattiva
      if(r.eq.0.d0)r=sigma/1000.d0
      rm1=1.d0/r
      f6=(sigma*rm1)**6
      f12=f6**2
      f=4*epsilon*(f12-f6)
      df=4*epsilon*(-12*f12+6*f6)*rm1
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine reset(iblk)
      implicit none
      include 'worm.h'
      integer i,k,l,iblk
      real*8 value,c0
      character*80 string
      character*48 word(20)

      flush(6)
      open(2,file='reset')
      do k=1,1000
       call readwords(2,20,word,i,string,0)
       if(i.eq.1)go to 1
       if(word(1).eq.'mu')then
        do l=1,mty
         if(word(2).eq.typename(l))then
          read(word(3),*)value
          if(abs(value-mu(l)).gt.1.d-10)then
           write(*,*)'reset mu '//typename(l),mu(l),' ---> ',value
           nblkeq=iblk
           mu(l)=value
          endif
         endif
        enddo
       elseif(word(1).eq.'c0')then
        read(word(2),*)value
        if(abs(value-c0).gt.1.d-10)then
         write(*,*)'reset c0 ',c0,' ---> ',value
         nblkeq=iblk
         c0=value
         c=c0/(vol*ntau)
        endif
       elseif(word(1).eq.'stop')then
        rewind(2)
        write(2,*)
        close(2)
        write(*,*)'stop'
        stop
       elseif(word(1).eq.'reset')then
        nblkeq=iblk
c       call nofr(1)
        call sofk(1)
        write(*,*)'reset cumulative averages'
       elseif(word(1).eq.'temperature')then
        read(word(2),*)value
        value=1.d0/value
        if(abs(value-beta).gt.1.d-10)then
         write(*,*)'reset T ',1.0/beta,' ---> ',1.d0/value
         nblkeq=iblk
         beta=value
         tau=beta/ntau
         do i=1,ntypes
          sigma(i)=sqrt(2*lambda(i)*tau)
         enddo
        endif
       endif
      enddo
    1 rewind(2)
      write(2,*)
      close(2)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine input
      implicit none
      include 'worm.h'
      integer i,k,l,seed,itau,nta,it,jt
      real*8 c0,kt,p(10,mta)
      character*80 string
      character*48 word(20),tname(mta),routinename(mta)
c initialize
      boltzmann=.false.
      c_n=0
      verbose=.false.
      irestart=0
      ntypes=0
      ngrid=1000
      ngridgratio=5
      rcut=0.5d0
      c0=1.d0
      m_bar=16
      nta=0
      iunit_nsp=70
      ncycles=0
c     nofr_stride=0
      nsofk=0
      call r_set(mty,mu,0.d0)
      call r_set(mdim,el,1.d0)
      call r_set(mdim,eli,1.d0)
      call i_set(nmove,emme_bar,m_bar)
      call r_set(nmove,pmove,1.d0)
      call move(0)
      pmove(idisplace)=0.02d0
c     pmove(iswty)=0.d0

      open(2,file='worm.in')
      do k=1,1000
       call readwords(2,20,word,i,string,1)
       if(i.eq.1)go to 1
       if(word(1).eq.'ndim')then
        read(word(2),*)ndim
       elseif(word(1).eq.'type')then
        ntypes=ntypes+1
        if(ntypes.gt.mty)stop'ntypes.gt.mty in input'
        read(word(2),'(a)')typename(ntypes)
        read(word(3),*)lambda(ntypes)       
       elseif(word(1).eq.'pbc')then
        rcut=1.d10
        vol=1.d0
        do l=1,ndim
         read(word(1+l),*)el(l)
         eli(l)=1.d0/el(l)
         rcut=min(rcut,0.4999d0*el(l))
         vol=vol*el(l)
        enddo
       elseif(word(1).eq.'vol')then
        read(word(2),*)vol
        do l=1,ndim
         el(l)=vol**(1.d0/ndim)
         eli(l)=1.d0/el(l)
        enddo
        rcut=0.4999d0*el(1)
       elseif(word(1).eq.'mu')then
        do l=1,mty
         if(word(2).eq.typename(l))then
          read(word(3),*)mu(l)
         endif
        enddo
       elseif(word(1).eq.'ntau')then
        read(word(2),*)ntau
        if(ntau.gt.mtau)then
         write(*,*)'reset ntau to ',mtau
         ntau=mtau
        endif
       elseif(word(1).eq.'verbose')then
        verbose=.true.
       elseif(word(1).eq.'canonico')then
        read(word(2),*)c_n
       elseif(word(1).eq.'restart')then
        read(word(2),*)irestart
       elseif(word(1).eq.'seed')then
        read(word(2),*)seed
        call init_random_sequence(seed)
       elseif(word(1).eq.'c0')then
        read(word(2),*)c0
       elseif(word(1).eq.'verme')then
        read(word(2),*)nblk     ! # blocks
        read(word(3),*)nblkeq   ! # equilibration blocks
        read(word(4),*)nstp     ! # steps
        read(word(5),*)nmoves   ! # moves
       elseif(word(1).eq.'temperature')then
        read(word(2),*)beta
        beta=1.d0/beta
       elseif(word(1).eq.'beta')then
        read(word(2),*)beta
       elseif(word(1).eq.'rcut')then
        read(word(2),*)rcut
       elseif(word(1).eq.'ngrid')then
        read(word(2),*)ngrid
       elseif(word(1).eq.'ngridgratio')then
        read(word(2),*)ngridgratio
       elseif(word(1).eq.'sofk')then
        read(word(2),*)nsofk
       elseif(word(1).eq.'cycles')then
        read(word(2),*)ncycles
c      elseif(word(1).eq.'nofr')then
c       read(word(2),*)nofr_stride
c       read(word(3),*)nofr_rcut
       elseif(word(1).eq.'v2')then
        do l=1,mty
         if(word(2).eq.typename(l))it=l
         if(word(3).eq.typename(l))jt=l
        enddo
        call add_tablename(word(4),nta,mta,tname)
        iv2(it,jt)=nta
        iv2(jt,it)=nta
        read(word(5),'(a)')routinename(nta)
        if(word(5).eq.'lj')then
         read(word(6),*)p(1,nta)
         read(word(7),*)p(2,nta)
        endif
       elseif(word(1).eq.'m_bar')then
        read(word(2),*)m_bar
        call i_set(nmove,emme_bar,m_bar)
       elseif(word(1).eq.'boltzmann')then
        pmove(iswap)=0.d0
        boltzmann=.true.
       elseif(word(1).eq.'move')then
        do i=1,nmove
         if(word(2).eq.name(i))then
          read(word(3),*)emme_bar(i)
          read(word(4),*)pmove(i)
         endif
        enddo
       else
        write(*,*)'no match for ',word(1)
        stop
       endif
      enddo
    1 close(2)

      drt=rcut/ngrid
      drti=1.d0/drt
      do i=1,nta
       call maketable(t(0,1,i),drt,ngrid,mgrid,tname(i),routinename(i)
     &               ,vshift(i),vtail(i),p(1,i))
      enddo
      nprops=0
      c=c0/(vol*ntau)
      tau=beta/ntau
      do i=1,ntypes
       sigma(i)=sqrt(2*lambda(i)*tau)
      enddo
      prev_tau(1)=ntau
      do itau=1,ntau-1
       next_tau(itau)=itau+1
       prev_tau(itau+1)=itau
      enddo
      next_tau(ntau)=1
      call set_pmove(nmove,pmove)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine init_random_sequence(seed)
      implicit none
      integer i,seed
      real*8 rnd
      do i=1,seed
       call random_number(rnd)
      enddo
      return
      end
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine insert
      implicit none
      include 'worm.h'
      integer itau,itau0,it,emme,j,ip,irandom_number
      real*8 rnd,quantov(0:m),x0(mdim),xnew(mdim,0:m),delta_azione,v,dv
      real*8 p
      if(gsector)return
      itau=irandom_number(ntau)
      itau0=itau
      it=irandom_number(ntypes)
      emme=irandom_number(emme_bar(iinsert))
      call randompoint(x0,el,ndim)
      call diffuse(mdim,ndim,emme,x0,xnew,sigma(it))
      call calcola_quantov(emme,quantov)
      delta_azione=0
      do j=0,emme
       call potential(xnew(1,j),v,dv,it,itau,0,quantov(j),.false.)
       delta_azione=delta_azione-tau*v
       itau=next_tau(itau)
      enddo
      p=c*vol*ntau*emme_bar(iinsert)*ntypes
     &   *exp(delta_azione+emme*mu(it)*tau)
      att(iinsert)=att(iinsert)+1
      call random_number(rnd)
      if(p.lt.rnd)return
      acc(iinsert)=acc(iinsert)+1
      call metti0(it,itau0,xnew(1,0),ip)
      call set_masha(ip,itau0,it)
      call metti(emme,ip,itau0,it,xnew)
      call set_ira(ip,itau0,it)
      if(np(itau0).gt.mnp)stop'np.gt.mnp in insert'
      zsector=.false.
      gsector=.true.
      if(verbose)write(*,*)'insert np(tau_ira) = ',np(itau0)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function irandom_number(n)
      implicit none
      integer irandom_number,n
      real*8 rnd
      call random_number(rnd)
      irandom_number=1+int(rnd*n)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine maketable(t,drt,ngrid,mgrid,tname,rname,vshift,vtail,p)
      implicit none
      integer i,j,ngrid,mgrid
      real*8 drt,t(0:mgrid,4),vshift,vtail,p(*)
      character*48 tname,rname,filename
      if(rname.eq.'hfdhe2')then
       call hfdhe2(0.d0,drt,0,ngrid,mgrid,t,vshift,vtail)
      elseif(rname.eq.'lj')then
       call lj(0.d0,drt,0,ngrid,mgrid,t,vshift,vtail,p)
      elseif(rname.eq.'zero')then
       call zero(0.d0,drt,0,ngrid,mgrid,t,vshift,vtail)
      elseif(rname.eq.'hspotential')then
       call hspotential(0.d0,drt,0,ngrid,mgrid,t,vshift,vtail)
      elseif(rname.eq.'dmix')then
       call dmix(0.d0,drt,0,ngrid,mgrid,t,vshift,vtail)
      elseif(rname.eq.'dipolar')then
       call dipolar(0.d0,drt,0,ngrid,mgrid,t,vshift,vtail)
      elseif(rname.eq.'dipolar_smoothed')then
       call dipolar_smoothed(0.d0,drt,0,ngrid,mgrid,t,vshift,vtail)
      endif
      filename=tname(1:index(tname,' ')-1)
      open(3,file=filename)
      write(3,*)'# ',ngrid,drt
      do i=0,ngrid
       write(3,'(4e20.12)')(t(i,j),j=1,4)
      enddo
      write(3,*)'# ',rname
      close(3)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine misure(i)
      implicit none
      integer i
      include 'worm.h'
      call z_sector(i)
c     call nofr(i)
      if(zsector.or.i.ne.2)then
       call energy(i)
       call numero(i)
       call rhos(i)
       call cycles(i)
       call sofk(i)
      endif
      call numerog(i)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cycles(i)
      implicit none
      include 'worm.h'
      character*13 string
      integer i,j(mty),k,l,n,iesimo,it,ndone,idone(mnp),iunit(mty),nppt
      real*8 cycles_it(mnp,mty)
      save j,string,iunit
      if(ncycles.eq.0)return
      if(i.eq.4)then; return
      elseif(i.eq.2)then
       if(gsector)return
       call r_set(mnp*mty,cycles_it(1,1),0.d0)
       ndone=0
       do n=1,np(1)
        call cerca_iesimo(iesimo,ndone,idone,jp,mnp,ntau,np(1))
        if(iesimo.eq.0)then
         do it=1,ntypes
          do l=1,ncycles
           k=j(it)+l-1
           call add_misura(k,l*1.d0/nppt(it,1),cycles_it(l,it))
          enddo
         enddo
         return
        endif
        call cycle_length(iesimo,it,l,idone,ndone)
        if(l.le.ncycles)cycles_it(l,it)=cycles_it(l,it)+1
       enddo
      elseif(i.eq.1)then
       do it=1,ntypes
        string=typename(it)
        string='cycles_'//string(1:6)
        call initialize_prop(j(it),ncycles,string,iunit(it))
       enddo
      elseif(i.eq.3)then
       do it=1,ntypes
c       call write_nonscalar_prop(j(it),ncycles,iunit(it))
        call write_cycles(j(it),ncycles,iunit(it))
       enddo
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cerca_iesimo(iesimo,ndone,idone,jp,mnp,ntau,np)
c ordinal number of first particle not in (idone(i),i=1,ndone)
      implicit none
      integer iesimo,ndone,idone(*),mnp,ntau,np,jp(mnp,ntau)
      integer primo,ip
      logical noxe
      save primo
      if(ndone.eq.0)primo=1
      do iesimo=primo,np
       ip=jp(iesimo,1)
       if(noxe(ip,idone,ndone))then
        primo=iesimo+1
        return
       endif
      enddo
      iesimo=0
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function noxe(ip,idone,ndone)
c true if ip is not in (idone(i),i=1,ndone)
      implicit none
      integer i,ip,ndone,idone(*)
      logical noxe
      noxe=.true.
      do i=1,ndone
       if(ip.eq.idone(i))then
        noxe=.false.
        return
       endif
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cycle_length(iesimo,it,l,idone,ndone)
      implicit none
      include 'worm.h'
      integer iesimo,l,ndone,idone(*),ip,itau,jtau,it
      ip=jp(iesimo,1)
      itau=1
      it=itype(ip,itau)
      do l=1,np(1)
       ndone=ndone+1
       idone(ndone)=ip
       do jtau=1,ntau
        ip=next(ip,itau)
        itau=next_tau(itau)
       enddo
       if(ip.eq.jp(iesimo,1))return
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sofk(i)
      implicit none
      include 'worm.h'
      integer i,k,l,iunit(mty),itau,iesimo,ip,mk,jjnn
      parameter(mk=14)
      real*8 sk(-mk:mk,-mk:mk),kount,tpieli(2),kr
      complex*16 rhok
      common/c_jjnn/jjnn
      save tpieli,sk,kount
      if(nsofk.eq.0)return
      if(i.eq.2)then
       do itau=1,ntau,10
        do k=-nsofk,nsofk
         do l=-nsofk,nsofk
          rhok=dcmplx(0.d0,0.d0)
          do iesimo=1,np(1)
           ip=jp(iesimo,itau)
           kr=k*tpieli(1)*x(1,ip,itau)+l*tpieli(2)*x(2,ip,itau)
           rhok=rhok+exp(dcmplx(0.d0,kr))
          enddo
          sk(k,l)=sk(k,l)+rhok*dconjg(rhok)
         enddo
        enddo
        kount=kount+1.d0
       enddo
      elseif(i.eq.1)then
       if(nsofk.gt.mk)then
        print*,'SOFK: nsofk reset to ',mk
        nsofk=mk
       endif
       tpieli(1)=2*acos(-1.d0)*eli(1)
       tpieli(2)=2*acos(-1.d0)*eli(2)
       kount=0.d0
       do k=-nsofk,nsofk
        do l=-nsofk,nsofk
         sk(l,k)=0.d0
        enddo
       enddo
      elseif(i.eq.3)then
       open(99,file='sofk')
c okkio
       sk(0,0)=0.25*(sk(1,0)+sk(0,1)+sk(-1,0)+sk(0,-1)) ! 0.d0
       do l=-nsofk,nsofk
        do k=-nsofk,nsofk
c        write(99,*)l*tpieli(1),k*tpieli(2),sk(l,k)/kount/np(1)
         write(99,*)l*tpieli(1),k*tpieli(2),sk(l,k)/kount/cml(jjnn)
        enddo
        write(99,*)
       enddo
       close(99)
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine nofr(i)
      implicit none
      include 'worm.h'
      character*13 string
      integer i,j(mty),k,l,iunit(mty),nofr_ngrid

      real*8 d,a,histo(0:mgrid),kount,r,omega,pi,dr,deltaz,deltag
      save j,string,iunit,nofr_ngrid,histo,kount,pi,deltaz,deltag
      if(nofr_stride.eq.0)return
      if(i.eq.2)then
       kount=kount+1.d0
       if(zsector)then
        deltaz=deltaz+1.d0
        return
       else
        deltag=deltag+1.d0
       endif
       if(itau_ira.ne.itau_masha)return
       d=0.d0
       do l=1,ndim
        a=x(l,ip_ira,itau_ira)-x(l,ip_masha,itau_masha)
        a=a-el(l)*nint(a*eli(l))
        d=d+a**2
       enddo
       d=sqrt(d)
c      if(d.gt.nofr_rcut)return
       k=int(d/(drt*nofr_stride))
       histo(k)=histo(k)+1.d0
      elseif(i.eq.1)then
       deltaz=0.d0
       deltag=0.d0
       pi=acos(-1.d0)
       kount=0.d0
       call r_set(mgrid+1,histo(0),0.d0)
       nofr_ngrid=ngrid/nofr_stride
      elseif(i.eq.3)then
       open(99,file='nofr')
       dr=0.d0
       do l=0,ngrid/nofr_stride
        omega=pi*dr**2
        dr=dr+drt*nofr_stride
        omega=pi*dr**2-omega
        a=(l+0.5)*drt*nofr_stride
        write(99,*)a,histo(l)/(omega*c*ntau*deltaz)/np(1)
c       write(99,*)a
c    &            ,histo(l)/(kount*omega*c*ntau*deltaz/(deltaz+deltag))
c    &            /np(1)
       enddo
       close(99)
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine energy(i)
      implicit none
      include 'worm.h'
      integer i,je,jk,jv,jpr,jg,jkt(mty),jvt(mty),jgt(mty*(mty+1)/2)
     &       ,k,l,ip,it,jt,itau,iesimo,npty(mty),nppt,kount
     &       ,iunit,g_iunit,gt_iunit(mty*(mty+1)/2),ngridg,jjnn
      real*8 v,xprev(mdim),d2(mty),dx,etot,ekin,epot,ekinty(mty),drtg
      real*8 dv,pres
      character*13 string,ztring
      common/c_jjnn/jjnn
      save je,jk,jv,jpr,jkt,jvt,jg,jgt,g_iunit,gt_iunit,ngridg,drtg
      if(i.eq.4)then; return
      elseif(i.eq.2)then
              if(verbose)write(*,*)'misura energia ',zsector
       if(gsector)return
       epot=0.d0
       pres=0.d0
       call r_set(mgridg+1,gofr(0),0.d0)
       if(ntypes.gt.1)then
        call r_set((mgridg+1)*mty*(mty+1)/2,gofrty(0,1),0.d0)
        call r_set(mty,vty,0.d0)
       endif
       kount=0
       do itau=1,ntau,10
        kount=kount+1
        do iesimo=1,np(1)
         ip=jp(iesimo,itau)
         it=itype(ip,itau)
         call potential(x(1,ip,itau),v,dv,it,itau,ip,1.d0,.true.)
         epot=epot+v
         pres=pres+dv
        enddo
       enddo
       epot=epot/kount
       pres=pres/kount
       call r_set(mdim*ntypes,wn,0.d0)
       call r_set(ntypes     ,d2,0.d0)
       itau=1
       do iesimo=1,np(itau)
        ip=jp(iesimo,itau)
        it=itype(ip,itau)
        do l=1,ndim
         xprev(l)=x(l,ip,itau)
        enddo
        do k=1,ntau
         ip=next(ip,itau)
         itau=next_tau(itau)
         do l=1,ndim
          dx=x(l,ip,itau)-xprev(l)
          dx=dx-el(l)*nint(dx*eli(l))
          wn(l,it)=wn(l,it)+dx
          d2(it)=d2(it)-dx**2
          xprev(l)=x(l,ip,itau)
         enddo
        enddo
       enddo
       do it=1,ntypes
        d2(it)=d2(it)/(ntau*4*lambda(it)*tau**2)
       enddo
       ekin=np(1)*ndim*0.5d0/tau
       do it=1,ntypes
        ekin=ekin+d2(it)
       enddo
       etot=ekin+epot
       pres=(2*ekin-pres)/(ndim*vol)
       if(np(1).ne.0)then
        etot=etot/np(1)
        ekin=ekin/np(1)
        epot=epot/np(1)
       endif
       if(ntypes.gt.1)then
        do it=1,ntypes
         npty(it)=nppt(it,1)
         ekinty(it)=npty(it)*ndim*0.5d0/tau+d2(it)
         if(npty(it).ne.0)ekinty(it)=ekinty(it)/npty(it)
         if(npty(it).ne.0)vty(it)   =vty(it)   /npty(it)/kount
        enddo
       endif
              if(verbose)write(*,*)'energia ',etot
       call add_tail(etot,epot)
       call add_misura(je,1.d0,etot)
       call add_misura(jk,1.d0,ekin)
       call add_misura(jv,1.d0,epot)
       call add_misura(jpr,1.d0,pres)
       call ngofr(ndim,np(1),1,1,1,vol,drtg,ngridg,gofr,cml(jjnn))
       do l=0,ngridg
        call add_misura(jg+l,1.d0/kount,gofr(l))
       enddo
       if(ntypes.gt.1)then
        do it=1,ntypes

         call add_misura(jkt(it),1.d0,ekinty(it))
         call add_misura(jvt(it),1.d0,vty(it)   )
         do jt=it,ntypes
          k=ikt(jt,it)
          call ngofr(ndim,npty,ntypes,it,jt,vol,drtg,ngridg,gofrty(0,k)
     &              ,0.d0)
          do l=0,ngridg
           call add_misura(jgt(ikt(jt,it))+l,1.d0/kount,gofrty(l,k))
          enddo
         enddo
        enddo
       endif
       
      elseif(i.eq.3)then
       call write_scalar_prop(je)
       call write_scalar_prop(jk)
       call write_scalar_prop(jv)
       call write_scalar_prop(jpr)
c      call write_nonscalar_prop(jg,ngridg,g_iunit)
       call write_gofr(jg,ngridg,g_iunit)
       if(ntypes.gt.1)then
        k=0
        do it=1,ntypes
         call write_scalar_prop(jkt(it))
         call write_scalar_prop(jvt(it))
         do jt=it,ntypes
          k=k+1
          call write_nonscalar_prop(jgt(k),ngridg,gt_iunit(k))
         enddo
        enddo
       endif
      elseif(i.eq.1)then
       ngridg=ngrid/ngridgratio
       drtg=drt*ngridgratio
       call initialize_prop(je,       1,'etot         ',iunit)
       call initialize_prop(jk,       1,'ekin         ',iunit)
       call initialize_prop(jv,       1,'epot         ',iunit)
       call initialize_prop(jpr,      1,'pres         ',iunit)
       call initialize_prop(jg,ngridg+1,'gofr         ',g_iunit)
       if(ntypes.gt.1)then
        do it=1,ntypes
         string=typename(it)
         string='ekin_'//string(1:8)
         call initialize_prop(jkt(it),1,string,iunit)
        enddo
        do it=1,ntypes
         string=typename(it)
         string='epot_'//string(1:8)
         call initialize_prop(jvt(it),1,string,iunit)
        enddo
        k=0
        do it=1,ntypes
         string=typename(it)
         string='gofr_'//string(1:8)
         do jt=it,ntypes
          ztring=typename(jt)
          ztring=string(1:index(string,' ')-1)//'_'//ztring(1:8)
          k=k+1
          ikt(it,jt)=k
          ikt(jt,it)=k
          call initialize_prop(jgt(k),ngridg+1,ztring,gt_iunit(k))
         enddo
        enddo
       endif
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine add_tail(etot,epot)
       implicit none
       include'worm.h'
       integer i,it,jt,nppt
       real*8 etot,epot
       if(.true.)return
       do it=1,ntypes
        do jt=1,ntypes
         i=iv2(it,jt)
         if(i.ne.0)then
          epot=epot+vtail(i)*nppt(jt,1)/vol
          etot=etot+vtail(i)*nppt(jt,1)/vol
          if(ntypes.gt.1)vty(it)=vty(it)+vtail(i)*nppt(jt,1)/vol
         endif
        enddo
       enddo
c for 1 type only: add vshift
       if(ntypes.eq.1)then
        epot=epot+vshift(i)*(np(1)-1)*0.5
        etot=etot+vshift(i)*(np(1)-1)*0.5
       endif
       return
       end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ngofr(ndim,np,ntypes,it,jt,vol,drt,ngrid,gofr,aux)
      implicit none
      integer ndim,ntypes,np(ntypes),it,jt,ngrid,i,icall
      real*8 vol,drt,gofr(0:ngrid),sd,f,den,aux
      save sd,icall
      data icall/0/
      if(icall.eq.0)then
       icall=1
       if(ndim.eq.1)sd=1.d0
       if(ndim.eq.2)sd=acos(-1.d0)
       if(ndim.eq.3)sd=4*acos(-1.d0)/3
      endif
      f=1.d0
      if(it.eq.jt)f=0.5d0
      if(aux.eq.0.d0)then
      do i=0,ngrid-1
       den=(np(it)*np(jt)*f*sd*((i+1)**ndim-i**ndim)*drt**ndim)
       if(den.ne.0.d0)then
        gofr(i)=gofr(i)*vol/(np(it)*np(jt)*f
     &                      *sd*((i+1)**ndim-i**ndim)*drt**ndim)
       else
        gofr(i)=0.d0
       endif
      enddo
      else
      do i=0,ngrid-1
       den=(aux*aux*f*sd*((i+1)**ndim-i**ndim)*drt**ndim)
       if(den.ne.0.d0)then
        gofr(i)=gofr(i)*vol/(aux*aux*f
     &                      *sd*((i+1)**ndim-i**ndim)*drt**ndim)
       else
        gofr(i)=0.d0
       endif
      enddo
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhos(i)
      implicit none
      include 'worm.h'
      integer i,j(mty),it,iunit,nppt
      real*8 rhos_it,ddot
      character*13 string
      data rhos_it/0/
      save j
      if(i.eq.4)then; return
      elseif(i.eq.2)then
       if(gsector)return
       do it=1,ntypes
        if(nppt(it,1).ne.0)rhos_it=ddot(ndim,wn(1,it),1,wn(1,it),1)
     &                     /(2*lambda(it)*ndim*ntau*tau*nppt(it,1))
        call add_misura(j(it),1.d0,rhos_it)
       enddo
      elseif(i.eq.1)then
       do it=1,ntypes
        string=typename(it)
        string='rhos_'//string(1:8)
        call initialize_prop(j(it),1,string,iunit)
       enddo
      elseif(i.eq.3)then
       do it=1,ntypes
        call write_scalar_prop(j(it))
       enddo
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine numerog(i)
      implicit none
      include 'worm.h'
      integer i,jng,iunit,it
      save jng
      if(c_n.eq.0)return
      if(i.eq.2)then
       it=prev_tau(itau_ira)
       if(gsector)call add_misura(jng ,1.d0,dble(np(it)))
      endif
      if(i.eq.1)then
       call initialize_prop(jng ,1,'ng           ',iunit)
      endif
      if(i.eq.3)then
       call write_scalar_prop(jng )
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine numero(i)
      implicit none
      include 'worm.h'
      integer i,jn,jrho,jnty(mty),jdnty(mty,mty),jrhoty(mty)
     &       ,iunit,it,jt,nppt,jjnn
      real*8 dnty
      character*13 string,ztring
      common/c_jjnn/jjnn
      save jn,jrho,jnty,jrhoty,jdnty
c     if(i.eq.4)then; return
c     elseif(i.eq.2)then
      if(i.eq.2)then
       if(gsector)return
       call add_misura(jn  ,1.d0,dble(np(1))    )
       call add_misura(jrho,1.d0,dble(np(1))/vol)
       if(ntypes.gt.1)then
        do it=1,ntypes
         call add_misura(jnty(it)  ,1.d0,dble(nppt(it,1))    )
         call add_misura(jrhoty(it),1.d0,dble(nppt(it,1))/vol)
        enddo
        do it=1,ntypes
         do jt=it+1,ntypes
          dnty=dble(nppt(jt,1)-nppt(it,1))**2
     &        /dble(nppt(jt,1)+nppt(it,1))**2
          call add_misura(jdnty(jt,it)  ,1.d0,dnty)
         enddo
        enddo
       endif
      elseif(i.eq.1)then
       call initialize_prop(jn  ,1,'np           ',iunit)
       call initialize_prop(jrho,1,'rho          ',iunit)
       jjnn=jn
       if(ntypes.gt.1)then
        do it=1,ntypes
         string=typename(it)
         string='np_'//string(1:10)
         call initialize_prop(jnty(it)  ,1,string,iunit)
         string=typename(it)
         string='rho_'//string(1:9)
         call initialize_prop(jrhoty(it),1,string,iunit)
        enddo
        do it=1,ntypes
         do jt=it+1,ntypes
          string=typename(it)
          ztring=typename(jt)
          string=ztring(1:index(ztring,' ')-1)//'_'//string
          string='dnp_'//string
          call initialize_prop(jdnty(jt,it),1,string,iunit)
         enddo
        enddo
       endif
      elseif(i.eq.3)then 
       call write_scalar_prop(jn  )
       call write_scalar_prop(jrho)
       if(ntypes.gt.1)then
        do it=1,ntypes
         call write_scalar_prop(jnty(it)  )
         call write_scalar_prop(jrhoty(it))
        enddo
        do it=1,ntypes
         do jt=it+1,ntypes
          call write_scalar_prop(jdnty(jt,it))
         enddo
        enddo
       endif
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initialize_prop(j,n,string,iunit)
      implicit none
      include 'worm.h'
      integer j,n,iunit
      character*13 string
      j=nprops+1
      write(pname(j),'(a13)')string
      nprops=nprops+n
      if(n.gt.1)then
       iunit_nsp=iunit_nsp+1
       iunit=iunit_nsp
       open(iunit,file=string)
      endif
c     print*,'j n nprops ',j,n,nprops,string
      if(nprops.gt.mprops)stop'increase mprops in worm.h'
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine z_sector(i)
      implicit none
      include 'worm.h'
      integer i,j,k
      real*8 p
      save j
      if(i.eq.2)then
       p=0.d0
       if(zsector)p=1.d0
       call add_misura(j,1.d0,p)
      elseif(i.eq.3)then; call write_scalar_prop(j)
c     elseif(i.eq.2)then; return
      elseif(i.eq.1)then; call initialize_prop(j,1,'zsector      ',k)
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine write_cycles(j,n,iunit)
      implicit none
      include 'worm.h'
      integer j,i,k,n,iunit
      rewind(iunit)
      do i=0,n-1
       k=j+i
       write(iunit,'(i4,2e19.11,e9.2,e19.11      )')i+1,blk(k),cml(k)
     &                                             ,err(k),blk_norm(k)
      enddo
      flush(iunit)
      return
      end

      subroutine write_gofr(j,n,iunit)
      implicit none
      include 'worm.h'
      integer j,i,k,n,iunit
      real*8 r
      rewind(iunit)
      do i=0,n-1
       r=(i+0.5)*drt*ngridgratio
       k=j+i
       write(iunit,'(f10.5,2e13.5,e9.2,e13.5      )')r,blk(k),cml(k)
     &                                          ,err(k),blk_norm(k)
      enddo
      flush(iunit)
      return
      end

      subroutine write_nonscalar_prop(j,n,iunit)
      implicit none
      include 'worm.h'
      integer j,i,k,n,iunit
      rewind(iunit)
      do i=0,n-1
       k=j+i
       write(iunit,'(2e19.11,e9.2,e19.11      )')blk(k),cml(k)
     &                                          ,err(k),blk_norm(k)
      enddo
      flush(iunit)
      return
      end

      subroutine write_scalar_prop(j)
      implicit none
      include 'worm.h'
      integer j
c     write(6,'(2e19.11,e9.2,e19.11,x,a13)')blk(j),cml(j)
c    &                                     ,err(j),blk_norm(j),pname(j)
      write(6,'(3e19.11,e9.2,x,a13)')blk(j),blk_norm(j)
     &                              ,cml(j),err(j),pname(j)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine masha_correzione(xnew,v,it,itau)
      implicit none
      include'worm.h'
      integer it,itau,itable,ind
      real*8 xnew(mdim),v,dv,correzione
      if(itau.ne.itau_masha)return
      correzione=0
      itable=iv2(it,it_masha)
      call add_v(ndim,xnew,x(1,ip_masha,itau_masha)
     &          ,el,eli,drt,drti,mgrid,ngrid,t(0,1,itable),correzione
     &          ,dv,ind) ! cc ,a2,rcut2)
      v=v-0.5d0*correzione
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine medie(j,iblk)
      implicit none
      include 'worm.h'
      integer i,j,iblk,iblk0
      real*8 nm1,cml_sum(mprops),cml2(mprops)
     &      ,v1(mprops),v2(mprops),blk_sum(mprops) 
      save
      if(j.eq.1)then
       do i=1,nmove
        acc(i)=0
        att(i)=0
       enddo
       iblk0=0
       if(iblk.eq.1.or.iblk.eq.nblkeq+1)then
        iblk0=1
c reset cml sum
        call r_set(nprops,cml_sum,0.d0)
        call r_set(nprops,cml2,0.d0)
        call r_set(nprops,v1,0.d0)
        call r_set(nprops,v2,0.d0)
       endif
c reset block sum
       call r_set(nprops,blk_sum,0.d0)
       call r_set(nprops,blk_norm,0.d0)
c in alternativa a initialize_prop
       call r_set(nprops,wt,0.d0)
       call r_set(nprops,prop,0.d0)
      elseif(j.eq.2)then
c update block sum
       do i=1,nprops
        blk_sum(i)=blk_sum(i)+prop(i)*wt(i)
        blk_norm(i)=blk_norm(i)+wt(i)
       enddo
       call r_set(nprops,wt,0.d0)
       call r_set(nprops,prop,0.d0)
      elseif(j.eq.3)then
c cumulative sum, error
       write(6,*)'===>> block ',iblk
       do i=1,nmove
        if(att(i).ne.0.d0)acc(i)=acc(i)/att(i)
        write(6,'(a9,f12.0,f6.3)')name(i),att(i),acc(i)
       enddo
       do i=1,nprops
        cml_sum(i)=cml_sum(i)+blk_sum(i)
        if(blk_norm(i).ne.0.d0)cml2(i)
     &                        =cml2(i)+blk_sum(i)**2/blk_norm(i)
        v2(i)=v2(i)+blk_norm(i)**2
        v1(i)=v1(i)+blk_norm(i)
        err(i)=0.d0
        if(v2(i).ne.0.d0)nm1=v1(i)**2/v2(i)-1.d0
        if(iblk0.ne.1.and.v1(i).ne.0.d0.and.nm1.ne.0.d0)
     &  err(i)=sqrt(abs(cml2(i)/v1(i)-(cml_sum(i)/v1(i))**2)/(nm1))
        if(blk_norm(i).ne.0.d0.and.v1(i).ne.0.d0)then
         blk(i)=blk_sum(i)/blk_norm(i)
         cml(i)=cml_sum(i)/v1(i)
        elseif(v1(i).ne.0.d0)then
         blk(i)=0.d0
         cml(i)=cml_sum(i)/v1(i)
        else
         blk(i)=0.d0
         cml(i)=0.d0
        endif
       enddo
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine metti0(it,itau,x0,ip)
c aggiunge un coso di tipo it al tempo itau nel punto x0. ip in output
      implicit none
      include 'worm.h'
      integer itau,it,ip,i
      real*8 x0(mdim)
      np(itau)=np(itau)+1
      ip=first_hole(itau)
      first_hole(itau)=next_hole(ip,itau)
      next_hole(ip,itau)=0
      jp(np(itau),itau)=ip
      do i=1,ndim
       x(i,ip,itau)=x0(i)
      enddo
      itype(ip,itau)=it
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine metti(emme,ip,itau,it,xnew)
c sposta ira di emme time steps lungo xnew
      implicit none
      include 'worm.h'
      integer emme,ip,itau,it,ip0,itau0,j
      real*8 xnew(mdim,0:emme)
      ip0=ip
      itau0=itau
      do j=1,emme
       itau=next_tau(itau0)
       call metti0(it,itau,xnew(1,j),ip)
       next(ip0,itau0)=ip
       prev(ip,itau)  =ip0
       itau0=itau
       ip0=ip
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine move(l)
      implicit none
      include 'worm.h'
      integer i,n,l
      real*8 s,rnd
      if(l.eq.0)then
       name(iadvance) ='advance '
       name(irecede)  ='recede  '
       name(iinsert)  ='insert  '
       name(iremove)  ='remove  '
       name(iopen)    ='open    '
       name(iclose)   ='close   '
       name(iswap)    ='swap    '
       name(iwiggle)  ='wiggle  '
       name(idisplace)='displace'
c      name(iswty)    ='swtype  '
       return
      endif
      call random_number(rnd)
      s=0.d0
      do i=1,nmove
       s=s+pmove(i)
       if(s.gt.rnd)then
        rnd=2.d0
        if(i.eq.iadvance )call advance
        if(i.eq.irecede  )call recede
        if(i.eq.iinsert  )call insert
        if(i.eq.iremove  )call remove
        if(i.eq.iopen    )call open
        if(i.eq.iclose   )call close
        if(i.eq.iswap    )call swap
        if(i.eq.iwiggle  )call wiggle
        if(i.eq.idisplace)call displace
c       if(i.eq.iswty    )call switch_type
       endif
      enddo
      if(c_n.ne.0.and.zsector.and.np(1).ne.c_n)stop'C:np'
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine next_emme_p_tau(emme,ip0,ip,itau0,itau) ! ,chi)
      implicit none
      include 'worm.h'
      integer emme,ip0,ip,itau0,itau,i ! ,chi
      itau=itau0
      ip=ip0
      do i=1,emme
       ip=next(ip,itau)
       itau=next_tau(itau)
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function next_emme_tau(emme,itau)
      implicit none
      include 'worm.h'
      integer next_emme_tau,emme,itau,i
      next_emme_tau=itau
      do i=1,emme
       next_emme_tau=next_tau(next_emme_tau)
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function prev_emme_tau(emme,itau)
      implicit none
      include 'worm.h'
      integer prev_emme_tau,emme,itau,i
      prev_emme_tau=itau
      do i=1,emme
       prev_emme_tau=prev_tau(prev_emme_tau)
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function nppt(jt,itau)
      implicit none
      include 'worm.h'
      integer nppt,it,jt,itau,i,ip
      nppt=0.d0
      do i=1,np(itau)
       ip=jp(i,itau)
       it=itype(ip,itau)
       if(it.eq.jt)nppt=nppt+1
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine open
      implicit none
      include 'worm.h'
      integer itau,ip,itau0,ip0,itau1,ip1,it,emme,j,irandom_number
      real*8 twosigma2,quantov(0:m),v,delta_azione,r2,rho,p,rnd
      real*8 dv
      if(gsector)return
      itau0=irandom_number(ntau)
      if(np(itau0).eq.0)return
      j=irandom_number(np(itau0))
      ip0=jp(j,itau0)
      it=itype(ip0,itau0)
      emme=irandom_number(emme_bar(iopen))
      twosigma2=2*emme*sigma(it)**2
      call next_emme_p_tau(emme,ip0,ip1,itau0,itau1) ! ,1)
      call gaussian(ndim,x(1,ip0,itau0),x(1,ip1,itau1)
     &             ,el,eli,twosigma2,r2,rho)
      if(r2.gt.4)return ! ekko
      ip=ip0
      itau=itau0
      call calcola_quantov(emme,quantov)
      delta_azione=0
      do j=0,emme
       call potential(x(1,ip,itau),v,dv,it,itau,ip,quantov(j),.false.)
       delta_azione=delta_azione+tau*v
       ip=next(ip,itau)
       itau=next_tau(itau)
      enddo
      p=c*emme_bar(iopen)*ntau*np(itau0)/rho
     &   *exp(delta_azione-emme*mu(it)*tau)
      att(iopen)=att(iopen)+1
      call random_number(rnd)
      if(p.lt.rnd)return
      acc(iopen)=acc(iopen)+1
      ip=prev(ip1,itau1)
      itau=prev_tau(itau1)
      call togli(emme-1,ip,itau)
      call set_masha(ip1,itau1,it)
      call set_ira(ip0,itau0,it)
      gsector=.true.
      zsector=.false.
      if(verbose)write(*,*)'open np(itau_ira) = ',np(itau0)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine potential(xp,v,dv,it,itau,salvo_ip,quanto,iall)
      implicit none
      include 'worm.h'
      integer it,itau,iesimo,kp,kt,itable,salvo_ip,ind,indg
      real*8 xp(mdim),v,dv,vij,dvij,quanto
      logical iall
      v=0.d0
      dv=0.d0
      if(iall)then ! sum_{i<j}, assume quanto=1, aggiunge gofr
       do iesimo=1,np(itau)
        kp=jp(iesimo,itau)
        if(kp.gt.salvo_ip)then
         kt=itype(kp,itau)
         itable=iv2(it,kt)
         vij=0.d0
         dvij=0.d0
         call add_v(ndim,xp,x(1,kp,itau),el,eli,drt,drti
     &            ,mgrid,ngrid,t(0,1,itable),vij,dvij,ind) ! cc ,a2,rcut2)
         indg=ind/ngridgratio
         gofr(indg)=gofr(indg)+1
         if(ind.lt.ngrid)vij=vij+vshift(itable)
          v=v+vij
          dv=dv+dvij
          if(ntypes.gt.1)then
           vty(it)=vty(it)+0.5d0*vij
           vty(kt)=vty(kt)+0.5d0*vij
           gofrty(indg,ikt(it,kt))=gofrty(indg,ikt(it,kt))+1
          endif
        endif
       enddo
      else
       do iesimo=1,np(itau)
        kp=jp(iesimo,itau)
        if(kp.ne.salvo_ip)then
         kt=itype(kp,itau)
         itable=iv2(it,kt)
         call add_v(ndim,xp,x(1,kp,itau),el,eli,drt,drti
     &            ,mgrid,ngrid,t(0,1,itable),v,dv,ind) ! cc ,a2,rcut2)
        endif
       enddo
       v=quanto*v
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine randompoint(x0,el,ndim)
      implicit none
      integer i,ndim
      real*8 x0(ndim),el(ndim),rnd
      do i=1,ndim
       call random_number(rnd)
       x0(i)=el(i)*(rnd-0.5d0)
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine readwords(iunit,mword,words,eof_flag,record,echo)
c read words from record. if iunit.ne.0 record is read from unit
      implicit none
      integer iunit,mword,iscan,n_items,i,lrec,j,eof_flag,echo
      character*80 record
      character*48 word,words(mword)
      character*1  char
      do i=1,mword
       words(i)=' '
      enddo
      do i=1,48
       word(i:i)=' '
      enddo
      n_items=0
      eof_flag=0
    4 continue
      if(iunit.ne.0)then
       read(iunit,'(a)',end=3)record
       if(echo.eq.1)write(6,*)record
      endif
c the first item is the number of subsequent items
      lrec=80
      iscan=0
c find next item
      do i=1,10000
c read next char
        iscan=iscan+1
c end of record
        if(iscan.eq.lrec+1)go to 2
        read(record(iscan:iscan),'(a)')char
        if(char.eq.'\\')then
         go to 4
        endif
        if(char.ne.' ')then
c item found
         n_items=n_items+1
         do j=1,100
          word(j:j)=char
          iscan=iscan+1
          read(record(iscan:iscan),'(a)')char
          if(char.eq.' ')go to 1
         enddo
    1    read(word,'(a)')words(n_items)
        endif
c reset word
        do j=1,48
         word(j:j)=' '
        enddo
      enddo
    2 return
    3 eof_flag=1
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine recede
      implicit none
      include 'worm.h'
      logical cemasha
      integer emme,ip,ipp,itau,itaup,it,j,irandom_number,prev_emme_tau
      real*8 quantov(0:m),delta_azione,v,dv,p,rnd
      if(zsector)return
      emme=irandom_number(emme_bar(irecede))
      call get_ira(ip,itau,it)
      if(cemasha(emme,ip,ipp,itau,itaup))return
      if(c_n.ne.0.and.np(prev_emme_tau(emme,itau)).lt.c_n-1)return
      call calcola_quantov(emme,quantov)
      delta_azione=0
      do j=0,emme
       call potential(x(1,ip,itau),v,dv,it,itau,ip,quantov(j),.false.)
       call masha_correzione(x(1,ip,itau),v,it,itau)
       delta_azione=delta_azione+tau*v
       ip=prev(ip,itau)
       itau=prev_tau(itau)
      enddo
      p=exp(delta_azione-emme*mu(it)*tau) 
      att(irecede)=att(irecede)+1
      call random_number(rnd)
      if(p.lt.rnd)return
      acc(irecede)=acc(irecede)+1
      call get_ira(ip,itau,it)
      call togli(emme,ip,itau)
      call set_ira(ip,itau,it)
      if(verbose)write(*,*)'recede np(itau_ira) = ',np(itau)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine remove
      implicit none
      include 'worm.h'
      logical ceira
      integer emme,ip,itau,it,j
      real*8 quantov(0:m),v,dv,delta_azione,p,rnd
      if(zsector)return
      if(c_n.ne.0.and.np(itau_ira).ne.c_n+1)return
      emme=emme_bar(iremove)
      call get_masha(ip,itau,it)
      if(.not.ceira(emme,ip,itau))return
      call calcola_quantov(emme,quantov)
      delta_azione=0
      do j=0,emme
       call potential(x(1,ip,itau),v,dv,it,itau,ip,quantov(j),.false.)
       delta_azione=delta_azione+tau*v
       ip=prev(ip,itau)
       itau=prev_tau(itau)
      enddo
      p=exp(delta_azione-emme*mu(it)*tau)
     &   /(c*vol*ntau*emme_bar(iremove)*ntypes)
      att(iremove)=att(iremove)+1
      call random_number(rnd)
      if(p.lt.rnd)return
      acc(iremove)=acc(iremove)+1
      call get_ira(ip,itau,it)
      call togli(emme+1,ip,itau)
      zsector=.true.
      gsector=.false.
      ip_ira=0
      ip_masha=0
      itau_ira=0
      itau_masha=0
      it_ira=0
      it_masha=0
      if(verbose)write(*,*)'remove np(1) = ',np(1)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine restart(i)
      implicit none
      include 'worm.h'
      integer i,j,l,k,itau,itau0,jtau,ip,it0,np0
      real*8 x0(mdim),a,b
      if(i.lt.0)then
       open(10,file='restart.ind',status='unknown')
       call scrivi_indici(10)
       close(10)
       open(10,file='restart.coord',status='unknown')
       call scrivi_coord(10)
       close(10)
       return
      endif
      if(irestart.eq.0)then
       if(c_n.ne.0)stop'C irestart 1 o -1'
       do itau=1,ntau
        np(itau)=0
        first_hole(itau)=1
        do ip=1,mnp-1
         next_hole(ip,itau)=ip+1
        enddo
        next_hole(mnp,itau)=1
       enddo
       zsector=.true.
       gsector=.false.
      elseif(irestart.lt.0)then ! legge np posizioni
c      open(10,file='worm.init',status='unknown')
c      write(*,*)'coordinate iniziali'
c      do ip=1,mnp
c       read(10,*,end=1)(x0(l),l=1,ndim),it0
c       write(*,*)(x0(l),l=1,ndim),it0
      
       call triangolo16(el,x(1,1,1))
       np0=16
       it0=1
       do ip=1,np0
        do itau=1,ntau
         do l=1,ndim
c         x(l,ip,itau)=x0(l)
          x(l,ip,itau)=x(l,ip,1)
         enddo
         next(ip,itau)=ip
         prev(ip,itau)=ip
         itype(ip,itau)=it0
        enddo
       enddo
c   1  np0=ip-1
       write(*,*)'initial conf. with N=16 on a triangular lattice'
       do itau=1,ntau
        np(itau)=np0
        do ip=1,np0
         jp(ip,itau)=ip
        enddo
        first_hole(itau)=np0+1
        do ip=1,mnp-1
         next_hole(ip,itau)=ip+1
        enddo
        next_hole(mnp,itau)=1
       enddo
       zsector=.true.
       gsector=.false.
      elseif(irestart.gt.0)then ! legge la configurazione
       open(10,file='restart.ind',status='unknown')
       call leggi_indici(10)
       close(10)
       open(10,file='restart.coord',status='unknown')
       call leggi_coord(10)
       close(10)
      endif
      return
      end

      subroutine triangolo16(el,x)
      implicit none
      real*8 el(2),x(2,16),a,b
      a=el(1)/4
      b=el(2)/2

      x(1,1)=0.0  ; x(2,1)=0.0
      x(1,2)=0.5*a; x(2,2)=0.5*b
      x(1,3)=a    ; x(2,3)=0.0
      x(1,4)=1.5*a; x(2,4)=0.5*b
      x(1,5)=2*a  ; x(2,5)=0.0
      x(1,6)=2.5*a; x(2,6)=0.5*b
      x(1,7)=3*a  ; x(2,7)=0.0
      x(1,8)=3.5*a; x(2,8)=0.5*b
      
      x(1, 9)=0.0  ; x(2, 9)=b
      x(1,10)=0.5*a; x(2,10)=1.5*b
      x(1,11)=a    ; x(2,11)=b
      x(1,12)=1.5*a; x(2,12)=1.5*b
      x(1,13)=2*a  ; x(2,13)=b
      x(1,14)=2.5*a; x(2,14)=1.5*b
      x(1,15)=3*a  ; x(2,15)=b
      x(1,16)=3.5*a; x(2,16)=1.5*b

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine leggi_indici(i)
      implicit none
      include 'worm.h'
      integer i,j,k,ip,itau,jtau,itau0
       read(i,*)k,itau,gsector
       read(i,*)ip_ira,itau_ira,it_ira
       read(i,*)ip_masha,itau_masha,it_masha
       if(k.ne.mnp)stop'restart: k.ne.mnp'
       read(i,*)next_hole
       itau0=itau
       do jtau=1,ntau
        read(i,*)itau,np(itau),first_hole(itau)
        do j=1,np(itau)
         read(i,*)ip,next(ip,itau),prev(ip,itau),itype(ip,itau)
         jp(j,itau)=ip
        enddo
        itau=next_tau(itau)
       enddo
       zsector=.not.gsector
       read(i,*)
       read(i,*)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine leggi_coord(i)
      implicit none
      include 'worm.h'
      integer i,j,l,ip,itau,itau0,jtau
       itau0=1
       if(gsector)itau0=itau_masha
       do j=1,np(itau0)
        itau=itau0
        ip=jp(j,itau)
        do jtau=1,ntau
         if(ip.ne.0)then
          read(i,*)itau,(x(l,ip,itau),l=1,ndim)
          ip=next(ip,itau)
          itau=next_tau(itau)
         endif
        enddo
       enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sample(n,rho0,sigma,k)
      implicit none
      integer k,n
      real*8 rho0(n),sigma,somma,rnd
      somma=0
      call random_number(rnd)
      do k=1,n
       somma=somma+rho0(k)
       if(somma.ge.rnd*sigma)return
      enddo
      stop'sample'
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine scrivi_coord(i)
      implicit none
      include 'worm.h'
      integer i,itau,itau0,jtau,j,l,ip,it
      real*8 xx(mdim)
       itau0=itau_masha
       if(zsector)itau0=1
       do j=1,np(itau0)
        itau=itau0
        ip=jp(j,itau)
        it=itype(ip,itau)
        do jtau=1,ntau
         if(ip.ne.0)then
          do l=1,ndim
           xx(l)=x(l,ip,itau)-el(l)*nint(x(l,ip,itau)*eli(l))
          enddo
          write(i,*)itau,(xx(l),l=1,ndim),'it ',it
          ip=next(ip,itau)
          itau=next_tau(itau)
         endif
        enddo
       enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine scrivi_indici(i)
      implicit none
      include 'worm.h'
      integer i,itau,jtau,j,ip
       itau=itau_masha
       if(zsector)itau=1
       write(i,*)mnp,itau,gsector
       write(i,*)ip_ira,itau_ira,it_ira
       write(i,*)ip_masha,itau_masha,it_masha
       write(i,*)next_hole
       do jtau=1,ntau
        write(i,*)itau,np(itau),first_hole(itau)
        do j=1,np(itau)
         ip=jp(j,itau)
         write(i,*)ip,next(ip,itau),prev(ip,itau),itype(ip,itau)
        enddo
        itau=next_tau(itau)
       enddo
       write(i,*)
       write(i,*)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine set_ira(ip,itau,it)
      implicit none
      include'worm.h'
      integer ip,itau,it
      ip_ira=ip
      itau_ira=itau
      it_ira=it
      next(ip,itau)=0
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine set_masha(ip,itau,it)
      implicit none
      include 'worm.h'
      integer ip,itau,it
      ip_masha=ip
      itau_masha=itau
      it_masha=it
      prev(ip,itau)=0
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine set_pmove(n,pmove)
      implicit none
      integer i,n
      real*8 pmove(n),s
      s=0.d0
      do i=1,n
       s=s+pmove(i)
      enddo
      do i=1,n
       pmove(i)=pmove(i)/s
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sommarho(somma,ipart,ip,itau,it,itau1,twosigma2,rho,k)
      implicit none
      include 'worm.h'
      integer ipart(mnp),ip,itau,it,itau1,iesimo,k,kp
      real*8 somma,twosigma2,rho(mnp),r2,rho0
      somma=0
      k=0
      do iesimo=1,np(itau1)
       kp=jp(iesimo,itau1)
       if(itype(kp,itau1).eq.it)then
        call gaussian(ndim,x(1,ip,itau),x(1,kp,itau1)
     &               ,el,eli,twosigma2,r2,rho0)
        if(r2.lt.4)then ! ekko
         k=k+1
         ipart(k)=kp
         rho(k)=rho0
         somma=somma+rho(k)
        endif
       endif
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine spline(mgrid,nt,drt,t)
      implicit none
      integer mgrid,nt,i
      real*8 drt,t(0:mgrid,4)
      do i=0,nt
       t(i,2)=t(i,2)*drt
       t(i,3)=t(i,3)*drt*drt*0.5d0
      enddo
      do i=0,nt-1
       t(i,3)= 3.d0*(t(i+1,1)-t(i,1))  -(t(i+1,2)+2.d0*t(i,2))
       t(i,4)=-2.d0*(t(i+1,1)-t(i,1))  +(t(i+1,2)+     t(i,2))
      enddo
      t(nt,1)=0.d0
      t(nt,2)=0.d0
      t(nt,3)=0.d0
      t(nt,4)=0.d0
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine swap
      implicit none
      include'worm.h'
      logical cemasha
      integer ip,itau,it,ip_zeta,itau_zeta,ip_alpha,itau_alpha
     &       ,emme,k,i,j,ipart(mnp),nppt,nterms,next_emme_tau
      real*8 twosigma2,sigmai,sigmaz,rho0(mnp),xnew(mdim,0:m),p,rnd
     &      ,v,dv,delta_azione
      if(zsector)return
      call get_ira(ip,itau_zeta,it)
      emme=emme_bar(iswap)
      twosigma2=2*emme*sigma(it)**2
      if(nppt(it,itau_zeta).lt.2)return
      itau_alpha=next_emme_tau(emme,itau_zeta)
      if(nppt(it,itau_alpha).lt.1)return
c     if(nppt(it,itau_zeta).lt.1)return
      call sommarho(sigmai,ipart,ip_ira,itau_ira,it,itau_alpha
     &             ,twosigma2,rho0,nterms)
      if(nterms.eq.0)return
      call sample(nterms,rho0,sigmai,k)
      ip_alpha=ipart(k)
      if(cemasha(emme,ip_alpha,ip_zeta,itau_alpha,itau))return
      if(itau_zeta.ne.itau)stop'swap merda'
      call sommarho(sigmaz,ipart,ip_zeta,itau_zeta,it,itau_alpha
     &             ,twosigma2,rho0,nterms)
      if(nterms.eq.0)return
      call bridge(x(1,ip_ira,itau_ira),x(1,ip_alpha,itau_alpha)
     &           ,emme,sigma(it),el,eli,ndim,mdim,xnew)
      delta_azione=0
      j=0
      call potential(x(1,ip_ira,itau_ira),v,dv,it,itau_ira,ip_ira,0.5d0
     &              ,.false.)
      delta_azione=delta_azione-tau*v
      call potential(x(1,ip_zeta,itau_zeta),v,dv,it,itau_zeta,ip_zeta
     &              ,0.5d0,.false.)
      delta_azione=delta_azione+tau*v
      ip=next(ip_zeta,itau_zeta)
      itau=next_tau(itau_zeta)
      do j=1,emme-1
       call potential(x(1,ip,itau),v,dv,it,itau,ip,1.0d0,.false.)
       call masha_correzione(x(1,ip,itau),v,it,itau)
       delta_azione=delta_azione+tau*v
       call potential(xnew(1,j),v,dv,it,itau,ip,1.0d0,.false.)
       call masha_correzione(xnew(1,j),v,it,itau)
       delta_azione=delta_azione-tau*v
       ip=next(ip,itau)
       itau=next_tau(itau)
      enddo
      p=exp(delta_azione)*sigmai/sigmaz
      att(iswap)=att(iswap)+1
      call random_number(rnd)
      if(p.lt.rnd)return
      acc(iswap)=acc(iswap)+1
      ip=next(ip_zeta,itau_zeta)
      itau=next_tau(itau_zeta)
      next(ip_ira,itau_ira)=ip
      prev(ip,itau)=ip_ira
      do j=1,emme-1
       do i=1,ndim
        x(i,ip,itau)=xnew(i,j)
       enddo
       ip=next(ip,itau)
       itau=next_tau(itau)
      enddo
      call set_ira(ip_zeta,itau_zeta,it)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine togli0(ip,itau)
c toglie il coso ip al tempo itau. 
      implicit none
      include 'worm.h'
      integer itau,ip,iesimo
      do iesimo=1,np(itau)
       if(jp(iesimo,itau).eq.ip)then
        jp(iesimo,itau)=jp(np(itau),itau)
        jp(np(itau),itau)=0 ! ip
        next_hole(ip,itau)=first_hole(itau)
        first_hole(itau)=ip
        np(itau)=np(itau)-1
        return
       endif
      enddo
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine togli(n,ip,itau)
c toglie n cosi a partire da (ip,itau)
      implicit none
      include'worm.h'
      integer n,itau,ip,itau0,ip0,j,i
      do j=1,n
       ip0=ip
       itau0=itau
       call togli0(ip0,itau0)
       ip=prev(ip0,itau0)
       itau=prev_tau(itau0)
       prev(ip0,itau0)=0
       next(ip0,itau0)=0
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     subroutine switch_type
c     implicit none
c     include 'worm.h'
c     integer ip,itau,it,jt,ip0,k,j,irandom_number,ip1,itau1
c     if(gsector)return
c     itau=1
c     if(np(itau).eq.0)return
c     j=irandom_number(np(itau))
c     ip=jp(j,itau)
c     it=itype(ip,itau)
c     jt=irandom_number(ntypes)
c     if(it.eq.jt)return
c     att(iswty)=att(iswty)+1
c     acc(iswty)=acc(iswty)+1
c     ip0=ip
c     do k=1,np(itau)
c      do j=1,ntau
c       itype(ip,itau)=jt
c       ip=next(ip,itau)
c       itau=next_tau(itau)
c      enddo
c      if(ip.eq.ip0)return
c     enddo
c     stop'switch_type failure'
c     end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine displace
      implicit none
      include 'worm.h'
      integer ip,itau,it,ip0,itau0,i,j,irandom_number
      real*8 rnd,dr(mdim),xnew(mdim),delta_azione,v,dv,p,delta
      if(gsector)return
      itau0=1
      if(np(itau0).eq.0)return
      delta=el(1)*0.05 ! /emme_bar(idisplace)
      j=irandom_number(np(itau0))
      ip0=jp(j,itau0)
      it=itype(ip0,itau0)
      ip=ip0
      itau=itau0
      do j=1,ntau
       ip=next(ip,itau)
       itau=next_tau(itau)
      enddo
      if(ip.ne.ip0)return
      do i=1,ndim
       call random_number(rnd)
       dr(i)=(rnd-0.5d0)*delta
      enddo
      delta_azione=0
      ip=ip0
      itau=itau0
      do j=1,ntau
       call potential(x(1,ip,itau),v,dv,it,itau,ip,1.0d0,.false.)
       delta_azione=delta_azione+tau*v
       do i=1,ndim
        xnew(i)=x(i,ip,itau)+dr(i)
       enddo
       call potential(xnew        ,v,dv,it,itau,ip,1.0d0,.false.)
       delta_azione=delta_azione-tau*v
       ip=next(ip,itau)
       itau=next_tau(itau)
      enddo
      p=exp(delta_azione)
      att(idisplace)=att(idisplace)+1
      call random_number(rnd)
      if(p.lt.rnd)return
      acc(idisplace)=acc(idisplace)+1
      ip=ip0
      itau=itau0
      do j=1,ntau
       do i=1,ndim
        x(i,ip,itau)=x(i,ip,itau)+dr(i)
       enddo
       ip=next(ip,itau)
       itau=next_tau(itau)
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function xi(sigma)
      implicit none
      real*8 xi,sigma,r1,r2,pi
      data pi/3.14159265358979d0/
      call random_number(r1)
      call random_number(r2)
      xi=cos(pi*r1)*sigma*sqrt(-2.d0*log(r2))
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine r_set(n,a,b)
      implicit none
      integer i,n
      real*8 a(n),b
      do i=1,n
       a(i)=b
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine i_set(n,a,b)
      implicit none
      integer i,n,a(n),b
      do i=1,n
       a(i)=b
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function ddot(n,x,ix,y,iy)
      implicit none
      integer i,n,ix,iy
      real*8 ddot,x(n),y(n)
      ddot=0.d0
      do i=1,n
       ddot=ddot+x(i)*y(i)
      enddo
      return
      end
