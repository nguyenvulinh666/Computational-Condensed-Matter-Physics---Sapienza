      integer mnp,m,mtau,mdim,mty,mgrid,mgridg,mta,mprops
      parameter(mnp=50,m=201,mtau=1001,mdim=2,mty=2,mgrid=1001
     &         ,mgridg=201,mta=6,mprops=5205)

      integer ngrid,iv2(mty,mty),iv1(mty)
      real*8 t(0:mgrid,4,mta),drt,drti,vshift(mta),vtail(mta)
      common /c_tables/t,drt,drti,vshift,vtail,ngrid,iv2,iv1

      integer ip_ira,ip_masha,itau_ira,itau_masha,it_ira,it_masha
      common /c_ira_masha/ip_ira,ip_masha,itau_ira,itau_masha
     &                   ,it_ira,it_masha

      integer jp(mnp,mtau),np(mtau),next(mnp,mtau),prev(mnp,mtau)
     &       ,ntypes,itype(mnp,mtau)
     &       ,first_hole(mtau),last_hole(mtau),next_hole(mnp,mtau)
      real*8 x(mdim,mnp,mtau)
      common /c_particles/x,jp,np,next,prev,ntypes,itype
     &                   ,first_hole,last_hole,next_hole

      integer next_tau(mtau),prev_tau(mtau),ntau
      real*8 tau
      common /c_tau/tau,ntau,next_tau,prev_tau

      integer ndim
      real*8 sigma(mty),lambda(mty)
      real*8 el(mdim),eli(mdim),rcut,rcut2,a2,vol
      real*8 mu(mty),beta
      character*48 typename(mty)
      common /c_system/mu,el,eli,rcut,rcut2,a2,vol,sigma,lambda
     &                ,beta,ndim,typename

      integer nblk,nblkeq,nstp,nmoves,irestart
      common /c_steps/nblk,nblkeq,nstp,nmoves,irestart

      logical zsector,gsector
      integer nmove,iadvance,irecede,iinsert,iremove,iopen,iclose
     &       ,iwiggle,iswap,idisplace ! ,iswty
      parameter(nmove=9,iadvance=1,irecede=2,iinsert=3,iremove=4
     &         ,iopen=5,iclose=6,iwiggle=7,iswap=8,idisplace=9) ! ,iswty=10)
      integer m_bar,emme_bar(nmove)
      real*8 c,pmove(nmove)
      real*8 att(nmove),acc(nmove)
      real*8 dr_tr
      character*9 name(nmove)
      common /c_verme/att,acc
     &               ,dr_tr
     &               ,c,pmove
     &               ,m_bar,emme_bar
     &               ,zsector,gsector
     &               ,name

      integer nprops
      character*13 pname(mprops)
      real*8 prop(mprops),wt(mprops)
     &      ,err(mprops),blk(mprops),cml(mprops),blk_norm(mprops)
      common /c_medie/prop,wt,err,blk,cml,blk_norm,nprops,pname

      real*8 nofr_rcut
      integer iunit_nsp,ncycles,nofr_stride,nsofk
      common /c_nsp/nofr_rcut,iunit_nsp,ncycles,nofr_stride,nsofk

      real*8 wn(mdim,mty)
      real*8 vty(mty),gofr(0:mgridg),gofrty(0:mgridg,mty*(mty+1)/2)
      integer ikt(mty,mty),ngridgratio
      common /c_pot/wn,vty,gofr,gofrty,ikt,ngridgratio

      logical verbose
      integer c_n
      common /c_ensemble/c_n,verbose

      logical boltzmann
      common /c_statistica/boltzmann
