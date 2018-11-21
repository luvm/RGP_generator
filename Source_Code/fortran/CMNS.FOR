      subroutine SETINFOdll (icomp,wmm,ttrp,tnbpt,tc,pc,Dc,Zc,
     &   acf,dip,Rgas)
cDEC$ ATTRIBUTES DLLEXPORT, Alias: "_SETINFOdll"::SETINFOdll
C  cDEC$ ATTRIBUTES STDCALL, REFERENCE::SETINFOdll
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c     dll_export SETINFOdll
      parameter (ncmax=20)        !max number of components in mixture
      parameter (refmax=10)        !max number of fluids for transport ECS
      parameter (n0=-ncmax-refmax,nx=ncmax)
      common /NCOMP/ nc,ic
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
      if (abs(icomp).le.nc) then
        wm(icomp)=wmm
        ttp(icomp)=ttrp
        tnbp(icomp)=tnbpt
        tcrit(icomp)=tc
        pcrit(icomp)=pc
        Dcrit(icomp)=Dc
        Zcrit(icomp)=Zc
        accen(icomp)=acf
        dipole(icomp)=dip
        R=Rgas
      else
        wmm=0.0d0
        ttrp=0.0d0
        tnbpt=0.0d0
        tc=0.0d0
        pc=0.0d0
        Dc=0.0d0
        Zc=0.0d0
        acf=0.0d0
        dip=0.0d0
        Rgas=R
      end if
c

      end
c ======================================================================
      subroutine GETINFOdll (icomp,wmm,ttrp,tnbpt,
     &                       tc,pc,Dc,Zc,acf,dip,Rgas)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
cDEC$ ATTRIBUTES DLLEXPORT, Alias: "_GETINFOdll"::GETINFOdll
C  cDEC$ ATTRIBUTES STDCALL, REFERENCE::GETINFOdll
c     dll_export GETINFOdll
      call INFO (icomp,wmm,ttrp,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas)
      end
c ======================================================================
      subroutine GETCMNdll (tc,rhoc,pc,tred,rhored,
     &                      n1,n2,coefhmx,n3,n4,coefcp0,coefxk)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (refmax=10)        !max number of fluids for transport ECS
      parameter (ncppmx=20)      !max number of Cp0 terms
      parameter (n0=-ncmax-refmax,nx=ncmax)
      parameter (mxtrm=72)
      parameter (mxcoef=72)
      dimension coefhmx(mxcoef),coefcp0(mxcoef),coefxk(mxcoef)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /WNTFEQ/ ntermf(n0:nx),ncoeff(n0:nx),
     &                ntp(n0:nx),ndp(n0:nx),nlp(n0:nx),
     &                ncrt(n0:nx),ncfcrt(n0:nx),
     &                nspare(n0:nx),ncfsp(n0:nx)
      common /WCFFEQ/ a(n0:nx,mxtrm),ti(n0:nx,mxtrm),di(n0:nx,mxtrm),
     &                gi(n0:nx,mxtrm),gi2(n0:nx,mxtrm),
     &                dli(n0:nx,mxtrm),tli(n0:nx,mxtrm),
     &                tpower(n0:nx,mxtrm),dpower(n0:nx,mxtrm),
     &                dlpower(n0:nx,mxtrm),
     &                rho0feq(n0:nx),t0feq(n0:nx),
     &                pcfeq(n0:nx),rhocfeq(n0:nx),tcfeq(n0:nx),
     &                wmf(n0:nx),Rfeq(n0:nx),
     &                pminfeq(n0:nx),rhotp(n0:nx),tminfeq(n0:nx),
     &                tmaxfeq(n0:nx),pmaxfeq(n0:nx)
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),ncosh(n0:nx),
     &                nsinh(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
      common /CPPSAV/ cp0sav(n0:nx),cpisav(n0:nx),cptsav(n0:nx),
     &                tsav(n0:nx)
cDEC$ ATTRIBUTES DLLEXPORT, Alias: "_GETCMNdll"::GETCMNdll
C  cDEC$ ATTRIBUTES DLLEXPORT, Alias: "_GETCMNdll@48"::GETCMNdll
C  cDEC$ ATTRIBUTES STDCALL, REFERENCE::GETCMNdll
c
c     dll_export GETCMNdll
      do i=1,mxcoef
        coefhmx(i)=0.d0
        coefcp0(i)=0.d0
        coefxk(i)=0.d0
      enddo
      tred=tz(1)
      rhored=rhoz(1)
      tc=tcfeq(1)
      rhoc=rhocfeq(1)
      pc=pcfeq(1)
      n1=ntermf(1)
      n2=ncrt(1)
      do i=1,n1+n2
        coefhmx(i)=a(1,i)
      enddo
      n3=ntermc(1)
      n4=nterme(1)
      do j=1,n3+n4
        coefcp0(j)=cpc(1,j)
        coefxk(j)=xk(1,j)
      enddo
      end
c ======================================================================
      subroutine SETCMNdll (tc,rhoc,pc,tred,rhored,
     &                      n1,n2,coefhmx,n3,n4,coefcp0,coefxk)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (refmax=10)        !max number of fluids for transport ECS
      parameter (ncppmx=20)      !max number of Cp0 terms
      parameter (n0=-ncmax-refmax,nx=ncmax)
      parameter (mxtrm=72)
      parameter (mxcoef=72)
      dimension coefhmx(mxcoef),coefcp0(mxcoef),coefxk(mxcoef)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /WNTFEQ/ ntermf(n0:nx),ncoeff(n0:nx),
     &                ntp(n0:nx),ndp(n0:nx),nlp(n0:nx),
     &                ncrt(n0:nx),ncfcrt(n0:nx),
     &                nspare(n0:nx),ncfsp(n0:nx)
      common /WCFFEQ/ a(n0:nx,mxtrm),ti(n0:nx,mxtrm),di(n0:nx,mxtrm),
     &                gi(n0:nx,mxtrm),gi2(n0:nx,mxtrm),
     &                dli(n0:nx,mxtrm),tli(n0:nx,mxtrm),
     &                tpower(n0:nx,mxtrm),dpower(n0:nx,mxtrm),
     &                dlpower(n0:nx,mxtrm),
     &                rho0feq(n0:nx),t0feq(n0:nx),
     &                pcfeq(n0:nx),rhocfeq(n0:nx),tcfeq(n0:nx),
     &                wmf(n0:nx),Rfeq(n0:nx),
     &                pminfeq(n0:nx),rhotp(n0:nx),tminfeq(n0:nx),
     &                tmaxfeq(n0:nx),pmaxfeq(n0:nx)
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),ncosh(n0:nx),
     &                nsinh(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
      common /CPPSAV/ cp0sav(n0:nx),cpisav(n0:nx),cptsav(n0:nx),
     &                tsav(n0:nx)
cDEC$ ATTRIBUTES DLLEXPORT, Alias: "_SETCMNdll"::SETCMNdll
C  cDEC$ ATTRIBUTES STDCALL, REFERENCE::SETCMNdll
c     dll_export SETCMNdll
c
      if (pc.gt.0.d0) pcfeq(1)=pc
      if (tc.gt.0.d0) tcfeq(1)=tc
      if (rhoc.gt.0.d0) rhocfeq(1)=rhoc
      if (tred.gt.0.d0) tz(1)=tred
      if (tred.gt.0.d0) t0feq(1)=tc
      if (rhored.gt.0.d0) rhoz(1)=rhored
      if (rhored.gt.0.d0) rho0feq(1)=rhoc
      if (n1.ne.0) then
        ntermf(1)=n1
        ncrt(1)=n2
        do i=1,n1+n2
          a(1,i)=coefhmx(i)
        enddo
        call SETEXP (1)
      endif
      if (n3.ne.0) then
        tsav(1)=0.d0
        ntermc(1)=n3
        nterme(1)=n4
        do j=1,n3+n4
          cpc(1,j)=coefcp0(j)
          xk(1,j)=coefxk(j)
        enddo
      endif
      end
c ======================================================================
      subroutine GETCMN (tc,rhoc,pc,tred,rhored,acf,
     &                   n1,n2,coefhmx,n3,n4,coefcp0,coefxk)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
cDEC$ ATTRIBUTES DLLEXPORT::GETCMN
c     dll_export GETCMN
      parameter (ncmax=20)        !max number of components in mixture
      parameter (refmax=10)        !max number of fluids for transport ECS
      parameter (ncppmx=20)      !max number of Cp0 terms
      parameter (n0=-ncmax-refmax,nx=ncmax)
      parameter (mxtrm=72)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /WNTFEQ/ ntermf(n0:nx),ncoeff(n0:nx),
     &                ntp(n0:nx),ndp(n0:nx),nlp(n0:nx),
     &                ncrt(n0:nx),ncfcrt(n0:nx),
     &                nspare(n0:nx),ncfsp(n0:nx)
      common /WCFFEQ/ a(n0:nx,mxtrm),ti(n0:nx,mxtrm),di(n0:nx,mxtrm),
     &                gi(n0:nx,mxtrm),gi2(n0:nx,mxtrm),
     &                dli(n0:nx,mxtrm),tli(n0:nx,mxtrm),
     &                tpower(n0:nx,mxtrm),dpower(n0:nx,mxtrm),
     &                dlpower(n0:nx,mxtrm),
     &                rho0feq(n0:nx),t0feq(n0:nx),
     &                pcfeq(n0:nx),rhocfeq(n0:nx),tcfeq(n0:nx),
     &                wmf(n0:nx),Rfeq(n0:nx),
     &                pminfeq(n0:nx),rhotp(n0:nx),tminfeq(n0:nx),
     &                tmaxfeq(n0:nx),pmaxfeq(n0:nx)
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),ncosh(n0:nx),
     &                nsinh(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
      common /CPPSAV/ cp0sav(n0:nx),cpisav(n0:nx),cptsav(n0:nx),
     &                tsav(n0:nx)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)

c
      parameter (mxcoef=72)
      dimension coefhmx(mxcoef),coefcp0(mxcoef),coefxk(mxcoef)
      do i=1,mxcoef
        coefhmx(i)=0.d0
        coefcp0(i)=0.d0
        coefxk(i)=0.d0
      enddo
      tred=tz(1)
      rhored=rhoz(1)
      tc=tcfeq(1)
      rhoc=rhocfeq(1)
      pc=pcfeq(1)
      n1=ntermf(1)
      n2=ncrt(1)
      acf = accen(1)
      do i=1,n1+n2
        coefhmx(i)=a(1,i)
      enddo
      n3=ntermc(1)
      n4=nterme(1)
      do j=1,n3+n4
        coefcp0(j)=cpc(1,j)
        coefxk(j)=xk(1,j)
      enddo
      end
c ======================================================================
      subroutine SETCMN (tc,rhoc,pc,tred,rhored,acf,
     &                   n1,n2,coefhmx,n3,n4,coefcp0,coefxk)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
cDEC$ ATTRIBUTES DLLEXPORT::SETCMN
c     dll_export SETCMN

      parameter (ncmax=20)        !max number of components in mixture
      parameter (refmax=10)        !max number of fluids for transport ECS
      parameter (ncppmx=20)      !max number of Cp0 terms
      parameter (n0=-ncmax-refmax,nx=ncmax)
      parameter (mxtrm=72)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /WNTFEQ/ ntermf(n0:nx),ncoeff(n0:nx),
     &                ntp(n0:nx),ndp(n0:nx),nlp(n0:nx),
     &                ncrt(n0:nx),ncfcrt(n0:nx),
     &                nspare(n0:nx),ncfsp(n0:nx)
      common /WCFFEQ/ a(n0:nx,mxtrm),ti(n0:nx,mxtrm),di(n0:nx,mxtrm),
     &                gi(n0:nx,mxtrm),gi2(n0:nx,mxtrm),
     &                dli(n0:nx,mxtrm),tli(n0:nx,mxtrm),
     &                tpower(n0:nx,mxtrm),dpower(n0:nx,mxtrm),
     &                dlpower(n0:nx,mxtrm),
     &                rho0feq(n0:nx),t0feq(n0:nx),
     &                pcfeq(n0:nx),rhocfeq(n0:nx),tcfeq(n0:nx),
     &                wmf(n0:nx),Rfeq(n0:nx),
     &                pminfeq(n0:nx),rhotp(n0:nx),tminfeq(n0:nx),
     &                tmaxfeq(n0:nx),pmaxfeq(n0:nx)
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),ncosh(n0:nx),
     &                nsinh(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
      common /CPPSAV/ cp0sav(n0:nx),cpisav(n0:nx),cptsav(n0:nx),
     &                tsav(n0:nx)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /FEQSAV/ phisav(n0:nx,mxtrm),delsav(n0:nx),tausav(n0:nx),
     &                taup(n0:nx,mxtrm),delp(n0:nx,mxtrm),
     &                delli(n0:nx,mxtrm),drvsav(n0:nx,16)

c
      parameter (mxcoef=72)
      dimension coefhmx(mxcoef),coefcp0(mxcoef),coefxk(mxcoef)

c  (re)initialize contents of /FEQSAV/ and /CPPSAV/ when a fluid's parameter are reset
      do 120 i=n0,nx
      delsav(i)=0.0d0
      tausav(i)=0.0d0
      cp0sav(i)=0.0d0
      cpisav(i)=0.0d0
      cptsav(i)=0.0d0
      tsav(i)=0.0d0
      do 100 j=1,mxtrm
      phisav(i,j)=0.0d0
      taup(i,j)=0.0d0
      delp(i,j)=0.0d0
      delli(i,j)=0.0d0
  100 continue
  120 continue
c
      call RESETA
      if (pc.gt.0.d0) pcfeq(1)=pc
      if (tc.gt.0.d0) tcfeq(1)=tc
      if (rhoc.gt.0.d0) rhocfeq(1)=rhoc
      if (tred.gt.0.d0) tz(1)=tred
      if (tred.gt.0.d0) t0feq(1)=tc
      if (rhored.gt.0.d0) rhoz(1)=rhored
      if (rhored.gt.0.d0) rho0feq(1)=rhoc
      if (abs(acf).gt.1.d-20) accen(1)=acf
      if (n1.ne.0) then
        ntermf(1)=n1
        ncrt(1)=n2
        do i=1,n1+n2
          a(1,i)=coefhmx(i)
        enddo
        call SETEXP (1)
      endif
      if (n3.ne.0) then
        tsav(1)=0.d0
        ntermc(1)=n3
        nterme(1)=n4
        do j=1,n3+n4
          cpc(1,j)=coefcp0(j)
          xk(1,j)=coefxk(j)
        enddo
      endif
      end
c ======================================================================
      subroutine GETALLCMN (tc,rhoc,pc,tred,rhored,acf,
     &                   n1,n2,coefhmx,tauexp,denexp,denlnexp
     &                      ,n3,n4,coefcp0,coefxk)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
cDEC$ ATTRIBUTES DLLEXPORT::GETALLCMN
c     dll_export GETALLCMN

      parameter (ncmax=20)        !max number of components in mixture
      parameter (refmax=10)        !max number of fluids for transport ECS
      parameter (ncppmx=20)      !max number of Cp0 terms
      parameter (n0=-ncmax-refmax,nx=ncmax)
      parameter (mxtrm=72)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /WNTFEQ/ ntermf(n0:nx),ncoeff(n0:nx),
     &                ntp(n0:nx),ndp(n0:nx),nlp(n0:nx),
     &                ncrt(n0:nx),ncfcrt(n0:nx),
     &                nspare(n0:nx),ncfsp(n0:nx)
      common /WCFFEQ/ a(n0:nx,mxtrm),ti(n0:nx,mxtrm),di(n0:nx,mxtrm),
     &                gi(n0:nx,mxtrm),gi2(n0:nx,mxtrm),
     &                dli(n0:nx,mxtrm),tli(n0:nx,mxtrm),
     &                tpower(n0:nx,mxtrm),dpower(n0:nx,mxtrm),
     &                dlpower(n0:nx,mxtrm),
     &                rho0feq(n0:nx),t0feq(n0:nx),
     &                pcfeq(n0:nx),rhocfeq(n0:nx),tcfeq(n0:nx),
     &                wmf(n0:nx),Rfeq(n0:nx),
     &                pminfeq(n0:nx),rhotp(n0:nx),tminfeq(n0:nx),
     &                tmaxfeq(n0:nx),pmaxfeq(n0:nx)
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),ncosh(n0:nx),
     &                nsinh(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
      common /CPPSAV/ cp0sav(n0:nx),cpisav(n0:nx),cptsav(n0:nx),
     &                tsav(n0:nx)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
      parameter (mxcoef=72)
      dimension coefhmx(mxcoef),coefcp0(mxcoef),coefxk(mxcoef)
      dimension tauexp(mxcoef),denexp(mxcoef),denlnexp(mxcoef)
      do i=1,mxcoef
        coefhmx(i)=0.d0
        tauexp(i)=0.d0
        denexp(i)=0.d0
        denlnexp(i)=0.d0
        coefcp0(i)=0.d0
        coefxk(i)=0.d0
      enddo
      tred=tz(1)
      rhored=rhoz(1)
      tc=tcfeq(1)
      rhoc=rhocfeq(1)
      pc=pcfeq(1)
      acf = accen(1)
      n1=ntermf(1)
      n2=ncrt(1)
      do i=1,n1+n2
        coefhmx(i)=a(1,i)
        tauexp(i)=ti(1,i)
        denexp(i)=di(1,i)
        denlnexp(i)=dli(1,i)
      enddo
      n3=ntermc(1)
      n4=nterme(1)
      do j=1,n3+n4
        coefcp0(j)=cpc(1,j)
        coefxk(j)=xk(1,j)
      enddo
      end

c ======================================================================
      subroutine GETTYPE(htype)
cDEC$ ATTRIBUTES DLLEXPORT::GETTYPE
c     dll_export GETTYPE


      character*3 htype
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (refmax=10)        !max number of fluids for transport ECS
      parameter (n0=-ncmax-refmax,nx=ncmax)
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)

        htype=heos
      end

c ======================================================================
      subroutine LIMITSET (tmin,tmax,Dmax,pmax)
cDEC$ ATTRIBUTES DLLEXPORT::LIMITSET

      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c     dll_export LIMITSET
      parameter (ncmax=20)       !max number of components in mixture
      parameter (refmax=10)        !max number of fluids for transport ECS
      parameter (n0=-ncmax-refmax,nx=ncmax)

      common /NCOMP/ nc,ic
      common /EOSLIM/ tmeos(n0:nx),txeos(n0:nx),peos(n0:nx),Deos(n0:nx)

c     reset limit only for the case of a pure component
      if (nc.eq.1) then
c  special case--pure component
        tmeos(1)=tmin
        txeos(1)=tmax
        Deos(1)=Dmax
        peos(1)=pmax
      end if
      end
C END LIMITSet
