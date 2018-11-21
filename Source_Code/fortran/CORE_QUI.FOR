c  begin file core_QUI.f
c
c  This file contains the routines implementing the Helmholtz form of
c  the pure fluid equation of state in the Quintic form.
c
c  contained here are:
c     function PHIQUI (icomp,itau,idel,tau,del)
c     subroutine CRTQUI (icomp,tcrit,pcrit,Dcrit)
c     subroutine REDQUI (icomp,tred,Dred)
c     subroutine SETQUI (nread,icomp,hcasno,ierr,herr)
c
c ======================================================================
c ======================================================================
c
      function PHIQUI (icomp,itau,idel,tau,del)
c
c  compute reduced Helmholtz energy or a derivative as functions
c  of dimensionless temperature and density for the Helmholtz-explicit
c  equation of state
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c     itau--flag specifying order of temperature derivative to calc
c     idel--flag specifying order of density derivative to calculate
c           when itau = 0 and idel = 0, compute A/RT
c           when itau = 0 and idel = 1, compute 1st density derivative
c           when itau = 1 and idel = 1, compute cross derivative
c           etc.
c      tau--dimensionless temperature (To/T)
c      del--dimensionless density (D/Do)
c  output (as function value):
c      phi--residual (real-gas) part of the Helmholtz energy, or one
c           of its derivatives (as specified by itau and idel),
c           in reduced form (A/RT)
c           itau  idel    output (dimensionless for all cases)
c             0    0      A/RT
c             1    0      tau*[d(A/RT)/d(tau)]
c             2    0      tau**2*[d**2(A/RT)/d(tau)**2]
c             0    1      del*[d(A/RT)/d(del)]
c             0    2      del**2*[d**2(A/RT)/d(del)**2]
c             1    1      tau*del*[d**2(A/RT)/d(tau)d(del)]
c                         etc.
c
c  written by D.E. Cristancho, NIST Thermophysics Division, Boulder, Colorado
c  06-18-09 DEC, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (mxtrm=72)
      character*16 drvflg(n0:nx)
c  numbers of terms associated with the "normal" Helmholtz function
c  (plus numbers of unique powers of temperature, density, etc.),
      common /WNTQUI/ neta(n0:nx),neps(n0:nx),nbb(n0:nx),ngam(n0:nx),
     &                nbet(n0:nx)
      common /WCFQUI/ aq(n0:nx,mxtrm),tiq(n0:nx,mxtrm),diq(n0:nx,mxtrm),
     &                rho0q(n0:nx),t0q(n0:nx),
     &                pcq(n0:nx),rhocq(n0:nx),tcq(n0:nx),
     &                wmfq(n0:nx),Rqui(n0:nx),
     &                pminq(n0:nx),rhotpq(n0:nx),tminq(n0:nx),
     &                tmaxq(n0:nx),pmaxq(n0:nx)
      common /QUISV2/ drvflg
      common /QUISAV/ phisav(n0:nx,mxtrm),delsav(n0:nx),tausav(n0:nx),
     &                taup(n0:nx,mxtrm),delp(n0:nx,mxtrm),
     &                drvsav(n0:nx,16)
      common /QUITERMS/ eta,eta10,eta20,eps,eps10,eps20,
     &                 bb,bb10,bb20,gam,gam10,gam20
c  common storing the fluid constants
      common /CCON/ tcc(n0:nx),pcc(n0:nx),rhocc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
      phiqui=0.d0
      if (del.le.1.0d-10) RETURN  !trivial solution at zero density for
      if (tau.le.0.d0) RETURN    !  any and all derivatives
c
      ncode=idel*4+itau+1
      nterm=neta(icomp)+neps(icomp)+nbb(icomp)+ngam(icomp)+nbet(icomp)
      if (abs(tau-tausav(icomp)).lt.1.0d-12 .and.
     &    abs(del-delsav(icomp)).lt.1.0d-16) then
c  retrieve value from previous call
        if (drvflg(icomp)(ncode:ncode).eq.'1') then
          phiqui=drvsav(icomp,ncode)
          RETURN
        endif
      else
c  otherwise, compute new values and save for possible future use
c  first compute needed powers of tau and del (and save for future use)
        drvflg(icomp)='0000000000000000'
        if (abs(tau-tausav(icomp)).gt.1.0d-12) then
          elntau=log(tau)
          tausav(icomp)=tau
          do j=1,nterm
            taup(icomp,j)=tiq(icomp,j)*elntau
          enddo
        end if
        if (abs(del-delsav(icomp)).gt.1.0d-16) then
          elndel=log(del)
          delsav(icomp)=del
          do j=1,nterm
            delp(icomp,j)=diq(icomp,j)*elndel
          enddo
        end if
c
        phisum=0.d0
        eta=0.d0
        eta10=0.d0
        eta20=0.d0
        eps=0.d0
        eps10=0.d0
        eps20=0.d0
        bb=0.d0
        bb10=0.d0
        bb20=0.d0
        gam=0.d0
        gam10=0.d0
        gam20=0.d0
        bet=0.d0
        bet10=0.d0
        bet20=0.d0
        j=0
        do k=1,neta(icomp)
          j=j+1
          eta=eta+aq(icomp,j)*tau**tiq(icomp,j)
          eta10=eta10+tiq(icomp,j)*aq(icomp,j)*tau**(tiq(icomp,j)-1.d0)
          eta20=eta20+tiq(icomp,j)*(tiq(icomp,j)-1.d0)*
     &           aq(icomp,j)*tau**(tiq(icomp,j)-2.d0)
        enddo
        do k=1,neps(icomp)
          j=j+1
          eps=eps+aq(icomp,j)*tau**tiq(icomp,j)
          eps10=eps10+tiq(icomp,j)*aq(icomp,j)*tau**(tiq(icomp,j)-1.d0)
          eps20=eps20+tiq(icomp,j)*(tiq(icomp,j)-1.d0)*
     &           aq(icomp,j)*tau**(tiq(icomp,j)-2.d0)
        enddo
        do k=1,nbb(icomp)
          j=j+1
          bb=bb+aq(icomp,j)*tau**tiq(icomp,j)
          bb10=bb10+tiq(icomp,j)*aq(icomp,j)*tau**(tiq(icomp,j)-1.d0)
          bb20=bb20+tiq(icomp,j)*(tiq(icomp,j)-1.d0)*
     &           aq(icomp,j)*tau**(tiq(icomp,j)-2.d0)
        enddo
        do k=1,ngam(icomp)
          j=j+1
          gam=gam+aq(icomp,j)*tau**tiq(icomp,j)
          gam10=gam10+tiq(icomp,j)*aq(icomp,j)*tau**(tiq(icomp,j)-1.d0)
          gam20=gam20+tiq(icomp,j)*(tiq(icomp,j)-1.d0)*
     &           aq(icomp,j)*tau**(tiq(icomp,j)-2.d0)
        enddo
        do k=1,nbet(icomp)
          j=j+1
          bet=bet+aq(icomp,j)*tau**tiq(icomp,j)
          bet10=bet10+tiq(icomp,j)*aq(icomp,j)*tau**(tiq(icomp,j)-1.d0)
          bet20=bet20+tiq(icomp,j)*(tiq(icomp,j)-1.d0)*
     &           aq(icomp,j)*tau**(tiq(icomp,j)-2.d0)
        enddo
        phisum=eta*(del-log(1.d0-bb*del)/bb)
     &        -eps*log((1.d0+gam*del)/(1.d0-bet*del))
c         ex=taup(icomp,k)+delp(icomp,k)
c         if (ex.lt.100.d0 .and. ex.gt.-200.d0) then
c           phisav(icomp,k)=aq(icomp,k)*EXP(ex)
c         else
c           phisav(icomp,k)=0.d0
c         endif
c         phisum=phisum+phisav(icomp,k)
        phiqui=phisum
        drvflg(icomp)(1:1)='1'
        drvsav(icomp,1)=phiqui
      end if
c
c  check if derivatives are requested, calculations make use of fact
c  that terms in derivative summations are very similar to A/RT terms
c
      if (idel.eq.1) then
c  compute derivative w.r.t. del (dimensionless density)
c  save individual terms for possible use in cross derivative
        phisum=0.d0
        phisum=eta*(1.d0+1.d0/(1.d0-bb*del))
     &        -eps*(gam/(1.d0+gam*del)+bet/(1.d0-bet*del))
        phiqui=phisum*del
c
      elseif (idel.eq.2) then
c  compute 2nd derivative w.r.t. del (dimensionless density)
c  save individual terms for possible use in cross derivative
        phisum=0.d0
c        do k=1,nterm
c          dik=diq(icomp,k)
c          phi02(k)=phisav(icomp,k)*(dik**2-dik)
c           phisum+phi02(k)
c          phisum=eta*(1.d0+1.d0/(1.d0-bb*del))
c              eps*(gam/(1.d0+gam*del)+bet/(1.d0-bet*del))
c         enddo
        phisum=eta*bb/(1.d0-bb*del)**2
     &        -eps*(-gam**2/(1.d0+gam*del)**2
     &        +bet**2/(1.d0-bet*del)**2)
        phiqui=phisum*del**2
c
      elseif (idel.eq.3) then
c  compute 3rd derivative w.r.t. del (dimensionless density)
        phisum=0.d0
c        do k=1,nterm
c          dik=diq(icomp,k)
c          phi03(k)=phisav(icomp,k)
c    &     +6.0d0*dik-3.0d0*dik-3.0d0*dik**2
c          phisum=phisum+phi03(k)
c        enddo
        phisum=2.d0*eta*bb**2/(1.d0-bb*del)**3
     &        -2.d0*eps*(gam**3/(1.d0+gam*del)**3
     &        +bet**3/(1.d0-bet*del)**3)
        phiqui=phisum*del**3
      end if
c
c
c   epsi0,bbio,gami0,beti0,etai0 are the i tau derivatives of fitting
c   parameters
c
      if (itau.eq.1) then
c  compute derivative w.r.t. tau (dimensionless temperature)
c  save individual terms for possible use in cross derivative
        phisum=0.d0
          phisum=eta10*(del-log(1.d0-bb*del)/bb)
     &        +eta*bb10/bb*(log(1.d0-bb*del)/bb+del/(1.d0-bb*del))
     &        -eps10*log((1.d0+gam*del)/(1.d0-bet*del))
     &        -eps*del*(gam10/(1.d0+gam*del)
     &        +bet10/(1.d0-bet*del))
        phiqui=phisum*tau
c
      elseif (itau.eq.2) then
c  compute 2nd derivative w.r.t. tau (dimensionless temperature)
c  save individual terms for possible use in cross derivative
        phisum=0.d0
c        do k=1,nterm
c          tik=tiq(icomp,k)
c          phi20(k)=phisav(icomp,k)*tik*(tik-1.0d0)
c          phisum=phisum+phi20(k)
c        enddo
          phisum=eta20*(del-log(1.d0-bb*del)/bb)
     &        +1.d0/bb*(2.d0*eta10*bb10+eta*bb20-eta*bb10**2/bb)
     &        *(log(1.d0-bb*del)/bb+del/(1.d0-bb*del))
     &        +eta*bb10**2/bb**2*(del*(2.d0*bb*del-1)/
     &        (1.d0-bb*del)**2-log(1.d0-bb*del)/bb)
     &        -eps20*log((1.d0+gam*del)/(1.d0-bet*del))
     &        -2.d0*eps10*del*(gam10/(1.d0+gam*del)
     &        +bet10/(1.d0-bet*del))
     &        -eps*del*(gam20/(1.d0+gam*del)
     &        +bet20/(1.d0-bet*del)-del*(gam10**2/
     &        (1.d0+gam*del)**2-bet10**2/(1.d0-bet*del)**2))

        phiqui=phisum*tau**2
c
C    not third derivative!!!!!
      elseif (itau.eq.3) then
c  compute 3rd derivative w.r.t. tau (dimensionless temperature)
        phisum=0.d0
        do k=1,nterm
          tik=tiq(icomp,k)
          phisum=phisum+phisav(icomp,k)*tik*(tik-1.d0)*(tik-2.d0)
        enddo
        phiqui=phisum
      end if
c
c
      if (itau.eq.1 .and. idel.eq.1) then
c  compute cross derivative using terms from 1st derivatives
        phisum=0.d0
          phisum=eta10*(1.d0+1.d0/(1.d0-bb*del))
     &        +eta*bb10*del/(1.d0-bb*del)**2
     &        -eps10*(gam/(1.d0+gam*del)+bet/(1.d0-bet*del))
     &        -eps*(gam10/(1.d0+gam*del)
     &        +bet10/(1.d0-bet*del)-del*(gam*gam10/
     &        (1.d0+gam*del)**2-bet*bet10/(1.d0-bet*del)**2))
        phiqui=phisum*del*tau
c
      elseif (itau.eq.2 .and. idel.eq.1) then
c  compute cross derivative using term from 1st derivative
        phisum=0.d0
C       do k=1,nterm
C         tik=tiq(icomp,k)
C         phisum=phisum+(tik*tik-tik)*phi01(k)
C       enddo
        phiqui=phisum
c
      elseif (itau.eq.1 .and. idel.eq.2) then
c  compute cross derivative using term from 2nd derivative
        phisum=0.d0
C       do k=1,nterm
C         phisum=phisum+tiq(icomp,k)*phi02(k)
C       enddo
        phiqui=phisum
c
      elseif (itau.eq.2 .and. idel.eq.2) then
c  compute cross derivative using terms from 2nd derivative
        phisum=0.d0
C       do k=1,nterm
C         tik=tiq(icomp,k)
C         phisum=phisum+(tik*tik-tik)*phi02(k)
C       enddo
        phiqui=phisum
c
      end if
c
      drvsav(icomp,ncode)=phiqui
      drvflg(icomp)(ncode:ncode)='1'
c
      RETURN
      end                                               !function PHIQUI
c
c ======================================================================
c
      subroutine CRTQUI (icomp,tcrit,pcrit,Dcrit)
c
c  returns critical parameters associated with Fundamental EOS
c
c  input:
c    icomp--pointer specifying component (1..nc)
c  outputs:
c    tcrit--critical temperature (K)
c    pcrit--critical pressure (kPa)
c    Dcrit--molar density (mol/L) at critical point
c
c  written by D.E. Cristancho, NIST Thermophysics Division, Boulder, Colorado
c  06-18-09 DEC, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
c
      parameter (mxtrm=72)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFQUI/ aq(n0:nx,mxtrm),tiq(n0:nx,mxtrm),diq(n0:nx,mxtrm),
     &                rho0q(n0:nx),t0q(n0:nx),
     &                pcq(n0:nx),rhocq(n0:nx),tcq(n0:nx),
     &                wmfq(n0:nx),Rqui(n0:nx),
     &                pminq(n0:nx),rhotpq(n0:nx),tminq(n0:nx),
     &                tmaxq(n0:nx),pmaxq(n0:nx)
c
      tcrit=tcq(icomp)
      pcrit=pcq(icomp)
      Dcrit=rhocq(icomp)
c
      RETURN
      end                                             !subroutine CRTQUI
c
c ======================================================================
c
      subroutine REDQUI (icomp,tred,Dred)
c
c  returns reducing parameters associated with Fundamental EOS;
c  used to calculate the 'tau' and 'del' which are the independent
c  variables in the EOS
c
c  input:
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c  outputs:
c     tred--reducing temperature (K)
c     Dred--reducing molar density (mol/L)
c
c  written by D.E. Cristancho, NIST Thermophysics Division, Boulder, Colorado
c  06-18-09 DEC, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
c
      parameter (mxtrm=72)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFQUI/ aq(n0:nx,mxtrm),tiq(n0:nx,mxtrm),diq(n0:nx,mxtrm),
     &                rho0q(n0:nx),t0q(n0:nx),
     &                pcq(n0:nx),rhocq(n0:nx),tcq(n0:nx),
     &                wmfq(n0:nx),Rqui(n0:nx),
     &                pminq(n0:nx),rhotpq(n0:nx),tminq(n0:nx),
     &                tmaxq(n0:nx),pmaxq(n0:nx)
c
      tred=t0q(icomp)
      Dred=rho0q(icomp)
c
      RETURN
      end                                             !subroutine REDQUI
c
c ======================================================================
c
      subroutine SETQUI (nread,icomp,hcasno,ierr,herr)
c
c  set up working arrays for use with Fundamental equation of state
c
c  inputs:
c    nread--file to read data from
c           <= 0 get data from block data
c           >0 read from logical unit nread (file should have already
c              been opened and pointer set by subroutine SETUP)
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c   hcasno--CAS number of component icomp (not req'd if reading from file)
c
c  outputs:
c     ierr--error flag:  0 = successful
c                        1 = error (e.g. fluid not found)
c     herr--error string (character*255 variable if ierr<>0)
c     other quantities returned via arrays in common /WCFQUI/
c
c  written by D.E. Cristancho, NIST Thermophysics Division, Boulder, Colorado
c  06-18-09 DEC, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (mxtrm=72)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*3 hpheq,heos,hmxeos,hmodcp
      character*12 hcasno
      character*255 herr
c     character*1 dummy
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c  commons associated with the nc components of current interest
c  ("working" commons and arrays)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c     common /CPMOD/ hmodcp(n0:nx)
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /WNTQUI/ neta(n0:nx),neps(n0:nx),nbb(n0:nx),ngam(n0:nx),
     &                nbet(n0:nx)
      common /WCFQUI/ aq(n0:nx,mxtrm),tiq(n0:nx,mxtrm),diq(n0:nx,mxtrm),
     &                rho0q(n0:nx),t0q(n0:nx),
     &                pcq(n0:nx),rhocq(n0:nx),tcq(n0:nx),
     &                wmfq(n0:nx),Rqui(n0:nx),
     &                pminq(n0:nx),rhotpq(n0:nx),tminq(n0:nx),
     &                tmaxq(n0:nx),pmaxq(n0:nx)
c  limits associated with the equation of state
      common /EOSLIM/ tmn(n0:nx),tmx(n0:nx),pmx(n0:nx),rhomx(n0:nx)
      common /QUISAV/ phisav(n0:nx,mxtrm),delsav(n0:nx),tausav(n0:nx),
     &                taup(n0:nx,mxtrm),delp(n0:nx,mxtrm),
     &                drvsav(n0:nx,16)
c
      ierr=0
      herr=' '
c  (re)initialize contents of /QUISAV/ when a new fluid is read in
      do i=n0,nx
        delsav(i)=0.d0
        tausav(i)=0.d0
        do j=1,mxtrm
          phisav(i,j)=0.d0
          taup(i,j)=0.d0
          delp(i,j)=0.d0
        enddo
      enddo
c
      if (nread.le.0 .or. hcasno.eq.' ') then
c  get coefficients from block data
c  identify specified fluid with entries in database via match of CAS no
        ierr=1
        herr='[SETQUI error] fluid input to SETQUI not found'//hnull
      else
c  read data from file
        read (nread,*) tminq(icomp)          !lower temperature limit
        read (nread,*) tmaxq(icomp)          !upper temperature limit
        read (nread,*) pmaxq(icomp)          !upper pressure limit
        read (nread,*) rhomx(icomp)         !upper density limit
        read (nread,2003) hmodcp(icomp)     !pointer to Cp0 model
        read (nread,*) wm(icomp)     !molecular weight
        wmfq(icomp)=wm(icomp)
        read (nread,*) ttp(icomp)    !triple point temperature
        read (nread,*) pminq(icomp)   !pressure at triple point
        read (nread,*) rhotpq(icomp)  !density at triple point
        read (nread,*) tnbp(icomp)   !normal boiling point temperature
        read (nread,*) accen(icomp)  !acentric factor
        read (nread,*) tcq(icomp),pcq(icomp),rhocq(icomp) !critical par
        tcrit(icomp)=tcq(icomp)
        pcrit(icomp)=pcq(icomp)
        Dcrit(icomp)=rhocq(icomp)
        ptp(icomp)=pminq(icomp)
        dtp(icomp)=rhotpq(icomp)
        dtpv(icomp)=0.d0
        dnbpl(icomp)=0.d0
        dnbpv(icomp)=0.d0
        read (nread,*) t0q(icomp),rho0q(icomp) !reducing parameters
        tz(icomp)=t0q(icomp)
        rhoz(icomp)=rho0q(icomp)
        read (nread,*) Rqui(icomp)    !gas constant used in fit
        if (nc.eq.1 .and. icomp.eq.1) R=Rqui(icomp)
        Reos(icomp)=Rqui(icomp)
        Zcrit(icomp)=pcq(icomp)/(Rqui(icomp)*tcq(icomp)*rhocq(icomp))
        read (nread,*) neta(icomp),neps(icomp),nbb(icomp),ngam(icomp),
     &                 nbet(icomp)
        nterm=neta(icomp)+neps(icomp)+nbb(icomp)+ngam(icomp)+nbet(icomp)
        do j=1,nterm
          read (nread,*) aq(icomp,j),tiq(icomp,j),diq(icomp,j)
        enddo
      end if
c
c  copy limits into /EOSLIM/ arrays
      tmn(icomp)=tminq(icomp)
      tmx(icomp)=tmaxq(icomp)
      pmx(icomp)=pmaxq(icomp)
c     rhomx(icomp)=rhomaxq(icomp)
c
      RETURN
 2003 format (a3)
      end                                             !subroutine SETQUI
c
c ======================================================================
c
c ======================================================================
c                                                    end file core_QUI.f
c ======================================================================
