c  begin file core_FEQ.f
c
c  This file contains the routines implementing the Helmholtz form of
c  the pure fluid equation of state (the so-called "Fundamental Eqn").
c
c  contained here are:
c     function PHIFEQ (icomp,itau,idel,tau,del)
c     subroutine CRTFEQ (icomp,tcrit,pcrit,Dcrit)
c     subroutine REDFEQ (icomp,tred,Dred)
c     subroutine SETFEQ (nread,icomp,hcasno,ierr,herr)
c     subroutine SETEXP (icomp)
c     block data BDFEQ
c
c ======================================================================
c ======================================================================
c
      function PHIFEQ (icomp,itau,idel,tau,del)
c
c  compute reduced Helmholtz energy or a derivative as functions
c  of dimensionless temperature and density for the Helmholtz-explicit
c  equation of state
c
c  based on Lemmon et al., J. Chem. Eng. Data, 2009.
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
c  N.B.  The reducing parameters To and Do are often, but not
c        necessarily, equal to the critical temperature and density.
c
c        The Helmholtz energy consists of ideal and residual (real-gas)
c        terms; this routine calculates only the residual part.
c
c        This function computes pure component properties only.
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  07-26-95  MM, original version
c  08-14-95  MM, put saved variables into common (rather than save stmt)
c  11-01-95  MM, increase parameter mxtrm to 52 (to accommodate steam)
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  01-16-96  MM, implicit integer (i-n); (include L)
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c  04-18-96  MM, apply tolerance to del=delsav + fix bug (phi not saved)
c  10-09-96  MM, fix potential underflow on exponential term
c  01-03-97  MM, use sorted powers of tau, del, stored in tpower, etc.
c                new arrays added to /WNTFEQ/,/WLFFEQ/,/WCFFEQ/
c  01-22-97 EWL, add critical-region terms of Wagner to Helmholtz equation.
c                increase parameter mxtrm to 56 (to accommodate steam)
c  01-31-97  MM, modify commons associated with critical-region terms
c  06-03-97 EWL, add third derivative of Helmholtz energy with respect to density.
c  07-02-97 EWL, add parameter for g in exp(-g*d^l) for fluids where g<>1
c                modify /WCFFEQ/ to accommodate g parameter
c  04-06-98 EWL, change mxcrt to 8
c  04-06-98 EWL, add term li2 to /WLFFQ2/, use li2 in exponential for tau power
c  04-06-98 EWL, increase parameter mxtrm to 72 (to accommodate krypton)
c  08-11-98  MM, add third temperature derivative and third-order cross derivatives
c                N.B. these not yet implemented for critical-region terms
c  12-28-98 EWL, change arrays in CRTSAV to include icomp dependence
c   8-25-99 EWL, change the calculation of DELB and TAUA at the critical point
c   4-24-00 EWL, change the exponent of del in the EXP term from integer to
c                real*8 and modify the common blocks WCFFEQ and WLFFEQ
c   7-11-00 EWL, remove krypton pieces
c  09-00-00 EWL, add critical region terms for idel=3
c                add additional logic to return if property has already been
c                calculated.  optimize calculations to increase speed.
c                return del*d(alpha)/d(del), etc. (put the del in the calculation)
c  12-00-00 EWL, changed dli2 to tli
c  12-00-00 EWL, make dli and tli integers if delb or taua is less than zero
c  07-16-01 EWL, change check on delsav from 10^-12 to 10^-16
c  11-04-04 OBD, add check for ex<-200
c  02-02-06 EWL, increase efficiency of code, calculation speed
c  10-19-06 EWL, remove INT in hxp(icomp,i)=alpha(icomp,i)*...**INT(dli(icomp,k)), and in txp()=...
c                it was causing calculated values for R125 to come out wrong
c  10-29-07 EWL, add 4th derivative of Helmholtz energy with respect to del^2 and tau^2
c  12-12-07 MLH, added general Helmholtz equation of Xiang and Deiters (2007)
c  09-08-09 EWL, optimize equations as given in propane publication
c  09-08-09 EWL, finished da^3/dt^2dd for the critical region terms
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (mxtrm=72,mxcrt=20)
      character*16 drvflg(n0:nx)
c  numbers of terms associated with the "normal" Helmholtz function
c  (plus numbers of unique powers of temperature, density, etc.),
c  the critical-region terms of Wagner, plus spare for future use
      common /WNTFEQ/ ntermf(n0:nx),ncoeff(n0:nx),
     &                ntp(n0:nx),ndp(n0:nx),nlp(n0:nx),
     &                ncrt(n0:nx),ncfcrt(n0:nx),
     &                nspare(n0:nx),ncfsp(n0:nx)
      common /WLFFEQ/ itp(n0:nx,mxtrm),
     &                idp(n0:nx,mxtrm),ilp(n0:nx,mxtrm)
      common /WCFFEQ/ a(n0:nx,mxtrm),ti(n0:nx,mxtrm),di(n0:nx,mxtrm),
     &                gi(n0:nx,mxtrm),gi2(n0:nx,mxtrm),
     &                dli(n0:nx,mxtrm),tli(n0:nx,mxtrm),
     &                tpower(n0:nx,mxtrm),dpower(n0:nx,mxtrm),
     &                dlpowr(n0:nx,mxtrm),
     &                rho0(n0:nx),t0(n0:nx),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),
     &                wmf(n0:nx),Rfeq(n0:nx),
     &                pmin(n0:nx),rhotp(n0:nx),tmin(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
c  parameters associated with the critical-region terms of Wagner
      common /WCFFQ2/ alpha(n0:nx,mxcrt),beta(n0:nx,mxcrt),
     &                gamma(n0:nx,mxcrt),delta(n0:nx,mxcrt),
     &                eta(n0:nx,mxcrt),eid(n0:nx,mxcrt),eit(n0:nx,mxcrt)
      common /CRTSAV/ delb(n0:nx,mxcrt),taua(n0:nx,mxcrt),
     &                txp(n0:nx,mxcrt),hxp(n0:nx,mxcrt),
     &                axp(n0:nx,mxcrt),dxp(n0:nx,mxcrt),
     &                ext(n0:nx,mxcrt),extd(n0:nx,mxcrt),
     &                extt(n0:nx,mxcrt),extdt(n0:nx,mxcrt),
     &                extt2(n0:nx,mxcrt),extd2(n0:nx,mxcrt)
c
      common /FEQSV2/ drvflg
      common /FEQSAV/ phisav(n0:nx,mxtrm),delsav(n0:nx),tausav(n0:nx),
     &                taup(n0:nx,mxtrm),delp(n0:nx,mxtrm),
     &                delli(n0:nx,mxtrm),drvsav(n0:nx,16)
      common /TERMS/ phi01(mxtrm),phi10(mxtrm),phi02(mxtrm),
     &               phi20(mxtrm),phi03(mxtrm),phi30(mxtrm)
      common /FEQSV3/ crtflg
      common /ALTFEQ/ igenfl(n0:nx)
c  common storing the fluid constants
      common /CCON/ tcc(n0:nx),pcc(n0:nx),rhocc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      logical crtflg(n0:nx)
c
c     write (*,*) ' PHIFEQ--enter w/ tau,del = ',tau,del
c
      phifeq=0.0d0
      if (del.le.1.0d-10) RETURN  !trivial solution at zero density for
      if (tau.le.0.0d0) RETURN    !  any and all derivatives
c
      ncode=idel*4+itau+1
      if (abs(tau-tausav(icomp)).lt.1.0d-12 .and.
     &    abs(del-delsav(icomp)).lt.1.0d-16) then
c  retrieve value from previous call
        if (drvflg(icomp)(ncode:ncode).eq.'1') then
          phifeq=drvsav(icomp,ncode)
          RETURN
        endif
c       write (*,*) ' PHIFEQ--using phisav itau,idel,tau:',itau,idel,tau
      else
c  otherwise, compute new values and save for possible future use
c  first compute needed powers of tau and del (and save for future use)
        drvflg(icomp)='0000000000000000'
        if (abs(tau-tausav(icomp)).gt.1.0d-12) then
          elntau=log(tau)
          tausav(icomp)=tau
          do j=1,ntp(icomp)
            taup(icomp,j)=tpower(icomp,j)*elntau
c         write (*,*) ' PHIFEQ--tau,tpower: ',tau,tpower(icomp,j)
          enddo
        end if
        if (abs(del-delsav(icomp)).gt.1.0d-16) then
          elndel=log(del)
          delsav(icomp)=del
          do j=1,ndp(icomp)
            delp(icomp,j)=dpower(icomp,j)*elndel
          enddo
          d2=del*del
          d3=del*d2
          d4=del*d3
          do j=1,nlp(icomp)
            if (ABS(dlpowr(icomp,j)).lt.1.d-20) then
              delli(icomp,j)=0.0d0
            else
              if (ABS(dlpowr(icomp,j)-1.0d0).lt.1.d-20) then
                delli(icomp,j)=gi2(icomp,j)*del
              elseif (ABS(dlpowr(icomp,j)-2.0d0).lt.1.d-20) then
                delli(icomp,j)=gi2(icomp,j)*d2
              elseif (ABS(dlpowr(icomp,j)-3).lt.1.d-20) then
                delli(icomp,j)=gi2(icomp,j)*d3
              elseif (ABS(dlpowr(icomp,j)-4).lt.1.d-20) then
                delli(icomp,j)=gi2(icomp,j)*d4
              else
                delli(icomp,j)=gi2(icomp,j)*del**dlpowr(icomp,j)
              endif
            end if
          enddo
c         write (*,*) ' PHIFEQ--del,dpower: ',del,dpower(icomp,j)
c         write (*,*) ' PHIFEQ--del,dlpowr: ',del,dlpowr(icomp,j)
        end if
c
c  check for presence of critical-region terms of Wagner (e.g. steam)
c
        if (ncrt(icomp).gt.0) then
          do i=1,ncrt(icomp)
          k=i+ntermf(icomp)
          ext(icomp,i)=1.0d0
          extd(icomp,i)=0.0d0
          extd2(icomp,i)=0.0d0
          extdt(icomp,i)=0.0d0
          extt(icomp,i)=0.0d0
          extt2(icomp,i)=0.0d0
          if (abs(eta(icomp,i)).lt.1.d-12) then
            delb(icomp,i)=del-delta(icomp,i)
            taua(icomp,i)=tau-gamma(icomp,i)
            if (dli(icomp,k).eq.2 .and. tli(icomp,i).eq.2) then
              dxp(icomp,i)=alpha(icomp,i)*delb(icomp,i)
     &                   **INT(dli(icomp,k)-1.d0)
              axp(icomp,i)= beta(icomp,i)*taua(icomp,i)
     &                   **INT(tli(icomp,i)-1.d0)
            else
              dxp(icomp,i)=alpha(icomp,i)*delb(icomp,i)
     &                   **(dli(icomp,k)-1.d0)
              axp(icomp,i)= beta(icomp,i)*taua(icomp,i)
     &                   **(tli(icomp,i)-1.d0)
            endif
            hxp(icomp,i)=dxp(icomp,i)*delb(icomp,i)
            txp(icomp,i)=axp(icomp,i)*taua(icomp,i)
          else
            del1=del-1.0d0
            tau1=1.0d0-tau
            if (ABS(del1).lt.1.0d-10) del1=SIGN(1.0d-10,del1)
            if (ABS(tau1).lt.1.0d-10) tau1=SIGN(1.0d-10,tau1)
            delb(icomp,i)=del1
            taua(icomp,i)=tau1
            tdg=tau1+gamma(icomp,i)*ABS(del1)**(1.0d0/beta(icomp,i))
            ext(icomp,i)=(tdg**2+eid(icomp,i)*(del1**2)**eit(icomp,i))
            s=SIGN(1.0d0,del1)
            extd(icomp,i)=2.0d0*tdg*gamma(icomp,i)/beta(icomp,i)
     &        *ABS(del1)**(1.0d0/beta(icomp,i)-1.0d0)*s
     &        +2.0d0*eid(icomp,i)*eit(icomp,i)
     &        *ABS(del1)**(2.0d0*eit(icomp,i)-1.0d0)*s
            extd2(icomp,i)=alpha(icomp,i)*((alpha(icomp,i)-1.0d0)
     &        *ext(icomp,i)**(alpha(icomp,i)-2.0d0)*extd(icomp,i)**2
     &        +ext(icomp,i)**(alpha(icomp,i)-1.0d0)
     &        *(2.0d0*tdg*gamma(icomp,i)/beta(icomp,i)
     &        *(1.0d0/beta(icomp,i)-1.0d0)
     &        *ABS(del1)**(1.0d0/beta(icomp,i)-2.0d0)
     &        +2.0d0*(gamma(icomp,i)/beta(icomp,i)
     &        *ABS(del1)**(1.0d0/beta(icomp,i)-1.0d0))**2
     &        +2.0d0*eid(icomp,i)*eit(icomp,i)
     &        *(2.0d0*eit(icomp,i)-1.0d0)
     &        *ABS(del1)**(2.0d0*eit(icomp,i)-2.0d0)))
            extdt(icomp,i)=-2.0d0*gamma(icomp,i)/beta(icomp,i)
     &        * ABS(del1)**(1.0d0/beta(icomp,i)-1.0d0)*s
            extt(icomp,i)=-2.0d0*tdg
            extt2(icomp,i)=alpha(icomp,i)*((alpha(icomp,i)-1.0d0)
     &        *ext(icomp,i)**(alpha(icomp,i)-2.0d0)*extt(icomp,i)**2
     &        +2.0d0*ext(icomp,i)**(alpha(icomp,i)-1.0d0))
            extdt(icomp,i)=alpha(icomp,i)*((alpha(icomp,i)-1.0d0)
     &        *ext(icomp,i)**(alpha(icomp,i)-2.0d0)
     &        *extt(icomp,i)*extd(icomp,i)
     &        +ext(icomp,i)**(alpha(icomp,i)-1.0d0)*extdt(icomp,i))
            extd(icomp,i)=alpha(icomp,i)
     &        *ext(icomp,i)**(alpha(icomp,i)-1.0d0)*extd(icomp,i)
            extt(icomp,i)=alpha(icomp,i)
     &        *ext(icomp,i)**(alpha(icomp,i)-1.0d0)*extt(icomp,i)
            ext(icomp,i)=ext(icomp,i)**alpha(icomp,i)
            hxp(icomp,i)=-delta(icomp,i)*del1**(INT(dli(icomp,k)))
            txp(icomp,i)=-eta(icomp,i)*tau1**(INT(dli(icomp,k)))
            taua(icomp,i)=-taua(icomp,i)
          endif
          enddo
        endif
c
c  end critical-region terms
c
        phisum=0.0d0
        do k=1,ntermf(icomp)
          ex=taup(icomp,itp(icomp,k))
     &      +delp(icomp,idp(icomp,k))-delli(icomp,ilp(icomp,k))
          if (ex.lt.100.d0 .and. ex.gt.-200.d0) then
            phisav(icomp,k)=a(icomp,k)*EXP(ex)
          else
            phisav(icomp,k)=0.d0
          endif
c         if (igenfl(icomp).eq.1) then
c           if (k.le.14) then
c             ! do nothing more to these terms
c           elseif (k.gt.14 .and. k.lt.29) then
c             ! multiply by acentric factor
c             phisav(icomp,k)=phisav(icomp,k)*accen(icomp)
c           else
c             ! multiply by theta factor
c             phisav(icomp,k)=phisav(icomp,k)*(zcrit(icomp)-0.29d0)**2
c           endif
c         endif
          phisum=phisum+phisav(icomp,k)
c
c       write (*,1010) k,phisav(icomp,k),phisum
c1010   format (1x,i3,2d30.20)   !write out each term for debugging
        enddo
        if (ncrt(icomp).gt.0) then
          do i=1,ncrt(icomp)
            k=i+ntermf(icomp)
            ex=taup(icomp,itp(icomp,k))
     &        +delp(icomp,idp(icomp,k))+hxp(icomp,i)+txp(icomp,i)
            if (.not.crtflg(icomp)) ex=ex+log(ext(icomp,i))
            if (ex.lt.100.d0 .and. ex.gt.-200.d0) then
              phisav(icomp,k)=a(icomp,k)*EXP(ex)
            else
              phisav(icomp,k)=0.d0
            endif
            phisum=phisum+phisav(icomp,k)
          enddo
        endif
        phifeq=phisum
        drvflg(icomp)(1:1)='1'
        drvsav(icomp,1)=phifeq
      end if
c
c     write (*,1012) tau,del,itau,idel,phifeq
c1012 format (1x,' PHIFEQ--tau,del,itau,idel,phi00: ',2f10.5,2i4,e22.12)
c
c  check if derivatives are requested, calculations make use of fact
c  that terms in derivative summations are very similar to A/RT terms
c
      if (idel.eq.1) then
c  compute derivative w.r.t. del (dimensionless density)
c  save individual terms for possible use in cross derivative
        phisum=0.0d0
        do k=1,ntermf(icomp)
          phi01(k)=phisav(icomp,k)
     &            *(di(icomp,k)-dli(icomp,k)*delli(icomp,ilp(icomp,k)))
          phisum=phisum+phi01(k)
        enddo
c  check for presence of critical-region terms
        if (ncrt(icomp).gt.0) then
          do i=1,ncrt(icomp)
            k=i+ntermf(icomp)
            if (abs(eta(icomp,i)).lt.1.d-12) then
              phi01(k)=phisav(icomp,k)
     &                *(di(icomp,k)+dli(icomp,k)*del*dxp(icomp,i))
            else
              phi01(k)=phisav(icomp,k)/ext(icomp,i)
     &                *(ext(icomp,i)*(di(icomp,k)
     &                +dli(icomp,k)*del*hxp(icomp,i)/delb(icomp,i))
     &                +extd(icomp,i)*del)
            endif
            phisum=phisum+phi01(k)
          enddo
        endif
        phifeq=phisum
c
      else if (idel.eq.2) then
c  compute 2nd derivative w.r.t. del (dimensionless density)
c  save individual terms for possible use in cross derivative
        phisum=0.0d0
        do k=1,ntermf(icomp)
          dik=di(icomp,k)
          dell=dli(icomp,k)*delli(icomp,ilp(icomp,k))
          phi02(k)=phisav(icomp,k)
     &            *((dik-dell)*(dik-1.d0-dell)-dli(icomp,k)*dell)
          phisum=phisum+phi02(k)
        enddo
c  check for presence of critical-region terms
        if (ncrt(icomp).gt.0) then
          do i=1,ncrt(icomp)
            k=i+ntermf(icomp)
            dik=di(icomp,k)
            if (abs(eta(icomp,i)).lt.1.d-12) then
              dl=dli(icomp,k)
              dlb=1.d0
              if (dl.ne.2.d0) dlb=delb(icomp,i)**(dl-2.d0)
              phi02(k)=phisav(icomp,k)*((dik+dl*del*dxp(icomp,i))**2
     &                -dik-(dl-dl**2)*del**2*alpha(icomp,i)*dlb)
            else
              d2=del**2
              dell=dli(icomp,k)*hxp(icomp,i)/delb(icomp,i)
              phi02(k)=phisav(icomp,k)/ext(icomp,i)
     &            *(d2*(extd2(icomp,i)+extd(icomp,i)*(dik/del+dell)
     &            +ext(icomp,i)*(-dik/d2+dell*(dli(icomp,k)-1.0d0)
     &            /delb(icomp,i)))
     &            +(extd(icomp,i)+ext(icomp,i)*(dik/del+dell))
     &            *(dik*del+d2*dell))
            endif
            phisum=phisum+phi02(k)
          enddo
        endif
        phifeq=phisum
c
      else if (idel.eq.3) then
c  compute 3rd derivative w.r.t. del (dimensionless density)
        phisum=0.0d0
        do k=1,ntermf(icomp)
          dik=di(icomp,k)
          dell=dli(icomp,k)*delli(icomp,ilp(icomp,k))
          dl=dli(icomp,k)
          phi03(k)=phisav(icomp,k)
     &     *(dik*(2.d0-3.d0*dik+dik**2)
     &     +dell*(-2.d0+6.d0*dik-3.d0*dik**2-3.d0*dik*dl+3.d0*dl-dl**2
     &     +dell*(3.d0*(dl+dik-1.d0)-dell)))
          phisum=phisum+phi03(k)
        enddo
c  check for presence of critical-region terms, note that the extended
c  terms (non-Gaussian) in the critical region are not included here.
        if (ncrt(icomp).gt.0) then
          do i=1,ncrt(icomp)
            if (abs(eta(icomp,i)).lt.1.d-12) then
              k=i+ntermf(icomp)
              dik=di(icomp,k)
              dl=dli(icomp,k)
              dlb=del**2*alpha(icomp,i)
              dlb2=0.d0
              if (dl.ne.2.d0) then
                dlb2=-dlb*del*(3.d0*dl**2-dl**3-2.d0*dl)
     &                *delb(icomp,i)**(dl-3.d0)
                dlb=dlb*delb(icomp,i)**(dl-2.d0)
              endif
              phi03(k)=phisav(icomp,k)*((dik+dl*del*dxp(icomp,i))**3
     &                -3.d0*dik**2+2.d0*dik+dlb2
     &                -3.d0*(dl-dl**2)*dik*dlb
     &                -3.d0*del*dxp(icomp,i)
     &              *(dik*dl-(dl**3-dl**2)*dlb))
            else
              phi03(k)=0.d0
            endif
            phisum=phisum+phi03(k)
          enddo
        endif
        phifeq=phisum
      end if
c
c
      if (itau.eq.1) then
c  compute derivative w.r.t. tau (dimensionless temperature)
c  save individual terms for possible use in cross derivative
        phisum=0.0d0
        do k=1,ntermf(icomp)
          phi10(k)=phisav(icomp,k)*ti(icomp,k)
          phisum=phisum+phi10(k)
        enddo
c  check for presence of critical-region terms
        if (ncrt(icomp).gt.0) then
          do i=1,ncrt(icomp)
            k=i+ntermf(icomp)
            if (abs(eta(icomp,i)).lt.1.d-12) then
              phi10(k)=phisav(icomp,k)
     &                *(ti(icomp,k)+tli(icomp,i)*tau*axp(icomp,i))
            else
              phi10(k)=phisav(icomp,k)*(ti(icomp,k)
     &              +tau*tli(icomp,i)*txp(icomp,i)/taua(icomp,i)
     &              +tau*extt(icomp,i)/ext(icomp,i))
            endif
            phisum=phisum+phi10(k)
          enddo
        endif
        phifeq=phisum
c
      else if (itau.eq.2) then
c  compute 2nd derivative w.r.t. tau (dimensionless temperature)
c  save individual terms for possible use in cross derivative
        phisum=0.0d0
        do k=1,ntermf(icomp)
          tik=ti(icomp,k)
          phi20(k)=phisav(icomp,k)*tik*(tik-1.0d0)
          phisum=phisum+phi20(k)
        enddo
c  check for presence of critical-region terms
        if (ncrt(icomp).gt.0) then
          do i=1,ncrt(icomp)
            k=i+ntermf(icomp)
            tik=ti(icomp,k)
            if (abs(eta(icomp,i)).lt.1.d-12) then
              tl=tli(icomp,i)
              tlb=1.d0
              if (tl.ne.2.d0) tlb=taua(icomp,i)**(tl-2.d0)
              phi20(k)=phisav(icomp,k)*((tik+tl*tau*axp(icomp,i))**2
     &                -tik-(tl-tl**2)*tau**2*beta(icomp,i)*tlb)
            else
              t2=tau**2
              tauu=tli(icomp,i)*txp(icomp,i)/taua(icomp,i)
              phi20(k)=phisav(icomp,k)*(t2
     &            *(extt2(icomp,i)+extt(icomp,i)*(tik/tau+tauu)
     &            +ext(icomp,i)*(-tik/t2+tauu*(tli(icomp,i)-1.0d0)
     &            /taua(icomp,i)))
     &            +(extt(icomp,i)+ext(icomp,i)*(tik/tau+tauu))
     &            *(tik*tau+t2*tauu))/ext(icomp,i)
            endif
            phisum=phisum+phi20(k)
          enddo
        endif
        phifeq=phisum
c
      else if (itau.eq.3) then
c  compute 3rd derivative w.r.t. tau (dimensionless temperature)
        phisum=0.0d0
        do k=1,ntermf(icomp)
          tik=ti(icomp,k)
          phisum=phisum+phisav(icomp,k)*tik*(tik-1.d0)*(tik-2.d0)
        enddo
        if (ncrt(icomp).gt.0) then
          do i=1,ncrt(icomp)
            if (abs(eta(icomp,i)).lt.1.d-12) then
              k=i+ntermf(icomp)
              tik=ti(icomp,k)
              tl=tli(icomp,i)
              tlb=tau**2*beta(icomp,i)
              tlb2=0.d0
              if (tl.ne.2.d0) then
                tlb2=-tlb*tau*(3.d0*tl**2-tl**3-2.d0*tl)
     &                *taua(icomp,i)**(tl-3.d0)
                tlb=tlb*taua(icomp,i)**(tl-2.d0)
              endif
              phi30(k)=phisav(icomp,k)*((tik+tl*tau*axp(icomp,i))**3
     &                -3.d0*tik**2+2.d0*tik+tlb2
     &                -3.d0*(tl-tl**2)*tik*tlb
     &                -3.d0*tau*axp(icomp,i)
     &              *(tik*tl-(tl**3-tl**2)*tlb))
              phisum=phisum+phi30(k)
            else
              phi30(k)=0.d0
            endif
          enddo
        endif
        phifeq=phisum
      end if
c
c
      if (itau.eq.1 .and. idel.eq.1) then
c  compute cross derivative using terms from 1st derivatives
        phisum=0.0d0
        do k=1,ntermf(icomp)
          if (ABS(phisav(icomp,k)).gt.1.0d-20)
     &      phisum=phisum+phi10(k)*phi01(k)/phisav(icomp,k)
        enddo
        if (ncrt(icomp).gt.0) then
          do i=1,ncrt(icomp)
            k=i+ntermf(icomp)
            if (abs(eta(icomp,i)).lt.1.d-12) then
              if (ABS(phisav(icomp,k)).gt.1.0d-20)
     &          phisum=phisum+phi10(k)*phi01(k)/phisav(icomp,k)
            else
              phisum=phisum+(ti(icomp,k)+tau*tli(icomp,i)
     &           *txp(icomp,i)/taua(icomp,i))*phi01(k)
     &           +phisav(icomp,k)/ext(icomp,i)*(del*tau*extdt(icomp,i)
     &           +tau*extt(icomp,i)*(di(icomp,k)
     &           +del*dli(icomp,k)*hxp(icomp,i)/delb(icomp,i)))
            endif
          enddo
        endif
        phifeq=phisum
c
      else if (itau.eq.2 .and. idel.eq.1) then
c  compute cross derivative using term from 1st derivative
        phisum=0.0d0
        do k=1,ntermf(icomp)+ncrt(icomp)
          if (ABS(phisav(icomp,k)).gt.1.0d-20)
     &      phisum=phisum+phi20(k)*phi01(k)/phisav(icomp,k)
        enddo
        phifeq=phisum
c
      else if (itau.eq.1 .and. idel.eq.2) then
c  compute cross derivative using term from 2nd derivative
        phisum=0.0d0
        do k=1,ntermf(icomp)+ncrt(icomp)
          if (ABS(phisav(icomp,k)).gt.1.0d-20)
     &      phisum=phisum+phi10(k)*phi02(k)/phisav(icomp,k)
        enddo
        phifeq=phisum
c
      else if (itau.eq.2 .and. idel.eq.2) then
c  compute cross derivative using terms from 2nd derivative
        phisum=0.0d0
        do k=1,ntermf(icomp)+ncrt(icomp)
          if (ABS(phisav(icomp,k)).gt.1.0d-20)
     &      phisum=phisum+phi20(k)*phi02(k)/phisav(icomp,k)
        enddo
        phifeq=phisum
      end if
c
c     write (*,1021) tau,del,itau,idel,phifeq
c1021 format (1x,' PHIFEQ--tau,del,itau,idel,phixx: ',2f10.5,2i4,e22.12)
c
      drvsav(icomp,ncode)=phifeq
      drvflg(icomp)(ncode:ncode)='1'
c
      RETURN
      end                                               !function PHIFEQ
c
c ======================================================================
c
      subroutine CRTFEQ (icomp,tcrit,pcrit,Dcrit)
c
c  returns critical parameters associated with Fundamental EOS
c
c  N.B.  these critical parameters may not necessarily be most
c        accurate values, but they are consistent with the EOS fit;
c        neither are they always equal to the reducing parameters
c
c  input:
c    icomp--pointer specifying component (1..nc)
c  outputs:
c    tcrit--critical temperature (K)
c    pcrit--critical pressure (kPa)
c    Dcrit--molar density (mol/L) at critical point
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  08-03-95  MM, original version
c  11-01-95  MM, increase parameter mxtrm to 52 (to accommodate steam)
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  01-16-96  MM, implicit integer (i-n); (include L)
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c  01-03-97  MM, new arrays added to /WCFFEQ/ (assoc w/ sorting of powers)
c  07-02-97 EWL, modify /WCFFEQ/ to accommodate g parameter
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
c
      parameter (mxtrm=72)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFFEQ/ a(n0:nx,mxtrm),ti(n0:nx,mxtrm),di(n0:nx,mxtrm),
     &                gi(n0:nx,mxtrm),gi2(n0:nx,mxtrm),
     &                dli(n0:nx,mxtrm),tli(n0:nx,mxtrm),
     &                tpower(n0:nx,mxtrm),dpower(n0:nx,mxtrm),
     &                dlpowr(n0:nx,mxtrm),
     &                rho0(n0:nx),t0(n0:nx),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),
     &                wmf(n0:nx),Rfeq(n0:nx),
     &                pmin(n0:nx),rhotp(n0:nx),tmin(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
c
      tcrit=tc(icomp)
      pcrit=pc(icomp)
      Dcrit=rhoc(icomp)
c
      RETURN
      end                                             !subroutine CRTFEQ
c
c ======================================================================
c
      subroutine REDFEQ (icomp,tred,Dred)
c
c  returns reducing parameters associated with Fundamental EOS;
c  used to calculate the 'tau' and 'del' which are the independent
c  variables in the EOS
c
c  N.B.  The reducing parameters are often, but not always, equal
c        to the critical temperature and density.
c
c  input:
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c  outputs:
c     tred--reducing temperature (K)
c     Dred--reducing molar density (mol/L)
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  08-03-95  MM, original version
c  11-01-95  MM, increase parameter mxtrm to 52 (to accommodate steam)
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  01-16-96  MM, implicit integer (i-n); (include L)
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c  01-03-97  MM, new arrays added to /WCFFEQ/ (assoc w/ sorting of powers)
c  07-02-97 EWL, modify /WCFFEQ/ to accommodate g parameter
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
c
      parameter (mxtrm=72)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFFEQ/ a(n0:nx,mxtrm),ti(n0:nx,mxtrm),di(n0:nx,mxtrm),
     &                gi(n0:nx,mxtrm),gi2(n0:nx,mxtrm),
     &                dli(n0:nx,mxtrm),tli(n0:nx,mxtrm),
     &                tpower(n0:nx,mxtrm),dpower(n0:nx,mxtrm),
     &                dlpowr(n0:nx,mxtrm),
     &                rho0(n0:nx),t0(n0:nx),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),
     &                wmf(n0:nx),Rfeq(n0:nx),
     &                pmin(n0:nx),rhotp(n0:nx),tmin(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
c
c     write (*,*) ' REDFEQ, i,t0,rho0: ',icomp,t0(icomp),rho0(icomp)
      tred=t0(icomp)
      Dred=rho0(icomp)
c
      RETURN
      end                                             !subroutine REDFEQ
c
c ======================================================================
c
      subroutine SETFEQ (nread,icomp,hcasno,ierr,herr)
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
c     other quantities returned via arrays in common /WCFFEQ/
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  07-20-95  MM, original version
c  09-13-95  MM, add ierr, herr to argument list
c  10-05-95  MM, adapt to file input, add nread to argument list
c  11-01-95  MM, increase parameter mxtrm to 52 (to accommodate steam)
c  11-02-95  MM, add common /CCAS/ to access CAS numbers
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  01-16-96  MM, implicit integer (i-n); (include L)
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                add Zcrit to common /CCON/
c  03-19-96  MM, add dipole moment to /CCON/ and /MSCFEQ/
c  06-03-96  MM, add limits to /EOSLIM/, reduce mxfeq from 20 to 2
c  01-03-97  MM, sort powers of tau, del, store in tpower, etc.
c                new arrays added to /WNTFEQ/,/WLFFEQ/,/WCFFEQ/
c  01-22-97 EWL, add variables for critical region for methane and water.
c  01-31-97  MM, modify commons associated with critical-region terms
c  02-06-97  MM, if steam, reset gas constant to fluid-specific value
c  02-07-97  MM, fix bug in reading of critical terms
c  05-27-97  MM, if nc = 1, set R to fluid-specific value
c  07-02-97 EWL, add parameter for g in exp(-g*d^l) for fluids where g<>1
c                modify /WCFFEQ/ to accommodate g parameter
c  11-13-97  MM, (re)initialize contents of /FEQSAV/ when a new fluid is read in
c  02-11-98  MM, store rho at triple point separate from rhomax
c  04-06-98 EWL, add check for change in exponential coefficient (gi)
c  04-06-98 EWL, change mxcrt to 8
c  04-06-98 EWL, read in powers of del and tau in exponential for critical terms
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-22-98 EWL, set Reos to Rfeq
c   4-24-00 EWL, change the exponent of del in the EXP term from integer to
c                real*8 and modify the common blocks WCFFEQ and WLFFEQ
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (mxfeq=2)         !max number of FEQ EOS in block data
      parameter (mxtrm=72,mxcrt=20)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*3 hcpfeq
      character*3 hpheq,heos,hmxeos,hmodcp
      character*12 hcas,hcasa,hcasno
      character*255 herr
c     character*1 dummy
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c  commons associated with the mxfeq fluids with FEQ equations stored
c  in block data BDFEQ
      common /CASFEQ/ hcasa(mxfeq)
      common /CPMFEQ/ hcpfeq(mxfeq)
      common /NTFEQ/ ntrmfa(mxfeq)
      common /LFFEQ/ dlia(mxfeq,mxtrm)
      common /CFFEQ/ aa(mxfeq,mxtrm),dia(mxfeq,mxtrm),tia(mxfeq,mxtrm),
     &               rho0a(mxfeq),t0a(mxfeq),
     &               pca(mxfeq),rhoca(mxfeq),tca(mxfeq),
     &               wmfa(mxfeq),Rfeqa(mxfeq),
     &               pmina(mxfeq),rhomxa(mxfeq),tmina(mxfeq),
     &               tmaxa(mxfeq),pmaxa(mxfeq)
      common /MSCFEQ/ ttpf(mxfeq),tnbpf(mxfeq),accenf(mxfeq),dipm(mxfeq)
c  commons associated with the nc components of current interest
c  ("working" commons and arrays)
      common /CCAS/ hcas(n0:nx)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c     common /CPMOD/ hmodcp(n0:nx)
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /WNTFEQ/ ntermf(n0:nx),ncoeff(n0:nx),
     &                ntp(n0:nx),ndp(n0:nx),nlp(n0:nx),
     &                ncrt(n0:nx),ncfcrt(n0:nx),
     &                nspare(n0:nx),ncfsp(n0:nx)
      common /WLFFEQ/ itp(n0:nx,mxtrm),
     &                idp(n0:nx,mxtrm),ilp(n0:nx,mxtrm)
      common /WCFFEQ/ a(n0:nx,mxtrm),ti(n0:nx,mxtrm),di(n0:nx,mxtrm),
     &                gi(n0:nx,mxtrm),gi2(n0:nx,mxtrm),
     &                dli(n0:nx,mxtrm),tli(n0:nx,mxtrm),
     &                tpower(n0:nx,mxtrm),dpower(n0:nx,mxtrm),
     &                dlpowr(n0:nx,mxtrm),
     &                rho0(n0:nx),t0(n0:nx),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),
     &                wmf(n0:nx),Rfeq(n0:nx),
     &                pmin(n0:nx),rhotp(n0:nx),tmin(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
c  parameters associated with the critical-region terms of Wagner
      common /WCFFQ2/ alpha(n0:nx,mxcrt),beta(n0:nx,mxcrt),
     &                gamma(n0:nx,mxcrt),delta(n0:nx,mxcrt),
     &                eta(n0:nx,mxcrt),eid(n0:nx,mxcrt),eit(n0:nx,mxcrt)
c  limits associated with the equation of state
      common /EOSLIM/ tmn(n0:nx),tmx(n0:nx),pmx(n0:nx),rhomx(n0:nx)
      common /FEQSAV/ phisav(n0:nx,mxtrm),delsav(n0:nx),tausav(n0:nx),
     &                taup(n0:nx,mxtrm),delp(n0:nx,mxtrm),
     &                delli(n0:nx,mxtrm),drvsav(n0:nx,16)
      common /FEQSV3/ crtflg
      logical crtflg(n0:nx)
c
      ierr=0
      crtflg(icomp)=.true.
      herr=' '
c  (re)initialize contents of /FEQSAV/ when a new fluid is read in
      do i=n0,nx
        delsav(i)=0.0d0
        tausav(i)=0.0d0
        do j=1,mxtrm
          phisav(i,j)=0.0d0
          taup(i,j)=0.0d0
          delp(i,j)=0.0d0
          delli(i,j)=0.0d0
        enddo
      enddo
c
      if (nread.le.0) then
c  get coefficients from block data
c  identify specified fluid with entries in database via match of CAS no
        do i=1,mxfeq
        if (hcasno.eq.hcasa(i)) then
C          hmodcp(icomp)=hcpfeq(i)   ! pointer to ideal gas model
C          ntermf(icomp)=ntrmfa(i)   ! number of terms
C          do j=1,ntermf(icomp)
C          dli(icomp,j)=dlia(i,j)      ! transfer coefficients
C          a(icomp,j)=aa(i,j)        ! into working arrays
C          di(icomp,j)=dia(i,j)
C          ti(icomp,j)=tia(i,j)
C          enddo
C          rho0(icomp)=rho0a(i)      ! reducing parameters
C          t0(icomp)=t0a(i)
C          rhoz(icomp)=rho0a(i)      ! reducing parameters
C          tz(icomp)=t0a(i)
C          pc(icomp)=pca(i)          ! critical parameters
C          rhoc(icomp)=rhoca(i)
C          tc(icomp)=tca(i)
C          wmf(icomp)=wmfa(i)        ! molecular weight
C          Rfeq(icomp)=Rfeqa(i)      ! gas constant
C          if (nc.eq.1 .and. icomp.eq.1) then
C            R=Rfeq(icomp)
C          end if
C          Reos(icomp)=Rfeq(icomp)
C          pmin(icomp)=pmina(i)      ! limits
C          rhomx(icomp)=rhomxa(i)
C          rhotp(icomp)=rhomxa(i)
C          tmin(icomp)=tmina(i)
C          tmax(icomp)=tmaxa(i)
C          pmax(icomp)=pmaxa(i)
C          ncoeff(icomp)=4           !# coeff per term
C          ncrt(icomp)=0             !no critical-region terms
C          nspare(icomp)=0
Cc  fill arrays in /CCON/
C          wm(icomp)=wmfa(i)
C          ttp(icomp)=ttpf(i)
C          tnbp(icomp)=tnbpf(i)
C          tcrit(icomp)=tc(icomp)
C          pcrit(icomp)=pc(icomp)
C          Dcrit(icomp)=rhoc(icomp)
C          Zcrit(icomp)=pc(icomp)/(Rfeq(icomp)*tc(icomp)*rhoc(icomp))
C          accen(icomp)=accenf(i)
C          dipole(icomp)=dipm(i)
C          goto 990
        end if
        enddo
        ierr=1
        herr='[SETFEQ error] fluid input to SETFEQ not found'//hnull
      else
c  read data from file
c       write (*,*) ' SETFEQ--read component',icomp,' from unit',nread
        read (nread,*) tmin(icomp)          !lower temperature limit
c       write (*,*) ' SETFEQ--first input: tmin: ',tmin(icomp)
        read (nread,*) tmax(icomp)          !upper temperature limit
        read (nread,*) pmax(icomp)          !upper pressure limit
        read (nread,*) rhomx(icomp)         !upper density limit
        read (nread,2003) hmodcp(icomp)     !pointer to Cp0 model
        read (nread,*) wm(icomp)     !molecular weight
        wmf(icomp)=wm(icomp)
        read (nread,*) ttp(icomp)    !triple point temperature
        read (nread,*) pmin(icomp)   !pressure at triple point
        read (nread,*) rhotp(icomp)  !density at triple point
        read (nread,*) tnbp(icomp)   !normal boiling point temperature
        read (nread,*) accen(icomp)  !acentric factor
        read (nread,*) tc(icomp),pc(icomp),rhoc(icomp) !critical par
        tcrit(icomp)=tc(icomp)
        pcrit(icomp)=pc(icomp)
        Dcrit(icomp)=rhoc(icomp)
        ptp(icomp)=pmin(icomp)
        dtp(icomp)=rhotp(icomp)
        dtpv(icomp)=0.0d0
        dnbpl(icomp)=0.0d0
        dnbpv(icomp)=0.0d0
        read (nread,*) t0(icomp),rho0(icomp) !reducing parameters
        tz(icomp)=t0(icomp)
        rhoz(icomp)=rho0(icomp)
        read (nread,*) Rfeq(icomp)    !gas constant used in fit
        if (nc.eq.1 .and. icomp.eq.1) then
          R=Rfeq(icomp)
c         write (*,*) ' SETFEQ--R set to ',R
        end if
        Reos(icomp)=Rfeq(icomp)
        Zcrit(icomp)=pc(icomp)/(Rfeq(icomp)*tc(icomp)*rhoc(icomp))
        read (nread,*) nterm,ncoeff(icomp),ncrt(icomp),ncfcrt(icomp),
     &                 nspare(icomp),ncfsp(icomp)
        ntermf(icomp)=nterm
c       write (*,*) ' SETFEQ--about to read ',nterm,' coefficients'
c  the gi term is a multiplier for the (rho or del) in only the exponential
c  terms; it is needed for e.g. the Bender EOS; set to 1.0 if not present
        do j=1,nterm
          if (ncoeff(icomp).eq.5) then
          read (nread,*)a(icomp,j),ti(icomp,j),di(icomp,j),dli(icomp,j),
     &                  gi(icomp,j)
          else
          read (nread,*) a(icomp,j),ti(icomp,j),di(icomp,j),dli(icomp,j)
          gi(icomp,j)=1.d0
          endif
        enddo
        if (ncrt(icomp).gt.0) then
        do j=1,ncrt(icomp)
          gi(icomp,j+nterm)=0.d0
          read (nread,*) a(icomp,j+nterm),ti(icomp,j+nterm),
     &       di(icomp,j+nterm),dli(icomp,j+nterm),tli(icomp,j),
     &       alpha(icomp,j),beta(icomp,j),gamma(icomp,j),delta(icomp,j),
     &       eta(icomp,j),eid(icomp,j),eit(icomp,j)
          if (abs(eta(icomp,j)).gt.1.d-20 .and.
     &        abs(eid(icomp,j)).gt.1.d-20 .and.
     &        abs(eit(icomp,j)).gt.1.d-20) crtflg(icomp)=.false.
        enddo
        end if
c       write (*,*) ' SETFEQ--final coefficient: ',a(icomp,nterm)
      end if
c
c 990 continue
c
      call SETEXP(icomp)
c
c  write out all coefficients for debugging
c
c     i=icomp
c     write (*,*) ' SETFEQ--coefficients for comp',i,', CAS # ',hcas(i)
c     write (*,*) ' Cp0 model: ',hmodcp(i)
c     write (*,*) ' Nterms :   ',ntermf(i)
c     write (*,2020) ((a(i,j),ti(i,j),di(i,j),dli(i,j)),j=1,ntermf(i))
c2020 format (1x,d24.14,2f9.3,i4)
c
c     write (*,*) ' SETFEQ--sorted coefs for comp',i,', CAS # ',hcas(i)
c     write (*,2020) ((a(i,j),tpower(i,itp(i,j)),dpower(i,idp(i,j)),
c    &               dlpowr(i,ilp(i,j))),j=1,ntermf(i))
c
c  copy limits into /EOSLIM/ arrays
      tmn(icomp)=tmin(icomp)
      tmx(icomp)=tmax(icomp)
      pmx(icomp)=pmax(icomp)
c     rhomx(icomp)=rhomax(icomp)
c
      RETURN
 2003 format (a3)
      end                                             !subroutine SETFEQ
c
c ======================================================================
c
      subroutine SETEXP (icomp)
c
c  scan and sort exponents to minimize computations;
c  the array tpower --stores the discrete powers of temperature
c                itp--is a pointer to the appropriate element of tpower
c        dpower, dtp--ditto for density
c        dlpowr, ltp--ditto for the powers of del in the exp terms
c  ntp, ndp, and nlp store the number of entries in the above arrays
c
c  inputs:
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  07-20-95  MM, original version
c   9-18-00 EWL, remove section of code from SETFEQ and place in new
c                subroutine SETEXP so that fitting routines can call it.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (mxtrm=72)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /WNTFEQ/ ntermf(n0:nx),ncoeff(n0:nx),
     &                ntp(n0:nx),ndp(n0:nx),nlp(n0:nx),
     &                ncrt(n0:nx),ncfcrt(n0:nx),
     &                nspare(n0:nx),ncfsp(n0:nx)
      common /WLFFEQ/ itp(n0:nx,mxtrm),
     &                idp(n0:nx,mxtrm),ilp(n0:nx,mxtrm)
      common /WCFFEQ/ a(n0:nx,mxtrm),ti(n0:nx,mxtrm),di(n0:nx,mxtrm),
     &                gi(n0:nx,mxtrm),gi2(n0:nx,mxtrm),
     &                dli(n0:nx,mxtrm),tli(n0:nx,mxtrm),
     &                tpower(n0:nx,mxtrm),dpower(n0:nx,mxtrm),
     &                dlpowr(n0:nx,mxtrm),
     &                rho0(n0:nx),t0(n0:nx),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),
     &                wmf(n0:nx),Rfeq(n0:nx),
     &                pmin(n0:nx),rhotp(n0:nx),tmin(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
c
      tpower(icomp,1)=ti(icomp,1)
      dpower(icomp,1)=di(icomp,1)
      dlpowr(icomp,1)=dli(icomp,1)
      gi2(icomp,1)=gi(icomp,1)
      ntp(icomp)=1
      ndp(icomp)=1
      nlp(icomp)=1
c
      do i=1,ntermf(icomp)+ncrt(icomp)
c  compare the power of temperature for each term with elements
c  already in the tpower array
        do j=1,ntp(icomp)
          if (abs(ti(icomp,i)-tpower(icomp,j)).lt.1.0d-10) then
            itp(icomp,i)=j
            goto 100
          end if
        enddo
c  power of temperature for term i does not match any previous term
        ntp(icomp)=ntp(icomp)+1
        tpower(icomp,ntp(icomp))=ti(icomp,i)
        itp(icomp,i)=ntp(icomp)
 100    continue
c  compare the power of density for each term with elements
c  already in the dpower array
        do j=1,ndp(icomp)
          if (abs(di(icomp,i)-dpower(icomp,j)).lt.1.0d-10) then
            idp(icomp,i)=j
            goto 110
          end if
        enddo
c  power of density for term i does not match any previous term
        ndp(icomp)=ndp(icomp)+1
        dpower(icomp,ndp(icomp))=di(icomp,i)
        idp(icomp,i)=ndp(icomp)
 110    continue
c  compare the power of density in exponential terms with elements
c  already in the dlpowr array
        do j=1,nlp(icomp)
          if (ABS(dli(icomp,i)-dlpowr(icomp,j)).lt.1.0d-10 .and.
     &        abs(gi(icomp,i)-gi2(icomp,j)).lt.1.0d-10) then
            ilp(icomp,i)=j
            goto 120
          end if
        enddo
c  exponential of density for term i does not match any previous term
        nlp(icomp)=nlp(icomp)+1
        dlpowr(icomp,nlp(icomp))=dli(icomp,i)
        ilp(icomp,i)=nlp(icomp)
        gi2(icomp,nlp(icomp))=gi(icomp,i)
 120    continue
      enddo
      end                                             !subroutine SETEXP
cc
cc ======================================================================
cc
c      block data BDFEQ
cc
cc  data for Helmholtz-explicit equations of state
cc
cc  explanation of parameters
cc     mxfeq:       maximum number of FEQ fits, used to dimension arrays
cc     mxtrm:       max number of terms per fit, used to dimension arrays
cc
cc  explanation of commons and constituent arrays
cc    /CASFEQ/      Chem Abstract number; used as unambiguous identifier
cc        hcas(i):  CAS number for fluid corresponding to equation "i"
cc
cc    /CPMFEQ/
cc      hmodcp(i):  pointer to ideal gas model to use with fluid "i"
cc
cc    /NTFEQ/
cc      ntermf(i):  number of terms in fit corresponding to fluid "i"
cc
cc    /LFFEQ/   parameters to FEQ fits for each of mxfeq fluids
cc      dli(i,1..N): power of del in the exponential multiplier
cc                  (if dli(j)=0, then multiplier = 1)
cc
cc    /CFFEQ/   parameters to FEQ fits for each of mxfeq fluids
cc      a(i,1..N):  N coefficients to "fundamental" equation of state
cc      di(i,1..N): exponents for the N density terms
cc      ti(i,1..N): exponents for the N temperature terms
cc      rho0(i):    reducing parameter for density (mol/L)
cc      t0(i):      reducing parameter for temperature (K)
cc      pc(i):      critical pressure (kPa)
cc      rhoc(i):    critical density (L/mol)
cc      tc(i):      critical temperature (K)
cc      wmf(i):     molecular weight (g/mol)
cc      Rfeq(i):    gas constant used in fit (J/(mol-K))
cc      pmin(i):    pressure at tmin(i) (kPa)
cc      rhotp(i):   density at tmin(i), e.g. triple point (L/mol)
cc      tmin(i):    low temperature limit of fit--often triple point
cc                  pmin(i),rhotp(i) are used for initial guesses, etc.
cc                  and are often approximate values only
cc      tmax(i):    upper temperature limit of fit (K)
cc      pmax(i):    upper pressure limit of fit (kPa)
cc
cc    /MSCFEQ/  miscellaneous fluid constants
cc      ttpf(i):    triple point temperature (K)
cc      tnbpf(i):   normal boiling point temperature (K)
cc      accenf(i):  acentric factor for fluid represented by eqn "i"
cc      dipm(i):    dipole moment [debye] (at Tnbp if t-dependent)
cc
cc    /FEQSAV/  used to save information between calls to PHIFEQ
cc      phisav(i,j):individual terms in summation
cc      delsav(i):  reduced density on last call to PHIFEQ
cc      tausav(i):  reduced temperature on last call to PHIFEQ
cc      taup(i,j):  reduced temperature raised to power
cc      delp(i,j):  reduced density raised to power
cc
cc    where "i" is the equation number
cc    N.B.--the "i" are, in general, not the same as fluid code numbers
cc
cc
cc  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
cc  07-26-95  MM, original version
cc  10-03-95  MM, /CPMFEQ/ changed to character variable
cc  11-01-95  MM, increase parameter mxtrm to 52 (to accommodate steam)
cc  11-29-95  MM, variable lower limit on coefficient/constant arrays
cc                to accommodate ECS reference fluid
cc  01-16-96  MM, implicit integer (i-n); (include L)
cc  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
cc  01-03-97  MM, arrays added to /FEQSAV/ (assoc with sorting of powers)
cc
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
c      parameter (ncmax=20)       !max number of components in mixture
c      parameter (nrefmx=10)      !max number of fluids for transport ECS
c      parameter (n0=-ncmax-nrefmx,nx=ncmax)
c      parameter (mxfeq=2)        !max number of FEQ EOS in block data
c      parameter (mxtrm=72)
c      parameter (nsave=(nx-n0+1)*(mxtrm+1))
c      parameter (nsave1=(nx-n0+1)*mxtrm)
c      parameter (nxsav=nx-n0+1)
c      character*3 hcpfeq
c      character*12 hcas
c      common /CASFEQ/ hcas(mxfeq)
c      common /CPMFEQ/ hcpfeq(mxfeq)
c      common /NTFEQ/ ntermf(mxfeq)
c      common /LFFEQ/ dli(mxfeq,mxtrm)
c      common /CFFEQ/ a(mxfeq,mxtrm),di(mxfeq,mxtrm),ti(mxfeq,mxtrm),
c     &               rho0(mxfeq),t0(mxfeq),
c     &               pc(mxfeq),rhoc(mxfeq),tc(mxfeq),
c     &               wmf(mxfeq),Rfeq(mxfeq),
c     &               pmin(mxfeq),rhotp(mxfeq),tmin(mxfeq),
c     &               tmax(mxfeq),pmax(mxfeq)
c      common /MSCFEQ/ ttpf(mxfeq),tnbpf(mxfeq),accenf(mxfeq),dipm(mxfeq)
c      common /FEQSAV/ phisav(n0:nx,mxtrm),delsav(n0:nx),tausav(n0:nx),
c     &                taup(n0:nx,mxtrm),delp(n0:nx,mxtrm),
c     &                delli(n0:nx,mxtrm)
cc
c      data phisav /nsave*0.0d0/
c      data delsav /nxsav*0.0d0/
c      data tausav /nxsav*0.0d0/
c      data taup /nsave1*0.0d0/
c      data delp /nsave1*0.0d0/
c      data delli /nsave1*0.0d0/
cc
cc
cc  R134a  1,1,1,2-tetrafluoroethane
c      data hcas(1) /'811-97-2'/
cc  use polynomial Cp0 model (at least for now)
c      data hcpfeq(1) /'CPP'/
c      data ntermf(1) /21/
c      data (a(1,j),j=1,21)/
c     &  0.5586817d-01,  0.4982230d+00,  0.2458698d-01,  0.8570145d-03,
c     &  0.4788584d-03, -0.1800808d+01,  0.2671641d+00, -0.4781652d-01,
c     &  0.1423987d-01,  0.3324062d+00, -0.7485907d-02,  0.1017263d-03,
c     & -0.5184567d+00, -0.8692288d-01,  0.2057144d+00, -0.5000457d-02,
c     &  0.4603262d-03, -0.3497836d-02,  0.6995038d-02, -0.1452184d-01,
c     & -0.1285458d-03/
cc  exponents for density terms
c      data (di(1,j),j=1,21)/
c     &       2.00d+00,       1.00d+00,       3.00d+00,       6.00d+00,
c     &       6.00d+00,       1.00d+00,       1.00d+00,       2.00d+00,
c     &       5.00d+00,       2.00d+00,       2.00d+00,       4.00d+00,
c     &       1.00d+00,       4.00d+00,       1.00d+00,       2.00d+00,
c     &       4.00d+00,       1.00d+00,       5.00d+00,       3.00d+00,
c     &      10.00d+00/
cc  exponents for temperature terms
c      data (ti(1,j),j=1,21)/
c     &      -0.50d+00,       0.00d+00,       0.00d+00,       0.00d+00,
c     &       1.50d+00,       1.50d+00,       2.00d+00,       2.00d+00,
c     &       1.00d+00,       3.00d+00,       5.00d+00,       1.00d+00,
c     &       5.00d+00,       5.00d+00,       6.00d+00,      10.00d+00,
c     &      10.00d+00,      10.00d+00,      18.00d+00,      22.00d+00,
c     &      50.00d+00/
cc  power of (-del) in the exponential multiplier
c      data (dli(1,j),j=1,21)/
c     &       0,              0,              0,              0,
c     &       0,              0,              0,              0,
c     &       1,              1,              1,              2,
c     &       2,              2,              2,              2,
c     &       2,              3,              3,              3,
c     &       4/
cc
c      data rho0(1),t0(1)
c     &  /4.978830171d0,374.18d0/          ! reducing parameters
c      data pc(1),rhoc(1),tc(1)
c     &  /4059.28d0,5.017053d0,374.21d0/   ! critical parameters
c      data wmf(1),Rfeq(1)
c     &  /102.032d0,8.314471d0/         ! mol weight and gas constant
c      data pmin(1),rhotp(1),tmin(1)
c     &  /0.391d0,15.594d0,169.85d0/    ! lower limits (= triple point)
c      data tmax(1),pmax(1)
c     &  /453.15d0,70000.0d0/           ! upper limits
cc
c      data ttpf(1) /169.85d0/
c      data tnbpf(1) /247.07d0/
c      data accenf(1) /0.32684d0/
c      data dipm(1) /2.058d0/  !dipole moment [Debye]; Meyer, (1991)
cc
c      end                                              !block data BDFEQ
c
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c ======================================================================
c                                                    end file core_FEQ.f
c ======================================================================
