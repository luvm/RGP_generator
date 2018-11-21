c  begin file core_PR.f
c
c  This file contains the routines implementing the Peng-Robinson and other
c  similar cubic equations of state.
c
c  contained here are:
c     function PHIPR (icomp,itau,idel,tau,del)
c     function PHIPRX (itau,idel,tau,del,x)
c     subroutine CRTPR (icomp,tc,pc,rhoc)
c     subroutine CRTPRX (x,tc,pc,rhoc)
c     subroutine REDPR (x,tc,rhoc)
c     subroutine SETPR (nread,icomp,hcasno,ierr,herr)
c     subroutine SETPRCO(i)
c     subroutine PREOS (i)
c     subroutine TPRHOPR (t,p,x,rho1,rho2)
c     subroutine CUBIC(a,z)
c     subroutine ESTPR(i,j,ibin,ierr,herr)
c
c ======================================================================
c ======================================================================
c
      function PHIPR (icomp,itau,idel,tau,del)
c
c  compute reduced Helmholtz energy or a derivative as functions
c  of dimensionless temperature and density for the Peng-Robinson or SRK
c  equations of state
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c     itau--flag specifying order of temperature derivative to calc
c     idel--flag specifying order of density derivative to calculate
c           when itau = 0 and idel = 0, compute A/RT
c           when itau = 0 and idel = 1, compute 1st density derivative
c           when itau = 1 and idel = 1, compute cross derivative
c           etc.
c      tau--dimensionless temperature (Tc/T)
c      del--dimensionless density (D/Dc)
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
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  07-21-03 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      parameter (ncmax=20)        !max number of components in mixture
      common /NCOMP/ nc,ic
      dimension x(ncmax)
c
      do i=1,nc
        x(i)=0.d0
      enddo
      x(icomp)=1.d0
      phipr=PHIPRX (itau,idel,tau,del,x)
      return
      end                                                !function PHIPR
c
c ======================================================================
c
      function PHIPRX (itau,idel,tau,del,x)
c
c  compute reduced Helmholtz energy or a derivative as functions
c  of dimensionless temperature and density for a volume-translated
c  modification of the Peng-Robinson equations of state
c  EOS: P=RT/(v+t+b)-a/((v+t)*(v+t+b)+b(v+t+b))
c  original translation concept given in Peneloux, A. and Rauzy, E.,
c  FPE 8 (1982) 7-23  (although it was applied to SRK)
c  Translation is considered as a constant (see Pfohl, FPE 163 (1999) 157-159)
c  to avoid crossing of isotherms at high density, near critical conditions
c  If translation is set to zero, original Peng Robinson is recovered.
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c     itau--flag specifying order of temperature derivative to calc
c     idel--flag specifying order of density derivative to calculate
c           when itau = 0 and idel = 0, compute A/RT
c           when itau = 0 and idel = 1, compute 1st density derivative
c           when itau = 1 and idel = 1, compute cross derivative
c           etc.
c      tau--dimensionless temperature (Tc/T)
c      del--dimensionless density (D/Dc)
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
c             0    3      del**3*[d**3(A/RT)/d(del)**3]
c             1    1      tau*del*[d**2(A/RT)/d(tau)d(del)]
c                         etc.
c
c  written by NIST Thermophysics Division, Boulder, Colorado
c  10-11-04 original version
c  09-11-09 EWL, remove x dependence for pure fluids
c  06-08-10 EWL, increase fkij array size
c
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (nmxprb=6)      !number of coeff for PR eos
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*16 drvflg(n0:nx)
      character*12 hcas
      COMMON /PRSV2/ drvflg
c     common /PRSAV/ delsav(n0:nx),tausav(n0:nx),drvsav(n0:nx,16)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /CCAS/ hcas(n0:nx)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /NCOMP/ nc,ic
      common /CFXPR/ fkij(nx,nx,nmxprb)
      common /CONPR/ prcoef(0:nx,20),iprflag(0:nx)
      double precision x(ncmax),ai(ncmax),bi(ncmax),ki(ncmax)
      dimension atci(ncmax)
      dimension daidt(ncmax),d2aidt2(ncmax)
      common /WCFPR/ pcpr(n0:nx),rhocpr(n0:nx),tcpr(n0:nx),
     &               wmfpr(n0:nx),Rpr(n0:nx),acnpr(n0:nx)
c
      call ISPURE (x,icomp)
c
      phiprx=0.d0
      if (del.le.1.d-10) RETURN  !trivial solution at zero density
      if (tau.le.0.d0) RETURN    !for any and all derivatives
c
c     ncode=idel*4+itau+1
c     if (abs(tau-tausav(icomp)).lt.1.0d-12 .and.
c    &    abs(del-delsav(icomp)).lt.1.0d-16) then
c  retrieve value from previous call
c       if (drvflg(icomp)(ncode:ncode).eq.'1') then
c         phiprx=drvsav(icomp,ncode)
c         RETURN
c       endif
c     else
c  otherwise, compute new values and save for possible future use
c  first compute needed powers of tau and del (and save for future use)
c       drvflg(icomp)='0000000000000000'
c
      call REDX (x,tred,Dred)
      if (icomp.ne.0 .and. hcas(icomp).eq.'811-97-2') then     !R-134a
        tred=374.18d0
        Dred=4.978830171d0
      endif
      t=tred/tau
      d=Dred*del
      v=1.d0/d
c
c  set up volume translation factor
      tshift=0.d0
      if (icomp.ne.0) then
        tshift=prcoef(icomp,1)
      else
        do i=1,nc
          tshift=tshift+x(i)*prcoef(i,1)
        enddo
      endif
c     tshift=0.d0      !For testing purposes only
c
      rt=r*t
c  parameters for the SRK equation
      u=1.d0
      w=0.d0
c  parameters for the Peng-Robinson equation
      u=2.d0
      w=-1.d0
c
      do i=1,nc
        bi(i)=0.0778d0*r*tcpr(i)/pcpr(i)
        ki(i)=0.37464d0+1.54226d0*acnpr(i)-0.26992d0*acnpr(i)**2
        atci(i)=0.45724d0*r**2*tcpr(i)**2/pcpr(i)
        ai(i)=atci(i)*(1.d0+ki(i)*(1.d0-SQRT(t/tcpr(i))))**2
      enddo
c
      if (icomp.ne.0) then
        a=ai(icomp)
        b=bi(icomp)
      else
        a=0.d0
        b=0.d0
        do i=1,nc
          do j=1,nc
            a=a+x(i)*x(j)*SQRT(ai(i)*ai(j))*(1.d0-fkij(i,j,1))
          enddo
          b=b+x(i)*bi(i)
        enddo
      endif
c
      vb=v-b
      q=b*SQRT(u**2-4.d0*w)
      v1=2.d0*v+u*b+q
      v2=2.d0*v+u*b-q
      v1n=v1+2.d0*tshift
      v2n=v2+2.d0*tshift
      dts=d*tshift
      b2=b**2
      bu=b*u
      b2dw=b2*d*w
      arg1=vb+tshift
      arg2=v2n/v1n
      if (arg1.lt.1.d-9) arg1=1.d-9 !prevent neg log arguments
      if (arg2.lt.1.d-9) arg2=1.d-9
      if (   v.lt.1.d-9)    v=1.d-9
      if (   q.lt.1.d-9)    q=1.d-9
      t1=d*(b-tshift)/(1.d0-b*d+dts)
      t2=((1.d0+dts)*(1.d0+d*(tshift+bu))+b2dw*d)/d
c
c  check if derivatives are requested
c
      if (itau.eq.0 .and. idel.eq.0) then
        phiprx=LOG(v)-LOG(arg1)+a/rt*LOG(arg2)/q
c
      elseif (itau.eq.1 .and. idel.eq.1) then
c  compute cross derivative using terms from 1st derivatives
        do i=1,nc
          daidt(i)=-atci(i)*ki(i)/SQRT(t*tcpr(i))
     &            *(1.d0+ki(i)*(1.d0-SQRT(t/tcpr(i))))
        enddo
        if (icomp.ne.0) then
          dadt=daidt(icomp)
        else
          dadt=0.d0
          do i=1,nc
          do j=1,nc
            da=SQRT(ai(i)/ai(j))*daidt(j)/2.d0+
     &         SQRT(ai(j)/ai(i))*daidt(i)/2.d0
            dadt=dadt+x(i)*x(j)*(1.d0-fkij(i,j,1))*da
          enddo
          enddo
        endif
        trm=2.d0*(1.d0+b*d+d*tshift)
        phiprx=4.d0*d/r*(dadt-a/t)/(trm-d*q)/(trm+d*q)
c
      elseif (itau.eq.2 .and. idel.eq.1) then
c  compute cross derivative using term from 1st derivative
        phiprx=0.d0      !not yet implemented
c
      elseif (itau.eq.1 .and. idel.eq.2) then
c  compute cross derivative using term from 2nd derivative
        phiprx=0.d0      !not yet implemented
c
      elseif (itau.eq.2 .and. idel.eq.2) then
        phiprx=0.d0      !not yet implemented
c
      elseif (idel.eq.1) then
c  compute derivative w.r.t. del (dimensionless density)
        phiprx=t1-a/rt/t2
c-------alternative method:
c       phiprx=-1.d0+1.d0/(1.d0-b*d+d*tshift)
c    &         +a/rt/q*(2.d0/v1n-2.d0/v2n)/d
c
      elseif (idel.eq.2) then
c  compute 2nd derivative w.r.t. del (dimensionless density)
        phiprx=t1**2+a/rt/t2**2
     &          *(bu+2.d0*(tshift+dts*tshift+b*dts*u+b2dw))
c-------alternative method:
c       bdt=1.d0-b*d+d*tshift
c       phiprx=1.d0-1.d0/bdt-(tshift*dred-b*dred)*del/bdt**2
c    &        +4.d0*a/rt/q/d*((1.d0/v2n-1.d0/v1n)
c    &        +(1.d0/v1n**2-1.d0/v2n**2)/d)
c
      elseif (idel.eq.3) then
c  compute 3rd derivative w.r.t. del (dimensionless density)
        phiprx=2.d0*t1**3+2.d0*a/rt*d/t2**2
     &        *(tshift**2+b*tshift*u+b2*w)-2.d0*a/rt/t2**3
     &        *(bu+2.d0*(tshift+dts*tshift+b*dts*u+b2dw))**2
c-------alternative method:
c       dts=d*tshift
c       db=d*b
c       bt2=2.d0*(b+tshift)
c       bqtp=bt2+q
c       bqtm=bt2-q
c       term2n=-2.d0+(2.d0+6.d0*db**2+6.d0*dts+6.d0*dts**2
c    &                    -6.d0*db-12.d0*db*dts)/(1.d0-db+dts)**3
c       terma=2.d0+d*bqtm
c       termb=2.d0+d*bqtp
c       termc=4.d0*(db*b+b+tshift*(1.d0+dts+2.d0*db))-d*q**2
c       term3n=-8.d0*a/rt*(q**2
c    &        -4.d0*b**2-8.d0*b*tshift-4.d0*tshift**2
c    &        +2.d0*bqtp*termc/termb
c    &        +2.d0*bqtm*termc/terma)/terma**2/termb**2
c       phiprx=term2n+term3n*d**3
c
      elseif (itau.eq.1) then
c  compute derivative w.r.t. tau (dimensionless temperature)
        do i=1,nc
          daidt(i)=-atci(i)*ki(i)/SQRT(t*tcpr(i))
     &            *(1.d0+ki(i)*(1.d0-SQRT(t/tcpr(i))))
        enddo
        if (icomp.ne.0) then
          dadt=daidt(icomp)
        else
          dadt=0.d0
          do i=1,nc
          do j=1,nc
            da=SQRT(ai(i)/ai(j))*daidt(j)/2.d0+
     &         SQRT(ai(j)/ai(i))*daidt(i)/2.d0
            dadt=dadt+x(i)*x(j)*(1.d0-fkij(i,j,1))*da
          enddo
          enddo
        endif
        phiprx=(a/t-dadt)*LOG(arg2)/r/q
c
      elseif (itau.eq.2) then
c  compute 2nd derivative w.r.t. tau (dimensionless temperature)
        do i=1,nc
          daidt(i)=-atci(i)*ki(i)/SQRT(t*tcpr(i))*
     &     (1.d0+ki(i)*(1.d0-SQRT(t/tcpr(i))))
          d2aidt2(i)=0.5d0*atci(i)*ki(i)/SQRT(tcpr(i))
     &              *(1.d0+ki(i))/t**1.5d0
        enddo
        if (icomp.ne.0) then
          dadt=daidt(icomp)
          d2adt2=d2aidt2(icomp)
        else
          dadt=0.d0
          d2adt2=0.d0
          do i=1,nc
          do j=1,nc
            da=SQRT(ai(i)/ai(j))*daidt(j)/2.d0+
     &         SQRT(ai(j)/ai(i))*daidt(i)/2.d0
            daij=0.5d0*(daidt(i)/SQRT(ai(i)*ai(j))-
     &                  daidt(j)*SQRT(ai(i))/(ai(j)**1.5d0))
            daji=0.5d0*(daidt(j)/SQRT(ai(j)*ai(i))-
     &                  daidt(i)*SQRT(ai(j))/(ai(i)**1.5d0))
            aij=SQRT(ai(i)/ai(j))
            aji=SQRT(ai(j)/ai(i))
            d2a=(aij*d2aidt2(j)+daidt(j)*daij+
     &          aji*d2aidt2(i)+daidt(i)*daji)*0.5d0
            dadt=dadt+x(i)*x(j)*(1.d0-fkij(i,j,1))*da
            d2adt2=d2adt2+x(i)*x(j)*(1.d0-fkij(i,j,1))*d2a
          enddo
          enddo
        endif
        phiprx=d2adt2*t*LOG(arg2)/r/q
c
      elseif (itau.eq.3) then
c  compute 3rd derivative w.r.t. tau (dimensionless temperature)
        phiprx=0.d0      !not yet implemented
      end if
c
c     drvsav(icomp,ncode)=phiprx
c     drvflg(icomp)(ncode:ncode)='1'
c
      RETURN
      end                                               !function PHIPRX
c ======================================================================
c
      subroutine CRTPR (icomp,tc,pc,rhoc)
c
c  returns critical parameters associated with PR model
c
c  input:
c    icomp--pointer specifying component (1..nc)
c  outputs:
c       tc--critical temperature [K]
c       pc--critical pressure [kPa]
c     rhoc--molar density [mol/L] at critical point
c
c  written by NIST Thermophysics Division, Boulder, Colorado
c  09-24-04 original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
      tc=tcrit(icomp)
      pc=pcrit(icomp)
      rhoc=Dcrit(icomp)
c
      RETURN
      end                                              !subroutine CRTPR
c
c ======================================================================
c
      subroutine CRTPRX (x,tc,pc,rhoc)
c
c  returns critical parameters associated with PR model
c
c  input:
c        x--composition array [mol frac]
c  outputs:
c       tc--critical temperature [K]
c       pc--critical pressure [kPa]
c     rhoc--molar density [mol/L] at critical point
c
c  written by NIST Thermophysics Division, Boulder, Colorado
c  10-21-04 original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      DIMENSION x(ncmax)
c
      call ISPURE (x,icomp)
      if (icomp.ne.0) then
        tc=tcrit(icomp)
        pc=pcrit(icomp)
        rhoc=Dcrit(icomp)
      else
        tc=0.d0
        pc=0.d0
        rhoc=0.d0
        do k=1,nc
          tc=tc+x(k)*tcrit(k)
          pc=pc+x(k)*pcrit(k)
          rhoc=rhoc+x(k)*Dcrit(k)
        enddo
      endif
c
      RETURN
      end                                             !subroutine CRTPRX
c
c ======================================================================
c
      subroutine REDPR (x,tred,Dred)
c
c  returns reducing parameters associated with PR model
c
c  input:
c        x--composition array [mol frac]
c  outputs:
c     tred--reducing temperature [K]
c     Dred--reducing molar density [mol/L]
c
c  written by NIST Thermophysics Division, Boulder, Colorado
c  09-24-04 original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
c     common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
c    &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
c    &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
c    &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /NCOMP/ nc,ic
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
c
      call ISPURE (x,icomp)
      if (icomp.ne.0) then
        tred=tz(icomp)
        Dred=rhoz(icomp)
      else
        tsum=0.d0
        Vsum=0.d0
        do i=1,nc
          tsum=tsum+x(i)*tz(i)
          Vsum=Vsum+x(i)/rhoz(i)
        enddo
        if (ABS(tsum).gt.0.5d-12 .and. ABS(Vsum).gt.0.5d-12) then
          tred=tsum
          Dred=1.d0/Vsum
        else
          tred=100
          Dred=10
        endif
c
        RETURN
      endif
      end                                              !subroutine REDPR
c
c ======================================================================
c
      subroutine SETPR (nread,icomp,hcasno,ierr,herr)
c
c  set up working arrays for use with cubic equations of state
c
c  inputs:
c    nread--file to read data from
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c   hcasno--CAS number of component icomp (not req'd if reading from file)
c
c  outputs:
c     ierr--error flag:  0 = successful
c                        1 = error (e.g. fluid not found)
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  11-10-04 EWL, original version
c  09-11-09 EWL, initialize prcoef
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*3 hpheq,heos,hmxeos,hmodcp,hmodcp2
      character*12 hcasno
      character*255 herr
c     character*1 dummy
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c     common /CPMOD/ hmodcp(n0:nx)
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
c  limits associated with the equation of state
      common /EOSLIM/ tmn(n0:nx),tmx(n0:nx),pmx(n0:nx),rhomx(n0:nx)
      common /WCFPR/ pcpr(n0:nx),rhocpr(n0:nx),tcpr(n0:nx),
     &               wmf(n0:nx),Rpr(n0:nx),acnpr(n0:nx)
      common /CONPR/ prcoef(0:nx,20),iprflag(0:nx)
c
      ierr=0
      herr=' '
c
      if (nread.eq.999) then
        acnpr(icomp)=accen(icomp)
        wmf(icomp)=wm(icomp)
        tcpr(icomp)=tcrit(icomp)
        pcpr(icomp)=pcrit(icomp)
        rhocpr(icomp)=Dcrit(icomp)
        Rpr(icomp)=Reos(icomp)
        RETURN
      endif

c
c  read data from file
      read (abs(nread),*) tmin                 !lower temperature limit
      read (abs(nread),*) tmax                 !upper temperature limit
      read (abs(nread),*) pmax                 !upper pressure limit
      read (abs(nread),*) rhomx2               !upper density limit
      read (abs(nread),2003) hmodcp2           !pointer to Cp0 model
      read (abs(nread),*) wmf(icomp)           !molecular weight
      read (abs(nread),*) acnpr(icomp)           !acentric factor
      read (abs(nread),*) tcpr(icomp)
      read (abs(nread),*) pcpr(icomp)
      read (abs(nread),*) rhocpr(icomp)
      read (abs(nread),*) RPr(icomp)           !gas constant
      if (nread.gt.0) then
        rhomx(icomp)=rhomx2
        hmodcp(icomp)=hmodcp2
        accen(icomp)=acnpr(icomp)
        wm(icomp)=wmf(icomp)
        tcrit(icomp)=tcpr(icomp)
        pcrit(icomp)=pcpr(icomp)
        Dcrit(icomp)=rhocpr(icomp)
        tz(icomp)=tcpr(icomp)
        rhoz(icomp)=rhocpr(icomp)
        ptp(icomp)=0.d0
        dtp(icomp)=0.d0
        dtpv(icomp)=0.d0
        dnbpl(icomp)=0.d0
        dnbpv(icomp)=0.d0
        if (nc.eq.1 .and. icomp.eq.1) R=RPr(icomp)
        Reos(icomp)=RPr(icomp)
        Zcrit(icomp)=pcpr(icomp)/(RPr(icomp)*tcpr(icomp)*rhocpr(icomp))
        tmn(icomp)=tmin
        tmx(icomp)=tmax
        pmx(icomp)=pmax
c       rhomx(icomp)=rhomax(icomp)
      endif
      read (abs(nread),*) nj
      if (icomp.ge.0) then
        prcoef(icomp,1)=0.d0
        if (nj.gt.0) then
          iprflag(icomp)=1  !values read in
          do j=1,nj
            read (abs(nread),*) prcoef(icomp,j)
          enddo
        endif
      endif
c
      RETURN
 2003 format (a3)
      end                                              !subroutine SETPR
c
c ======================================================================
c
      subroutine SETPRCO(i)
c
c     compute a value for the translation factor in translated Peng Robinson EOS
c     for component i
c
c  written by NIST Thermophysics Division, Boulder, Colorado
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /NCOMP/ nc,ic
      common /CONPR/ prcoef(0:nx,20),iprflag(0:nx)
      common /WCFPR/ pcpr(n0:nx),rhocpr(n0:nx),tcpr(n0:nx),
     &               wmfpr(n0:nx),Rpr(n0:nx),acnpr(n0:nx)
c
c     compute translation factor if not pre-loaded
c     default prescription is constant portion of translation in
c     Magoulas, K and Tassios, D. FPE 56 (1990) 119-140
c     do not use a t-dependent translation- see Pfohl, O. FPE 163 (1999) 157-159.
      if (iprflag(i).lt.1 .and. pcpr(i).gt.0.d0) then
         tr0=-0.014471d0+0.067498d0*acnpr(i)-0.084852d0*acnpr(i)**2
     &      +0.067298d0*acnpr(i)**3-0.017366d0*acnpr(i)**4
         prcoef(i,1)=(R*tcpr(i)/pcpr(i))*tr0
      endif
      RETURN
      end                                            !subroutine SETPRCO
c
c ======================================================================
c
      subroutine PREOS (i)
c
c  turn on or off the use of the PR cubic equation
c
c  inputs:
c        i--flag specifying use of PR:
c           0 - Use full equation of state (Peng-Robinson off)
c           1 - Use full equation of state with Peng-Robinson for sat. conditions
c               (not currently working)
c           2 - Use Peng-Robinson equation for all calculations
c           if i=-1, then i is returned with current usage of PR:  0, 1, or 2.
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  07-21-03 EWL, original version
c  08-11-06 MLH, use variable hrf instead of 'DEF'
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PREOS
c     dll_export PREOS
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*3 hpheq,heos,hmxeos,hmodcp
      character*3 hsav,hmxsav
      character*3 hrefst,hrefdf
      character*255 herr
      character hrf*3
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /PRHSAV/ hsav,hmxsav(n0:nx)
      common /FSHSAV/ lsatt,lsatp
      common /REFST/ hrefst,hrefdf(n0:nx)
      dimension x0(ncmax)
      lsatt=.false.
      lsatp=.false.
      k=1
      hrf='DEF'
      if (i.eq.-1) then
        i=0
        if (heos.eq.'PR') i=1
      elseif (i.eq.1) then
        if (heos.ne.'PR') hsav=heos
        do j=1,nx
          if (hmxeos(j).ne.'PR') hmxsav(j)=hmxeos(j)
          hmxeos(j)='PR'
        enddo
        heos='PR'
        hrefst=' '
        call SETREF (hrf,k,x0,h0,s0,t0,p0,ierr,herr)
      elseif (i.eq.2) then       ! this should be for full prs
        if (heos.ne.'PR') hsav=heos
        do j=1,nx
          if (hmxeos(j).ne.'PR') hmxsav(j)=hmxeos(j)
          hmxeos(j)='PR'
        enddo
        heos='PR'
        hrefst=' '
        call SETREF (hrf,k,x0,h0,s0,t0,p0,ierr,herr)
      elseif (hsav.ne.' ') then
        do j=1,nx
          if (hmxsav(j).ne.' ') hmxeos(j)=hmxsav(j)
        enddo
        if (hsav.ne.' ') heos=hsav
        hrefst=' '
        call SETREF (hrf,k,x0,h0,s0,t0,p0,ierr,herr)
      endif
c
      RETURN
      end                                              !subroutine PREOS
c
c ======================================================================
c
      subroutine TPRHOPR (t,p,x,rho1,rho2)
c  compute density using a volume-translated modification of the Peng-Robinson equation of state
c  EOS: P=RT/(v+t+b)-a/((v+t)*(v+t+b)+b(v+t+b)), t is a translation constant
c  Translation concept given in Peneloux, A. and Rauzy, E., FPE 8 (1982) 7-23
c
c  written by NIST Thermophysics Division, Boulder, Colorado
c  10-11-04 original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      character*12 hcas
c     character*255 herr
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (nmxprb=6)      !number of coeff for PR eos
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /CCAS/ hcas(n0:nx)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /NCOMP/ nc,ic
      common /CFXPR/ fkij(nx,nx,nmxprb)
      common /CONPR/ prcoef(0:nx,20),iprflag(0:nx)
      common /WCFPR/ pcpr(n0:nx),rhocpr(n0:nx),tcpr(n0:nx),
     &               wmfpr(n0:nx),Rpr(n0:nx),acnpr(n0:nx)
      double precision x(ncmax),ai(ncmax),bi(ncmax),ki(ncmax)
      dimension atci(ncmax),az(3),zz(3)
      icomp=0
      call ISPURE (x,icomp)
c     ierr=0
c     herr=' '
c
      call REDX (x,tred,Dred)
      if (icomp.ne.0 .and. hcas(icomp).eq.'811-97-2') then     !R-134a
        tred=374.18d0
        Dred=4.978830171d0
      endif
      rt=r*t
c     set volume translation factor
      tshift=0.d0
      if (icomp.ne.0) then
        tshift=prcoef(icomp,1)
      else
        do i=1,nc
          tshift=tshift+x(i)*prcoef(i,1)
        enddo
      endif
c     tshift=0.d0      !For testing purposes only

      do i=1,nc
        bi(i)=0.0778d0*r*tcpr(i)/pcpr(i)
        ki(i)=0.37464d0+1.54226d0*acnpr(i)-0.26992d0*acnpr(i)**2
        atci(i)=0.45724d0*r**2*tcpr(i)**2/pcpr(i)
        ai(i)=atci(i)*(1.d0+ki(i)*(1.d0-SQRT(t/tcpr(i))))**2
      enddo
      if (icomp.ne.0) then
        a=ai(icomp)
        b=bi(icomp)
      else
        a=0.d0
        b=0.d0
        do i=1,nc
          do j=1,nc
            a=a+x(i)*x(j)*SQRT(ai(i)*ai(j))*(1.d0-fkij(i,j,1))
          enddo
          b=b+x(i)*bi(i)
        enddo
      endif
c
c
c
      aa=a*p/r**2/T**2
      bb=b*p/r/T
      az(1)=(-aa+bb+bb**2)*bb
      az(2)=aa-bb*(3.d0*bb+2.d0)
      az(3)=bb-1.d0
c     revised translated eos
c     solve the cubic equation
c     only use largest and smallest roots; middle one does not have a physical meaning
      call CUBIC(az,zz)
c     compute the translation term
      rho1=p/(zz(1)*rt) !largest density root
      rho2=p/(zz(2)*rt) !smallest density root
      v1=1.d0/rho1
      v2=1.d0/rho2
      rho1=1.d0/(v1-tshift)
      rho2=1.d0/(v2-tshift)
c
c     the following lines can be uncommented to demonstrate that the shifted equation
c     is equivalent to the non-shifted version
c     if the shift is zero, it is identical to the PR EOS
c
c      CC=tshift*p/r/T
c      az(1)=(-AA*BB +BB*BB +BB*BB*BB) +
c     :      (-3.*BB*BB*CC + AA*CC + BB*CC*CC + CC*CC*CC -2.*BB*CC-CC*CC)
c      az(2) = (AA-3.*BB*BB-2.*BB) +
c     :      (2.*BB*CC + 3.*CC*CC - 2.*CC)
c      az(3) = (BB-1.d0) +3.*CC
c      call CUBIC(az,zz)
c       compute the translation term
c      rho1 = p / (zz(1) * rt) !largest density root
c      rho2 = p / (zz(2) * rt) !smallest density root
c
      RETURN
      end                                            !subroutine TPRHOPR
c
c ======================================================================
c
      subroutine CUBIC (a,z)
c
c     General purpose routine to solve a cubic equation of the form
c       z**3 + a(3)*z**2 +a(2)*z + a(1) = 0
c
c     The smallest root is returned in z(1) and the largest in z(2).
c     The middle root, if it exists, is discarded.
c
c  written by NIST Thermophysics Division, Boulder, Colorado
c  09-24-04 original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      dimension a(3),z(3)
c  All of the digits are necessary for liquid densities with pressures very close to zero
      pi=3.14159265358979323846d0
      add1=2.d0*pi/3.d0
      add2=2.d0*add1
      a0=a(1)
      a1=a(2)
      a2=a(3)/3.d0
      q=a1/3.d0-a2*a2
      r=(a1*a2-a0)/2.d0-a2*a2*a2
      test=q*q*q+r*r
      if (test.lt.0.d0) then       ! three distinct roots
        test=SQRT(-test)
        q=2.d0*CBRTX(SQRT(r*r+test*test))
        theta=pi
        if (abs(r).gt.1.d-20) theta=ATAN(test/r)
        if (theta.lt.0.d0) theta=theta+pi
        theta=theta/3.d0
        z(1)=q*COS(theta)-a2
        z(2)=q*COS(theta+add1)-a2
        z(3)=q*COS(theta+add2)-a2
      else                          ! either one real, or three real w/ 2 identical
        test=SQRT(test)
        s1=CBRTX(r+test)
        s2=CBRTX(r-test)
        z(1)=s1+s2-a2
        z(2)=z(1)
        z(3)=z(1)
        if (test.gt.0.d0) RETURN   ! only one real
        z(2)=-0.5d0*(s1+s2)-a2
        z(3)=z(2)                   ! three real, but 2 identical
      endif
      zmin=MIN(z(1),z(2),z(3))
      z(2)=MAX(z(1),z(2),z(3))
      z(1)=zmin
      if (zmin.lt.0.d0) z(1)=z(2)
      RETURN
      end                                              !subroutine CUBIC
c ======================================================================
c
      subroutine ESTPR(i,j,ibin,ierr,herr)
c
c  estimate binary parameters for the Peng Robinson mixture model
c  estimates based on the following references
c     Nishiumi, H., Arai, T. and Takeuchi, K., FPE 42 (1988) 43-62
c     Valderrama and Reyes, FPE 13 (1983) 195-202
c
c
c  inputs:
c        i--component i
c        j--component j
c     ibin--counter for binary pair (1,2 = 1; 1,3 = 2; 2,3 = 3; etc.)
c  outputs:
c     ierr--error flag:  0 = successful
c                     -117 = no binary parameters found (warning)
c                      117 = no binary parameters found (critical error)
c
c     Binary interaction parameters fkij are written to common block CFXPR
c
c     original version 10.14.04
c  10-04-07  HH, change MXINF2 to MXINF1 as is the case in other code
c  05-21-10 EWL, do not set hbin in this routine, overwrites HMX file
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (nmxprb=6)      !number of coeff for PR eos
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (nbin=ncmax*(ncmax-1)/2)   !# possible binary pairs
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
      character*1 htab,hnull
      character*255 family,UNNumb,herr,hbin
      character*12 hcas
      common /HCHAR/ htab,hnull
      common /NCOMP/ nc,ic
      common /FAML/ family(n0:nx),UNNumb(n0:nx)  !fluid family type
      common /CFXPR/ fkij(nx,nx,nmxprb)
      common /CCAS/ hcas(n0:nx)
c     common /MXINF1/ hbin(nbin)
      common /WCFPR/ pcpr(n0:nx),rhocpr(n0:nx),tcpr(n0:nx),
     &               wmfpr(n0:nx),Rpr(n0:nx),acnpr(n0:nx)
      dimension ifam(n0:nx),IJ(2)
c
c     legend
c     igroup    type
c
c     0         not covered in the correlation
c     1         c1-c16 alkane, branched alkanes
c     2         cycloalkane
c     3         alkene
c     4         aromatic
c     5         c18-c20 alkanes
c     6         CO2
c     7         N2
c     8         H2S
c     9         acetylene
c
c               mixtures with hydrogen are covered separately using Valderram correlation
c
      dimension igroup(0:11,0:11),binco(5,14)
c
      DATA (igroup(k, 0),k=0,11)/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA (igroup(k, 1),k=0,11)/0, 1, 1, 2, 3, 4, 5, 7, 8,10, 0, 0/
      DATA (igroup(k, 2),k=0,11)/0, 1, 1, 0, 3, 0, 5, 0, 0, 0, 0, 0/
      DATA (igroup(k, 3),k=0,11)/0, 2, 0, 2, 3, 4, 6, 7 ,8, 9, 0, 0/
      DATA (igroup(k, 4),k=0,11)/0, 3, 3, 3, 3, 0, 6, 7,12,11, 0, 0/
      DATA (igroup(k, 5),k=0,11)/0, 4, 0, 4, 0, 0, 4, 0, 0, 0, 0, 0/
      DATA (igroup(k, 6),k=0,11)/0, 5, 5, 6, 6, 4, 0,12,13, 0, 0, 0/
      DATA (igroup(k, 7),k=0,11)/0, 7, 0, 7, 7, 0,12, 0,14, 0, 0, 0/
      DATA (igroup(k, 8),k=0,11)/0, 8, 0, 8,12, 0,13,14, 0, 0, 0, 0/
      DATA (igroup(k, 9),k=0,11)/0, 10,0, 9,11, 0, 0, 0, 0, 0, 0, 0/
      DATA (igroup(k,10),k=0,11)/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA (igroup(k,11),k=0,11)/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA (binco(k,1), k=1,5)/1.041, 0.110,-0.0403, 0.0367, 0.0000/
      DATA (binco(k,2), k=1,5)/1.017,-0.417,-0.0124, 0.0852, 0.0000/
      DATA (binco(k,3), k=1,5)/1.025, 0.317,-0.0385,-0.0258, 0.0000/
      DATA (binco(k,4), k=1,5)/0.823, 0.000, 0.0673, 0.0000,-0.0051/
      DATA (binco(k,5), k=1,5)/0.883, 0.000, 0.0023, 0.0000, 0.0000/
      DATA (binco(k,6), k=1,5)/0.948, 0.000,-0.0084, 0.0000, 0.0000/
      DATA (binco(k,7), k=1,5)/0.982, 0.000,-0.0241, 0.0000, 0.0000/
      DATA (binco(k,8), k=1,5)/0.907, 0.000, 0.0109, 0.0000, 0.0000/
      DATA (binco(k,9), k=1,5)/1.090, 0.000,-0.1435, 0.0000, 0.0000/
      DATA (binco(k,10),k=1,5)/0.855, 0.000, 0.0000, 0.0000, 0.0000/
      DATA (binco(k,11),k=1,5)/0.965, 0.000, 0.0000, 0.0000, 0.0000/
      DATA (binco(k,12),k=1,5)/1.016, 0.000, 0.0000, 0.0000, 0.0000/
      DATA (binco(k,13),k=1,5)/0.894, 0.000, 0.0000, 0.0000, 0.0000/
      DATA (binco(k,14),k=1,5)/0.848, 0.000, 0.0000, 0.0000, 0.0000/
c
c     assign family code to components i and j
c     initialize to zero
      do k=1,nmxprb
        fkij(i,j,k)=0.0d0
        fkij(j,i,k)=0.0d0
      enddo
      ierr=0
      ij(1)=i
      ij(2)=j
      if (i.ne.j) then
        do k=1,2
          L=ij(k)
          if ((family(L).eq.'br-alkane').or.
     &        (family(L).eq.'n-alkane')) then
            if (wm(L).lt.227) then                          !c1 to c16
              ifam(L)=1
            elseif (wm(L).gt.227 .and. wm(L).lt.283) then !c17 to c20
              ifam(L)=5
            else
              ifam(L)=5                                   !treat like c17-c20 :untested
            endif
          elseif (family(L).eq.'napthene') then
            ifam(L)=2
          elseif (family(L).eq.'br-alkene') then
            ifam(L)=3
          elseif (family(L).eq.'n-alkene') then
            ifam(L)=3
          elseif (family(L).eq.'aromatic') then
            ifam(L)=4
          elseif (hcas(L).eq.'124-38-9')   then            !co2
            ifam(L)=6
          elseif (hcas(L).eq.'7727-37-9')  then            !n2
            ifam(L)=7
          elseif (hcas(L).eq.'7783-06-4')  then            !h2s
            ifam(L)=8
          elseif (hcas(L).eq.'74-86-2')    then            !acetylene
            ifam(L)=9
          elseif (family(L).eq.'other')    then            !other fluids
            ifam(L)=10
          else
            ifam(L)=0                                     !not covered
          endif
        enddo
c
c    The following is the correlation of Nishiumi and Arai
        if (igroup(ifam(i),ifam(j)).eq.0) then
          fkij(i,j,1)=0.d0
          ierr=-117
        else
          wdif=ABS(acnpr(i)-acnpr(j))
          kk=igroup(ifam(i),ifam(j))
          bigc=binco(1,kk)+binco(2,kk)*wdif
          bigd=binco(3,kk)+binco(4,kk)*wdif
          bige=binco(5,kk)
          vcrat=dcrit(i)/dcrit(j)
          if (vcrat.lt.1.d0) vcrat=1.d0/vcrat
          fkij(i,j,1)=1.d0-(bigc+bigd*vcrat+bige*vcrat**2)
          fkij(j,i,1)=fkij(i,j,1)
          ierr=0
        endif
c
c     The following is the correlation of Valderrama, J.O. and Reyes, L.R.
c     FPE 13, 195-202.  for mixtures with hydrogen,
c     but evaluated at the normal boiling point of the
c     non-hydrogen component instead of using full temperature dependence.
        if (hcas(i).eq.'1333-74-0' .or. hcas(j).eq.'1333-74-0') then
          if (hcas(i).eq.'1333-74-0') then
            wi2=acnpr(j)
            tr2=tnbp(j)/tcrit(j)
          else
            wi2=acnpr(i)
            tr2=tnbp(i)/tcrit(i)
          endif
          alpha2=0.1805d0+3.21d0*wi2+2.437d0*wi2**2
          beta2=0.1323d0+0.5507d0*wi2+3.5994d0*wi2**2
          fkij(i,j,1)=alpha2-beta2/tr2
          fkij(j,i,1)=fkij(i,j,1)
          ierr=0
        endif
      else
        fkij(i,j,1)=0.d0
        ierr=0
      endif
c
c     The following is based loosely upon correlation of water+alkane given in
c     Daridon, J.L., Lagourette, B., Saint-Guirons, H. and Xans, P., FPE 91(1993),31-54
c     does not implement temp or composition dependent contribution; thus just a rough estimate
      if (hcas(i).eq.'7732-18-5') then
         if (family(j).eq.'n-alkane') then
           nmet=(INT(wm(j))-2)/14 -2  !number of -ch2- groups
           fkij(i,j,1)=0.5d0 -0.00834d0*float(nmet)
           fkij(j,i,1)=fkij(i,j,1)
           ierr=0
         endif
      elseif (hcas(j).eq.'7732-18-5') then
         if (family(i).eq.'n-alkane') then
           nmet=(INT(wm(i))-2)/14 -2  !number of -ch2- groups
           fkij(j,i,1)=0.5d0 -0.00834d0*float(nmet)
           fkij(i,j,1)=fkij(j,i,1)
           ierr=0
         endif
      endif
c
      fkij(j,i,1)=fkij(i,j,1)
c
c     set warning messages for GUI
      if (ierr.eq.-117) then
        herr='[SETUP error -117] Binary interaction parameters are not '
     &     //'presently available for this mixture; calculations '
     &     //'will be made assuming ideal solution behavior. '//hnull
c       hbin(ibin)='Estimation of mixing parameters is not available '
c    &     //'for this mixture.'//hnull
c  do not call this yet until this routine is fully functional
c       call ERRMSG (ierr,herr)
      endif
c
      RETURN
      END                                             ! subroutine ESTPR
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c ======================================================================
c                                                     end file core_PR.f
c ======================================================================
