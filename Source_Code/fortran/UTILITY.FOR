c  begin file utility.f
c
c  This file contains various utility subroutines to retrieve information
c  about the components.
c
c  contained here are:
c     subroutine INFO (icomp,wmm,ttrp,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas)
c     subroutine NAME (icomp,hname,hn80,hcas)
c     function WMOL (x)
c     subroutine XMASS (xmol,xkg,wmix)
c     subroutine XMOLE (xkg,xmol,wmix)
c     subroutine LIMITX (htyp,t,D,p,x,tmin,tmax,Dmax,pmax,ierr,herr)
c     subroutine LIMITK (htyp,icomp,t,D,p,tmin,tmax,Dmax,pmax,ierr,herr)
c     subroutine LIMITS (htyp,x,tmin,tmax,Dmax,pmax)
c     subroutine ERRMSG (ierr,herr)
c     subroutine QMASS (qmol,xl,xv,qkg,xlkg,xvkg,wliq,wvap,ierr,herr)
c     subroutine QMOLE (qkg,xlkg,xvkg,qmol,xl,xv,wliq,wvap,ierr,herr)
c     subroutine GOLD (x0i,x1i,nc,lmax,z,z2,z3,z4,GEVAL,xopt,yopt,ierr)
c     subroutine DOTFILL (nc,x,ptest,filrat,ierr,herr)
c     function CBRTX (x)
c
c  various arrays are dimensioned with parameter statements
c     parameter (ncmax=20)        !max number of components in mixture
c     parameter (nrefmx=10)       !max number of fluids for transport ECS
c     parameter (n0=-ncmax-nrefmx,nx=ncmax)
c
c ======================================================================
c ======================================================================
c
      subroutine INFO (icomp,wmm,ttrp,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas)
c
c  provides fluid constants for specified component
c
c  input:
c    icomp--component number in mixture; 1 for pure fluid
c  outputs:
c      wmm--molecular weight [g/mol]
c     ttrp--triple point temperature [K]
c    tnbpt--normal boiling point temperature [K]
c       tc--critical temperature [K]
c       pc--critical pressure [kPa]
c       Dc--critical density [mol/L]
c       Zc--compressibility at critical point [pc/(Rgas*Tc*Dc)]
c      acf--acentric factor [-]
c      dip--dipole moment [debye]
c     Rgas--gas constant [J/mol-K]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  01-31-97  MM, original version
c  02-19-97  MM, add check that input icomp is within bounds
c  10-02-97  MM, add compiler switch to allow access by DLL
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  03-08-00 EWL, change names of inputs to avoid conflicts with CCON
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: INFO
c     dll_export INFO
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
      if (abs(icomp).le.nc) then
        wmm=wm(icomp)
        ttrp=ttp(icomp)
        tnbpt=tnbp(icomp)
        tc=tcrit(icomp)
        pc=pcrit(icomp)
        Dc=Dcrit(icomp)
        Zc=Zcrit(icomp)
        acf=accen(icomp)
        dip=dipole(icomp)
        Rgas=R
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
      RETURN
      end                                               !subroutine INFO
c
c ======================================================================
c
      subroutine NAME (icomp,hname,hn80,hcas)
c
c  provides name information for specified component
c
c  input:
c    icomp--component number in mixture; 1 for pure fluid
c  outputs:
c    hname--component name [character*12]
c     hn80--component name--long form [character*80]
c     hcas--CAS (Chemical Abstracts Service) number [character*12]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  02-06-97  MM, original version
c  02-19-97  MM, add check that input icomp is within bounds
c  10-02-97  MM, add compiler switch to allow access by DLL
c  12-01-97  MM, add synonyms to /CNAM80/
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: NAME
c     dll_export NAME
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*12 hn,hcasn,hcas,hname
      character*80 hn80
      character*255 hnam80,hsyn1,hsyn2
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCAS/ hcasn(n0:nx)
      common /CNAM/ hn(n0:nx)
      common /CNAM80/ hnam80(n0:nx),hsyn1(n0:nx),hsyn2(n0:nx)
c
      if (abs(icomp).le.nc) then
        hname=hn(icomp)
        hn80=hnam80(icomp)(1:80)
        hcas=hcasn(icomp)
      else
        hname='not defined'
        hn80='not defined                             '//
     &       '                                        '
        hcas='not defined'
      end if
c
      RETURN
      end                                               !subroutine NAME
c
c ======================================================================
c
      function WMOL (x)
c
c  molecular weight for a mixture of specified composition
c
c  input:
c        x--composition array [array of mol frac]
c
c  output (as function value):
c     WMOL--molar mass [g/mol], a.k.a. "molecular weight"
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  01-10-96  MM, original version
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                add Zcrit to common /CCON/
c  03-19-96  MM, add dipole moment to /CCON/
c  10-02-97  MM, add compiler switch to allow access by DLL
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-02-98 EWL, remove compositional dependence for pure fluids
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: WMOL
c     dll_export WMOL
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      dimension x(ncmax)
      common /NCOMP/ nc,ic
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
      call ISPURE (x,icomp)
      if (icomp.ne.0) then
        WMOL=wm(icomp)
      else
        wsum=0.0d0
        do i=1,nc
          wsum=wsum+x(i)*wm(i)
        enddo
        WMOL=wsum
      endif
c
      RETURN
      end                                                 !function WMOL
c
c ======================================================================
c
      subroutine XMASS (xmol,xkg,wmix)
c
c  converts composition on a mole fraction basis to mass fraction
c
c  input:
c     xmol--composition array [array of mol frac]
c  outputs:
c      xkg--composition array [array of mass frac]
c     wmix--molar mass of the mixture [g/mol], a.k.a. "molecular weight"
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  04-08-96  MM, original version
c  10-02-97  MM, add compiler switch to allow access by DLL
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: XMASS
c     dll_export XMASS
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      dimension xmol(ncmax),xkg(ncmax),xsumi(ncmax)
      common /NCOMP/ nc,ic
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
      call ISPURE (xmol,icomp)
      if (icomp.ne.0) then
        wmix=wm(icomp)
        do i=1,nc
          xkg(i)=0.d0
        enddo
        xkg(icomp)=1.d0
        RETURN
      endif
      xsum=0.0d0
      do i=1,nc
        xsumi(i)=xmol(i)*wm(i)
        xsum=xsum+xsumi(i)
      enddo
      wmix=xsum
      xsinv=1.d0
      if (ABS(xsum).gt.1.0d-10) xsinv=1.0d0/xsum
      do i=1,nc
        xkg(i)=xsumi(i)*xsinv
      enddo
c
      RETURN
      end                                              !subroutine XMASS
c
c ======================================================================
c
      subroutine XMOLE (xkg,xmol,wmix)
c
c  converts composition on a mass fraction basis to mole fraction
c
c  input:
c      xkg--composition array [array of mass frac]
c  outputs:
c     xmol--composition array [array of mol frac]
c     wmix--molar mass of the mixture [g/mol], a.k.a. "molecular weight"
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  04-08-96  MM, original version
c  10-02-97  MM, add compiler switch to allow access by DLL
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  10-24-02 EWL, set xmol=1 for a pure fluid
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: XMOLE
c     dll_export XMOLE
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      dimension xmol(ncmax),xkg(ncmax),xsumi(ncmax)
      common /NCOMP/ nc,ic
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
      call ISPURE (xkg,icomp)
      if (icomp.ne.0) then
        wmix=wm(icomp)
        do i=1,nc
          xmol(i)=0.d0
        enddo
        xmol(icomp)=1.d0
        RETURN
      endif
      xsum=0.0d0
      wsum=0.0d0
      do i=1,nc
        xsumi(i)=xkg(i)/wm(i)
        xsum=xsum+xsumi(i)
      enddo
      wmix=xsum
      xsinv=1.d0
      if (ABS(xsum).gt.1.0d-10) xsinv=1.0d0/xsum
      do i=1,nc
        xmol(i)=xsumi(i)*xsinv
        wsum=wsum+xmol(i)*wm(i)
      enddo
      wmix=wsum
c
      RETURN
      end                                              !subroutine XMOLE
c
c ======================================================================
c
      subroutine LIMITX (htyp,t,D,p,x,tmin,tmax,Dmax,pmax,ierr,herr)
c
c  returns limits of a property model as a function of composition
c  and/or checks input t, D, p against those limits
c
c  Pure fluid limits are read in from the .fld files; for mixtures, a
c  simple mole fraction weighting in reduced variables is used.
c
c  Attempting calculations below the minimum temperature and/or above
c  the maximum density will result in an error.  These will often
c  correspond to a physically unreasonable state; also many equations of
c  state do not extrapolate reliably to lower T's and higher D's.
c
c  A warning is issued if the temperature is above the maximum but below
c  1.5 times the maximum; similarly pressures up to twice the maximum
c  result in only a warning. Most equations of state may be
c  extrapolated to higher T's and P's.  Temperatures and/or pressures
c  outside these extended limits will result in an error.
c
c  When calling with an unknown temperature, set t to -1 to avoid performing
c  the melting line check
c
c  inputs:
c     htyp--flag indicating which models are to be checked [character*3]
c           'EOS':  equation of state for thermodynamic properties
c           'ETA':  viscosity
c           'TCX':  thermal conductivity
c           'STN':  surface tension
c        t--temperature [K]
c        D--molar density [mol/L]
c        p--pressure [kPa]
c        x--composition array [mol frac]
c     N.B.--all inputs must be specified, if one or more are not
c           available, (or not applicable as in case of surface tension)
c           use reasonable values, such as:
c           t = tnbp
c           D = 0
c           p = 0
c  outputs:
c     tmin--minimum temperature for model specified by htyp [K]
c     tmax--maximum temperature [K]
c     Dmax--maximum density [mol/L]
c     pmax--maximum pressure [kPa]
c     ierr--error flag:  0 = all inputs within limits
c                      <>0 = one or more inputs outside limits:
c                       -1 = 1.5*tmax > t > tmax
c                        1 = t < tmin or t > 1.5*tmax
c                        2 = D > Dmax or D < 0
c                       -4 = 2*pmax > p > pmax
c                        4 = p < 0 or p > 2*pmax
c                        8 = component composition < 0 or > 1
c                            and/or composition sum < 0 or > 1
c                       16 = p>pmelt
c                      -16 = t<ttrp (important for water)
c           if multiple inputs are outside limits, ierr = abs[sum(ierr)]
c           with the sign determined by the most severe excursion
c           (ierr > 0 indicate an error--calculations not possible,
c            ierr < 0 indicate a warning--results may be questionable)
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  05-30-96  MM, original version
c  06-03-96  MM, p > 2*pmax and t > 1.5*tmax result in error
c                add htyp to argument list
c  02-21-97  MM, add checks for viscosity and thermal conductivity
c  06-03-97  MM, initialize ierr = 0 and herr = hnull
c  06-04-97 EWL, zero delsum and xsum before do loop
c  06-17-97  MM, change format on xsum error
c  10-02-97  MM, add compiler switch to allow access by DLL
c  12-05-97  MM, x(i) missing in summation for max density
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-01-98 EWL, reduce Tmin for mixtures by 5 degrees
c  12-01-98 EWL, reduce Tmin by 1.0d-10 to avoid machine precision problems
c  12-22-98 EWL, add ammonia-water triple-point line
c  03-05-99 EWL, add checks for melting line pressure and sublimation pressure
c  09-02-99 EWL, increase check for p>pmelt by pmelt+.0001
c  10-12-99 EWL, add transport limits for the ECS model
c  10-20-99 EWL, change herrx to 140 to accommodate full string
c  10-20-99 EWL, increased the tolerance in sum(x)<>1 for setup with single prec.
c  10-22-99 EWL, call new subroutine LIMITS to get tmin,tmax,dmax, and pmax
c  01-25-00 EWL, do not allow p>pmax for parahydrogen (bad Younglove EOS)
c  01-25-00 EWL, do not allow t>tmax for krypton (bad Juza EOS)
c  02-25-00 EWL, do not check p>pmax when t=-1
c  05-25-00 EWL, reorganize ltemp logic
c  05-25-00 EWL, remove old code for calculating pmlt
c  07-11-00 EWL, remove check on t>tmax for krypton, equation was replaced
c  11-20-01 EWL, check for t<ttrp
c  11-20-01 EWL, allow mixtures to go below ttrp by 25 degrees
c  06-30-04 EWL, allow mixtures to go below ttrp by 25 degrees only if Tc's
c                differ by more than 50 (thus, air is not included)
c  03-28-06 EWL, allow p>pmax for parahydrogen with new Leachman equation
c  10-05-07  HH, remove common block WLMETA
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: LIMITX
c     dll_export LIMITX
c
      logical lerr,lwarn,ltemp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (nrf0=n0)    !lower limit for transport ref fluid arrays
      dimension x(ncmax)
      character*1 htab,hnull
      character*3 htyp
      character*75 herrt,herrd,herrp
      character*140 herrx
      character*255 herr
      character*3 heta,hetak,htcx,htcxk
      character*12 hcasn
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /EOSLIM/ tmeos(n0:nx),txeos(n0:nx),peos(n0:nx),Deos(n0:nx)
      common /WLMSTN/ tmstn(n0:nx),txstn(n0:nx)
      common /WLMTCX/ tmtcx(nrf0:nx),txtcx(nrf0:nx),ptcx(nrf0:nx),
     &                Dtcx(nrf0:nx)
      common /WLMTRN/ tmecs(nrf0:nx),txecs(nrf0:nx),pecs(nrf0:nx),
     &                Decs(nrf0:nx)
      common /TRNMOD/ heta,hetak(nrf0:nx),htcx,htcxk(nrf0:nx)
      common /CCAS/ hcasn(n0:nx)
c
c  initialize flags and strings
c
c     write (*,*) ' LIMITX--entering with htyp,t = ',htyp,t
      ierr=0
      ierrt=0
      ierrd=0
      ierrp=0
      ierrx=0
      xsum=0.0d0
      herr=' '
      herrt=' '
      herrd=' '
      herrp=' '
      herrx=' '
      nchart=1
      nchard=1
      ncharp=1
      ncharx=1
      lerr=.false.
      lwarn=.false.
      pmlt=1.0d15
      call ISPURE (x,icomp)
c
      call LIMITS(htyp,x,tmin,tmax,Dmax,pmax)
      if (ABS(p).gt.1.0d-10 .and. t.gt.0) call MELTT(t,x,pmlt,ierr,herr)
c
c  Set EOS and transport routines so that they cannot be extrapolated
c  to lower temps
      ltemp=.true.
      if (htyp.eq.'STN' .or. htyp.eq.'stn') then
        ltemp=.false. !surface tension may be extrapolated to lower T's
      end if
c
      if (icomp.eq.0) then
c  general mixture case
        xsum=0.0d0
        do i=1,nc
          xsum=xsum+x(i)
          if (x(i).lt.-1.0d-10 .or. x(i).gt.1.0000000001d0) then
            lerr=.true.
            ierrx=8
          end if
        enddo
        if (xsum.lt.0.999999d0 .or. xsum.gt.1.000001d0) then
          lerr=.true.
          ierrx=8
        end if
      end if
c
c  check inputs against limits
c
      if (t.lt.tmin-1.0d-10 .and. abs(t+1.0d0).gt.1.d-15) then
        pm=0
        ierr=0
        tsub=0.d0
        if (icomp.ne.0) then
          if (t.lt.ttp(icomp)) call SUBLT (t,x,pm,ierr,herr)
        else
c  set error number if p>pmelt.  If p=0, also set error number
          do i=1,nc-1
            do j=i+1,nc
              if (abs(tcrit(i)-tcrit(j)).gt.50) tsub=25.d0
            enddo
          enddo
        endif
        if (tsub.gt.0.d0) then
          call REDX (x,tred,Dred)
          if (tred.lt.200) tsub=5.d0
        endif
        if (ierr.ne.0 .or. p.gt.pm+1.0d-10 .or. p.lt.1.d-20) then
          if (icomp.ne.0 .or. t.lt.tmin-tsub) then
            if (ltemp) then
              lerr=.true.
              ierrt=1
            else
              lwarn=.true.
              ierrt=-1
            end if
            write (herrt,1010) t,tmin
            nchart=70
 1010       format (' temperature below lower limit, T =',g11.5,
     &              ' K, Tmin =',g11.5,' K;')
          endif
        endif
      else if (t.gt.1.5d0*tmax) then
        lerr=.true.
        ierrt=1
        write (herrt,1011) t,tmax
        nchart=72
 1011   format (' temperature > 1.5 x upper limit, T =',g11.5,
     &          ' K, Tmax =',g11.5,' K;')
      else if (t.gt.tmax) then
c       if (hcasn(icomp).eq.'7439-90-9' .and. icomp.ne.0) then
c         lerr=.true.
c         ierrt=1
c       else
          lwarn=.true.
          ierrt=-1
c       endif
        write (herrt,1012) t,tmax
        nchart=70
 1012   format (' temperature above upper limit, T =',g11.5,
     &          ' K, Tmax =',g11.5,' K;')
      end if
      if (D.lt.0.0d0) then
        lerr=.true.
        ierrd=2
        write (herrd,1020) D
        nchard=36
 1020   format (' density < 0, D =',g11.5,' mol/L;')
      else if (D.gt.Dmax) then
        lerr=.true.
        ierrd=2
        write (herrd,1021) D,Dmax
        nchard=74
 1021   format (' density above upper limit, D =',g11.5,
     &          ' mol/L, Dmax =',g11.5,' mol/L;')
      end if
      if (p.lt.0.0d0) then
        lerr=.true.
        ierrp=4
        write (herrp,1040) p/1000.0d0
        ncharp=34
 1040   format (' pressure < 0, P =',g11.5,' MPa;')
      else if (p.gt.2.0d0*pmax) then
        lerr=.true.
        ierrp=4
        write (herrp,1042) p/1000.0d0,pmax/1000.0d0
        ncharp=71
 1042   format (' pressure > 2 x upper limit, P =',g11.5,
     &          ' MPa, Pmax =',g11.5,' MPa;')
      else if (p.gt.pmax+1.0d-8) then
        lwarn=.true.
        ierrp=-4
        write (herrp,1044) p/1000.0d0,pmax/1000.0d0
        ncharp=71
 1044   format (' pressure above upper limit, P =',g11.5,
     &          ' MPa, Pmax =',g11.5,' MPa;')
      end if
      if (p*0.9995d0.gt.pmlt .and. abs(pmlt).gt.1.d-20) then
        lerr=.true.
        ierrp=16
        write (herrp,1050) p/1000.0d0,pmlt/1000.0d0
        ncharp=72
 1050   format (' pressure > melting pressure, P =',g11.5,
     &          ' MPa, Pmelt =',g11.5,' MPa;')
      end if
      if (ierrx.ne.0) then
c       write (*,1080) xsum,(x(i),i=1,nc)
        write (herrx,1080) xsum,(x(i),i=1,MIN(nc,5))
        ncharx=131
 1080   format (' composition(s) out of range, Xsum =',f13.10,
     &          ' mol frac, X(i) =',5f13.10)
      end if
      if (icomp.ne.0) then
        if (t.lt.ttp(icomp)-1.d-10 .and. ierrt.eq.0.and.t.gt.0.d0) then
c  check for cases where the temperature is less than the triple point
c  temperature, but still in a valid liquid region (like water between
c  251.165 and 273.16 K.)
          lwarn=.true.
          ierrt=-16
          write (herrt,1090) t,ttp(icomp)
          nchart=75
 1090     format (' temperature less than the triple point,',
     &          ' T=',g10.5,' K, Ttrp=',g10.5,' K;')
          if (hcasn(icomp).eq.'7732-18-5' .and. abs(p).gt.1.d-20) then
            call MLTH2O (t,p1,p2)
            if (p.lt.p2*0.9999d0 .or. p.gt.p1/0.9999d0) then
              lerr=.true.
              ierrt=16
              write (herrt,1095)
              nchart=34
 1095   format (' inputs are within the solid phase;')
            endif
          endif
        endif
      endif
c
c  compose error string and compute overall value of ierr
c
      if (lerr .or. lwarn) then
c       write (*,*) ' LIMITX--nchart,d,p,x:',nchart,nchard,nchard,ncharx
c       write (*,1999) ' LIMITX--herrt: ',herrt(1:nchart)
c       write (*,1999) ' LIMITX--herrd: ',herrd(1:nchard)
c       write (*,1999) ' LIMITX--herrp: ',herrp(1:ncharp)
c       write (*,1999) ' LIMITX--herrx: ',herrx(1:ncharx)
c1999   format (1x,a16,a80)
c       write (*,*) ' LIMITX--#char: t,d,p,x,sum:  ',nchart,nchard,
c    &              ncharp,ncharx,nchart+nchard+ncharp+ncharx
        herr='one or more inputs are out of range: '//herrt(1:nchart)
     &      //herrd(1:nchard)//herrp(1:ncharp)//herrx(1:ncharx)
     &      //hnull
        if (lerr) then
          ierr=abs(ierrt)+abs(ierrd)+abs(ierrp)+abs(ierrx)
        else if (lwarn) then
          ierr=-abs(ierrt)-abs(ierrd)-abs(ierrp)-abs(ierrx)
        end if
        call ERRMSG (ierr,herr)
      end if
c
      RETURN
      end                                             !subroutine LIMITX
c
c ======================================================================
c
      subroutine LIMITK (htyp,icomp,t,D,p,tmin,tmax,Dmax,pmax,ierr,herr)
c
c  returns limits of a property model (read in from the .fld files) for
c  a mixture component and/or checks input t, D, p against those limits
c
c  This routine functions in the same manner as LIMITX except that the
c  composition x is replaced by the component number icomp.
c
c  Attempting calculations below the minimum temperature and/or above
c  the maximum density will result in an error.  These will often
c  correspond to a physically unreasonable state; also many equations of
c  state do not extrapolate reliably to lower T's and higher D's.
c
c  A warning is issued if the temperature is above the maximum but below
c  1.5 times the maximum; similarly pressures up to twice the maximum
c  result in only a warning. Most equations of state may be
c  extrapolated to higher T's and P's.  Temperatures and/or pressures
c  outside these extended limits will result in an error.
c
c  inputs:
c     htyp--flag indicating which models are to be checked [character*3]
c           'EOS':  equation of state for thermodynamic properties
c           'ETA':  viscosity
c           'TCX':  thermal conductivity
c           'STN':  surface tension
c    icomp--component number in mixture; 1 for pure fluid
c        t--temperature [K]
c        D--molar density [mol/L]
c        p--pressure [kPa]
c     N.B.--all inputs must be specified, if one or more are not
c           available, (or not applicable as in case of surface tension)
c           use reasonable values, such as:
c           t = tnbp
c           D = 0
c           p = 0
c  outputs:
c     tmin--minimum temperature for model specified by htyp [K]
c     tmax--maximum temperature [K]
c     Dmax--maximum density [mol/L]
c     pmax--maximum pressure [kPa]
c     ierr--error flag:  0 = all inputs within limits
c                      <>0 = one or more inputs outside limits:
c                       -1 = 1.5*tmax > t > tmax
c                        1 = t < tmin or t > 1.5*tmax
c                        2 = D > Dmax or D < 0
c                       -4 = 2*pmax > p > pmax
c                        4 = p < 0 or p > 2*pmax
c                       16 = p>pmelt
c                      -16 = t<ttrp (important for water)
c           if multiple inputs are outside limits, ierr = abs[sum(ierr)]
c           with the sign determined by the most severe excursion
c           (ierr > 0 indicate an error--calculations not possible,
c            ierr < 0 indicate a warning--results may be questionable)
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by M. McLinden, NIST Phys & Chem Properties Div, Boulder, CO
c  03-18-97  MM, original version; based on LIMITX
c  06-03-97  MM, initialize ierr = 0 and herr = hnull
c  10-02-97  MM, add compiler switch to allow access by DLL
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-01-98 EWL, reduce Tmin by 1.0d-10 to avoid machine precision problems
c  10-12-99 EWL, add transport limits for the ECS model
c  10-22-99 EWL, add checks for melting line pressure and sublimation pressure
c  08-16-00 EWL, add dimension x(ncmax)
c  11-20-01 EWL, check for t<ttrp
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: LIMITK
c     dll_export LIMITK
c
      logical lerr,lwarn,ltemp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (nrf0=n0)    !lower limit for transport ref fluid arrays
      character*1 htab,hnull
      character*3 htyp
      character*75 herrt,herrd,herrp
      character*120 herrx
      character*255 herr
      character*3 heta,hetak,htcx,htcxk
      character*12 hcasn
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /EOSLIM/ tmeos(n0:nx),txeos(n0:nx),peos(n0:nx),Deos(n0:nx)
      common /WLMSTN/ tmstn(n0:nx),txstn(n0:nx)
      common /WLMTCX/ tmtcx(nrf0:nx),txtcx(nrf0:nx),ptcx(nrf0:nx),
     &                Dtcx(nrf0:nx)
      common /WLMETA/ tmeta(nrf0:nx),txeta(nrf0:nx),peta(nrf0:nx),
     &                Deta(nrf0:nx)
      common /WLMTRN/ tmecs(nrf0:nx),txecs(nrf0:nx),pecs(nrf0:nx),
     &                Decs(nrf0:nx)
      common /TRNMOD/ heta,hetak(nrf0:nx),htcx,htcxk(nrf0:nx)
      common /CCAS/ hcasn(n0:nx)
c
c  initialize flags and strings
c
c     write (*,*) ' LIMITK--entering with htyp,t = ',htyp,t
      ierr=0
      ierrt=0
      ierrd=0
      ierrp=0
      herr=' '
      herrt=' '
      herrd=' '
      herrp=' '
      herrx=' '
      nchart=1
      nchard=1
      ncharp=1
      ncharx=1
      lerr=.false.
      lwarn=.false.
      i=icomp
c
      pmlt=0.0d0
      if (htyp.eq.'EOS' .or. htyp.eq.'eos') then
c  equation of state
        ltemp=.true.       !EOS cannot be extrapolated to lower temps
        tmin=tmeos(i)
        tmax=txeos(i)
        Dmax=Deos(i)
        pmax=peos(i)
        if (d.gt.dtp(icomp) .or. ABS(d).lt.1.0d-10) then
          if (t.gt.tmeos(icomp).and.ABS(p).gt.1.0d-10) then
            call MELTK (icomp,t,pmlt,ierr,herr)
          end if
        endif
      else if (htyp.eq.'STN' .or. htyp.eq.'stn') then
c  surface tension model
        ltemp=.false.     !STN may be extrapolated to lower temps
        tmin=tmstn(i)
        tmax=txstn(i)
        Dmax=-9.99d99     !density and pressure limits not applicable
        pmax=-9.99d99     !for surface tension--set to large number
      else if (htyp.eq.'TCX' .or. htyp.eq.'tcx') then
c  thermal conductivity
        ltemp=.true.    !transport cannot be extrapolated to lower temps
        tmin=tmtcx(i)
        tmax=txtcx(i)
        Dmax=Dtcx(i)
        pmax=ptcx(i)
      else if (htyp.eq.'ETA' .or. htyp.eq.'eta') then
c  viscosity
        ltemp=.true.    !transport cannot be extrapolated to lower temps
        tmin=tmeta(i)
        tmax=txeta(i)
        Dmax=Deta(i)
        pmax=peta(i)
      else
c  unknown model specification--use EOS limits
        ltemp=.true.
        tmin=tmeos(i)
        tmax=txeos(i)
        Dmax=Deos(i)
        pmax=peos(i)
      end if
c
c  check inputs against limits
c
      if (t.lt.tmin-1.0d-10 .and. abs(t).gt.1.d-20) then
        ierr=0
        if (i.ne.0) call SUBLK (icomp,t,pm,ierr,herr)
c  set error number if p>pmelt.  If p=0, also set error number
        if (ierr.ne.0 .or. p.gt.pm+1.0d-10 .or. p.lt.1.d-20) then
          if (ltemp) then
            lerr=.true.
            ierrt=1
          else
            lwarn=.true.
            ierrt=-1
          end if
          write (herrt,1010) t,tmin
          nchart=70
 1010     format (' temperature below lower limit, T =',g11.5,
     &            ' K, Tmin =',g11.5,' K;')
        endif
      else if (t.gt.1.5d0*tmax) then
        lerr=.true.
        ierrt=1
        write (herrt,1011) t,tmax
        nchart=72
 1011   format (' temperature > 1.5 x upper limit, T =',g11.5,
     &          ' K, Tmax =',g11.5,' K;')
      else if (t.gt.tmax) then
c       if (hcasn(icomp).eq.'7439-90-9') then
c         lerr=.true.
c         ierrt=1
c       else
          lwarn=.true.
          ierrt=-1
c       endif
        write (herrt,1012) t,tmax
        nchart=70
 1012   format (' temperature above upper limit, T =',g11.5,
     &          ' K, Tmax =',g11.5,' K;')
      end if
      if (D.lt.0.0d0) then
        lerr=.true.
        ierrd=2
        write (herrd,1020) D
        nchard=36
 1020   format (' density < 0, D =',g11.5,' mol/L;')
      else if (D.gt.Dmax) then
        lerr=.true.
        ierrd=2
        write (herrd,1021) D,Dmax
        nchard=74
 1021   format (' density above upper limit, D =',g11.5,
     &          ' mol/L, Dmax =',g11.5,' mol/L;')
      end if
      if (p.lt.0.0d0) then
        lerr=.true.
        ierrp=4
        write (herrp,1040) p/1000.0d0
        ncharp=34
 1040   format (' pressure < 0, P =',g11.5,' MPa;')
      else if (p.gt.2.0d0*pmax) then
        lerr=.true.
        ierrp=4
        write (herrp,1042) p/1000.0d0,pmax/1000.0d0
        ncharp=71
 1042   format (' pressure > 2 x upper limit, P =',g11.5,
     &          ' MPa, Pmax =',g11.5,' MPa;')
      else if (p.gt.pmax+1.0d-8) then
        if (hcasn(icomp).eq.'1333-74-0p') then
          lerr=.true.
          ierrp=4
        else
          lwarn=.true.
          ierrp=-4
        endif
        write (herrp,1044) p/1000.0d0,pmax/1000.0d0
        ncharp=71
 1044   format (' pressure above upper limit, P =',g11.5,
     &          ' MPa, Pmax =',g11.5,' MPa;')
      end if
      if (p*0.9995d0.gt.pmlt .and. abs(pmlt).gt.1.d-20) then
        lerr=.true.
        ierrp=16
        write (herrp,1050) p/1000.0d0,pmlt/1000.0d0
        ncharp=72
 1050   format (' pressure > melting pressure, P =',g11.5,
     &          ' MPa, Pmelt =',g11.5,' MPa;')
      end if
c
      if (icomp.ne.0) then
        if (t.lt.ttp(icomp)-1.d-10 .and. ierrt.eq.0.and.t.gt.0.d0) then
c  check for cases where the temperature is less than the triple point
c  temperature, but still in a valid liquid region (like water between
c  251.165 and 273.16 K.)
          lwarn=.true.
          ierrt=-16
          write (herrt,1090) t,ttp(icomp)
          nchart=75
 1090     format (' temperature less than the triple point,',
     &          ' T=',g10.5,' K, Ttrp=',g10.5,' K;')
          if (hcasn(icomp).eq.'7732-18-5' .and. abs(p).gt.1.d-20) then
            call MLTH2O (t,p1,p2)
            if (p.lt.p2 .or. p.gt.p1) then
              lerr=.true.
              ierrt=16
              write (herrt,1095)
              nchart=34
 1095   format (' inputs are within the solid phase;')
            endif
          endif
        endif
      endif
c
c  compose error string and compute overall value of ierr
c
      if (lerr .or. lwarn) then
c       write (*,*) ' LIMITK--nchart,d,p,x:',nchart,nchard,nchard,ncharx
c       write (*,1999) ' LIMITK--herrt: ',herrt(1:nchart)
c       write (*,1999) ' LIMITK--herrd: ',herrd(1:nchard)
c       write (*,1999) ' LIMITK--herrp: ',herrp(1:ncharp)
c       write (*,1999) ' LIMITK--herrx: ',herrx(1:ncharx)
c1999   format (1x,a16,a80)
c       write (*,*) ' LIMITK--#char: t,d,p,x,sum:  ',nchart,nchard,
c    &              ncharp,ncharx,nchart+nchard+ncharp+ncharx
        herr='one or more inputs are out of range: '//herrt(1:nchart)
     &      //herrd(1:nchard)//herrp(1:ncharp)//herrx(1:ncharx)
     &      //hnull
        if (lerr) then
          ierr=abs(ierrt)+abs(ierrd)+abs(ierrp)
        else if (lwarn) then
          ierr=-abs(ierrt)-abs(ierrd)-abs(ierrp)
        end if
        call ERRMSG (ierr,herr)
      end if
c
      RETURN
      end                                             !subroutine LIMITK
c
c ======================================================================
c
      subroutine LIMITS (htyp,x,tmin,tmax,Dmax,pmax)
c
c  returns limits of a property model as a function of composition
c
c  Pure fluid limits are read in from the .fld files; for mixtures, a
c  simple mole fraction weighting in reduced variables is used.
c
c  inputs:
c     htyp--flag indicating which models are to be checked [character*3]
c           'EOS':  equation of state for thermodynamic properties
c           'ETA':  viscosity
c           'TCX':  thermal conductivity
c           'STN':  surface tension
c        x--composition array [mol frac]
c  outputs:
c     tmin--minimum temperature for model specified by htyp [K]
c     tmax--maximum temperature [K]
c     Dmax--maximum density [mol/L]
c     pmax--maximum pressure [kPa]
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  10-22-99 EWL, original version
c  04-21-08 EWL, do not subtract 5 K if tmin<10 K
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: LIMITS
c     dll_export LIMITS
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (nrf0=n0)    !lower limit for transport ref fluid arrays
      dimension x(ncmax)
      dimension tred(ncmax),Dred(ncmax)
      dimension tmn(ncmax),tmx(ncmax),Dmx(ncmax),pmx(ncmax)
      character*3 htyp
      character*3 heta,hetak,htcx,htcxk
      character*12 hcasn
      common /NCOMP/ nc,ic
      common /CCAS/ hcasn(n0:nx)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /EOSLIM/ tmeos(n0:nx),txeos(n0:nx),peos(n0:nx),Deos(n0:nx)
      common /WLMTCX/ tmtcx(nrf0:nx),txtcx(nrf0:nx),ptcx(nrf0:nx),
     &                Dtcx(nrf0:nx)
      common /WLMETA/ tmeta(nrf0:nx),txeta(nrf0:nx),peta(nrf0:nx),
     &                Deta(nrf0:nx)
      common /WLMTRN/ tmecs(nrf0:nx),txecs(nrf0:nx),pecs(nrf0:nx),
     &                Decs(nrf0:nx)
      common /WLMSTN/ tmstn(n0:nx),txstn(n0:nx)
      common /TRNMOD/ heta,hetak(nrf0:nx),htcx,htcxk(nrf0:nx)
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
      if (htyp.eq.'EOS' .or. htyp.eq.'eos') then
c  equation of state
        do i=1,nc
          tmn(i)=tmeos(i)
          tmx(i)=txeos(i)
          Dmx(i)=Deos(i)
          pmx(i)=peos(i)
        enddo
      else if (htyp.eq.'STN' .or. htyp.eq.'stn') then
c  surface tension model
        do i=1,nc
          tmn(i)=tmstn(i)
          tmx(i)=txstn(i)
          Dmx(i)=9.99d99     !density and pressure limits not applicable
          pmx(i)=9.99d99     !for surface tension--set to large number
c         write (*,*) ' LIMITX-STN; i,Tmin,Tmax: ',i,tmn(i),tmx(i)
        enddo
      else if (htyp.eq.'TCX' .or. htyp.eq.'tcx') then
c  thermal conductivity
        do i=1,nc
          tmn(i)=tmtcx(i)
          tmx(i)=txtcx(i)
          Dmx(i)=Dtcx(i)
          pmx(i)=ptcx(i)
          if (htcxk(i).eq.'ECS') then
            tmn(i)=tmecs(i)
            tmx(i)=txecs(i)
            Dmx(i)=Decs(i)
            pmx(i)=pecs(i)
          endif
        enddo
      else if (htyp.eq.'ETA' .or. htyp.eq.'eta') then
c  viscosity
        do i=1,nc
          tmn(i)=tmeta(i)
          tmx(i)=txeta(i)
          Dmx(i)=Deta(i)
          pmx(i)=peta(i)
          if (hetak(i).eq.'ECS') then
            tmn(i)=tmecs(i)
            tmx(i)=txecs(i)
            Dmx(i)=Decs(i)
            pmx(i)=pecs(i)
          endif
        enddo
      else
c  unknown model specification--use EOS limits
        do i=1,nc
          tmn(i)=tmeos(i)
          tmx(i)=txeos(i)
          Dmx(i)=Deos(i)
          pmx(i)=peos(i)
        enddo
      end if
c
      call ISPURE (x,icomp)
      if (icomp.ne.0) then
c  special case--pure component
        tmin=tmn(icomp)
        tmax=tmx(icomp)
        Dmax=Dmx(icomp)
        pmax=pmx(icomp)
        if (hcasn(icomp).eq.'7732-18-5' .and. tmin.gt.251.165d0)
     &      tmin=251.165d0
      else
c  general mixture case
        taumin=0.0d0
        taumax=0.0d0
        delmax=0.0d0
        pmax=0.0d0
        do i=1,nc
          tred(i)=tz(i)
          Dred(i)=rhoz(i)
          if (tmn(i).gt.0.d0) taumin=taumin+x(i)*tred(i)/tmn(i)
          if (tmx(i).gt.0.d0) taumax=taumax+x(i)*tred(i)/tmx(i)
          if (Dred(i).gt.0.d0) delmax=delmax+x(i)*Dmx(i)/Dred(i)
          pmax=pmax+x(i)*pmx(i)
        enddo
        call REDX (x,trmix,Drmix)
        Dmax=Drmix*delmax
        if (taumin.le.0.d0) taumin=100.d0
        if (taumax.le.0.d0) taumax=300.d0
        if (taumax.gt.0.d0) tmax=trmix/taumax
        if (taumin.gt.0.d0) then
          tmin=trmix/taumin
          if (tmin.gt.10.d0) tmin=tmin-5.d0
        endif
        if (tmin.lt.0.d0) tmin=100.d0
        if (tmax.lt.0.d0) tmax=300.d0
        if (iamwat.gt.0) then
          i=iamwat
          if (x(i).lt.0.33367d0) then
            tmin=273.16d0*(1.0d0-0.3439823d0*x(i)
     &          -1.3274271d0*x(i)**2-274.973d0*x(i)**7)
          elseif (x(i).lt.0.58396d0) then
            tmin=193.549d0*(1.0d0-4.987368d0*(x(i)-0.5d0)**2)
          elseif (x(i).lt.0.81473d0) then
            tmin=194.380d0*(1.0d0-4.886151d0*(x(i)-2.0d0/3.0d0)**2
     &          +10.37298d0*(x(i)-2.0d0/3.0d0)**3)
          else
            tmin=195.495d0*(1.0d0-0.323998d0*(1.0d0-x(i))
     &          -15.87560d0*(1.0d0-x(i))**4)
          endif
        endif
      end if
c
      RETURN
      end                                             !subroutine LIMITS
c ======================================================================
c
      subroutine ERRMSG (ierr,herr)
c
c  write error messages to default output; this subroutine should be
c  called immediately after any call to a subroutine which potentially
c  can error out
c
c  inputs:
c     ierr--error flag:  0 = successful (no message will be written)
c                       <0 = warning
c                       >0 = error
c     herr--error string (character*255 variable)
c
c  outputs:
c     if iprnterr in common block prnterr is equal to zero:
c     error string written to default output
c
c     if iprnterr is equal to 1:
c     error string written to screen if ierr is positive
c
c     if iprnterr is equal to -1:
c     error string written to screen
c
c     if iprnterr is equal to 3, -3:
c     same as 1, -1, but program also pauses
c
c     Note:  no information should be written to the screen
c     when compiling the DLL.
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  09-11-95  MM, original version
c  07-30-96  MM, also write out the value of ierr
c  03-26-97 EWL, change name ERROR --> ERRMSG as 'error' is a standard
c                routine in Lahey Fortran90
c  10-02-97  MM, add compiler switch to allow access by DLL
c  12-01-98 EWL, add common block prterr
c  12-01-98 EWL, don't print extra spaces at end of herr
c  10-21-99 EWL, remove use of i outside of do loop
c  10-22-99 EWL, remove extra spaces before "K", "MPa", etc.
c  11-17-99  MM, new variable for printed message
c
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: ERRMSG
c     dll_export ERRMSG
c
      character*1 htab,hnull
      common /HCHAR/ htab,hnull
      common /prnterr/ iprnterr
      character*255 herr,hout
      character*4 hstrg
c
      j=0
 100  continue
      i=index(herr,'   K')
      if (i.eq.0) i=index(herr,'   MPa')
      if (i.eq.0) i=index(herr,'   mol')
      if (i.eq.0) i=index(herr,'   J/mol')
      if (i.ne.0) then
        j=j+1
        herr=herr(1:i)//herr(i+3:255)
        if (j.lt.10) goto 100
      endif
      i=index(herr,'=0')
      if (i.ne.0) then
        j=j+1
        herr=herr(1:i)//' '//herr(i+1:255)
        if (j.lt.10) goto 100
      endif
      i=index(herr,';   ')
      if (i.ne.0) herr=herr(1:i-1)//herr(i+1:255)
      i=index(herr,'; '//hnull)
      if (i.ne.0) herr=herr(1:i-1)//herr(i+1:255)
      i=index(herr,';  '//hnull)
      if (i.ne.0) herr=herr(1:i-1)//herr(i+1:255)
c
      if (ierr.ne.0 .and. iprnterr.ne.0) then
        do i=255,1,-1
          if (herr(i:i).ne.' ' .and. herr(i:i).ne.hnull) then
            ilast=i
            goto 110
          end if
        enddo
        ilast=1
 110    continue
        write (hstrg,'(i4)') ierr
        if (ierr.gt.0 .or. iprnterr.lt.0) then
          hout=herr(1:ilast)//' (ierr='//hstrg//')'
          write (*,1000) hout
c         if (ABS(iprnterr).ge.2)  call error('Called')
          if (ABS(iprnterr).ge.3)  pause
        end if
      end if
 1000 format (1x,a255)
c
      RETURN
      end                                             !subroutine ERRMSG
c
c ======================================================================
c
      subroutine QMASS (qmol,xl,xv,qkg,xlkg,xvkg,wliq,wvap,ierr,herr)
c
c  converts quality and composition on a mole basis to a mass basis
c
c  inputs:
c     qmol--molar quality [moles vapor/total moles]
c           qmol = 0 indicates saturated liquid
c           qmol = 1 indicates saturated vapor
c           0 < qmol < 1 indicates a two-phase state
c           qmol < 0 or qmol > 1 are not allowed and will result in warning
c       xl--composition of liquid phase [array of mol frac]
c       xv--composition of vapor phase [array of mol frac]
c  outputs:
c      qkg--quality on mass basis [mass of vapor/total mass]
c     xlkg--mass composition of liquid phase [array of mass frac]
c     xvkg--mass composition of vapor phase [array of mass frac]
c     wliq--molecular weight of liquid phase [g/mol]
c     wvap--molecular weight of vapor phase [g/mol]
c     ierr--error flag:  0 = all inputs within limits
c           -19:  input q < 0 or > 1
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by M. McLinden, NIST Phys & Chem Properties Div, Boulder, CO
c  02-23-98  MM, original version, based on XMASS
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: QMASS
c     dll_export QMASS
c
      parameter (ncmax=20)        !max number of components in mixture
      dimension xl(ncmax),xv(ncmax),xlkg(ncmax),xvkg(ncmax)
      character*1 htab,hnull
      character*255 herr
      common /HCHAR/ htab,hnull
      data eps /1.0d-8/
c
      call XMASS (xl,xlkg,wliq)
      call XMASS (xv,xvkg,wvap)
      if (qmol.lt.-eps .or. qmol.gt.1.0d0+eps) then
        ierr=-19
        herr='[QMASS warning 19] input quality out of range'//hnull
        call ERRMSG (ierr,herr)
        qkg=qmol
      else
        ierr=0
        herr=' '
        qkg=qmol*wvap/(qmol*wvap+(1.0d0-qmol)*wliq)
      end if
c
      RETURN
      end                                              !subroutine QMASS
c
c ======================================================================
c
      subroutine QMOLE (qkg,xlkg,xvkg,qmol,xl,xv,wliq,wvap,ierr,herr)
c
c  converts quality and composition on a mass basis to a molar basis
c
c  inputs:
c      qkg--quality on mass basis [mass of vapor/total mass]
c           qkg = 0 indicates saturated liquid
c           qkg = 1 indicates saturated vapor
c           0 < qkg < 1 indicates a two-phase state
c           qkg < 0 or qkg > 1 are not allowed and will result in warning
c     xlkg--mass composition of liquid phase [array of mass frac]
c     xvkg--mass composition of vapor phase [array of mass frac]
c  outputs:
c     qmol--quality on mass basis [mass of vapor/total mass]
c       xl--molar composition of liquid phase [array of mol frac]
c       xv--molar composition of vapor phase [array of mol frac]
c     wliq--molecular weight of liquid phase [g/mol]
c     wvap--molecular weight of vapor phase [g/mol]
c     ierr--error flag:  0 = all inputs within limits
c           -19:  input q < 0 or > 1
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by M. McLinden, NIST Phys & Chem Properties Div, Boulder, CO
c  02-23-98  MM, original version, based on XMOLE
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: QMOLE
c     dll_export QMOLE
c
      parameter (ncmax=20)        !max number of components in mixture
      dimension xl(ncmax),xv(ncmax),xlkg(ncmax),xvkg(ncmax)
      character*1 htab,hnull
      character*255 herr
      common /HCHAR/ htab,hnull
      data eps /1.0d-8/
c
      call XMOLE (xlkg,xl,wliq)
      call XMOLE (xvkg,xv,wvap)
      if (qkg.lt.-eps .or. qkg.gt.1.0d0+eps) then
        ierr=-19
        herr='[QMOLE warning 19] input quality out of range'//hnull
        call ERRMSG (ierr,herr)
        qmol=qkg
      else
        ierr=0
        herr=' '
        qmol=qkg/wvap/(qkg/wvap+(1.0d0-qkg)/wliq)
      end if
c
      RETURN
      end                                              !subroutine QMOLE
c
c ======================================================================
c
      subroutine GOLD (x0i,x1i,nc,lmax,z,z2,z3,z4,bt,xopt,yopt,ierr)
c
c  This subroutine carries out a Fibonnaci (golden) search
c  technique to locate an extremum (maximum or minimum) in a
c  function within a specified interval.
c
c  inputs:
c      x0i--lower bound of interval containing extremum
c      x1i--upper bound of interval containing extremum
c       nc--number of function evaluations carried out;
c           the extremum is located with an interval of size:
c           (x1i - x0i)*0.618**nc
c           (e.g. 29 evaluations required to reduce interval to 10**-6
c           of its original size)
c     lmax--logical flag; if lmax =
c           .true. -  locate maximum value
c           .false. - locate minimum value
c        z--composition array
c       z2--additional independent variable
c       z3--additional independent variable
c       z4--additional independent variable
c       bt--function type, 'H'-call ENTHAL, etc.
c  outputs:
c     xopt--location of extremum
c     yopt--value of function at extremum
c     ierr--integer output from external function (e.g. error flag)
c
c  taken from:
c    McLinden, M.O. (1988). Working fluid selection for space-based
c    two-phase heat transport systems. National Bureau of Standards,
c    NBSIR 88-3812.
c  by M. McLinden, NIST Chemical & Physical Properties Div, Boulder, CO
c  09-05-97  MM, transcribed from original reference, add additional comments
c                add z,z2,z3,z4 for passing to GEVAL
c  03-26-98  MM, make z an array (e.g. for passing composition)
c  01-19-01 EWL, remove external dependence
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      logical lmax
      parameter (ncmax=20)        !max number of components in mixture
      dimension z(ncmax)
      character*1 bt
      data GR /0.61803398875/
c
c  routine always finds a maximum, so to locate a minimum multiply by -1
      ierr=0
      y38=0.d0
      if (lmax) then
        xmax=1.0d0
      else
        xmax=-1.0d0
      end if
      x0=x0i  !x0 stores the lower bound of interval containing extremum
      x1=x1i  !x0 stores the upper bound
      x62=x0+GR*(x1-x0)          !x62 is 62% of way across interval
      if (bt.eq.'H') then
        call ENTHAL (z2,x62,z,b)
      elseif (bt.eq.'E') then
        call ENERGY (z2,x62,z,b)
      elseif (bt.eq.'S') then
        call ENTRO (z2,x62,z,b)
      endif
      y62=xmax*b
c
      do ig=1,nc
        x38=x0+(1.0d0-GR)*(x1-x0)  !x38 is 38% of way across interval
        if (bt.eq.'H') then
          call ENTHAL (z2,x38,z,b)
        elseif (bt.eq.'E') then
          call ENERGY (z2,x38,z,b)
        elseif (bt.eq.'s') then
          call ENTRO (z2,x38,z,b)
        endif
        y38=xmax*b
        if (y62.lt.y38) then
c  keep sub-interval containing maximum value of GEVAL
          x1=x62
          x62=x38
          y62=y38
        else
          x0=x1
          x1=x38
        end if
      enddo
c
      xopt=0.5d0*(x0+x1)
      yopt=xmax*MAX(y38,y62)
c
      RETURN
      end                                               !subroutine GOLD
c
c ======================================================================
c
      function UCASE (b,l)
c
c  Make all the characters in the string b uppercase, from the first
c  character to character l
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      character ucase*255, a*255, b*(*)
      a=b
      do i=1,l
        j=ichar(a(i:i))
        if (j.gt.96 .and. j.le.122) j=j-32
        a(i:i)=char(j)
      enddo
      UCASE=a
      RETURN
      end                                                !function UCASE
c
c ======================================================================
c
      subroutine DOTFILL (x,ptest,filrat,ierr,herr)
c
c  Calculate filling ratio according to UN P200 document
c  Packing instructions for hazardous substances
c
c  Input
c    x         composition, mol fraction
c    ptest     test pressure, absolute kPa  (may also be output variable)
C
c  Output
c    filrat    filling ratio
c    ierr      error flag
c    herr      error message
c
c  Version 0.0   10.07.05
c
      implicit DOUBLE PRECISION (a-h,o-z)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DOTFILL
c     dll_export DOTFILL
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),y(ncmax),z(ncmax)
      character herr*255
      character*1 htab,hnull
      INTEGER unflag
      common /NCOMP/ ncomp,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /CONUN/ tminUN(n0:nx),tmaxUN(n0:nx),pmaxUN(n0:nx),
     &               rhomUN(n0:nx),prmUN(n0:nx,10),ntrmUN(n0:nx),
     &               unflag(n0:nx)
c
c
      herr=' '
      ierr=0
      rho0=999.10262d0 !kg/m3 at 15C, 0.101325 MPa (water)
      filrat=0.0d0
      ione=1
      ttest65=65.0d0+273.15d0
      ttest60=60.0d0+273.15d0
      ttest50=50.0d0+273.15d0

      if (unflag(1).ne.1) then !do not allow calculations
        filrat=0.0d0
        ptest= 0.0d0
        ierr=1
        herr='User defined calculations not permitted for this fluid'
      else
c  obtain critical temperature to determine if fluid is high or low pressure
c  and get molecular weight of mixture
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
          tcc=tcrit(icomp)
          wmx=wm(icomp)
        else
          call CRTHMX (x,tcc,pcmix,Dcmix,ierr,herr)
          wmx=wmol(x)
        endif
        tcc=tcc-273.15d0 !convert to C
c  get the saturation pressure at 65 C
        call SATT (ttest65,x,ione,psat65,dl,dv,xliq,xvap,ierr,herr)
c
c  begin procedure for low pressure fluids (tc<65C)
c  for these fluids, there is only one choice for test pressure- Psat at 65C
        if (tcc.gt.65.d0) then
c  obtain the density at 50C and saturation
          call SATT (ttest50,x,ione,p,dl,dv,xliq,xvap,ierr,herr)
          rhotest=0.95d0*dl !mol/L
c  perform first check; get sat liquid rho at 60C
          call SATT (ttest60,x,ione,p,dl,dv,xliq,xvap,ierr,herr)
          rhof60=dl !mol/L
c  select the lower of these two values as the fill density
          rhofill=MIN(rhof60,rhotest)
          rhofill=rhofill*wmx  !convert to kg/m3
          if (rho0.gt.0.d0) filrat=rhofill/rho0
          ptest=psat65 !kPa,abs remember to change to gauge for final results
        else
c  begin procedure for High pressure fluids (tc>65C)
          call TPFLSH (ttest65,ptest,z,D,Dl,Dv,x,y,q,e,h,s,cv,cp,w,
     &                 ierr,herr)
          rhof=d*wmx    !convert to kg/m3
          if (rho0.gt.0.d0) filrat=rhof/rho0
        endif
c
c     adjust to match accepted un values
        if (ntrmun(1).gt.0) then
          if (icomp.ne.0) filrat=filrat*prmUN(icomp,icomp)
        endif  !won't work for mixtures, must revise
c
c  test pressure must exceed saturation pressure at 65C
        if (ptest.lt.psat65) then
          ierr=-1
          herr='Test pressure must exceed saturation pressure at 65 C'
          ptest=psat65
          filrat=-999
        endif
      endif
      RETURN
      end                                            !subroutine DOTFILL
c
c ======================================================================
c
      function CBRTX (x)
c
c  cube root function--allows negative arguments
c
c  input:
c        x--value to be acted upon
c  output (as function value):
c     CBRTX--cube root of x
c
c  written by M. McLinden, NIST Phys & Chem Properties Div, Boulder, CO
c  08-25-97 MM, original version
c  11-06-01 MLH, changed name to CBRTX to be standard conforming.
c  09-29-04 MLH, fixed bug for neg arguments
c
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (thrd=1.0d0/3.0d0)
c
      xa=ABS(x)
      xathrd=xa**thrd
      if (x.ge.0.0d0) then
        cbrtx=xathrd
      else
        cbrtx=-xathrd
      end if
c
      RETURN
      end                                                !function CBRTX

c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c ======================================================================
c                                                     end file utility.f
c ======================================================================
