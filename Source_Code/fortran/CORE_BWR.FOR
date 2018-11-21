c  begin file core_BWR.f
c
c  This file contains the functions implementing the MBWR equation of
c  state for pure fluids.
c
c  contained here are:
c     function PBWR (icomp,t,rho)
c     function ABWR (icomp,t,rho)
c     function DABWR (icomp,t,rho)
c     function D2ABWR (icomp,t,rho)
c     function DPDBWR (icomp,t,rho)
c     function D2PBWR (icomp,t,rho)
c     function DPTBWR (icomp,t,rho)
c     function PHIBWR (icomp,itau,idel,tau,del)
c     subroutine CRTBWR (icomp,tcrit,pcrit,Dcrit)
c     subroutine SETBWR (nread,icomp,hcasno,ierr,herr)
c     block data BDBWR
c
c  these routines set and/or use values in the following common blocks
c
c  commons associated with the MBWR equations stored in block data
c     common /CASBWR/ hcas(mxbwr)
c     common /CPMBWR/ icpbwr(mxbwr)
c     common /CFBWR/ ba(mxbwr,32),
c    &               pca(mxbwr),rhoca(mxbwr),tca(mxbwr),
c    &               Rbwra(mxbwr),ptra(mxbwr),rhotra(mxbwr),ttra(mxbwr),
c    &               gammaa(mxbwr),tmaxa(mxbwr),pmaxa(mxbwr)
c     common /MSCBWR/ wmb(mxbwr),ttpb(mxbwr),tnbpb(mxbwr),accenb(mxbwr),
c    &                dipm(mxbwr)
c  commons associated with the nc components of current interest
c  ("working" commons and arrays)
c     common /WCFBWR/ b(n0:nx,32),
c    &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),
c    &                Rbwr(n0:nx),pmin(n0:nx),rhomin(n0:nx),tmin(n0:nx),
c    &                gamma(n0:nx),tmax(n0:nx),pmax(n0:nx)
c     common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
c    &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
c    &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
c    &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
c  these routines use the following common blocks from other files
c     common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c     common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
c
c  various arrays are dimensioned with parameter statements
c     parameter (mxbwr=2)         !max number of MBWR EOS in block data
c     parameter (ncmax=20)        !max number of components in mixture
c     parameter (nrefmx=10)       !max number of fluids for transport ECS
c     parameter (n0=-ncmax-nrefmx,nx=ncmax)
c
c ======================================================================
c ======================================================================
c
      function PBWR (icomp,t,rho)
c
c  compute pressure as a function of temperature and density
c  for MBWR equation of state, this is the most basic form of MBWR
c
c  based on Younglove & McLinden (1994), JPCRD 23:731-779
c  equations B1 and B2
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c        t--temperature (K)
c      rho--molar density (mol/L)
c  output: (as function value):
c        p--pressure (kPa)
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-06-94  MM, original version
c  07-20-95  MM, separate polynomial Cp0 routines from core_BWR
c                and restructure coefficient arrays
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFBWR/ b(n0:nx,32),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),Rbwr(n0:nx),
     &                ptr(n0:nx),rhotr(n0:nx),ttr(n0:nx),gamma(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension an(15)
c
      j=icomp              !pointer to appropriate fluid
      PBWR=0.d0
      if (t.le.0.d0) RETURN
c
c  calculate a(n) terms in MBWR (Eq B2)
c
      tinv=1.0d0/t
      tinv2=tinv*tinv
      tinv3=tinv2*tinv
      tinv4=tinv3*tinv
c
      an(1)=Rbwr(j)*t
      an(2)= b(j,1)*t+b(j,2)*SQRT(t)+b(j,3)+b(j,4)*tinv+b(j,5)*tinv2
      an(3)= b(j,6)*t+b(j,7)+b(j,8)*tinv+b(j,9)*tinv2
      an(4)=b(j,10)*t+b(j,11)+b(j,12)*tinv
      an(5)=b(j,13)
      an(6)=b(j,14)*tinv+b(j,15)*tinv2
      an(7)=b(j,16)*tinv
      an(8)=b(j,17)*tinv+b(j,18)*tinv2
      an(9)=b(j,19)*tinv2
      an(10)=b(j,20)*tinv2+b(j,21)*tinv3
      an(11)=b(j,22)*tinv2+b(j,23)*tinv4
      an(12)=b(j,24)*tinv2+b(j,25)*tinv3
      an(13)=b(j,26)*tinv2+b(j,27)*tinv4
      an(14)=b(j,28)*tinv2+b(j,29)*tinv3
      an(15)=b(j,30)*tinv2+b(j,31)*tinv3+b(j,32)*tinv4
c
c  summation of terms 1-9 in Eq B1
      psum=0.0d0
      rhon=1.0d0
      do n=1,9
        rhon=rhon*rho
        psum=psum+an(n)*rhon
      enddo
c  summation of terms 10-15 (exponential terms) in Eq B1
      expsum=0.0d0
      rhon=rho
      rho2=rho*rho
      do n=10,15
        rhon=rhon*rho2
        expsum=expsum+an(n)*rhon
      enddo
c
c  sum terms and convert from bar to kPa
      PBWR=(psum+exp(-(rho/gamma(j))**2)*expsum)*100.0d0
c     PBWR=(psum+exp(-(rho/gamma(j))**2)*expsum)*(R/Rbwr(j))
c
      RETURN
      end                                                 !function PBWR
c
c ======================================================================
c
      function ABWR (icomp,t,rho)
c
c  compute residual Helmholtz free energy for MBWR equation of state
c
c  based on Younglove & McLinden (1994), JPCRD 23:731-779
c  equations A4 and B6
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c        t--temperature (K)
c      rho--molar density (mol/L)
c  output (as function value):
c       Ar--(A-A0) (J/mol)
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-06-94  MM, original version
c  07-20-95  MM, separate polynomial Cp0 routines from core_BWR
c                and restructure coefficient arrays
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFBWR/ b(n0:nx,32),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),Rbwr(n0:nx),
     &                ptr(n0:nx),rhotr(n0:nx),ttr(n0:nx),gamma(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension an(15)
c
      j=icomp              !pointer to appropriate fluid
      ABWR=0.d0
      if (t.le.0.d0) RETURN
c
c  calculate a(n) terms in MBWR (Eq B2)
c
      tinv=1.0d0/t
      tinv2=tinv*tinv
      tinv3=tinv2*tinv
      tinv4=tinv3*tinv
c
      an(1)=Rbwr(j)*t
      an(2)= b(j,1)*t+b(j,2)*SQRT(t)+b(j,3)+b(j,4)*tinv+b(j,5)*tinv2
      an(3)= b(j,6)*t+b(j,7)+b(j,8)*tinv+b(j,9)*tinv2
      an(4)=b(j,10)*t+b(j,11)+b(j,12)*tinv
      an(5)=b(j,13)
      an(6)=b(j,14)*tinv+b(j,15)*tinv2
      an(7)=b(j,16)*tinv
      an(8)=b(j,17)*tinv+b(j,18)*tinv2
      an(9)=b(j,19)*tinv2
      an(10)=b(j,20)*tinv2+b(j,21)*tinv3
      an(11)=b(j,22)*tinv2+b(j,23)*tinv4
      an(12)=b(j,24)*tinv2+b(j,25)*tinv3
      an(13)=b(j,26)*tinv2+b(j,27)*tinv4
      an(14)=b(j,28)*tinv2+b(j,29)*tinv3
      an(15)=b(j,30)*tinv2+b(j,31)*tinv3+b(j,32)*tinv4
c
c  summation of terms 2-9 in Eq B6
      ar=0.0d0
      rhon=1.0d0
      do n=2,9
        rhon=rhon*rho
        ar=ar+an(n)*rhon/real(n-1)
      enddo
c  summation of term 10 (first exponential term)
      delsq=(rho/gamma(j))**2
      rhoc2=gamma(j)**2
      rhocn=rhoc2
      expdel=exp(-delsq)
      ar=ar-0.5d0*an(10)*rhocn*(expdel-1.0d0)
c  summation of terms 11-15 (remaining exponential terms) in Eq B7
      expmul=1.0d0
      expsub=1.0d0
      delnew=1.0d0
      do n=11,15
        rhocn=rhocn*rhoc2
        delnew=delnew*delsq
        expmul=delnew+expmul*real(n-10)
        expsub=expsub*real(n-10)
        ar=ar-0.5d0*an(n)*rhocn*(expdel*expmul-expsub)
      enddo
c
c  convert from L-bar/mol to J/mol
      ABWR=ar*(100.0d0)
c     ABWR=ar*(R/Rbwr(j))
c
      RETURN
      end                                                 !function ABWR
c
c ======================================================================
c
      function DABWR (icomp,t,rho)
c
c  compute temperature derivative of residual Helmholtz free energy
c  for MBWR equation of state
c
c  based on derivations in Younglove & McLinden (1994), JPCRD 23:731-779
c  equations B4 and B7
c
c  inputs:
c    icomp--pointer specifying fluid to calculate, refers to code numbers
c        t--temperature (K)
c      rho--molar density (mol/L)
c  output (as function value):
c    dardt--dAr/dT (J/(mol-K))
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-05-94  MM, original version
c  07-20-95  MM, separate polynomial Cp0 routines from core_BWR
c                and restructure coefficient arrays
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFBWR/ b(n0:nx,32),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),Rbwr(n0:nx),
     &                ptr(n0:nx),rhotr(n0:nx),ttr(n0:nx),gamma(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension dandt(15)
c
      j=icomp              !pointer to appropriate fluid
      DABWR=0.d0
      if (t.le.0.d0) RETURN
c
c  calculate temperature derivatives of a(n) terms in MBWR (Eq B4)
c
      tinv=1.0d0/t
      tinv2=tinv*tinv
      tinv3=tinv2*tinv
      tinv4=tinv3*tinv
      tinv5=tinv4*tinv
c
      dandt(1)=Rbwr(j)
      dandt(2)=b(j,1)+0.5d0*b(j,2)/SQRT(t)-b(j,4)*tinv2
     &         -2.0d0*b(j,5)*tinv3
      dandt(3)=b(j,6)-b(j,8)*tinv2-2.0d0*b(j,9)*tinv3
      dandt(4)=b(j,10)-b(j,12)*tinv2
      dandt(5)=0.0d0
      dandt(6)=-b(j,14)*tinv2-2.0d0*b(j,15)*tinv3
      dandt(7)=-b(j,16)*tinv2
      dandt(8)=-b(j,17)*tinv2-2.0d0*b(j,18)*tinv3
      dandt(9)=-2.0d0*b(j,19)*tinv3
      dandt(10)=-2.0d0*b(j,20)*tinv3-3.0d0*b(j,21)*tinv4
      dandt(11)=-2.0d0*b(j,22)*tinv3-4.0d0*b(j,23)*tinv5
      dandt(12)=-2.0d0*b(j,24)*tinv3-3.0d0*b(j,25)*tinv4
      dandt(13)=-2.0d0*b(j,26)*tinv3-4.0d0*b(j,27)*tinv5
      dandt(14)=-2.0d0*b(j,28)*tinv3-3.0d0*b(j,29)*tinv4
      dandt(15)=-2.0d0*b(j,30)*tinv3-3.0d0*b(j,31)*tinv4
     &          -4.0d0*b(j,32)*tinv5
c
c  summation of terms 2-9 in Eq B7
      dardt=0.0d0
      rhon=1.0d0
      do n=2,9
        rhon=rhon*rho
        dardt=dardt+dandt(n)*rhon/real(n-1)
      enddo
c  summation of term 10 (first exponential term)
      delsq=(rho/gamma(j))**2
      rhoc2=gamma(j)**2
      rhocn=rhoc2
      expdel=exp(-delsq)
      dardt=dardt-0.5d0*dandt(10)*rhocn*(expdel-1.0d0)
c  summation of terms 11-15 (remaining exponential terms) in Eq B7
      expmul=1.0d0
      expsub=1.0d0
      delnew=1.0d0
      do n=11,15
        rhocn=rhocn*rhoc2
        delnew=delnew*delsq
        expmul=delnew+expmul*real(n-10)
        expsub=expsub*real(n-10)
        dardt=dardt-0.5d0*dandt(n)*rhocn*(expdel*expmul-expsub)
      enddo
c
c  convert from L-bar/(mol-K) to J/(mol-K)
      DABWR=dardt*100.0d0
c     DABWR=dardt*(R/Rbwr(j))
c
      RETURN
      end                                                !function DABWR
c
c ======================================================================
c
      function D2ABWR (icomp,t,rho)
c
c  compute second temperature derivative of residual Helmholtz free
c  energy for MBWR equation of state
c
c  based on derivations in Younglove & McLinden (1994), JPCRD 23:731-779
c  equations B8 and B9
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c        t--temperature (K)
c      rho--molar density (mol/L)
c  output (as function value):
c   d2ardt--d2Ar/dT2 (J/(mol-K**2))
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-06-94  MM, original version
c  07-20-95  MM, separate polynomial Cp0 routines from core_BWR
c                and restructure coefficient arrays
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFBWR/ b(n0:nx,32),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),Rbwr(n0:nx),
     &                ptr(n0:nx),rhotr(n0:nx),ttr(n0:nx),gamma(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension d2andt(15)
c
      j=icomp              !pointer to appropriate fluid
c
c  calculate temperature derivatives of a(n) terms in MBWR (Eq B4)
c
      tinv=1.0d0/t
      tinv2=tinv*tinv
      tinv3=tinv2*tinv
      tinv4=tinv3*tinv
      tinv5=tinv4*tinv
      tinv6=tinv5*tinv
c
      d2andt(1)=0.0d0
      d2andt(2)=-0.25d0*b(j,2)/(t**1.5)+2.0d0*b(j,4)*tinv3
     &          +6.0d0*b(j,5)*tinv4
      d2andt(3)=2.0d0*b(j,8)*tinv3+6.0d0*b(j,9)*tinv4
      d2andt(4)=2.0d0*b(j,12)*tinv3
      d2andt(5)=0.0d0
      d2andt(6)=2.0d0*b(j,14)*tinv3+6.0d0*b(j,15)*tinv4
      d2andt(7)=2.0d0*b(j,16)*tinv3
      d2andt(8)=2.0d0*b(j,17)*tinv3+6.0d0*b(j,18)*tinv4
      d2andt(9)=6.0d0*b(j,19)*tinv4
      d2andt(10)=6.0d0*b(j,20)*tinv4+12.0d0*b(j,21)*tinv5
      d2andt(11)=6.0d0*b(j,22)*tinv4+20.0d0*b(j,23)*tinv6
      d2andt(12)=6.0d0*b(j,24)*tinv4+12.0d0*b(j,25)*tinv5
      d2andt(13)=6.0d0*b(j,26)*tinv4+20.0d0*b(j,27)*tinv6
      d2andt(14)=6.0d0*b(j,28)*tinv4+12.0d0*b(j,29)*tinv5
      d2andt(15)=6.0d0*b(j,30)*tinv4+12.0d0*b(j,31)*tinv5
     &           +20.0d0*b(j,32)*tinv6
c
c  summation of terms 2-9 in Eq B8
      d2ardt=0.0d0
      rhon=1.0d0
      do n=2,9
        rhon=rhon*rho
        d2ardt=d2ardt+d2andt(n)*rhon/real(n-1)
      enddo
c  summation of term 10 (first exponential term)
      delsq=(rho/gamma(j))**2
      rhoc2=gamma(j)**2
      rhocn=rhoc2
      expdel=exp(-delsq)
      d2ardt=d2ardt-0.5d0*d2andt(10)*rhocn*(expdel-1.0d0)
c  summation of terms 11-15 (remaining exponential terms) in Eq B7
      expmul=1.0d0
      expsub=1.0d0
      delnew=1.0d0
      do n=11,15
        rhocn=rhocn*rhoc2
        delnew=delnew*delsq
        expmul=delnew+expmul*real(n-10)
        expsub=expsub*real(n-10)
        d2ardt=d2ardt-0.5d0*d2andt(n)*rhocn*(expdel*expmul-expsub)
      enddo
c
c  convert from L-bar/(mol-K**2) to J/(mol-K**2)
      D2ABWR=d2ardt*(100.0d0)
c     D2ABWR=d2ardt*(R/Rbwr(j))
c
      RETURN
      end                                               !function D2ABWR
c
c ======================================================================
c
      function DPDBWR (icomp,t,rho)
c
c  compute partial derivative of pressure with respect to density
c  for MBWR equation of state
c
c  based on derivations in Younglove & McLinden (1994), JPCRD 23:731-779
c  equation B5
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c        t--temperature (K)
c      rho--molar density (mol/L)
c  output (as function value):
c   DPDBWR--dPdD (kPa-L/mol)
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-07-94  MM, original version
c  07-20-95  MM, separate polynomial Cp0 routines from core_BWR
c                and restructure coefficient arrays
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c  06-06-97 EWL, check if rho < 1.0d-10, avoid divide by zero
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFBWR/ b(n0:nx,32),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),Rbwr(n0:nx),
     &                ptr(n0:nx),rhotr(n0:nx),ttr(n0:nx),gamma(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension an(15)
c
      DPDBWR=R*t           !ideal-gas
      if (rho.lt.1.0d-10 .or. t.le.0.d0) then
        RETURN             !use ideal-gas solution
      end if
      j=icomp              !pointer to appropriate fluid
c
c  calculate a(n) terms in MBWR (Eq B2)
c
      tinv=1.0d0/t
      tinv2=tinv*tinv
      tinv3=tinv2*tinv
      tinv4=tinv3*tinv
c
      an(1)=Rbwr(j)*t
      an(2)= b(j,1)*t+b(j,2)*SQRT(t)+b(j,3)+b(j,4)*tinv+b(j,5)*tinv2
      an(3)= b(j,6)*t+b(j,7)+b(j,8)*tinv+b(j,9)*tinv2
      an(4)=b(j,10)*t+b(j,11)+b(j,12)*tinv
      an(5)=b(j,13)
      an(6)=b(j,14)*tinv+b(j,15)*tinv2
      an(7)=b(j,16)*tinv
      an(8)=b(j,17)*tinv+b(j,18)*tinv2
      an(9)=b(j,19)*tinv2
      an(10)=b(j,20)*tinv2+b(j,21)*tinv3
      an(11)=b(j,22)*tinv2+b(j,23)*tinv4
      an(12)=b(j,24)*tinv2+b(j,25)*tinv3
      an(13)=b(j,26)*tinv2+b(j,27)*tinv4
      an(14)=b(j,28)*tinv2+b(j,29)*tinv3
      an(15)=b(j,30)*tinv2+b(j,31)*tinv3+b(j,32)*tinv4
c
c  summation of terms 1-9 in Eq B1
      psum=0.0d0
      rhon=1.0d0/rho
      do n=1,9
        rhon=rhon*rho
        psum=psum+real(n)*an(n)*rhon
      enddo
c  summation of terms 10-15 (exponential terms) in Eq B1
      expsum=0.0d0
      rhon=1.0d0
      rho2=rho*rho
      delsq=(rho/gamma(j))**2
      do n=10,15
        rhon=rhon*rho2
        expsum=expsum+an(n)*rhon*(2.0d0*real(n)-17.0d0-2.0d0*delsq)
      enddo
c
c  collect terms and convert from L-bar/mol to kPa-L/mol
      DPDBWR=(psum+exp(-delsq)*expsum)*(100.0d0)
c     DPDBWR=(psum+exp(-delsq)*expsum)*(R/Rbwr(j))
c
      RETURN
      end                                               !function DPDBWR
c
c ======================================================================
c
      function D2PBWR (icomp,t,rho)
c
c  compute the second partial derivative of pressure with respect to density
c  for the MBWR equation of state
c
c  based on derivations in Younglove & McLinden (1994), JPCRD 23:731-779
c  equation B5
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c        t--temperature (K)
c      rho--molar density (mol/L)
c  output (as function value):
c   D2PBWR--d^2P/dD^2 (kPa-L^2/mol^2)
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  06-04-97 EWL, original version
c  04-09-98 EWL, split expression for expsum (Lahey compiler choked on it)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFBWR/ b(n0:nx,32),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),Rbwr(n0:nx),
     &                ptr(n0:nx),rhotr(n0:nx),ttr(n0:nx),gamma(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension an(15)
c
      D2PBWR=0.0d0
      if (rho.lt.1d-10 .or. t.le.0.d0) then
        RETURN             !dP2dD2 zero for rho = 0, avoid divide by 0
      end if
      j=icomp              !pointer to appropriate fluid
c
c  calculate a(n) terms in MBWR (Eq B2)
c
      tinv=1.0d0/t
      tinv2=tinv*tinv
      tinv3=tinv2*tinv
      tinv4=tinv3*tinv
c
      an(1)=Rbwr(j)*t
      an(2)=b(j,1)*t+b(j,2)*SQRT(t)+b(j,3)+b(j,4)*tinv+b(j,5)*tinv2
      an(3)=b(j,6)*t+b(j,7)+b(j,8)*tinv+b(j,9)*tinv2
      an(4)=b(j,10)*t+b(j,11)+b(j,12)*tinv
      an(5)=b(j,13)
      an(6)=b(j,14)*tinv+b(j,15)*tinv2
      an(7)=b(j,16)*tinv
      an(8)=b(j,17)*tinv+b(j,18)*tinv2
      an(9)=b(j,19)*tinv2
      an(10)=b(j,20)*tinv2+b(j,21)*tinv3
      an(11)=b(j,22)*tinv2+b(j,23)*tinv4
      an(12)=b(j,24)*tinv2+b(j,25)*tinv3
      an(13)=b(j,26)*tinv2+b(j,27)*tinv4
      an(14)=b(j,28)*tinv2+b(j,29)*tinv3
      an(15)=b(j,30)*tinv2+b(j,31)*tinv3+b(j,32)*tinv4
c
c  summation of terms 1-9 in Eq B1
      psum=0.0d0
      rhon=1.0d0/rho**2
      do n=1,9
        rhon=rhon*rho
        psum=psum+real(n)*real(n-1)*an(n)*rhon
      enddo
c  summation of terms 10-15 (exponential terms) in Eq B1
      expsum=0.0d0
      rhon=1.0d0/rho
      rho2=rho*rho
      delsq=(rho/gamma(j))**2
      do n=10,15
        rhon=rhon*rho2
        sum=-35.0d0*REAL(n)+2.0d0*REAL(n)**2
        sum=sum+153.0d0+33.0d0*delsq+2.0d0*delsq**2-4.0d0*REAL(n)*delsq
        expsum=expsum+an(n)*rhon*2.0d0*sum
c       expsum=expsum+an(n)*rhon*2.0d0*(-35.0d0*REAL(n)+2.0d0*REAL(n)**2
c    &      +153.0d0+33.0d0*delsq+2.0d0*delsq**2-4.0d0*REAL(n)*delsq)
      enddo
c
c  collect terms and convert from L-bar/mol to kPa-L/mol
      D2PBWR=(psum+exp(-delsq)*expsum)*(100.0d0)
c     D2PBWR=(psum+exp(-delsq)*expsum)*(R/Rbwr(j))
c
      RETURN
      end                                               !function D2PBWR
c
c ======================================================================
c
      function DPTBWR (icomp,t,rho)
c
c  compute partial derivative of pressure with respect to temperature
c  for MBWR equation of state
c
c  based on derivations in Younglove & McLinden (1994), JPCRD 23:731-779
c  equations B3 and B4
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c        t--temperature (K)
c      rho--molar density (mol/L)
c  output (as function value):
c   DPTBWR--dPdT (kPa/K)
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-07-94  MM, original version
c  07-20-95  MM, separate polynomial Cp0 routines from core_BWR
c                and restructure coefficient arrays
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFBWR/ b(n0:nx,32),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),Rbwr(n0:nx),
     &                ptr(n0:nx),rhotr(n0:nx),ttr(n0:nx),gamma(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension dandt(15)
c
      j=icomp              !pointer to appropriate fluid
c
c  calculate temperature derivatives of a(n) terms in MBWR (Eq B4)
c
      tinv=1.0d0/t
      tinv2=tinv*tinv
      tinv3=tinv2*tinv
      tinv4=tinv3*tinv
      tinv5=tinv4*tinv
c
      dandt(1)=Rbwr(j)
      dandt(2)=b(j,1)+0.5d0*b(j,2)/SQRT(t)-b(j,4)*tinv2
     &         -2.0d0*b(j,5)*tinv3
      dandt(3)=b(j,6)-b(j,8)*tinv2-2.0d0*b(j,9)*tinv3
      dandt(4)=b(j,10)-b(j,12)*tinv2
      dandt(5)=0.0d0
      dandt(6)=-b(j,14)*tinv2-2.0d0*b(j,15)*tinv3
      dandt(7)=-b(j,16)*tinv2
      dandt(8)=-b(j,17)*tinv2-2.0d0*b(j,18)*tinv3
      dandt(9)=-2.0d0*b(j,19)*tinv3
      dandt(10)=-2.0d0*b(j,20)*tinv3-3.0d0*b(j,21)*tinv4
      dandt(11)=-2.0d0*b(j,22)*tinv3-4.0d0*b(j,23)*tinv5
      dandt(12)=-2.0d0*b(j,24)*tinv3-3.0d0*b(j,25)*tinv4
      dandt(13)=-2.0d0*b(j,26)*tinv3-4.0d0*b(j,27)*tinv5
      dandt(14)=-2.0d0*b(j,28)*tinv3-3.0d0*b(j,29)*tinv4
      dandt(15)=-2.0d0*b(j,30)*tinv3-3.0d0*b(j,31)*tinv4
     &          -4.0d0*b(j,32)*tinv5
c
c  summation of terms 2-9 in Eq B3
      psum=0.0d0
      rhon=1.0d0
      do n=1,9
        rhon=rhon*rho
        psum=psum+dandt(n)*rhon
      enddo
c  summation of terms 10-15 (exponential terms) in Eq B3
      expsum=0.0d0
      rhon=rho
      rho2=rho*rho
      do n=10,15
        rhon=rhon*rho2
        expsum=expsum+dandt(n)*rhon
      enddo
c
c  collect terms and convert from L-bar/mol to kPa-L/mol
      DPTBWR=(psum+exp(-(rho/gamma(j))**2)*expsum)*(100.0d0)
c     DPTBWR=(psum+exp(-(rho/gamma(j))**2)*expsum)*(R/Rbwr(j))
c
      RETURN
      end                                               !function DPTBWR
c
c ======================================================================
c
      function PHIBWR (icomp,itau,idel,tau,del)
c
c  compute reduced Helmholtz energy or a derivative as functions
c  of dimensionless temperature and density for the MBWR equation of
c  state; uses the dimensioned functions ABWR, PBWR, etc. in this file
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
c  11-03-95  MM, original version
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c  06-04-97 EWL, added third derivative with respect to density
c  09-00-00 EWL, remove del**idel*tau**itau in calculation of phibwr
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
c
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
      phibwr=0.0d0
      if (del.le.1.0d-10) then    !trivial solution at zero density
        RETURN                    !for any and all derivatives
      end if
c
      t=tz(icomp)/tau
      rho=rhoz(icomp)*del
c
      if (itau.eq.0 .and. idel.eq.0) then
c  compute reduced Helmholtz energy
        phibwr=ABWR(icomp,t,rho)/(R*t)
      else if (itau.eq.1 .and. idel.eq.1) then
c  compute cross derivative
        phibwr=(PBWR(icomp,t,rho)/t-DPTBWR(icomp,t,rho))/(R*rho)
      else if (itau.eq.1 .and. idel.eq.0) then
c  compute first temperature derivative
        phibwr=(ABWR(icomp,t,rho)/t-DABWR(icomp,t,rho))/R
      else if (itau.eq.2 .and. idel.eq.0) then
c  compute second derivative w.r.t. tau
        phibwr=D2ABWR(icomp,t,rho)*t/R
      else if (idel.eq.1 .and. itau.eq.0) then
c  compute first density derivative
        phibwr=PBWR(icomp,t,rho)/(R*t*rho)-1.0d0
      else if (idel.eq.2 .and. itau.eq.0) then
c  compute second derivative w.r.t. del
        phibwr=(DPDBWR(icomp,t,rho)-2.0d0*PBWR(icomp,t,rho)/rho)/
     &         (R*t)+1.0d0
      else if (idel.eq.3 .and. itau.eq.0) then
c  compute third derivative w.r.t. del
        phibwr=(D2PBWR(icomp,t,rho)*rho-4.0d0*DPDBWR(icomp,t,rho)
     &        +6.0d0*PBWR(icomp,t,rho)/rho)/(R*t)-2.0d0
      else
        phibwr=0
      end if
c
      RETURN
      end                                               !function PHIBWR
c
c ======================================================================
c
      subroutine CRTBWR (icomp,tcrit,pcrit,Dcrit)
c
c  returns critical parameters associated with MBWR EOS
c
c  N.B.  these critical parameters may not necessarily be most
c        accurate values, but they are consistent with MBWR fit
c
c  input:
c    icomp--pointer specifying component (1..nc)
c  outputs:
c    tcrit--critical temperature (K)
c    pcrit--critical pressure (kPa)
c    Dcrit--molar density (mol/L) at critical point
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  07-20-95  MM, original version
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /WCFBWR/ b(n0:nx,32),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),Rbwr(n0:nx),
     &                ptr(n0:nx),rhotr(n0:nx),ttr(n0:nx),gamma(n0:nx),
     &                tmax(n0:nx),pmax(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
      tcrit=tc(icomp)
      pcrit=pc(icomp)*100.0d0    ! convert bar -> kPa
      Dcrit=rhoc(icomp)
c
      RETURN
      end                                             !subroutine CRTBWR
c
c ======================================================================
c
      subroutine SETBWR (nread,icomp,hcasno,ierr,herr)
c
c  set up working arrays for use with MBWR equation of state
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
c     other quantities returned via arrays in common /WCFBWR/
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  07-20-95  MM, original version
c  09-13-95  MM, add ierr, herr to argument list
c  10-03-95  MM, adapt to file input, add nread to argument list
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                add Zcrit to common /CCON/
c  03-19-19  MM, add dipole moment to /CCON/ and /MSCBWR/
c  03-22-96  MM, replace /CPMOD/ with /EOSMOD/
c  06-03-96  MM, add limits to /EOSLIM/, reduce mxbwr from 20 to 2
c  05-27-97  MM, if nc = 1, set R to fluid-specific value
c  02-11-98  MM, store rho at triple point separate from rhomax
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-22-98 EWL, set Reos to Rbwr
c  02-15-00 EWL, change read statement to aid UNIX machines
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (mxbwr=2)         !max number of MBWR EOS in block data
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*3 hcpbwr
      character*3 hpheq,heos,hmxeos,hmodcp
      character*12 hcasno,hcas
      character*255 herr
      common /NCOMP/ nc,ic
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /HCHAR/ htab,hnull
c  commons associated with the mxbwr fluids with MBWR equations stored
c  in block data BDBWR
      common /CASBWR/ hcas(mxbwr)
      common /CPMBWR/ hcpbwr(mxbwr)
      common /CFBWR/ ba(mxbwr,32),
     &               pca(mxbwr),rhoca(mxbwr),tca(mxbwr),
     &               Rbwra(mxbwr),ptra(mxbwr),rhotra(mxbwr),ttra(mxbwr),
     &               gammaa(mxbwr),tmaxa(mxbwr),pmaxa(mxbwr)
      common /MSCBWR/ wmb(mxbwr),ttpb(mxbwr),tnbpb(mxbwr),accenb(mxbwr),
     &                dipm(mxbwr)
c  commons associated with the nc components of current interest
c  ("working" commons and arrays)
      common /WCFBWR/ b(n0:nx,32),
     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),
     &                Rbwr(n0:nx),pmin(n0:nx),rhotp(n0:nx),tmin(n0:nx),
     &                gamma(n0:nx),tmax(n0:nx),pmax(n0:nx)
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c     common /CPMOD/ hmodcp(n0:nx)
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
c  limits associated with the equation of state
      common /EOSLIM/ tmn(n0:nx),tmx(n0:nx),pmx(n0:nx),rhomx(n0:nx)
c
      if (nread.le.0) then
c  get coefficients from block data
c  identify specified fluid with entries in database via match of CAS no
        do i=1,mxbwr
          if (hcasno.eq.hcas(i)) then
c           write (*,*) ' SETBWR--coeff from block data for CAS#',hcasno
            hmodcp(icomp)=hcpbwr(i)   !pointer to Cp0 model
            do j=1,32
              b(icomp,j)=ba(i,j)        !32 coefficients of MBWR EOS
            enddo
            pc(icomp)=pca(i)          !critical parameters, limits, etc.
            rhoc(icomp)=rhoca(i)
            tc(icomp)=tca(i)
            Rbwr(icomp)=Rbwra(i)
            if (nc.eq.1 .and. icomp.eq.1) then
              R=Rbwr(icomp)*100.0d0   !MBWR uses pressure in bar
            end if
            Reos(icomp)=Rbwr(icomp)*100.0d0
            pmin(icomp)=ptra(i)       !MBWR uses pressure in bar
            rhotp(icomp)=rhotra(i)
            rhomx(icomp)=rhotra(i)
            tmin(icomp)=ttra(i)
            gamma(icomp)=gammaa(i)
            tmax(icomp)=tmaxa(i)
            pmax(icomp)=pmaxa(i)      !MBWR uses bar
c  fill arrays in /CCON/
            wm(icomp)=wmb(i)
            ttp(icomp)=ttpb(i)
            tnbp(icomp)=tnbpb(i)
            tcrit(icomp)=tc(icomp)
            pcrit(icomp)=pc(icomp)*100.0d0    !MBWR uses bar
            Dcrit(icomp)=rhoc(icomp)
            Zcrit(icomp)=pc(icomp)/(Rbwr(icomp)*100.0d0  ! MBWR uses bar
     &                  *tc(icomp)*rhoc(icomp))
            accen(icomp)=accenb(i)
            dipole(icomp)=dipm(i)
            ierr=0
            herr=' '
            tz(icomp)=tcrit(icomp)
            rhoz(icomp)=Dcrit(icomp)
            RETURN
          end if
        enddo
        ierr=1
        herr='[SETBWR error] Fluid input to SETBWR not found'//hnull
      else
c  read data from file
c       write (*,*) ' SETBWR--read component',icomp,' from unit',nread
        read (nread,*) tmin(icomp)          !lower temperature limit
        read (nread,*) tmax(icomp)          !upper temperature limit
        read (nread,*) pmax(icomp)          !upper pressure limit
        read (nread,*) rhomx(icomp)         !upper density limit
        read (nread,2003) hmodcp(icomp)     !pointer to Cp0 model
        read (nread,*) wm(icomp)     !molecular weight
        read (nread,*) ttp(icomp)    !triple point temperature
        read (nread,*) pmin(icomp)   !pressure at triple point
        read (nread,*) rhotp(icomp)  !density at triple point
        read (nread,*) tnbp(icomp)   !normal boiling point temperature
        read (nread,*) accen(icomp)  !acentric factor
        read (nread,*) tc(icomp),pcrit(icomp),rhoc(icomp) !critical par
        tcrit(icomp)=tc(icomp)
        pc(icomp)=pcrit(icomp)*0.01d0!MBWR uses pressure in bar
        Dcrit(icomp)=rhoc(icomp)
        tz(icomp)=tcrit(icomp)
        rhoz(icomp)=Dcrit(icomp)
        ptp(icomp)=pmin(icomp)
        dtp(icomp)=rhotp(icomp)
        dtpv(icomp)=0.0d0
        dnbpl(icomp)=0.0d0
        dnbpv(icomp)=0.0d0
        read (nread,*) tred,Dred     !reducing parameters (same as crit)
        read (nread,*) gamma(icomp)  !gamma (usually equal to rhoc)
        read (nread,*) Rbwr(icomp)   !gas constant used in fit
        if (nc.eq.1 .and. icomp.eq.1) then
          R=Rbwr(icomp)*100.0d0      !MBWR uses pressure in bar
c         write (*,*) ' SETBWR--R set to ',R
        end if
        Reos(icomp)=Rbwr(icomp)*100.0d0
        Zcrit(icomp)=pcrit(icomp)/(Rbwr(icomp)*100.0d0  ! MBWR uses bar
     &              *tc(icomp)*rhoc(icomp))
        read (nread,*) nterm,ncoeff  !always equal to 32, 1
c  The following format is needed in some cases on Unix machines:
        do j=1,28,3
          read (nread,*) b(icomp,j),b(icomp,j+1),b(icomp,j+2)
        enddo
        read (nread,*) b(icomp,31),b(icomp,32)
c       read (nread,*) (b(icomp,j),j=1,32)  !the 32 coefficients
        if (ABS(tred-tc(icomp)).gt.1.0d-4) then
          ierr=-104
          write (herr,1104) tred,tc(icomp),hnull
          call ERRMSG (ierr,herr)
 1104     format ('[SETUP warning 104] error in specification of BWR ',
     &            'model:  reducing temperature not equal to critical ',
     &            'temperature; T_red = ',f10.3,' K; T_crit = ',f10.3,
     &            ' K.',a1)
        end if
        if (ABS(Dred-rhoc(icomp)).gt.1.0d-6) then
          ierr=-104
          write (herr,1114) Dred,rhoc(icomp),hnull
          call ERRMSG (ierr,herr)
 1114     format ('[SETUP warning 104] error in specification of BWR ',
     &            'model:  reducing density not equal to critical ',
     &            'density; D_red = ',f10.5,' mol/L; D_crit = ',f10.5,
     &            ' mol/L.',a1)
        end if
        if (nterm.ne.32 .or. ncoeff.ne.1) then
          ierr=-104
          write (herr,1124) nterm,ncoeff,hnull
          call ERRMSG (ierr,herr)
 1124     format ('[SETUP warning 104] error in specification of BWR ',
     &            'model:  Nterm = ',i3,' (must be 32); Ncoeff = ',i3,
     &            '(must be 1).',a1)
        end if
c       write (*,*) ' SETBWR--final coefficient: ',b(icomp,32)
        ierr=0
        herr=' '
      end if
c
c  copy limits into /EOSLIM/ arrays
      tmn(icomp)=tmin(icomp)
      tmx(icomp)=tmax(icomp)
      pmx(icomp)=pmax(icomp)
c     rhomx(icomp)=rhomax(icomp)
c
      RETURN
 2003 format (a3)
      end                                             !subroutine SETBWR
cc
cc ======================================================================
cc
c      block data BDBWR
cc
cc  data for MBWR equations of state
cc
c      implicit double precision (a-h,o-z)
c      parameter (mxbwr=2)
c      character*3 hcpbwr
c      character*12 hcas
c      common /CASBWR/ hcas(mxbwr)
c      common /CFBWR/ b(mxbwr,32),
c     &               pc(mxbwr),rhoc(mxbwr),tc(mxbwr),
c     &               Rbwr(mxbwr),pmin(mxbwr),rhomin(mxbwr),tmin(mxbwr),
c     &               gamma(mxbwr),tmax(mxbwr),pmax(mxbwr)
c      common /MSCBWR/ wmb(mxbwr),ttpb(mxbwr),tnbpb(mxbwr),accenb(mxbwr),
c     &                dipm(mxbwr)
c      common /CPMBWR/ hcpbwr(mxbwr)
cc
cc  explanation of parameter
cc     mxbwr:        maximum number of MBWR fits, used to dimension arrays
cc
cc  explanation of commons and constituent arrays
cc    /CASBWR/    Chem Abstract number; used as unambiguous identifier
cc      hcas(i):  CAS number for fluid corresponding to equation "i"
cc
cc    /CFBWR/   parameters to MBWR fits for each of mxbwr fluids
cc      b(i,1..32):  32 coefficients to MBWR equation of state
cc                   units are L, mol, bar, K
cc      pc(i):       critical pressure (bar)
cc      rhoc(i):     critical density (L/mol)
cc      tc(i):       critical temperature (K)
cc      Rbwr(i):     gas constant used in MBWR fit (L-bar/(mol-K))
cc      pmin(i):     pressure at tmin(i) (bar)
cc      rhotp(i):    density at tmin(i), e.g. triple point (L/mol)
cc      tmin(i):     low temperature limit of MBWR--often triple point
cc                   pmin(i),rhotp(i) are used for initial guesses, etc.
cc                   and are often approximate values only
cc      gamma(i):    term in exponential of MBWR--usually same as rhoc(i)
cc      tmax(i):     upper temperature limit of MBWR (K)
cc      pmax(i):     upper pressure limit of MBWR (bar)
cc
cc    /MSCBWR/  miscellaneous fluid constants
cc      wmb(i):      molecular mass (g/mol)
cc      ttpb(i):     triple point temperature (K)
cc      tnbpb(i):    normal boiling point temperature (K)
cc      accenb(i):   acentric factor for fluid represented by eqn "i"
cc      dipm(i):     dipole moment [debye] (at Tnbp if t-dependent)
cc
cc    /CPMBWR/
cc      hmodcp(i)    pointer to Cp0 model to use with fluid "i"
cc
cc      where "i" is the equation number
cc
cc
cc  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
cc  11-20-94  MM, original version
cc  07-21-95  MM, remove Cp0 data and add pointer to Cp0 model
cc  03-19-96  MM, add dipole moment to /MSCBWR/
cc
cc
cc  R134a  1,1,1,2-tetrafluoroethane
c      data hcas(1) /'811-97-2'/
cc  use polynomial Cp0 model
c      data hcpbwr(1) /'CPP'/
cc  Huber fit 48m (4-11-92)
cc  pub as Huber & McLinden (1992), Int Refrig Conf, Purdue, 453-462
cc
c      data (b(1,i),i=1,32)/
c     &   0.965209362217d-01,  -0.401824768889d+01,   0.395239532858d+02,
c     &   0.134532868960d+04,  -0.139439741347d+07,  -0.309281355175d-02,
c     &   0.292381512283d+01,  -0.165146613555d+04,   0.150706003118d+07,
c     &   0.534973948313d-04,   0.543933317622d+00,  -0.211326049762d+03,
c     &  -0.268191203847d-01,  -0.541067125950d+00,  -0.851731779398d+03,
c     &   0.205188253646d+00,  -0.733050188093d-02,   0.380655963862d+01,
c     &  -0.105832087589d+00,  -0.679243084424d+06,  -0.126998378601d+09,
c     &  -0.426234431829d+05,   0.101973338234d+10,  -0.186699526782d+03,
c     &  -0.933426323419d+05,  -0.571735208963d+01,  -0.176762738787d+06,
c     &  -0.397282752308d-01,   0.143016844796d+02,   0.803085294260d-04,
c     &  -0.171959073552d+00,   0.226238385661d+01/
cc
c      data pc(1),rhoc(1),tc(1),
c     &     Rbwr(1),pmin(1),rhomin(1),tmin(1),
c     &     gamma(1),tmax(1),pmax(1)/
c     &  40.56d0,5.0308d0,374.179d0,
c     &  0.08314471d0,0.003935d0,15.609d0,169.853d0,
c     &  5.0308d0,600.0d0,400.0d0/
cc
c      data wmb(1) /102.031d0/
c      data ttpb(1) /169.853d0/
c      data tnbpb(1) /247.082d0/
c      data accenb(1) /0.32705d0/
c      data dipm(1) /2.058d0/  !dipole moment [Debye]; Meyer, (1991)
cc
cc
cc  R123  2,2-dichloro-1,1,1-trifluoroethane
c      data hcas(2) /'306-83-2'/
cc  use polynomial Cp0 model
c      data hcpbwr(2) /'CPP'/
cc  fit of Younglove & McLinden (1994), JPCRD 23:731-779
c      data (b(2,i),i=1,32)/
c     &  -0.657453133659d-02,   0.293479845842d+01,  -0.989140469845d+02,
c     &   0.201029776013d+05,  -0.383566527886d+07,   0.227587641969d-02,
c     &  -0.908726819450d+01,   0.434181417995d+04,   0.354116464954d+07,
c     &  -0.635394849670d-03,   0.320786715274d+01,  -0.131276484299d+04,
c     &  -0.116360713718d+00,  -0.113354409016d+02,  -0.537543457327d+04,
c     &   0.258112416120d+01,  -0.106148632128d+00,   0.500026133667d+02,
c     &  -0.204326706346d+01,  -0.249438345685d+07,  -0.463962781113d+09,
c     &  -0.284903429588d+06,   0.974392239902d+10,  -0.637314379308d+04,
c     &   0.314121189813d+06,  -0.145747968225d+03,  -0.843830261449d+07,
c     &  -0.241138441593d+01,   0.108508031257d+04,  -0.106653193965d-01,
c     &  -0.121343571084d+02,  -0.257510383240d+03/
cc
c      data pc(2),rhoc(2),tc(2),
c     &     Rbwr(2),pmin(2),rhomin(2),tmin(2),
c     &     gamma(2),tmax(2),pmax(2)/
c     &  36.618d0, 3.596417d0,456.831d0,
c     &  0.08314510d0,0.00004d0,11.6d0,166.0d0,
c     &  3.596417d0,600.0d0,400.0d0/
cc
c      data wmb(2) /152.931d0/
c      data ttpb(2) /166.0d0/
c      data tnbpb(2) /300.973d0/
c      data accenb(2) /0.28192d0/
c      data dipm(2) /1.356d0/  !dipole moment [Debye]; Meyer, (1991)
cc
c      end                                              !block data BDBWR
c
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c ======================================================================
c                                                    end file core_BWR.f
c ======================================================================
