c  begin file core_CPP.f
c
c  This file contains the functions implementing a polynomial fit
c  for the ideal gas heat capacity.
c
c  routines for the polynomial fit of ideal gas heat capacity
c
c  contained here are:
c     subroutine SETCPP (nread,icomp,hcasno,ierr,herr)
c     function CP0CPP (icomp,t)
c     function CPICPP (icomp,t)
c     function CPTCPP (icomp,t)
c     function PH0CPP (icomp,itau,idel,t,rho)
c     block data SAVCPP
c     block data BDCPP
c
c  these routines use the following common blocks from other files
c     common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
c     common /HCHAR/ htab,hnull
c     common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
c  various arrays are dimensioned with parameter statements
c     parameter (mxcpp=2)         !max number of fluids in block data
c     parameter (ncmax=20)        !max number of components in mixture
c     parameter (nrefmx=10)       !max number of fluids for transport ECS
c     parameter (n0=-ncmax-nrefmx,nx=ncmax)
c     parameter (ncppmx=20)       !max number of Cp0 terms
c
c ======================================================================
c ======================================================================
c
      subroutine SETCPP (nread,icomp,hcasno,ierr,herr)
c
c  set up working arrays for polynomial form of ideal heat gas Cp
c  uses cubic polynomial correlation of Cp0 (dimensional form):
c     Cp0 = cpc(i,0) + cpc(i,1)*T + cpc(i,2)*T*T + cpc(i,3)*T**3
c
c  inputs:
c    nread--file to read data from
c           <= 0 get data from block data
c           >0 read from logical unit nread (file should have already
c              been opened and pointer set by subroutine SETUP)
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c           zero and negative numbers designate ECS reference fluids
c   hcasno--CAS number of component icomp (not req'd if reading from file)
c
c  outputs:
c     ierr--error flag:  0 = successful
c                        1 = error (e.g. fluid not found)
c     herr--error string (character*255 variable if ierr<>0)
c     coefficients, etc. returned via arrays in commons /xxxCPP/
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  07-20-95  MM, original version
c  09-13-95  MM, add ierr, herr to argument list
c  10-04-95  MM, adapt to file input, add nread to argument list
c  11-10-95  MM, convert from cubic to general polynomial
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c  06-17-96  MM, read in exponential terms (Einstein function)
c  11-13-97  MM, (re)initialize contents of /CPPSAV/ when a new fluid is read in
c  12-05-97  MM, change ncppmx from 15 to 20 (accommodate H2 in NIST14)
c  04-02-98 EWL, read in number of cosh and sinh terms
c  06-29-98 MEV, initialize ncosh and nsinh for version previous to 6.1
c  08-13-98  MM, delete obsolete (unused) format statement
c  03-24-10 EWL, set variables in common blocks to 0
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (mxcpp=2)         !max number of fluids in block data
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (ncppmx=20)       !max number of Cp0 terms
      character*1 htab,hnull
      character*12 hcasno,hcas
      character*255 herr
      common /HCHAR/ htab,hnull
c  commons associated with the mxcpp fluids with FEQ equations stored
c  in block data BDCPP
      common /NTCPP/ ntrmca(mxcpp),ntrmce(mxcpp)
      common /RDCPP/ treda(mxcpp),Creda(mxcpp)
      common /LMCPP/ tmna(mxcpp),tmxa(mxcpp),pmxa(mxcpp),rhomxa(mxcpp)
      common /CPCPP/ cpca(mxcpp,ncppmx),xka(mxcpp,ncppmx)
      common /CASCPP/ hcas(mxcpp)
c  commons associated with the nc components of current interest
c  ("working" commons and arrays)
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),ncosh(n0:nx),
     &                nsinh(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WRDCPP/ tred(n0:nx),Cred(n0:nx)
      common /WLMCPP/ tmin(n0:nx),tmax(n0:nx),pmax(n0:nx),rhomax(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
      common /CPPSAV/ cp0sav(n0:nx),cpisav(n0:nx),cptsav(n0:nx),
     &                tsav(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /VERS/ verfl(n0:nx),vermx    !fluid & mix file version nos.
c
c  (re)initialize contents of /CPPSAV/ when a new fluid is read in
      do i=n0,nx
        cp0sav(i)=0.d0
        cpisav(i)=0.d0
        cptsav(i)=0.d0
        tsav(i)=0.d0
      enddo
      do i=1,ncppmx
        cpc(icomp,i)=0.d0
        xk(icomp,i)=0.d0
        cph(icomp,i)=0.d0
        xth(icomp,i)=0.d0
        xh(icomp,i)=0.d0
      enddo
c
      if (nread.le.0) then
c  get coefficients from block data
c  identify specified fluid with entries in database via match of CAS no
        do k=1,mxcpp
          if (hcasno.eq.hcas(k)) then
c         write (*,*) ' SETCPP--coeff from block data for CAS # ',hcasno
            tmin(icomp)=tmna(k)
            tmax(icomp)=tmxa(k)
            pmax(icomp)=pmxa(k)
            rhomax(icomp)=rhomxa(k)
            tred(icomp)=treda(k)
            Cred(icomp)=Creda(k)
            ntermc(icomp)=ntrmca(k)
            nterme(icomp)=ntrmce(k)
            ncosh(icomp)=0
            nsinh(icomp)=0
            nsp1(icomp)=0         !spares for future use
            nsp2(icomp)=0
            nsp3(icomp)=0
            do j=1,ntrmca(k)
              cpc(icomp,j)=cpca(k,j)
              xk(icomp,j)=xka(k,j)
            enddo
c         write (*,*) ' SETCPP--final coeff: ',cpc(icomp,ntermc(icomp))
            ierr=0
            herr=' '
            RETURN
          end if
        enddo
        ierr=1
        herr='[SETCPP error] Fluid input not found'//hnull
        RETURN
      else
c  read data from file
c       write (*,*) ' SETCPP--read component',icomp,' from unit',nread
        read (nread,*) tmin(icomp)              !lower temperature limit
        read (nread,*) tmax(icomp)              !upper temperature limit
        read (nread,*) pmax(icomp)              !upper pressure limit
        read (nread,*) rhomax(icomp)            !upper density limit
        read (nread,*) tred(icomp),Cred(icomp)  !reducing parameters
c  read number of polynomial and exponential terms
        if (verfl(icomp).ge.6.099d0) then
c  hyperbolic terms (and spares for future) added for file version 6.100
          read (nread,*) ntermc(icomp),nterme(icomp),
     &                   ncosh(icomp),nsinh(icomp),
     &                   nsp1(icomp),nsp2(icomp),nsp3(icomp)  !spares
        else
          read (nread,*) ntermc(icomp),nterme(icomp)
          ncosh(icomp)=0  !these terms not used in files prior to v6.100
          nsinh(icomp)=0
          nsp1(icomp)=0
          nsp2(icomp)=0
          nsp3(icomp)=0
        endif
        jterm=0
        if (ntermc(icomp).ge.1) then
          do j=1,ntermc(icomp)        !read polynomial coefficients
            jterm=jterm+1
            read (nread,*) cpc(icomp,jterm),xk(icomp,jterm)
          enddo
        end if
        if (nterme(icomp).ge.1) then
          do j=1,nterme(icomp)        !read exponential coefficients
            jterm=jterm+1
            read (nread,*) cpc(icomp,jterm),xk(icomp,jterm)
          enddo
        end if
        if (ncosh(icomp).ge.1) then
          do j=1,ncosh(icomp)     !read cosh coefficients
            jterm=jterm+1
            read (nread,*) cpc(icomp,jterm),xk(icomp,jterm),
     &                     cph(icomp,jterm),xth(icomp,jterm),
     &                                      xh(icomp,jterm)
          enddo
        end if
        if (nsinh(icomp).ge.1) then
          do j=1,nsinh(icomp)     !read sinh coefficients
            jterm=jterm+1
            read (nread,*) cpc(icomp,jterm),xk(icomp,jterm),
     &                     cph(icomp,jterm),xth(icomp,jterm),
     &                                      xh(icomp,jterm)
          enddo
        end if
c       write (*,*) ' SETCPP--final coeff: ',cpc(icomp,jterm)
        ierr=0
        herr=' '
      end if
c
      RETURN
      end                                             !subroutine SETCPP
c
c ======================================================================
c
      function CP0CPP (icomp,t)
c
c  compute Cp0
c
c  uses polynomial correlation of Cp0:
c     Cp0/Cred = SUM [cpc(i,k)*(T/Tred(i))**xk(i,k)]
c              + SUM [cpc(i,k)*u(i,k)**2*exp{u(i,k)}/(1-exp{u(i,k)})**2]
c     where cpc(i,k), xk(i,k) are k-th coefficients for component i
c           u(i,k) = xk(i,k)/T
c           Cred is reducing parameter for Cp0 (e.g. gas constant, R)
c           Tred is reducing parameter for t (e.g. critical temperature)
c           (the reducing parameters are specified in the .fld file)
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c        t--temperature (K)
c  output (as function value):
c   CP0CPP--Cp0 (J/(mol-K))
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-06-94  MM, original version
c  07-20-95  MM, separate polynomial Cp0 routines from core_BWR
c  11-10-95  MM, convert from cubic to general polynomial
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c  06-17-96  MM, add provision for exponential terms
c  06-18-96  MM, actually add the exponential terms
c  04-02-98 EWL, add hyperbolic terms
c  06-29-98 MEV, break up expression for the exponential terms to avoid
c                a possible overflow error.
c  01-24-00 EWL, check for ui>1000 to avoid exp(ui)=overflow
c  05-02-01 EWL, change u1>1000 to ui>200
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (ncppmx=20)       !max number of Cp0 terms
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),ncosh(n0:nx),
     &                nsinh(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WRDCPP/ tred(n0:nx),Cred(n0:nx)
      common /WLMCPP/ tmin(n0:nx),tmax(n0:nx),pmax(n0:nx),rhomax(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
      CP0CPP=0.d0
      if (t.le.0.d0) return
      j=icomp                          !pointer to appropriate fluid
c
      cpsum=0.0d0
      iterm=0
c  polynomial terms
      if (ntermc(j).ge.1) then
        do i=1,ntermc(j)
          iterm=iterm+1
          cpsum=cpsum+cpc(j,iterm)*(t/tred(j))**xk(j,iterm)
        enddo
      end if
c  exponential terms
      if (nterme(j).ge.1) then
        do i=1,nterme(j)
          iterm=iterm+1
          ui=xk(j,iterm)/t
          if (ui.lt.200.0d0) then
            expui=EXP(ui)
            cp0cp1=expui/(1.0d0-expui)
            cp0cp1=cp0cp1/(1.0d0-expui)
            cpsum=cpsum+cpc(j,iterm)*ui*ui*cp0cp1
          endif
c       cpsum=cpsum+cpc(j,iterm)*ui*ui*expui/(1.0d0-expui)**2
        enddo
      end if
c  cosh terms
      if (ncosh(j).ge.1) then
        do i=1,ncosh(j)
          iterm=iterm+1
          ttau=cph(j,iterm)*(t/tred(j))**xth(j,iterm)
          if (ttau.lt.700.d0) then
            cpsum=cpsum+cpc(j,iterm)*(t/tred(j))**xk(j,iterm)*
     &          COSH(ttau)**xh(j,iterm)
          else
            cpsum=1.d6
          endif
        enddo
      end if
c  sinh terms
      if (nsinh(j).ge.1) then
        do i=1,nsinh(j)
          iterm=iterm+1
          ttau=cph(j,iterm)*(t/tred(j))**xth(j,iterm)
          if (ttau.lt.700.d0) then
            cpsum=cpsum+cpc(j,iterm)*(t/tred(j))**xk(j,iterm)*
     &          SINH(ttau)**xh(j,iterm)
          else
            cpsum=1.d6
          endif
        enddo
      end if
c
      CP0CPP=cpsum*Cred(j)
c     write (*,1020) j,t,tred(j),Cred(j),CP0CPP,(cpc(j,i),i=1,iterm)
c1020 format (1x,' CP0CPP--j,t,tred,Cred,CP0,cpc(i): ',i3,4f8.2,8d16.6)
c
      RETURN
      end                                               !function CP0CPP
c
c ======================================================================
c
      function CPICPP (icomp,t)
c
c  compute integral of Cp0 over limits of Tref to T
c  for use in enthalpy calculation
c
c  uses polynomial correlation of Cp0:
c     Cp0/Cred = SUM [cpc(i,k)*(T/Tred(i))**xk(i,k)]
c              + SUM [cpc(i,k)*u(i,k)**2*exp{u(i,k)}/(1-exp{u(i,k)})**2]
c     where cpc(i,k), xk(i,k) are k-th coefficients for component i
c           u(i,k) = xk(i,k)/T
c           Cred is reducing parameter for Cp0 (e.g. gas constant, R)
c           Tred is reducing parameter for t (e.g. critical temperature)
c           (the reducing parameters are specified in the .fld file)
c
c  based on derivations in Younglove & McLinden (1994), JPCRD 23:731-779
c  equation C11
c
c  inputs:
c    icomp--pointer specifying component (1..nc);
c        t--temperature (K)
c  output: (as function value):
c   CPICPP--int (Cp0 dT)|T-Tref (J/mol)
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-06-94  MM, original version
c  07-20-95  MM, separate polynomial Cp0 routines from core_BWR
c  11-10-95  MM, convert from cubic to general polynomial
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c  06-17-96  MM, add provision for exponential terms
c  06-18-96  MM, actually add the exponential terms
c  04-02-98 EWL, add hyperbolic terms
c  11-04-02 EWL, include contribution of hyp. terms from the ref. temperature
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (ncppmx=20)       !max number of Cp0 terms
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),ncosh(n0:nx),
     &                nsinh(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WRDCPP/ tred(n0:nx),Cred(n0:nx)
      common /WLMCPP/ tmin(n0:nx),tmax(n0:nx),pmax(n0:nx),rhomax(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
      j=icomp                          !pointer to appropriate fluid
c
      cpsum=0.0d0
      CPICPP=0.d0
      if (t.le.0) RETURN
      iterm=0
c  polynomial terms
      if (ntermc(j).ge.1) then
        do i=1,ntermc(j)
          iterm=iterm+1
          xki=xk(j,iterm)
          xk1=xki+1.0d0
          if (abs(xk1).lt.1.0d-6) then
c  any term with temperature exponent of -1 handled differently
            cpsum=cpsum+cpc(j,iterm)*tred(j)*LOG(t/tref(j))
          else
            cpsum=cpsum+cpc(j,iterm)*(t**xk1-tref(j)**xk1)
     &           /(xk1*tred(j)**xki)
          end if
        enddo
      end if
c  exponential terms
      if (nterme(j).ge.1) then
        do i=1,nterme(j)
          iterm=iterm+1
          ui=xk(j,iterm)/t
          if (ui.lt.200.0d0) then
            expui=EXP(ui)
            exptr=EXP(xk(j,iterm)/tref(j))
            cpsum=cpsum+cpc(j,iterm)*(-0.5d0*xk(j,iterm))
     &        *((1.0d0+expui)/(1.0d0-expui)-(1.0d0+exptr)/(1.0d0-exptr))
          endif
        enddo
      end if
c  cosh terms, only works for c1/t^2/COSH(c2/t)^2
      if (ncosh(j).ge.1) then
        do i=1,ncosh(j)
          iterm=iterm+1
          cpsum=cpsum-cpc(j,iterm)/cph(j,iterm)*TANH(cph(j,iterm)/t)
     &         +cpc(j,iterm)/cph(j,iterm)*TANH(cph(j,iterm)/tref(j))
        enddo
      end if
c  sinh terms, only works for c1/t^2/SINH(c2/t)^2
      if (nsinh(j).ge.1) then
        do i=1,nsinh(j)
          iterm=iterm+1
          cpsum=cpsum+cpc(j,iterm)/cph(j,iterm)/TANH(cph(j,iterm)/t)
     &             -cpc(j,iterm)/cph(j,iterm)/TANH(cph(j,iterm)/tref(j))
        enddo
      end if
c
      CPICPP=cpsum*Cred(j)
c     write (*,1020) j,t,tred(j),Cred(j),CPICPP,(cpc(j,i),i=1,iterm)
c1020 format (1x,' CPICPP--j,t,tred,Cred,CPI,cpc(i): ',i3,4f8.2,8d16.6)
c
      RETURN
      end                                               !function CPICPP
c
c ======================================================================
c
      function CPTCPP (icomp,t)
c
c  compute integral of Cp0/T over limits of Tref to T
c  for use in entropy calculation
c
c  uses polynomial correlation of Cp0:
c     Cp0/Cred = SUM [cpc(i,k)*(T/Tred(i))**xk(i,k)]
c              + SUM [cpc(i,k)*u(i,k)**2*exp{u(i,k)}/(1-exp{u(i,k)})**2]
c     where cpc(i,k), xk(i,k) are k-th coefficients for component i
c           u(i,k) = xk(i,k)/T
c           Cred is reducing parameter for Cp0 (e.g. gas constant, R)
c           Tred is reducing parameter for t (e.g. critical temperature)
c           (the reducing parameters are specified in the .fld file)
c
c  based on derivations in Younglove & McLinden (1994), JPCRD 23:731-779
c  equation C12
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c        t--temperature (K)
c  output (as function value):
c   CPTCPP--int (Cp0/T dT)|T-Tref (J/(mol-K))
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-06-94  MM, original version
c  07-20-95  MM, separate polynomial Cp0 routines from core_BWR
c  11-10-95  MM, convert from cubic to general polynomial
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c  06-17-96  MM, add provision for exponential terms
c  06-18-96  MM, actually add the exponential terms
c  04-02-98 EWL, add hyperbolic terms
c  11-04-02 EWL, include contribution of hyp. terms from the ref. temperature
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (ncppmx=20)       !max number of Cp0 terms
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),ncosh(n0:nx),
     &                nsinh(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WRDCPP/ tred(n0:nx),Cred(n0:nx)
      common /WLMCPP/ tmin(n0:nx),tmax(n0:nx),pmax(n0:nx),rhomax(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
      j=icomp                          !pointer to appropriate fluid
c
      cpsum=0.0d0
      CPTCPP=0.d0
      if (t.le.0.d0) RETURN
c
      iterm=0
c  polynomial terms
      if (ntermc(j).ge.1) then
        do i=1,ntermc(j)
          iterm=iterm+1
          xki=xk(j,iterm)
          if (abs(xki).lt.1.0d-6) then
c  any term with temperature exponent of 0 handled differently
            cpsum=cpsum+cpc(j,iterm)*LOG(t/tref(j))
          else
            cpsum=cpsum+cpc(j,iterm)*(t**xki-tref(j)**xki)
     &           /(xki*tred(j)**xki)
          end if
        enddo
c  exponential terms
      end if
      if (nterme(j).ge.1) then
        do i=1,nterme(j)
          iterm=iterm+1
          ui=xk(j,iterm)/t
          if (ui.lt.200.0d0) then
            uiref=xk(j,iterm)/tref(j)
            expui=EXP(ui)
            exptr=EXP(uiref)
            cpsum=cpsum+cpc(j,iterm)*(LOG((1.0d0-exptr)/(1.0d0-expui))
     &           +ui*expui/(expui-1.0d0)-uiref*exptr/(exptr-1.0d0))
          endif
        enddo
      end if
c  cosh terms, only works for c1/t^2/COSH(c2/t)^2
      if (ncosh(j).ge.1) then
        do i=1,ncosh(j)
          iterm=iterm+1
          ttau=cph(j,iterm)/t
          if (ttau.lt.700.d0) then
           cpsum=cpsum
     &     -cpc(j,iterm)/cph(j,iterm)/t      *TANH(ttau)
     &     +cpc(j,iterm)/cph(j,iterm)/tref(j)*TANH(cph(j,iterm)/tref(j))
     &     +cpc(j,iterm)/cph(j,iterm)**2*LOG(COSH(ttau))
     &     -cpc(j,iterm)/cph(j,iterm)**2*LOG(COSH(cph(j,iterm)/tref(j)))
          else
            cpsum=1.d6
          endif
        enddo
      end if
c  sinh terms, only works for c1/t^2/SINH(c2/t)^2
      if (nsinh(j).ge.1) then
        do i=1,nsinh(j)
          iterm=iterm+1
          ttau=cph(j,iterm)/t
          if (ttau.lt.700.d0) then
           cpsum=cpsum
     &     +cpc(j,iterm)/cph(j,iterm)/t      /TANH(ttau)
     &     -cpc(j,iterm)/cph(j,iterm)/tref(j)/TANH(cph(j,iterm)/tref(j))
     &     -cpc(j,iterm)/cph(j,iterm)**2*LOG(SINH(ttau))
     &     +cpc(j,iterm)/cph(j,iterm)**2*LOG(SINH(cph(j,iterm)/tref(j)))
          else
            cpsum=1.d6
          endif
        enddo
      endif
c
      CPTCPP=cpsum*Cred(j)
c     write (*,1020) j,t,tred(j),Cred(j),CPTCPP,(cpc(j,i),i=iterm)
c1020 format (1x,' CPTCPP--j,t,tred,Cred,CPT,cpc(i): ',i3,4f8.2,8d16.6)
c
      RETURN
      end                                               !function CPTCPP
c
c ======================================================================
c
      function PH0CPP (icomp,itau,idel,t,rho)
c
c  compute the ideal gas part of the reduced Helmholtz energy or a
c  derivative as functions of temperature and density; for
c  use with a Helmholtz-explicit equation of state
c
c  inputs:
c    icomp--pointer specifying component (1..nc)
c     itau--flag specifying order of temperature derivative to calc
c     idel--flag specifying order of density derivative to calculate
c           when itau = 0 and idel = 0, compute A0/RT
c           when itau = 1 and idel = 0, 1st temperature derivative
c           when itau = 2 and idel = 0, 2nd temperature derivative
c           when itau = 0 and idel = 1, 1st density derivative
c           (actually the derivatives are with respect to the dimensionless
c           quantities tau and del)
c        t--temperature (K)
c      rho--density (mol/L)
c  output (as function value):
c   ph0cpp--ideal-gas part of the Helmholtz energy in reduced form (A/RT);
c           the tau derivatives (as specified by itau) are multiplied by
c           the corresponding power of tau; i.e.
c           when itau = 1, the quantity returned is tau*d(ph0cpp)/d(tau)
c                itau = 2, tau**2*d2(ph0cpp)/d(tau)**2 is returned
c                itau = 3, tau**3*d3(ph0cpp)/d(tau)**3 is returned
c           similarly, the del derivatives (as specified by idel) are
c           multiplied by the corresponding power of del (the derivatives
c           usually appear with this factor and this approach neatly
c           avoids a possible divide by zero)
c
c  N.B.  While the real-gas part of the Helmholtz energy is calculated
c        in terms of dimensionless temperature and density, the ideal-
c        gas part is calculated in terms of absolute temperature and
c        density.  (This distinction is necessary for mixtures.)
c
c        The Helmholtz energy consists of ideal-gas and residual
c        (real-gas) terms; this routine calculates only the ideal part.
c
c        This function computes pure component properties only.
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  08-04-95  MM, original version
c  08-21-95  MM, put saved variables into common (rather than save stmt)
c  10-03-95  MM, change /MODEL/ + /CPMOD/:  models specified by strings
c  11-06-95  MM, insert HMX model
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c  03-13-96  MM, replace calls to REDFEQ, etc with generic REDK
c  03-21-96  MM, delete /MODEL/, not needed because of above change
c  04-18-96  MM, apply tolerance to t-tsav test
c  04-19-96  MM, fix bug in expression for PHI; change input p -> rho
c  05-10-96  MM, eliminate duplication of href, sref terms with THERM
c  07-05-96  MM, change derivative outputs:  tau*d(phi)/d(tau), etc
c  08-20-97  MM, call ERRMSG if itau out of range; drop idel=idel
c  08-12-98  MM, add third temperature derivative
c                N.B. third derivative not yet implemented for hyperbolic terms
c  08-13-98  MM, add cross derivatives (all zero) and 3rd density deriv
c  08-16-07 NRN, add 3rd derivatives for cosh and sinh (Ryan Nannan, Delft University)
c  03-21-10 EWL, add call to RMIX
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (ncppmx=20)       !max number of Cp0 terms
      character*1 htab,hnull
      character*255 herr
      common /HCHAR/ htab,hnull
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /CPPSAV/ cp0sav(n0:nx),cpisav(n0:nx),cptsav(n0:nx),
     &                tsav(n0:nx)
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),ncosh(n0:nx),
     &                nsinh(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WRDCPP/ tred(n0:nx),Cred(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
c
      call INFO (icomp,wmm,ttrp,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas)
      PH0CPP=0.0d0    !initialize in case of error
      if (t.le.0) RETURN
      if (abs(t-tsav(icomp)).lt.1.0d-8) then
c  use values from previous call (values already in cp0sav, etc.)
c       write (*,*) 'PH0CPP--using stored values for itau,t = ',itau,t
      else
c  otherwise, compute new values and save for possible future use
        cp0sav(icomp)=CP0CPP(icomp,t)
        cpisav(icomp)=CPICPP(icomp,t)
        cptsav(icomp)=CPTCPP(icomp,t)
        tsav(icomp)=t
      end if
c
      Rt=Rgas*t
c
      if (itau*idel.ne.0) then
c  all cross derivates are zero
        PH0CPP=0.0d0
      else if (itau.eq.0 .and. idel.eq.0) then
c  compute reduced Helmholtz
        if (rho.gt.0) PH0CPP=cpisav(icomp)/Rt-cptsav(icomp)/Rgas
     &        +LOG(t*rho/(tref(icomp)*rhoref(icomp)))
     &        -1.0d0
c
      else if (itau.eq.1) then
c  compute derivative w.r.t. tau (dimensionless temperature)
        PH0CPP=cpisav(icomp)/Rt-1.0d0      !return tau*d(ph0cpp)/d(tau)
c
      else if (itau.eq.2) then
c  compute 2nd derivative w.r.t. tau (dimensionless temperature)
c  return tau**2*d2(ph0cpp)/d(tau**2)
        PH0CPP=1.0d0-cp0sav(icomp)/Rgas
c
      else if (itau.eq.3) then
        j=icomp                          !pointer to appropriate fluid
        tau=tz(icomp)/t
c  note that t0 is reducing temperature for EOS (conversion between
c  absolute temperature and tau) while tred is a (possible) separate
c  reducing temperature appearing only in Cp0
        PH0sum=2.0d0*(-1.0d0+cp0sav(icomp)/Rgas)
        iterm=0
c       write (*,*) 'PH0CPP--initial PH0sum:        ',iterm,PH0sum
c  polynomial terms
        if (ntermc(j).ge.1) then
          do i=1,ntermc(j)
            iterm=iterm+1
            if (ABS(xk(j,iterm)).gt.1.0d-10) then
c  these terms arise from d(Cp0)/d(tau) and derivative of constant
c  terms is zero
              PH0sum=PH0sum+cpc(j,iterm)*xk(j,iterm)
     &             *(tz(icomp)/tred(j))**xk(j,iterm)*tau**(-xk(j,iterm))
     &             *Cred(j)/R

            end if
c         write (*,*) 'PH0CPP--i-polynomial,PH0sum:   ',i,PH0sum
          enddo
        end if
c  exponential terms
        if (nterme(j).ge.1) then
          do i=1,nterme(j)
            iterm=iterm+1
            ui=xk(j,iterm)/t
            if (ui.lt.200.0d0) then
              expui=EXP(ui)
              PH0sum=PH0sum-tau*cpc(j,iterm)*xk(j,iterm)*expui
     &              /(tz(icomp)*(1.0d0-expui)**2)
     &              *(2.0d0*ui+ui*ui*(1.0d0+2.0d0*expui/(1.0d0-expui)))
c           write (*,*) 'PH0CPP--i-exponential,PH0sum:  ',i,PH0sum
            endif
          enddo
        end if
c  cosh terms
        sum_term = 0.d0
        if (ncosh(j).ge.1) then
          do i=1,ncosh(j)
            iterm=iterm+1
            tt   = t/tred(j)
            ttau = cph(j,iterm) * tt**xth(j,iterm)
            term = cpc(j,iterm) * tt**xk(j,iterm) *
     &             cosh(ttau)**xh(j,iterm) * (xk(j,iterm) + xh(j,iterm)*
     &             xth(j,iterm) * cph(j,iterm) * tt**xth(j,iterm) *
     &             tanh(ttau))
            sum_term = sum_term + term
          enddo
        end if
c  sinh terms
        if (nsinh(j).ge.1) then
          do i=1,nsinh(j)
            iterm=iterm+1
            tt   = t/tred(j)
            ttau = cph(j,iterm) * tt**xth(j,iterm)
            term = cpc(j,iterm) * tt**xk(j,iterm) *
     &             sinh(ttau)**xh(j,iterm) * (xk(j,iterm) + xh(j,iterm)*
     &             xth(j,iterm) * cph(j,iterm) * tt**xth(j,iterm) /
     &             tanh(ttau))
            sum_term = sum_term + term
          enddo
        end if
        PH0CPP=PH0sum+sum_term
c
      else if (idel.eq.1) then
c  the three density derivatives are 1/del, -1/del**2, and 1/del**3
c  the quantity returned is derivative multiplied by corresponding power
c  of del, avoiding a possible divide by zero
        PH0CPP=1.0d0
      else if (idel.eq.2) then
        PH0CPP=-1.0d0
      else if (idel.eq.3) then
        PH0CPP=2.0d0
c
      else
c  invalid value of itau
        ierr=99
        write (herr,1099) itau,idel,hnull
 1099   format ('[PH0CPP warning] invalid input; itau =',i4,'; idel =',
     &          i4,a1)
        call ERRMSG (ierr,herr)
        PH0CPP=0.0d0
      end if
c
c     write (*,*) ' PH0CPP:  output phi: ',ph0cpp
c
      RETURN
      end                                               !function PH0CPP
cc
cc ======================================================================
cc
c      block data BDCPP
cc
cc  data for polynomial form of ideal gas heat capacity
cc     Cp0/Cred = SUM [cpc(i,k)*(T/Tred(i))**xk(i,k)]
cc     where Cred is reducing parameter for Cp0 (e.g. gas constant, R)
cc           Tred is reducing parameter for t (e.g. critical temperature)
cc
c      implicit double precision (a-h,o-z)
c      parameter (mxcpp=2)        !max number of fluids in block data
c      parameter (ncppmx=20)      !max number of Cp0 terms
c      character*12 hcas
c      common /NTCPP/ ntrmca(mxcpp),ntrmce(mxcpp)
c      common /RDCPP/ treda(mxcpp),Creda(mxcpp)
c      common /LMCPP/ tmna(mxcpp),tmxa(mxcpp),pmxa(mxcpp),rhomxa(mxcpp)
c      common /CPCPP/ cpca(mxcpp,ncppmx),xka(mxcpp,ncppmx)
c      common /CASCPP/ hcas(mxcpp)
cc
cc  explanation of commons and constituent arrays
cc    /NTCPP/  number of terms in polynomial fit
cc    /RDCPP/  reducing parameters for temperature, Cp0 (often =1)
cc    /LMCPP/  limits of fit: tmin, tmax, pmax (= 0), rhomax (= 0)
cc    /CPCPP/  parameters to polynomial Cp0 fit
cc    /CASCPP/    Chem Abstract number; used as unambiguous identifier
cc      hcas(i):  CAS number for fluid corresponding to equation "i"
cc
cc
cc  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
cc  07-20-95  MM, original version
cc  11-10-95  MM, convert from cubic to general polynomial
cc  06-17-96  MM, add provision for exponential terms
cc
cc      N.B.--the "i" have no particular significance
cc
cc
cc   R134a  1,1,1,2-tetrafluoroethane
c      data hcas(1) /'811-97-2'/
cc  fit of McLinden et al. (1989), ASHRAE Trans 95(pt 2):263-283
cc  also in Huber & McLinden (1992), Int Refrig Conf, Purdue, 453-462
c      data ntrmca(1),ntrmce(1) /3,0/
c      data treda(1),Creda(1) /1.0d0,1.0d0/
c      data tmna(1),tmxa(1),pmxa(1),rhomxa(1)
c     &  /150.0d0,500.0d0,0.0d0,0.0d0/
c      data (xka(1,i),i=1,3)
c     &  / 0.0d0,   1.0d0,      2.0d0/
c      data (cpca(1,i),i=1,3)
c     &  /19.4006d0,0.258531d0,-1.29665d-4/
cc
cc  R123  2,2-dichloro-1,1,1-trifluoroethane
c      data hcas(2) /'306-83-2'/
cc  fit of Younglove & McLinden (1994), JPCRD 23:731-779
c      data ntrmca(2),ntrmce(2) /4,0/
c      data treda(2),Creda(2) /456.831d0,8.31451d0/
c      data tmna(2),tmxa(2),pmxa(2),rhomxa(2)
c     &  /150.0d0,500.0d0,0.0d0,0.0d0/
c      data (xka(2,i),i=1,4)
c     &  /0.0d0,      1.0d0,       2.0d0,     3.0d0/
c      data (cpca(2,i),i=1,4)
c     &  /2.046009d0,22.231991d0,-11.658491d0,2.691665/
cc
c      end                                              !block data BDCPP
c
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c ======================================================================
c                                                    end file core_CPP.f
c ======================================================================
