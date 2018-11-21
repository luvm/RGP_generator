c  begin file prop_sub.f
c
c  This file contains the basic (non-iterative) routines to calculate
c  various properties of fluids and mixtures.  These routines must first
c  be initialized by a call to the subroutine SETUP.
c
c  contained here are:
c     subroutine CRITP (x,tcrit,pcrit,Dcrit,ierr,herr)
c     subroutine THERM (t,rho,x,p,e,h,s,cv,cp,w,hjt)
c     subroutine THERM2 (t,rho,x,p,e,h,s,cv,cp,w,Z,hjt,A,G,xkappa,beta,
c    &                   dPdrho,d2PdD2,dPT,drhodT,drhodP,
c    &                   d2PT2,d2PdTD,spare3,spare4)
c     subroutine THERM3 (t,rho,x,
c    &           xkappa,beta,xisenk,xkt,betas,bs,xkkt,thrott,pi,spht)
c     subroutine THERM0 (t,rho,x,p0,e0,h0,s0,cv0,cp0,w0,A0,G0)
c     subroutine RESIDUAL (t,rho,x,pr,er,hr,sr,cvr,cpr,Ar,Gr)
c     subroutine ENTRO (t,rho,x,s)
c     subroutine ENTHAL (t,rho,x,h)
c     subroutine ENERGY (t,rho,x,e)
c     subroutine CVCP (t,rho,x,cv,cp)
c     subroutine CVCPK (icomp,t,rho,cv,cp)
c     subroutine GIBBS (t,rho,x,Ar,Gr)
c     subroutine AG (t,rho,x,a,g)
c     subroutine PRESS (t,rho,x,p)
c     subroutine DPDD (t,rho,x,dpdrho)
c     subroutine DPDDK (icomp,t,rho,dpdrho)
c     subroutine DPDD2 (t,rho,x,d2PdD2)
c     subroutine DPDT (t,rho,x,dpt)
c     subroutine DPDTK (icomp,t,rho,dpt)
c     subroutine DDDP (t,rho,x,drhodp)
c     subroutine DDDT (t,rho,x,drhodt)
c     subroutine DHD1(t,rho,x,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
c     subroutine FGCTY (t,rho,x,f)
c     subroutine FGCTY2 (t,rho,x,f,ierr,herr)
c     subroutine FUGCOF (t,rho,x,phi,ierr,herr)
c     subroutine CHEMPOT (t,rho,x,u,ierr,herr)
c     subroutine ACTVY (t,rho,x,actv,gamma,ierr,herr)
c     subroutine PHIDERV (i,j,t,rho,x,dadn,dnadn,ierr,herr)
c     subroutine VIRB (t,x,b)
c     subroutine DBDT (t,x,dbt)
c     subroutine DBDT2 (t,x,dbt2)
c     subroutine VIRC (t,x,c)
c     subroutine DCDT (t,x,dct)
c     subroutine DCDT2 (t,x,dct2)
c     subroutine VIRD (t,x,d)
c     subroutine VIRBA (t,x,ba)
c     subroutine VIRCA (t,x,ca)
c     subroutine B12 (t,x,b)
c     subroutine EXCESS (t,p,x,kph,rho,vE,eE,hE,sE,aE,gE,ierr,herr)
c     subroutine FPV (t,rho,p,x,f)
c     subroutine RMIX (x)
c     subroutine RMIX2 (x,Rgas)
c     subroutine VIRBCD (icomp,idel,tau,vir)
c     subroutine SPECGR (t,rho,p,gr)
c     subroutine HEAT (t,rho,x,hg,hn,ierr,herr)
c     subroutine ENTHHC (icmb,t1,t2,h)
c     subroutine ISPURE (x,icomp)
c
c  these routines use the following common blocks from other files
c     common /MODEL/ hrefst,heos,hpheq,h2eos(n0:nx),hmixp,htran,hsten
c     common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
c     common /NCOMP/ nc,ic
c     common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
c  various arrays are dimensioned with parameter statements
c     parameter (ncmax=20)        !max number of components in mixture
c     parameter (nrefmx=10)       !max number of fluids for transport ECS
c     parameter (n0=-ncmax-nrefmx,nx=ncmax)
c
c ======================================================================
c ======================================================================
c
      subroutine CRITP (x,tcrit,pcrit,Dcrit,ierr,herr)
c
c  critical parameters as a function of composition
c
c  input:
c        x--composition [array of mol frac]
c  outputs:
c    tcrit--critical temperature [K]
c    pcrit--critical pressure [kPa]
c    Dcrit--critical density [mol/L]
c     ierr--error flag:  0 = successful
c                        1 = did not converge
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  11-20-94  MM, original version
c  07-21-95  MM, call CRTBWR instead of accessing arrays directly
c  08-07-95  MM, add call to Fundamental (Helmholtz) EOS
c  09-13-95  MM, add ierr, herr to argument list
c  10-03-95  MM, change /MODEL/--models specified by strings
c  11-02-95  MM, add call mixture Helmholtz model (HMX)
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                add Zcrit to common /CCON/
c                add call to ECS model
c  03-19-19  MM, add dipole moment to /CCON/
c  03-22-96  MM, replace /MODEL/ with /EOSMOD/
c  10-01-97  MM, add compiler switch to allow access by DLL
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  03-07-07 EWL, add check for tc, pc, dc less than 0
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: CRITP
c     dll_export CRITP
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*3 hpheq,heos,hmxeos,hmodcp
      character*255 herr
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      dimension x(ncmax)
c
      ierr=0
      herr=' '
      if (heos.eq.'FEQ') then
c  pure fluid Fundamental (Helmholtz) EOS
        icomp=1
        call CRTFEQ (icomp,tcrit,pcrit,Dcrit)
      else if (heos.eq.'QUI') then
c  pure fluid quintic equation of state
        icomp=1
        call CRTQUI (icomp,tcrit,pcrit,Dcrit)
      else if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state
        icomp=1
        call CRTBWR (icomp,tcrit,pcrit,Dcrit)
      else if (heos.eq.'ECS') then
c  pure fluid ECS-thermo model
        icomp=1
        call CRTECS (icomp,tcrit,pcrit,Dcrit)
      else if (heos.eq.'HMX') then
c  mixture Helmholtz model
c       write (*,1022) (x(i),i=1,nc)
c1022   format (1x,' CRITP--about to call CRTHMX w/ x = ',5f12.8)
        call CRTHMX (x,tcrit,pcrit,Dcrit,ierr,herr)
      else if (heos.eq.'AGA') then
        call CRTHMX (x,tcrit,pcrit,Dcrit,ierr,herr)
      else if (heos.eq.'PR') then
        call CRTHMX (x,tcrit,pcrit,Dcrit,ierr,herr)
      else
        tcrit=300
        pcrit=1000
        dcrit=10
        ierr=1
        herr='[CRITP error] Specified model not found'//hnull
        call ERRMSG (ierr,herr)
      end if
      if (tcrit.le.0) tcrit=300
      if (pcrit.le.0) pcrit=1000
      if (dcrit.le.0) dcrit=10
c
      RETURN
      end                                              !subroutine CRITP
c
c ======================================================================
c
      subroutine THERM (t,rho,x,p,e,h,s,cv,cp,w,hjt)
c
c  compute thermal quantities as a function of temperature, density,
c  and compositions using core functions (Helmholtz free energy, ideal
c  gas heat capacity and various derivatives and integrals)
c
c  Based on derivations in Younglove & McLinden, JPCRD 23 #5, 1994,
c  Appendix A for pressure-explicit equations (e.g. MBWR) and
c  Baehr & Tillner-Roth, Thermodynamic Properties of Environmentally
c  Acceptable Refrigerants, Berlin:  Springer-Verlag (1995) for
c  Helmholtz-explicit equations (e.g. FEQ).
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c        p--pressure [kPa]
c        e--internal energy [J/mol]
c        h--enthalpy [J/mol]
c        s--entropy [J/mol-K]
c       Cv--isochoric heat capacity [J/mol-K]
c       Cp--isobaric heat capacity [J/mol-K]
c        w--speed of sound [m/s]
c      hjt--isenthalpic Joule-Thompson coefficient [K/kPa]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-11-94  MM, original version
c  08-04-95  MM, add calls to Fundamental (Helmholtz) EOS
c  10-03-95  MM, change /MODEL/--models specified by strings
c  10-10-95  MM, compute ideal gas pressure and pass to PHI0
c  11-03-95  MM, add calls to mixture Helmholtz (HMX) model
c  11-06-95  MM, add calls to mixture ideal gas function
c  11-08-95  MM, convert calls to PHI0, CP0, CPI, CPT to mixture form
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  12-13-95  MM, compute entropy using Cp0, etc rather than PHI0
c  01-18-96  MM, fix s and h ref state for HMX: s = s - sum[x(i)*sref(i)]
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                add Zcrit to common /CCON/
c                replace calls to PHIHMX, PHIFEQ with general PHIX, PHIK
c  03-19-96  MM, add dipole moment to /CCON/
c  03-22-96  MM, replace /MODEL/ with /EOSMOD/
c  04-19-96  MM, change call to PHI0:  pass rho instead of pideal
c                calculate e,h,s using PHI0 rather than Cp0
c  07-05-96  MM, change e, Cv:  PHI0 returns tau*d(phi0)/d(tau), etc.
c  04-22-97  MM, lower bound on rho for s calc set to 1.0d-20
c  10-01-97  MM, add compiler switch to allow access by DLL
c  03-30-98 EWL, add Joule-Thompson coeff to MBWR case
c  03-31-98  MM, special case for Joule-Thompson coeff for rho = 0
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-02-98 EWL, remove compositional dependence for pure fluids
c  12-02-98 EWL, restructure to closely mimic THERM2
c  12-22-98 EWL, recalculate R for mixtures based on values for pure fluids
c  05-25-00 EWL, moved calculation of Z to bottom AFTER p is calculated!
c  09-00-00 EWL, removed the del from del*phi01, etc.  The del's and tau's
c                are now included in the core routines.  Put the reducing
c                variables tz and rhoz directly in the common blocks.
c  09-05-02 EWL, change check on rho from 1.d-10 to 1.0d-8 to avoid
c                zero's coming back from core_feq at low rhos.
c  10-04-06 EWL, change remaining checks on rho from 1.d-10 to 1.0d-8
c  06-18-07 EWL, change all checks for rho in prop_sub.for from 1.d-20 to 1.d-40
c                for ideal gas situations
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: THERM
c     dll_export THERM
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      dimension x(ncmax)
c
      call RMIX (x)
      RT=R*t
      w2=0.d0
      p=0.d0
      e=0.d0
      h=0.d0
      s=0.d0
      cv=0.d0
      cp=0.d0
      w=0.d0
      hjt=0.d0
      if (t.le.0.d0) return
      if (rho.lt.1.0d-40) then
c  entropy calc will crash if rho = 0
        rhos=1.0d-40
      else
        rhos=rho
      end if
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        p=PBWR(icomp,t,rho)
        a=ABWR(icomp,t,rho)
        dadt=DABWR(icomp,t,rho)
        cpiint=CPI(t,x)
        cptint=CPT(t,x)
        dPdrho=DPDBWR(icomp,t,rho)
        dPT=DPTBWR(icomp,t,rho)
        e=a-t*dadt
c    &    +cpiint-R*(t-tref(icomp))  ! R*tref is const, merge w/ href
     &    +cpiint-RT
     &    -href(icomp)
c
        s=-dadt+R*log(rhoref(icomp)/rhos)+cptint-R*log(t/tref(icomp))
     &    -sref(icomp)
        cv=-t*D2ABWR(icomp,t,rho)+CP0(t,x)-R
        if (rho.gt.1.0d-8) then
          cp=cv+t/rho**2*dPT**2/dPdrho
        else
          cp=CP0(t,x)
        end if
        if (rho.gt.1.0d-8) then
          h=e+p/rho
          hjt=(t/rho*dPT/dPdrho-1.0d0)/rho/cp
        else
          h=e+RT
          call VIRB (t,x,b)
          call DBDT (t,x,dbt)
          hjt=(dbt*t-b)/cp
        end if
c  if any of the factors in speed of sound are negative (e.g. in two-
c  phase region) return 0.0
c       w=SQRT(1.0d3/wm*cp/cv*dPdrho)
        w2=cp/cv*dPdrho
        if (w2.gt.0.0d0) then
          w=SQRT(1.0d3/wm(icomp)*w2)
        else
          w=0.0d0
        end if
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi00=PHIK(icomp,0,0,tau,del)  !real-gas terms
          phi01=PHIK(icomp,0,1,tau,del)
          phi10=PHIK(icomp,1,0,tau,del)
          phi11=PHIK(icomp,1,1,tau,del)
          phi02=PHIK(icomp,0,2,tau,del)
          phi20=PHIK(icomp,2,0,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi00=PHIX(0,0,tau,del,x)   !real-gas terms
          phi01=PHIX(0,1,tau,del,x)
          phi10=PHIX(1,0,tau,del,x)
          phi11=PHIX(1,1,tau,del,x)
          phi02=PHIX(0,2,tau,del,x)
          phi20=PHIX(2,0,tau,del,x)
        end if
c       write (*,1003) t,rho,phi00,phi01,phi10,phi11,phi02,phi20
c1003   format (1x,' THERM--t,rho,PHIs:   ',f8.2,f12.6,6d16.6)
c
        phig00=PHI0(0,0,t,rhos,x)      !ideal-gas terms
        phig10=PHI0(1,0,t,rho,x)
        phig20=PHI0(2,0,t,rho,x)
c       write (*,1005) (x(j),j=1,ncmax)
c1005   format (1x,' THERM--output x(i): ',5f14.8)
c       write (*,1024) phig00,phig10,phig20
c1024   format (1x,' THERM--phig-00/01/02:',20x,3d16.6)
        p=RT*rho*(1.0d0+phi01)
        e=RT*(phig10+phi10)
        if (icomp.ne.0) then
          e=e-href(icomp)
        else
          do i=1,nc
            e=e-x(i)*href(i)
          enddo
        endif
        h=e+RT*(1.0d0+phi01)
        s=R*(phig10+phi10-phig00-phi00)
c       write (*,*) ' THERM--t,rho,x,s,sref:  ',t,rho,x(1),s,sref(1)
        if (icomp.ne.0) then
          s=s-sref(icomp)
        else
          do i=1,nc
            s=s-x(i)*sref(i)
          enddo
        endif
        cv=-R*(phi20+phig20)
c       write (*,*) ' THERM--tau,del,Cv_resid:  ',tau,del,phi20
        cp=cv+R*(1.0d0+phi01-phi11)**2/(1.0d0+2.0d0*phi01+phi02)
        if (cv.gt.0.d0) w2=RT*cp/cv*(1.0d0+2.0d0*phi01+phi02)
c  if any of the factors in speed of sound are negative (e.g. in two-
c  phase region) return 0.0
        if (w2.gt.0.0d0) then
          w=SQRT(w2*1.0d3/WMOL(x))  !convert from molar to mass units
        else
          w=0.0d0
        end if
        if (rho.gt.1.0d-8) then
          hjt=-1.0d0/(cp*rho)*(phi01+phi02+phi11)/
     &        (1.0d0+2.0d0*phi01+phi02)
        else
          call VIRB (t,x,b)
          call DBDT (t,x,dbt)
          hjt=(dbt*t-b)/cp
        end if
      end if
c
      RETURN
      end                                              !subroutine THERM
c
c ======================================================================
c
      subroutine THERM2 (t,rho,x,p,e,h,s,cv,cp,w,Z,hjt,A,G,
     &                   xkappa,beta,dPdrho,d2PdD2,dPT,drhodT,drhodP,
     &                   d2PT2,d2PdTD,spare3,spare4)
c
c  compute thermal quantities as a function of temperature, density,
c  and compositions using core functions (Helmholtz free energy, ideal
c  gas heat capacity and various derivatives and integrals)
c
c  this routine is the same as THERM, except that additional properties
c  are calculated
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c        p--pressure [kPa]
c        e--internal energy [J/mol]
c        h--enthalpy [J/mol]
c        s--entropy [J/mol-K]
c       Cv--isochoric heat capacity [J/mol-K]
c       Cp--isobaric heat capacity [J/mol-K]
c        w--speed of sound [m/s]
c        Z--compressibility factor (= PV/RT) [dimensionless]
c      hjt--isenthalpic Joule-Thompson coefficient [K/kPa]
c        A--Helmholtz energy [J/mol]
c        G--Gibbs free energy [J/mol]
c   xkappa--isothermal compressibility (= -1/V dV/dP = 1/rho dD/dP) [1/kPa]
c     beta--volume expansivity (= 1/V dV/dT = -1/rho dD/dT) [1/K]
c   dPdrho--derivative dP/drho [kPa-L/mol]
c   d2PdD2--derivative d^2P/drho^2 [kPa-L^2/mol^2]
c      dPT--derivative dP/dT [kPa/K]
c   drhodT--derivative drho/dT [mol/(L-K)]
c   drhodP--derivative drho/dP [mol/(L-kPa)]
c    d2PT2--derivative d2P/dT2 [kPa/K^2]
c   d2PdTD--derivative d2P/dTd(rho) [J/mol-K]
c   sparei--2 space holders for possible future properties
c
c  written by M. McLinden, NIST Physical & Chem Properties Div, Boulder, CO
c  03-16-98  MM, original version; based on THERM
c  03-30-98 EWL, add Joule-Thompson coeff to MBWR case
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-02-98 EWL, remove compositional dependence for pure fluids
c  12-02-98 EWL, restructure to closely mimic THERM
c  12-02-98 EWL, add the reference states to the calculation of A and G
c  12-22-98 EWL, recalculate R for mixtures based on values for pure fluids
c  02-11-99 EWL, skip calculation of d2PdD2 if rho=0
c  04-27-01 EWL, change order of calculation of a and g in BWR section
c  04-27-01 DGF, change sign before sref in the calculation of a and g
c  09-05-02 EWL, add ideal gas isothermal compressibility and d2pdD2
c  09-05-02 EWL, change check on rho from 1.d-10 to 1.0d-8 to avoid
c                zero's coming back from core_feq at low rhos.
c  06-14-06 EWL, change dPdD to dPdrho, etc, to avoid compiler problems with subroutine dPdD, etc.
c  03-04-08 EWL, add checks for T=0
c  12-04-09 EWL, add d2P/dT2
c  03-26-10 EWL, add d2P/dTdrho
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: THERM2
c     dll_export THERM2
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c  common block containing flags to GUI
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      dimension x(ncmax)
c

      call ISPURE (x,icomp)
      d2PT2=xnotc   !flag indicating not calculated
      d2PdTD=xnotc
      spare3=xnotc
      spare4=xnotc
      p=0.d0
      e=0.d0
      h=0.d0
      s=0.d0
      cv=0.d0
      cp=0.d0
      w=0.d0
      hjt=0.d0
      Z=0.d0
      A=0.d0
      G=0.d0
      xkappa=0.d0
      beta=0.d0
      dPdrho=0.d0
      d2PdD2=0.d0
      dPT=0.d0
      drhodT=0.d0
      drhodP=0.d0
      if (t.le.0.d0) return
c
      call RMIX (x)
      RT=R*t
      call DBDT (t,x,dbt)
      call VIRB (t,x,b)
      rhos=rho
c  entropy calc will crash if rho = 0
      if (rho.lt.1.0d-40) rhos=1.0d-20
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        p=PBWR(icomp,t,rho)
        ar=ABWR(icomp,t,rho)
        dadt=DABWR(icomp,t,rho)
        cpiint=CPI(t,x)
        cptint=CPT(t,x)
        dPdrho=DPDBWR(icomp,t,rho)
        dPT=DPTBWR(icomp,t,rho)
        e=ar-t*dadt
c    &    +cpiint-R*(t-tref(icomp))  ! R*tref is const, merge w/ href
     &    +cpiint-RT
     &    -href(icomp)
c  additional properties added to THERM2 not in THERM
        d2PdD2=D2PBWR(icomp,t,rho)
        drhodT=-dPT/dPdrho
        drhodP=1.0d0/dPdrho
c
        s=-dadt+R*log(rhoref(icomp)/rhos)+cptint-R*log(t/tref(icomp))
     &    -sref(icomp)
        cv=-t*D2ABWR(icomp,t,rho)+CP0(t,x)-R
        A=e-t*s
        if (rho.gt.1.0d-8) then
          G=A+p/rho
          beta=-drhodT/rho
          xkappa=drhodP/rho
        else
          G=A+R*t
          beta=1.0d0/t  !if rho = 0, then ideal-gas behavior
          xkappa=xnotc
          if (p.gt.0.d0) xkappa=1.d0/p
        end if
        if (rho.gt.1.0d-8) then
          cp=cv+t/rho**2*dPT**2/dPdrho
        else
          cp=CP0(t,x)
        end if
        if (rho.gt.1.0d-8) then
          h=e+p/rho
          hjt=(t/rho*DPTBWR(icomp,t,rho)/dPdrho-1.0d0)/rho/cp
        else
          h=e+RT
          hjt=(dbt*t-b)/cp
        end if
c  if any of the factors in speed of sound are negative (e.g. in two-
c  phase region) return 0.0
c       w=SQRT(1.0d3/wm*cp/cv*dPdrho)
        w2=cp/cv*dPdrho
        if (w2.gt.0.0d0) then
          w=SQRT(1.0d3/wm(icomp)*w2)
        else
          w=0.0d0
        end if
c
      else
c  call general PHIK or PHIX routines for all other models
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi00=PHIK(icomp,0,0,tau,del)  !real-gas terms
          phi01=PHIK(icomp,0,1,tau,del)
          phi10=PHIK(icomp,1,0,tau,del)
          phi11=PHIK(icomp,1,1,tau,del)
          phi12=PHIK(icomp,1,2,tau,del)
          phi21=PHIK(icomp,2,1,tau,del)
          phi02=PHIK(icomp,0,2,tau,del)
          phi20=PHIK(icomp,2,0,tau,del)
          phi03=PHIK(icomp,0,3,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi00=PHIX(0,0,tau,del,x)   !real-gas terms
          phi01=PHIX(0,1,tau,del,x)
          phi10=PHIX(1,0,tau,del,x)
          phi11=PHIX(1,1,tau,del,x)
          phi12=PHIX(1,2,tau,del,x)
          phi21=PHIX(2,1,tau,del,x)
          phi02=PHIX(0,2,tau,del,x)
          phi20=PHIX(2,0,tau,del,x)
          phi03=PHIX(0,3,tau,del,x)
        end if
c       write (*,1003) t,rho,phi00,phi01,phi10,phi11,phi02,phi20
c1003   format (1x,' THERM--t,rho,PHIs:   ',f8.2,f12.6,6d16.6)
c
        phig00=PHI0(0,0,t,rhos,x)      !ideal-gas terms
        phig10=PHI0(1,0,t,rho,x)
        phig20=PHI0(2,0,t,rho,x)
c       write (*,1005) (x(j),j=1,ncmax)
c1005   format (1x,' THERM--output x(i): ',5f14.8)
c       write (*,1024) phig00,phig10,phig20
c1024   format (1x,' THERM--phig-00/01/02:',20x,3d16.6)
        p=RT*rho*(1.0d0+phi01)
        e=RT*(phig10+phi10)
        if (icomp.ne.0) then
          e=e-href(icomp)
        else
          do i=1,nc
            e=e-x(i)*href(i)
          enddo
        endif
        h=e+RT*(1.0d0+phi01)
        s=R*(phig10+phi10-phig00-phi00)
c       write (*,*) ' THERM--t,rho,x,s,sref:  ',t,rho,x(1),s,sref(1)
        if (icomp.ne.0) then
          s=s-sref(icomp)
        else
          do i=1,nc
            s=s-x(i)*sref(i)
          enddo
        endif
        cv=-R*(phi20+phig20)
c       write (*,*) ' THERM--tau,del,Cv_resid:  ',tau,del,phi20
        cp=cv+R*(1.0d0+phi01-phi11)**2/
     &     (1.0d0+2.0d0*phi01+phi02)
        w2=RT*cp/cv*(1.0d0+2.0d0*phi01+phi02)
c  if any of the factors in speed of sound are negative (e.g. in two-
c  phase region) return 0.0
        if (w2.gt.0.0d0) then
          w=SQRT(w2*1.0d3/WMOL(x))  !convert from molar to mass units
        else
          w=0.0d0
        end if
        if (rho.gt.1.0d-8) then
          hjt=-1.0d0/(cp*rho)*(phi01+phi02+phi11)/
     &        (1.0d0+2.0d0*phi01+phi02)
        else
          hjt=(dbt*t-b)/cp
        end if
c  additional properties added to THERM2 not in THERM
        A=RT*(phi00+phig00)
        G=A+RT*(1.0d0+phi01)
        if (icomp.ne.0) then
          a=a-href(icomp)+sref(icomp)*t
          g=g-href(icomp)+sref(icomp)*t
        else
          do i=1,nc
            a=a-x(i)*(href(i)-sref(i)*t)
            g=g-x(i)*(href(i)-sref(i)*t)
          enddo
        endif
        dPdrho=RT*(1.0d0+2.0d0*phi01+phi02)
        dPT=R*rho*(1.0d0+phi01-phi11)
        d2PT2=R*rho*phi21        !d2P/dT=d3A/dT2dV=1/T*dCv/dV
        d2PdTD=R*(1.d0+2.d0*phi01+phi02-2.d0*phi11-phi12)
        drhodP=1.0d0/(RT*(1.0d0+2.0d0*phi01+phi02))
        drhodT=-rho*(1.0d0+phi01-phi11)/(t*(1.0d0+2.0d0*phi01+phi02))
        if (rho.gt.1.0d-8) then
          d2PdD2=RT/rho*(2.0d0*phi01+4.0d0*phi02+phi03)
          beta=-drhodT/rho
          xkappa=drhodP/rho
        else
          d2PdD2=2.d0*b*R*t
          beta=xnotc
          beta=1.0d0/t  !if rho = 0, then ideal-gas behavior
          xkappa=xnotc
          if (p.gt.0.d0) xkappa=1.d0/p
        end if
      end if
      if (rho.lt.1.0d-40) then
        Z=1.0d0       !if rho = 0, then ideal-gas behavior
      else
        Z=p/(RT*rho)
      end if
c
      RETURN
      end                                             !subroutine THERM2
c
c ======================================================================
c
      subroutine THERM0 (t,rho,x,p0,e0,h0,s0,cv0,cp00,w0,A0,G0)
c
c  compute ideal gas thermal quantities as a function of temperature, density,
c  and compositions using core functions
c
c  this routine is the same as THERM, except it only calculates ideal gas
c  properties (Z=1) at any temperature and density
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c       p0--pressure [kPa]
c       e0--internal energy [J/mol]
c       h0--enthalpy [J/mol]
c       s0--entropy [J/mol-K]
c      Cv0--isochoric heat capacity [J/mol-K]
c     Cp00--isobaric heat capacity [J/mol-K]
c       w0--speed of sound [m/s]
c       A0--Helmholtz energy [J/mol]
c       G0--Gibbs free energy [J/mol]
c
c  11-26-02 EWL, original version; based on THERM
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: THERM0
c     dll_export THERM0
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      p0=0.d0
      e0=0.d0
      h0=0.d0
      s0=0.d0
      cv0=0.d0
      cp00=0.d0
      w0=0.d0
      A0=0.d0
      G0=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      RT=R*t
      rhos=rho
c  entropy calc will crash if rho = 0
      if (rho.lt.1.0d-40) rhos=1.0d-20
      phig00=PHI0(0,0,t,rhos,x)      !ideal-gas terms
      phig10=PHI0(1,0,t,rho,x)
      phig20=PHI0(2,0,t,rho,x)
      p0=RT*rho
      e0=RT*phig10
      s0=R*(phig10-phig00)
      A0=RT*phig00
      call ISPURE (x,icomp)
      if (icomp.ne.0) then
        e0=e0-href(icomp)
        s0=s0-sref(icomp)
        a0=a0-href(icomp)+sref(icomp)*t
      else
        do i=1,nc
          e0=e0-x(i)*href(i)
          s0=s0-x(i)*sref(i)
          a0=a0-x(i)*(href(i)-sref(i)*t)
        enddo
      endif
      cv0=-R*phig20
      cp00=cv0+R
      h0=e0+RT
      G0=A0+RT
      w2=RT*cp00/cv0
      if (w2.gt.0.0d0) then
        w0=SQRT(w2*1.0d3/WMOL(x))
      else
        w0=0.0d0
      end if
c
      RETURN
      end                                             !subroutine THERM0
c
c ======================================================================
c
      subroutine RESIDUAL (t,rho,x,pr,er,hr,sr,cvr,cpr,Ar,Gr)
c
c  compute the residual quantities as a function of temperature, density,
c  and compositions (where the residual is the property minus the ideal
c  gas portion).
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c       pr--residual pressure [kPa]  (p-rho*R*T)
c       er--residual internal energy [J/mol]
c       hr--residual enthalpy [J/mol]
c       sr--residual entropy [J/mol-K]
c      Cvr--residual isochoric heat capacity [J/mol-K]
c      Cpr--residual isobaric heat capacity [J/mol-K]
c       Ar--residual Helmholtz energy [J/mol]
c       Gr--residual Gibbs free energy [J/mol]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  07-07-10 EWL, original version; based on THERM2
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: RESIDUAL
c     dll_export RESIDUAL
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c

      pr=0.d0
      er=0.d0
      hr=0.d0
      sr=0.d0
      cvr=0.d0
      cpr=0.d0
      Ar=0.d0
      Gr=0.d0
      if (t.le.0.d0) return
      call ISPURE (x,icomp)
      call RMIX (x)
      RT=R*t
      if (icomp.ne.0) then
c  pure fluid
        tau=tz(icomp)/t
        del=rho/rhoz(icomp)
        phi00=PHIK(icomp,0,0,tau,del)  !real-gas terms
        phi01=PHIK(icomp,0,1,tau,del)
        phi10=PHIK(icomp,1,0,tau,del)
        phi11=PHIK(icomp,1,1,tau,del)
        phi12=PHIK(icomp,1,2,tau,del)
        phi21=PHIK(icomp,2,1,tau,del)
        phi02=PHIK(icomp,0,2,tau,del)
        phi20=PHIK(icomp,2,0,tau,del)
        phi03=PHIK(icomp,0,3,tau,del)
      else
c  mixture
        call REDX (x,t0,rho0)
        tau=t0/t
        del=rho/rho0
        phi00=PHIX(0,0,tau,del,x)   !real-gas terms
        phi01=PHIX(0,1,tau,del,x)
        phi10=PHIX(1,0,tau,del,x)
        phi11=PHIX(1,1,tau,del,x)
        phi12=PHIX(1,2,tau,del,x)
        phi21=PHIX(2,1,tau,del,x)
        phi02=PHIX(0,2,tau,del,x)
        phi20=PHIX(2,0,tau,del,x)
        phi03=PHIX(0,3,tau,del,x)
      end if
c
      pr=RT*rho*phi01
      er=RT*phi10
      hr=er+RT*phi01
      sr=R*(phi10-phi00)
      cvr=-R*phi20
      cpr=cvr+R*(1.0d0+phi01-phi11)**2/(1.0d0+2.0d0*phi01+phi02)-R
      Ar=RT*phi00
      Gr=Ar+RT*phi01
c
      RETURN
      end                                           !subroutine RESIDUAL
c
c ======================================================================
c
      subroutine ENTRO (t,rho,x,s)
c
c  compute entropy as a function of temperature, density and composition
c  using core functions (temperature derivative of Helmholtz free energy
c  and ideal gas integrals)
c
c  based on derivations in Younglove & McLinden, JPCRD 23 #5, 1994,
c  equations A5, A19 - A26
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  output:
c        s--entropy [J/mol-K]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-05-94  MM, original version
c  10-03-95  MM, change /MODEL/--models specified by strings
c  10-10-95  MM, compute ideal gas pressure and pass to PHI0
c  11-03-95  MM, add calls to mixture Helmholtz (HMX) model
c  11-08-95  MM, convert calls to PHI0, CP0, CPI, CPT to mixture form
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  01-19-96  MM, fix ref state for HMX: s = s - sum[x(i)*sref(i)]
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                replace calls to PHIHMX, PHIFEQ with general PHIX, PHIK
c  03-22-96  MM, replace /MODEL/ with /EOSMOD/
c  04-19-96  MM, change call to PHI0:  pass rho instead of pideal
c                calculate s using PHI0 rather than Cp0
c  07-19-96  MM, change general calls to PHI0 (same as THERM)
c  04-22-97  MM, lower bound on rho for s calc set to 1.0d-20
c  10-01-97  MM, add compiler switch to allow access by DLL
c  12-02-98 EWL, remove compositional dependence for pure fluids
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: ENTRO
c     dll_export ENTRO
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      s=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (rho.lt.1.0d-40) then
c  entropy calc will crash if rho = 0
        rhos=1.0d-40
      else
        rhos=rho
      end if
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        s=-DABWR(icomp,t,rho)+R*log(rhoref(icomp)/rhos)+CPT(t,x)
     &    -R*log(t/tref(icomp))-sref(icomp)
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi00=PHIK(icomp,0,0,tau,del)  !real-gas terms
          phi10=PHIK(icomp,1,0,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi00=PHIX(0,0,tau,del,x)  !real-gas terms
          phi10=PHIX(1,0,tau,del,x)
        end if
c
        phig00=PHI0(0,0,t,rhos,x)      !ideal-gas terms
        phig10=PHI0(1,0,t,rho,x)
        s=R*(phig10+phi10-phig00-phi00)
        if (icomp.ne.0) then
          s=s-sref(icomp)
        else
          do i=1,nc
            s=s-x(i)*sref(i)
          enddo
        endif
      end if
c
      RETURN
      end                                              !subroutine ENTRO
c
c ======================================================================
c
      subroutine ENTHAL (t,rho,x,h)
c
c  compute enthalpy as a function of temperature, density, and
c  composition using core functions (Helmholtz free energy and ideal
c  gas integrals)
c
c  based on derivations in Younglove & McLinden, JPCRD 23 #5, 1994,
c  equations A7, A18, A19
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  output:
c        h--enthalpy [J/mol]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-06-94  MM, original version
c  10-03-95  MM, change /MODEL/--models specified by strings
c  10-10-95  MM, compute ideal gas pressure and pass to PHI0
c  11-03-95  MM, add calls to mixture Helmholtz (HMX) model
c  11-08-95  MM, convert calls to PHI0, CP0, CPI, CPT to mixture form
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  01-19-96  MM, fix ref state for HMX: h = h - sum[x(i)*href(i)]
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                replace calls to PHIHMX, PHIFEQ with general PHIX, PHIK
c  03-22-96  MM, replace /MODEL/ with /EOSMOD/
c  07-19-96  MM, change general calls to PHI0 (same as THERM)
c  10-01-97  MM, add compiler switch to allow access by DLL
c  12-02-98 EWL, remove compositional dependence for pure fluids
c  04-05-00 EWL, check for rho=0 and avoid division by zero
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: ENTHAL
c     dll_export ENTHAL
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      h=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        if (rho.lt.1.d-8) then
          h=ABWR(icomp,t,rho)-t*DABWR(icomp,t,rho)+CPI(t,x)-href(icomp)
        else
          h=ABWR(icomp,t,rho)-t*DABWR(icomp,t,rho)+
     &     PBWR(icomp,t,rho)/rho-R*t+CPI(t,x)-href(icomp)
        endif
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi01=PHIK(icomp,0,1,tau,del)
          phi10=PHIK(icomp,1,0,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi01=PHIX(0,1,tau,del,x)  !real-gas terms
          phi10=PHIX(1,0,tau,del,x)
        end if
c
        RT=R*t
        phig10=PHI0(1,0,t,rho,x)
        e=RT*(phig10+phi10)
        if (icomp.ne.0) then
          e=e-href(icomp)
        else
          do i=1,nc
            e=e-x(i)*href(i)
          enddo
        endif
        h=e+RT*(1.0d0+phi01)
      end if
c
      RETURN
      end                                             !subroutine ENTHAL
c
c ======================================================================
c
      subroutine ENERGY (t,rho,x,e)
c
c  compute energy as a function of temperature, density, and
c  composition using core functions (Helmholtz free energy and ideal
c  gas integrals)
c
c  based on derivations in Younglove & McLinden, JPCRD 23 #5, 1994
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  output:
c        e--energy [J/mol]
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  12-13-00 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: ENERGY
c     dll_export ENERGY
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      e=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        e=ABWR(icomp,t,rho)-t*DABWR(icomp,t,rho)-R*t+CPI(t,x)
     &   -href(icomp)
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi10=PHIK(icomp,1,0,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi10=PHIX(1,0,tau,del,x)
        end if
c
        phig10=PHI0(1,0,t,rho,x)
        e=R*t*(phig10+phi10)
        if (icomp.ne.0) then
          e=e-href(icomp)
        else
          do i=1,nc
            e=e-x(i)*href(i)
          enddo
        endif
      end if
c
      RETURN
      end                                             !subroutine ENERGY
c
c ======================================================================
c
      subroutine CVCP (t,rho,x,cv,cp)
c
c  compute isochoric (constant volume) and isobaric (constant pressure)
c  heat capacity as functions of temperature, density, and composition
c  using core functions
c
c  based on derivations in Younglove & McLinden, JPCRD 23 #5, 1994,
c  equation A15, A16
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c       cv--isochoric heat capacity [J/mol-K]
c       cp--isobaric heat capacity [J/mol-K]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-06-94  MM, original version
c  10-03-95  MM, change /MODEL/--models specified by strings
c  10-10-95  MM, compute ideal gas pressure and pass to PHI0
c  11-03-95  MM, add calls to mixture Helmholtz (HMX) model
c  11-08-95  MM, convert calls to PHI0, CP0, CPI, CPT to mixture form
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                replace calls to PHIHMX, PHIFEQ with general PHIX, PHIK
c  03-22-96  MM, replace /MODEL/ with /EOSMOD/
c  07-19-96  MM, change general calls to PHI0 (same as THERM)
c  10-01-97  MM, add compiler switch to allow access by DLL
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: CVCP
c     dll_export CVCP
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      cv=0.d0
      cp=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        cv=-t*D2ABWR(icomp,t,rho)+CP0(t,x)-R
        if (rho.gt.1.0d-8) then
          cp=cv+t/rho**2*DPTBWR(icomp,t,rho)**2/DPDBWR(icomp,t,rho)
        else
          cp=cv+R
        end if
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi01=PHIK(icomp,0,1,tau,del)  !real-gas terms
          phi02=PHIK(icomp,0,2,tau,del)
          phi11=PHIK(icomp,1,1,tau,del)
          phi20=PHIK(icomp,2,0,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi01=PHIX(0,1,tau,del,x)  !real-gas terms
          phi02=PHIX(0,2,tau,del,x)
          phi11=PHIX(1,1,tau,del,x)
          phi20=PHIX(2,0,tau,del,x)
        end if
c
        phig20=PHI0(2,0,t,rho,x)         !ideal-gas term
        cv=-R*(phi20+phig20)
        cp=cv+R*(1.0d0+phi01-phi11)**2/
     &     (1.0d0+2.0d0*phi01+phi02)
      end if
c
      RETURN
      end                                               !subroutine CVCP
c
c ======================================================================
c
      subroutine CVCPK (icomp,t,rho,cv,cp)
c
c  compute isochoric (constant volume) and isobaric (constant pressure)
c  heat capacity as functions of temperature for a given component
c
c  analogous to CVCP, except for component icomp, this is used by transport
c  routines to calculate Cv & Cp for the reference fluid (component zero)
c
c  inputs:
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c        t--temperature [K]
c      rho--molar density [mol/L]
c  outputs:
c       cv--isochoric heat capacity [J/mol-K]
c       cp--isobaric heat capacity [J/mol-K]
c
c  written by M. McLinden, NIST Physical & Chem Properties Div, Boulder, CO
c  06-16-97  MM, original version; based on CVCP
c  10-01-97  MM, add compiler switch to allow access by DLL
c  03-06-98  MM, check hmxeos, not heos, for 'BWR' (crash for icomp = 0)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: CVCPK
c     dll_export CVCPK
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
      cv=0.d0
      cp=0.d0
      if (t.le.0.d0) return
      if (hmxeos(icomp).eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        if (rho.gt.1.0d-8) then
          cv=-t*D2ABWR(icomp,t,rho)+CP0K(icomp,t)-R
          cp=cv+t/rho**2*DPTBWR(icomp,t,rho)**2/DPDBWR(icomp,t,rho)
        else
          cp=CP0K(icomp,t)
          cv=cp-R
        end if
c
      else
c  call general PHIK or PHIX routines for all other models
        tau=tz(icomp)/t
        del=rho/rhoz(icomp)
        phi01=PHIK(icomp,0,1,tau,del)  !real-gas terms
        phi02=PHIK(icomp,0,2,tau,del)
        phi11=PHIK(icomp,1,1,tau,del)
        phi20=PHIK(icomp,2,0,tau,del)
c
        phig20=PHI0K(icomp,2,0,t,rho)  !ideal-gas term
        cv=-R*(phi20+phig20)
        cp=cv+R*(1.0d0+phi01-phi11)**2/(1.0d0+2.0d0*phi01+phi02)
      end if
c
      RETURN
      end                                              !subroutine CVCPK
c
c ======================================================================
c
      subroutine GIBBS (t,rho,x,Ar,Gr)
c
c  compute residual Helmholtz and Gibbs free energy as a function of
c  temperature, density, and composition using core functions
c
c  N.B.  The quantity calculated is
c
c        G(T,rho) - G0(T,P*) = G(T,rho) - G0(T,rho) + RTln(RTrho/P*)
c
c        where G0 is the ideal gas state and P* is a reference pressure
c        which is equal to the current pressure of interest.  Since Gr
c        is used only as a difference in phase equilibria calculations
c        where the temperature and pressure of the phases are equal, the
c        (RT/P*) part of the log term will cancel and is omitted.
c
c        "normal" (not residual) A and G are computed by subroutine AG
c
c  based on derivations in Younglove & McLinden, JPCRD 23 #5, 1994,
c  equations A8 - A12
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c       Ar--residual Helmholtz free energy [J/mol]
c       Gr--residual Gibbs free energy [J/mol]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-07-94  MM, original version
c  08-07-95  MM, add calls to Fundamental (Helmholtz) EOS
c  10-03-95  MM, change /MODEL/--models specified by strings
c  11-03-95  MM, add calls to mixture Helmholtz (HMX) model
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                replace calls to PHIHMX, PHIFEQ with general PHIX, PHIK
c  03-22-96  MM, replace /MODEL/ with /EOSMOD/
c  10-01-97  MM, add compiler switch to allow access by DLL
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: GIBBS
c     dll_export GIBBS
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      Ar=0.d0
      Gr=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        Ar=ABWR(icomp,t,rho)
        Gr=Ar+PBWR(icomp,t,rho)/rho+R*t*(-1.0d0+log(rho))
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi00=PHIK(icomp,0,0,tau,del)  !real-gas terms
          phi01=PHIK(icomp,0,1,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi00=PHIX(0,0,tau,del,x)  !real-gas terms
          phi01=PHIX(0,1,tau,del,x)
        end if
c
        RT=R*t
        Ar=RT*phi00
        Gr=Ar+RT*(1.0d0+phi01)+RT*(-1.0d0+log(rho))
      end if
c
      RETURN
      end                                              !subroutine GIBBS
c
c ======================================================================
c
      subroutine AG (t,rho,x,a,g)
c
c  compute Helmholtz and Gibbs energies as a function of temperature,
c  density, and composition.
c
c  N.B.  These are not residual values (those are calculated by GIBBS).
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c        a--Helmholtz energy [J/mol]
c        g--Gibbs free energy [J/mol]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c   3-27-98 EWL, original version
c  12-02-98 EWL, reorganize code so to eliminate x(i) in pure fluid calculation
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: AG
c     dll_export AG
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      a=0.d0
      g=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (rho.lt.1.0d-40) then
c  calc will crash if rho = 0
        rhos=1.0d-40
      else
        rhos=rho
      end if
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        p=PBWR(icomp,t,rho)
        a=ABWR(icomp,t,rho)
        cpiint=CPI(t,x)
        cptint=CPT(t,x)
        a=a+cpiint-R*t-href(icomp)
     &   -t*(R*log(rhoref(icomp)/rhos)+cptint-R*log(t/tref(icomp))
     &    -sref(icomp))
        if (rho.gt.0.d0) then
          g=a+p/rho
        else
          g=a+R*t
        endif
c
      else
        RT=R*t
        phig00=PHI0(0,0,t,rhos,x)      !ideal-gas terms
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi00=PHIK(icomp,0,0,tau,del)  !real-gas terms
          phi01=PHIK(icomp,0,1,tau,del)
          a=RT*(phig00+phi00)-(href(icomp)-t*sref(icomp))
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi00=PHIX(0,0,tau,del,x)   !real-gas terms
          phi01=PHIX(0,1,tau,del,x)
          a=RT*(phig00+phi00)
          do i=1,nc
            a=a-x(i)*(href(i)-t*sref(i))
          enddo
        end if
        g=a+RT*(1.0d0+phi01)
      end if
c
      RETURN
      end                                                 !subroutine AG
c
c ======================================================================
c
      subroutine PRESS (t,rho,x,p)
c
c  compute pressure as a function of temperature,
c  density, and composition using core functions
c
c  direct implementation of core function of corresponding model
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  output:
c        p--pressure [kPa]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  11-19-94  MM, original version
c  08-07-95  MM, add calls to Fundamental (Helmholtz) EOS
c  10-03-95  MM, change /MODEL/--models specified by strings
c  11-03-95  MM, add calls to mixture Helmholtz (HMX) model
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                replace calls to PHIHMX, PHIFEQ with general PHIX, PHIK
c  03-22-96  MM, replace /MODEL/ with /EOSMOD/
c  10-01-97  MM, add compiler switch to allow access by DLL
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PRESS
c     dll_export PRESS
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      p=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        p=PBWR(icomp,t,rho)
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          p=R*t*rho*(1.0d0+PHIK(icomp,0,1,tau,del))
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          p=R*t*rho*(1.0d0+PHIX(0,1,tau,del,x))
        end if
      end if
c
      RETURN
      end                                              !subroutine PRESS
c
c ======================================================================
c
      subroutine DPDD (t,rho,x,dpdrho)
c
c  compute partial derivative of pressure w.r.t. density at constant
c  temperature as a function of temperature, density, and composition
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  output:
c   dpdrho--dP/drho [kPa-L/mol]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  04-23-95  MM, original version
c  08-07-95  MM, add calls to Fundamental (Helmholtz) EOS
c  10-03-95  MM, change /MODEL/--models specified by strings
c  11-03-95  MM, add calls to mixture Helmholtz (HMX) model
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                replace calls to PHIHMX, PHIFEQ with general PHIX, PHIK
c  03-22-96  MM, replace /MODEL/ with /EOSMOD/
c  10-16-96  MM, change name from DPRHO to DPDD
c  10-01-97  MM, add compiler switch to allow access by DLL
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DPDD
c     dll_export DPDD
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      dpdrho=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        dpdrho=DPDBWR(icomp,t,rho)
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi01=PHIK(icomp,0,1,tau,del)  !real-gas terms
          phi02=PHIK(icomp,0,2,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi01=PHIX(0,1,tau,del,x)  !real-gas terms
          phi02=PHIX(0,2,tau,del,x)
        end if
        dpdrho=R*t*(1.0d0+2.0d0*phi01+phi02)
      end if
c
      RETURN
      end                                               !subroutine DPDD
c
c ======================================================================
c
      subroutine DPDDK (icomp,t,rho,dPdrho)
c
c  compute partial derivative of pressure w.r.t. density at constant
c  temperature as a function of temperature and density for a specified
c  component
c
c  analogous to DPDD, except for component icomp, this is used by transport
c  routines to calculate dP/dD for the reference fluid (component zero)
c
c  inputs:
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c        t--temperature [K]
c      rho--molar density [mol/L]
c  output:
c   dPdrho--dP/drho [kPa-L/mol]
c
c  written by M. McLinden, NIST Physical & Chem Properties Div, Boulder, CO
c  06-16-97  MM, original version; based on DPDD
c  09-29-97  MM, if component uses MBWR, call DPDBWR
c  10-01-97  MM, add compiler switch to allow access by DLL
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DPDDK
c     dll_export DPDDK
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
      dpdrho=0.d0
      if (t.le.0.d0) return
      if (hmxeos(icomp).eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        dpdrho=DPDBWR(icomp,t,rho)
c
      else
c  call general PHIK or PHIX routines for all other models
        tau=tz(icomp)/t
        del=rho/rhoz(icomp)
        phi01=PHIK(icomp,0,1,tau,del)  !real-gas terms
        phi02=PHIK(icomp,0,2,tau,del)
        dpdrho=R*t*(1.0d0+2.0d0*phi01+phi02)
      end if
c
      RETURN
      end                                              !subroutine DPDDK
c
c ======================================================================
c
      subroutine DPDD2 (t,rho,x,d2PdD2)
c
c  compute second partial derivative of pressure w.r.t. density at const
c  temperature as a function of temperature, density, and composition
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  output:
c   d2pdD2--d^2P/drho^2 [kPa-L^2/mol^2]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  06-03-97 EWL, original version
c  10-01-97  MM, add compiler switch to allow access by DLL
c  02-11-99 EWL, skip calculation of d2PdD2 if rho=0
c  09-05-02 EWL, add ideal gas d2PdD2
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DPDD2
c     dll_export DPDD2
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      d2PdD2=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        d2PdD2=D2PBWR(icomp,t,rho)
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi01=PHIK(icomp,0,1,tau,del)  !real-gas terms
          phi02=PHIK(icomp,0,2,tau,del)
          phi03=PHIK(icomp,0,3,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi01=PHIX(0,1,tau,del,x)  !real-gas terms
          phi02=PHIX(0,2,tau,del,x)
          phi03=PHIX(0,3,tau,del,x)
        end if
        if (rho.gt.1.0d-8) then
          d2PdD2=R*t/rho*(2.0d0*phi01+4.0d0*phi02+phi03)
        else
          call VIRB (t,x,b)
          d2PdD2=2.d0*b*R*t
        end if
      end if
c
      RETURN
      end                                              !subroutine DPDD2
c
c ======================================================================
c
      subroutine DPDT (t,rho,x,dpt)
c
c  compute partial derivative of pressure w.r.t. temperature at constant
c  density as a function of temperature, density, and composition
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  output:
c      dpt--dP/dT [kPa/K]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  10-16-96  MM, original version, based on DPDD
c  10-28-96  MM, insert missing rho into form using PHI's
c  10-01-97  MM, add compiler switch to allow access by DLL
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DPDT
c     dll_export DPDT
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      dpt=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        dpt=DPTBWR(icomp,t,rho)
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi01=PHIK(icomp,0,1,tau,del)  !real-gas terms
          phi11=PHIK(icomp,1,1,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi01=PHIX(0,1,tau,del,x)  !real-gas terms
          phi11=PHIX(1,1,tau,del,x)
        end if
        dpt=R*rho*(1.0d0+phi01-phi11)
      end if
c
      RETURN
      end                                               !subroutine DPDT
c
c ======================================================================
c
      subroutine DPDTK (icomp,t,rho,dpt)
c
c  compute partial derivative of pressure w.r.t. temperature at constant
c  density as a function of temperature and density for a specified component
c
c  analogous to DPDT, except for component icomp, this is used by transport
c  routines to calculate dP/dT
c
c  inputs:
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c        t--temperature [K]
c      rho--molar density [mol/L]
c  output:
c      dpt--dP/dT [kPa/K]
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  07-07-98 EWL, original version, based on DPDT
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DPDTK
c     dll_export DPDTK
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
      dpt=0.d0
      if (t.le.0.d0) return
      if (hmxeos(icomp).eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        dpt=DPTBWR(icomp,t,rho)
c
      else
c  call general PHIK or PHIX routines for all other models
        tau=tz(icomp)/t
        del=rho/rhoz(icomp)
        phi01=PHIK(icomp,0,1,tau,del)  !real-gas terms
        phi11=PHIK(icomp,1,1,tau,del)
        dpt=R*rho*(1.0d0+phi01-phi11)
      end if
c
      RETURN
      end                                              !subroutine DPDTK
c
c ======================================================================
c
      subroutine DDDP (t,rho,x,drhodp)
c
c  compute partial derivative of density w.r.t. pressure at constant
c  temperature as a function of temperature, density, and composition
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  output:
c   drhodp--drho/dP [mol/(L-kPa)]
c
c  written by M. McLinden, NIST Phys & Chem Properties Div, Boulder, CO
c  08-29-97  MM, original version, based on DPDD (just the inverse)
c  10-01-97  MM, add compiler switch to allow access by DLL
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DDDP
c     dll_export DDDP
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      drhodp=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        drhodp=1.0d0/DPDBWR(icomp,t,rho)
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi01=PHIK(icomp,0,1,tau,del)  !real-gas terms
          phi02=PHIK(icomp,0,2,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi01=PHIX(0,1,tau,del,x)  !real-gas terms
          phi02=PHIX(0,2,tau,del,x)
        end if
        drhodp=1.0d0/(R*t*(1.0d0+2.0d0*phi01+phi02))
      end if
c
      RETURN
      end                                               !subroutine DDDP
c
c ======================================================================
c
      subroutine DDDT (t,rho,x,drhodt)
c
c  compute partial derivative of density w.r.t. temperature at constant
c  pressure as a function of temperature, density, and composition
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  output:
c   drhodt--drho/dT [mol/(L-K)]
c
c   d(rho)/d(T) = -d(rho)/dP x dP/dT = -dP/dT / (dP/d(rho))
c
c  written by M. McLinden, NIST Phys & Chem Properties Div, Boulder, CO
c  08-29-97  MM, original version, based on DPDD and DPDT
c  10-01-97  MM, add compiler switch to allow access by DLL
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DDDT
c     dll_export DDDT
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      drhodt=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        drhodt=-DPTBWR(icomp,t,rho)/DPDBWR(icomp,t,rho)
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi01=PHIK(icomp,0,1,tau,del)  !real-gas terms
          phi11=PHIK(icomp,1,1,tau,del)
          phi02=PHIK(icomp,0,2,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi01=PHIX(0,1,tau,del,x)  !real-gas terms
          phi11=PHIX(1,1,tau,del,x)
          phi02=PHIX(0,2,tau,del,x)
        end if
        drhodt=-rho*(1.0d0+phi01-phi11)/(t*(1.0d0+2.0d0*phi01+phi02))
      end if
c
      RETURN
      end                                               !subroutine DDDT
c
c ======================================================================
c
      subroutine DHD1(t,rho,x,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
c
c  compute partial derivatives of enthalpy w.r.t. t, p, or rho at constant
c  t, p, or rho as a function of temperature, density, and composition
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c   dhdt_d--dH/dT at constant density [J/(mol-K)]
c   dhdt_p--dH/dT at constant pressure [J/(mol-K)]
c   dhdd_t--dH/drho at constant temperature [(J/mol)/(mol/L)]
c   dhdd_p--dH/drho at constant pressure [(J/mol)/(mol/L)]
c   dhdp_t--dH/dP at constant temperature [J/(mol-kPa)]
c   dhdp_d--dH/dP at constant density [J/(mol-kPa)]
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  10-02-00 EWL, original version
c  05-30-06 EWL, change subroutine name from DHDT to DHD1, and add other
c                derivates of h with respect to t, p, or rho
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DHD1
c     dll_export DHD1
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      dimension x(ncmax)
c
      dhdt_d=0.d0
      dhdt_p=0.d0
      dhdd_t=0.d0
      dhdd_p=0.d0
      dhdp_t=0.d0
      dhdp_d=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      rhos=rho
      if (rho.lt.1.0d-10) rhos=1.0d-10
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        tau=tz(icomp)/t
        del=rho/rhoz(icomp)
        phi01=PHIBWR(icomp,0,1,tau,del)  !real-gas terms
        phi11=PHIBWR(icomp,1,1,tau,del)
        phi20=PHIBWR(icomp,2,0,tau,del)
        phi02=PHIBWR(icomp,0,2,tau,del)
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi01=PHIK(icomp,0,1,tau,del)  !real-gas terms
          phi11=PHIK(icomp,1,1,tau,del)
          phi20=PHIK(icomp,2,0,tau,del)
          phi02=PHIK(icomp,0,2,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi01=PHIX(0,1,tau,del,x)  !real-gas terms
          phi11=PHIX(1,1,tau,del,x)
          phi20=PHIX(2,0,tau,del,x)
          phi02=PHIX(0,2,tau,del,x)
        end if
      end if
      phig20=PHI0(2,0,t,rho,x)
      phig11=PHI0(1,1,t,rho,x)
      call THERM2 (t,rhos,x,p,e,h,s,cv,cp,w,Z,hjt,A,G,xkappa,beta,
     &      dPdrho,d2PdD2,dPT,drhodT,drhodP,d2PT2,d2PdTD,spare3,spare4)
      dhdt_p=cp
      dhdt_d=R*(-phig20-phi20+phi01-phi11+1.d0)
      if (rho.gt.1.0d-8) then
        dhdp_t=1.d0/rho+t*drhodT/rho**2
        dhdd_t=R*T/rho*(phig11+phi11+phi01+phi02)
        dhdp_d=dhdp_t+dhdt_p/dPT
        dhdd_p=dhdd_t+dhdt_d/drhodT
      else
        call VIRB (t,x,b)
        call DBDT (t,x,dbt)
        dhdp_t=1.d0/rhos+t*drhodT/rhos**2
        dhdd_t=-r*t**2*dbt+r*t*b
        dhdp_d=xinf
        dhdd_p=xinf
      endif
c
      RETURN
      end                                               !subroutine DHD1
c
c ======================================================================
c
      subroutine FGCTY2 (t,rho,x,f,ierr,herr)
c
c  compute fugacity for each of the nc components of a mixture by
c  analytical differentiation of the dimensionless residual Helmholtz energy
c
c  based on derivations in the GERG-2004 document for natural gas
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c        f--array (1..nc) of fugacities [kPa]
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  02-22-10 EWL, original version
c
cDEC$ ATTRIBUTES DLLEXPORT :: FGCTY2
c     dll_export FGCTY2
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      dimension x(ncmax),f(ncmax)
      common /NCOMP/ nc,ic
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,xerr,xnotd,xnotc
      character*255 herr
c
      ierr=0
      herr=' '
      call ISPURE (x,icomp)
      call RMIX (x)
      do i=1,nc
        f(i)=0.0d0
      enddo
      if (t.le.0.d0 .or. rho.lt.1.0d-40) RETURN
      if (icomp.ne.0) goto 10  !Call old FGCTY routine for pure fluids
c
      RTrho=R*t*rho
      do i=1,nc
        if (x(i).gt.0.d0) then
          call PHIDERV (i,0,t,rho,x,dadn,dnadn,ierr,herr)
          if (ierr.ne.0) goto 10
          f(i)=xerr
          if (ABS(dnadn).lt.100.0d0) f(i)=x(i)*RTrho*exp(dnadn)
        endif
      enddo
      RETURN
c
 10   continue
      call FGCTY (t,rho,x,f)
      ierr=0
      herr=' '
      RETURN
      end                                             !subroutine FGCTY2
Cc
Cc ======================================================================
Cc
C      subroutine GIBDERV (i,t,rho,x,dgdxi,dgdxi2,dvdxar01,
C     &                    da0dxi,da0dxii,ierr,herr)
Cc
Cc  calculate derivatives of the Gibbs energy with respect to composition
Cc
Cc  based on derivations in the GERG-2004 document for natural gas
Cc
Cc  inputs:
Cc        i--component number of which to take derivative
Cc        t--temperature (K)
Cc      rho--density (mol/L)
Cc        x--composition [array of mol frac]
Cc  outputs: (where n is mole number)
Cc    dgdxi--partial of reduced Gibbs energy with respect to xi (Eq. x.xx in GERG)
Cc   dgdxi2--2nd partial of reduced Gibbs energy with respect to xi (Eq. x.xx in GERG)
Cc
Cc  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
Cc  03-09-10 EWL, original version
Cc
C      implicit double precision (a-h,o-z)
C      implicit integer (i-k,m,n)
C      parameter (ncmax=20)        !max number of components in mixture
C      parameter (nrefmx=10)       !max number of fluids for transport ECS
C      parameter (n0=-ncmax-nrefmx,nx=ncmax)
C      dimension x(ncmax),xliq(ncmax),xvap(ncmax)
C      common /NCOMP/ nc,ic
C      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
C      character*255 herr
Cc
C      call RMIX (x)
C      dgdxi=0
C      dgdxi2=0
C      call RDXHMX (0,i,0,x,t0,D0,ierr,herr)
C      call RDXHMX (1,i,0,x,dtdx,dvdx,ierr,herr)
C      if (ierr.ne.0) RETURN
C      tau=t0/t
C      del=rho/d0
C      call PRESS (t,rho,x,p)
C      call PHIDERV (i,0,t,rho,x,dadn,dnadn,ierr,herr)
C      if (ierr.ne.0) RETURN
C      call RESETA
C      icc=ic
C      ic=i
C      call TPFLSH (t,p,x,dd,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
C      call INFO (i,wm,ttp,tnbp,tc,pc,dc,zc,acf,dip,rgas)
C      ic=icc
C      a0i=PHI0K(i,0,0,t,rho)
C      dd=p/r/t
C      da0dxi=a0i+log(x(i))+1.d0
C      da0dxii=1.d0/x(i)
C      ar01=PHIX(0,1,tau,del,x)
C      dvdxar01=d0*dvdx*ar01
C      dii=0  !Fix this
C      dadxii=0  !Fix this
Cc     dgdxi=da0dxi+dadxi+dvdxar01+daddx(i)
C      dgdxi=da0dxi+dadxi         +daddx(i)
Cc     dgdxi=       dadxi         +daddx(i)
C      dgdxi2=da0dxii+dadxii+dii*ar01+2.d0*daddx(i)+daddxii
C      RETURN
C      end                                            !subroutine GIBDERV
c
c ======================================================================
c
      subroutine PHIDERV (i,j,t,rho,x,dadn,dnadn,ierr,herr)
c
c  calculate various derivatives needed for VLE determination
c
c  based on derivations in the GERG-2004 document for natural gas
c
c  inputs:
c        i--component number of which to take derivative
c        j--component number of which to take derivative (can be set to zero)
c        t--temperature (K)
c      rho--density (mol/L)
c        x--composition [array of mol frac]
c  outputs: (where n is mole number)
c     dadn--n*partial(alphar)/partial(ni)                 (Eq. 7.16 in GERG)
c    dnadn--partial(n*alphar)/partial(ni)                 (Eq. 7.15 in GERG)
c
c     dtdn--n*[partial(Tred)/partial(ni)]/Tred            (Eq. 7.19 in GERG)
c     dvdn--n*[partial(Vred)/partial(ni)]/Vred            (Eq. 7.18 in GERG)
c           (=-n*[partial(Dred)/partial(ni)]/Dred)
c    daddn--del*n*partial(darddel)/partial(ni)            (Eq. 7.17 in GERG)
c           where darddel=partial(alphar)/partial(del)
c   d2adnn--n*partial^2(n*alphar)/partial(ni)/partial(nj) (Eq. 7.46 in GERG)
c  the following are at constant tau and/or del
c    dadxi--partial(alphar)/partial(xi)                   (Eq. 7.21g in GERG)
c   sdadxi--sum[xi*partial(alphar)/partial(xi)]           (Eq. 7.21g in GERG)
c   dadxij--partial^2(alphar)/partial(xi)/partial(xj)     (Eq. 7.21i in GERG)
c    daddx--del*partial^2(alphar)/partial(xi)/partial(del) (Eq. 7.21j in GERG)
c  daddxii--del*partial^3(alphar)/partial(xi)/partial(xj)/partial(del)
c
c  other calculated variables:
c   d2addn--del*par.(n*(par.(alphar)/par.(n)/par.(del)    (Eq. 7.50 in GERG)
c   d2adtn--tau*par.(n*(par.(alphar)/par.(n)/par.(tau)    (Eq. 7.51 in GERG)
c   d2adxn--par.(n*(par.(alphar)/par.(n)/par.(xj)         (Eq. 7.52 in GERG)
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  02-22-10 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      dimension x(ncmax)
      common /NCOMP/ nc,ic
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      character*255 herr
      common /PHIDR/ dtdn,dvdn,dadxi,daddn,d2adnn,d2addn,d2adtn,d2adxn,
     &               dadxij,daddx,dadtx,ddrdxn,dtrdxn,dpdn
      dimension dtdn(ncmax),dvdn(ncmax),dadxij(ncmax,ncmax),
     &          dadxi(ncmax),daddx(ncmax),dadtx(ncmax),
     &          ddrdxn(ncmax,ncmax),dtrdxn(ncmax,ncmax),
     &          d2adxn(ncmax,ncmax),dpdn(ncmax)
      dimension aok(ncmax),aok01(ncmax),aok10(ncmax),
     &          sdadxi(ncmax),sdaddxi(ncmax),sdadtxi(ncmax),
     &          tr01(ncmax),dr01(ncmax),sdadxx(ncmax)
      common /PHIDSV/ phimxk(ncmax,ncmax),
     &                phimxk01(ncmax,ncmax),phimxk10(ncmax,ncmax),
     &                dr11(ncmax,ncmax),tr11(ncmax,ncmax)
     &                ,ierrdr
c
c     data ierrdr /0/
      ierr=0
      herr=' '
      call RMIX (x)
      call REDX (x,t0,rho0)
      tau=t0/t
      del=rho/rho0
      RT=R*t
      daddxii=0.d0

      daddn=0.d0
      d2adnn=0.d0
      d2addn=0.d0
      d2adtn=0.d0
      do m=1,nc
        daddx(m)=0.d0
        dadtx(m)=0.d0
        dtdn(m)=0.d0
        dvdn(m)=0.d0
        dadxi(m)=0.d0
        aok(m)=0.d0
        aok01(m)=0.d0
        aok10(m)=0.d0
        sdadxi(m)=0.d0
        sdaddxi(m)=0.d0
        sdadtxi(m)=0.d0
        tr01(m)=0.d0
        dr01(m)=0.d0
        dpdn(m)=0.d0
        sdadxx(m)=0.d0
      enddo
      do k=1,nc
      do m=1,nc
        ddrdxn(k,m)=0.d0
        dtrdxn(k,m)=0.d0
        d2adxn(k,m)=0.d0
        dadxij(k,m)=0.d0     !if i=j, derivative=0
      enddo
      enddo
c
      if (i.lt.1 .or. x(i).le.0.d0) then
        ierr=1
        RETURN
      endif
c
c  get Helmholtz energy of pure fluid
      ar  =PHIX(0,0,tau,del,x)
      ar01=PHIX(0,1,tau,del,x)
      ar10=PHIX(1,0,tau,del,x)

c *** check ierrdr for old error ***

      do k=1,nc
        phimxk(k,k)=0.d0
        if (x(k).gt.0.d0) then
          call RDXHMX (1,k,0,x,tr01(k),dr01(k),ierr,herr)
          ierrdr=ierr
          if (ierr.ne.0) RETURN
          aok(k)=PHIK(k,0,0,tau,del)
          if (k.ne.nc) then
            do m=k+1,nc
              if (x(m).gt.0.d0) then
                phimxk(m,k)=PHIMIX(m,k,0,0,tau,del,x)
                phimxk(k,m)=phimxk(m,k)
              endif
            enddo
          endif
        endif
      enddo
      if (j.ne.0) then
        do k=1,nc
          phimxk01(k,k)=0.d0
          phimxk10(k,k)=0.d0
          if (x(k).gt.0.d0) then
            aok01(k)=PHIK(k,0,1,tau,del)
            aok10(k)=PHIK(k,1,0,tau,del)
            if (k.ne.nc) then
              do m=k+1,nc
                if (x(m).gt.0.d0) then
                  phimxk01(m,k)=PHIMIX(m,k,0,1,tau,del,x)
                  phimxk10(m,k)=PHIMIX(m,k,1,0,tau,del,x)
                  phimxk01(k,m)=phimxk01(m,k)
                  phimxk10(k,m)=phimxk10(m,k)
                endif
              enddo
            endif
            do m=1,nc
              if (x(m).gt.0.d0) then
                ij=11
                if (m.eq.k) ij=2
                call RDXHMX (ij,k,m,x,tr11(k,m),dr11(k,m),ierr,herr)
                dr11(k,m)=2*rho0**3*dr01(m)*dr01(k)-rho0**2*dr11(k,m)
                ierrdr=ierr
                if (ierr.ne.0) RETURN
              endif
            enddo
          endif
        enddo
      endif
c
c  calculate first order derivatives only
      if (j.eq.0) then
        do k=1,nc
          if (x(k).gt.0.d0) then
c  get derivatives of reducing parameters
            dtdn(i)=dtdn(i)-x(k)*tr01(k)
            dvdn(i)=dvdn(i)-x(k)*dr01(k)
            sdadxi(i)=sdadxi(i)+x(k)*aok(k)
            if (i.eq.k) then
              dtdn(i)=dtdn(i)+tr01(k)
              dvdn(i)=dvdn(i)+dr01(k)
              dadxi(i)=dadxi(i)+aok(k)
            else
c  add excess Helmholtz energy of i-k interaction
              dadxi(i)=dadxi(i)+phimxk(i,k)/x(i)
            endif
            if (k.ne.nc) then
              do m=k+1,nc
c  subtract excess Helmholtz energy of m-k interaction twice (to include k-m)
                if (x(m).ne.0.d0) sdadxi(i)=sdadxi(i)+phimxk(m,k)*2.d0
              enddo
            endif
          endif
        enddo
c
c  *** do this only while testing numerical derivatives!!! ***
        call RDXHMX (-1,0,0,x,t0,rho0,ierr,herr)

        dtdn(i)=dtdn(i)/t0
        dvdn(i)=dvdn(i)*rho0
        dadn=ar01*(1.d0+dvdn(i))+ar10*dtdn(i)+dadxi(i)-sdadxi(i)     !Eq. 7.16 in GERG
        dnadn=ar+dadn                                    !Eq. 7.15
c
c  calculate first and second order derivatives
      else
        if (x(j).le.0.d0) then
          ierr=2
          RETURN
        endif
        ar02=PHIX(0,2,tau,del,x)
        ar20=PHIX(2,0,tau,del,x)
        ar11=PHIX(1,1,tau,del,x)

        do m=1,nc    !j loop
        do n=1,nc    !i loop
        ddrdxn(n,m)= dr01(m)*rho0**2                              !Eq. 7.55
        dtrdxn(n,m)=-tr01(m)                                      !Eq. 7.56
        do k=1,nc
          if (x(k).gt.0.d0 .and. x(n).gt.0.d0 .and.
     &        x(m).gt.0.d0) then
            if (ierr.ne.0) RETURN
            ddrdxn(n,m)=ddrdxn(n,m)-x(k)*dr11(k,m)
            dtrdxn(n,m)=dtrdxn(n,m)-x(k)*tr11(k,m)
            if (n.eq.k) then
              ddrdxn(n,m)=ddrdxn(n,m)+dr11(k,m)
              dtrdxn(n,m)=dtrdxn(n,m)+tr11(k,m)
            endif
          endif
        enddo
        enddo
        enddo
        do n=1,nc
        do k=1,nc
          if (x(k).gt.0.d0 .and. x(n).gt.0.d0) then
            dtdn(n)=dtdn(n)-x(k)*tr01(k)
            dvdn(n)=dvdn(n)-x(k)*dr01(k)
            sdadxi(n)=sdadxi(n)+x(k)*aok(k)
            sdaddxi(n)=sdaddxi(n)+x(k)*aok01(k)
            sdadtxi(n)=sdadtxi(n)+x(k)*aok10(k)
            if (n.eq.k) then
              daddx(n)=daddx(n)+aok01(k)
              dadtx(n)=dadtx(n)+aok10(k)
              dtdn(n)=dtdn(n)+tr01(k)
              dvdn(n)=dvdn(n)+dr01(k)
              dadxi(n)=dadxi(n)+aok(k)
            endif
c  add excess Helmholtz energy of i-k interaction
            dadxi(n)=dadxi(n)+phimxk(n,k)/x(n)
            daddx(n)=daddx(n)+phimxk01(n,k)/x(n)
            dadtx(n)=dadtx(n)+phimxk10(n,k)/x(n)
            if (k.ne.nc) then
              do m=k+1,nc
c  subtract excess Helmholtz energy of m-k interaction twice (to include k-m)
                sdadxi(n) =sdadxi(n) +phimxk(m,k)*2.d0
                sdaddxi(n)=sdaddxi(n)+phimxk01(m,k)*2.d0
                sdadtxi(n)=sdadtxi(n)+phimxk10(m,k)*2.d0
              enddo
            endif
          endif
        enddo
        dtdn(n)=dtdn(n)/t0
        dvdn(n)=dvdn(n)*rho0
        enddo
c
        dvdn1=1.d0+dvdn(i)
        dadn =ar01*dvdn1+ar10*dtdn(i)+dadxi(i)-sdadxi(i)     !Eq. 7.16 in GERG
        daddn=ar02*dvdn1+ar11*dtdn(i)+daddx(i)-sdaddxi(i)    !Eq. 7.17
        dnadn=ar+dadn                                        !Eq. 7.15
        d2addn=(ar01+ar02)*dvdn1+ar11*dtdn(i)+daddx(i)-sdaddxi(i)  !Eq. 7.50
        d2adtn=ar11*dvdn1+(ar10+ar20)*dtdn(i)+dadtx(i)-sdadtxi(i)  !Eq. 7.51


        sd2adxn=0.d0
        do k=1,nc
          if (x(k).gt.0.d0) then
          daddnk=ar02*(1.d0+dvdn(k))+ar11*dtdn(k)+daddx(k)-sdaddxi(k) !Eq. 7.17
          dpdn(k)=rho*R*t*(1.d0+ar01*(2.d0+dvdn(k))+daddnk)  !Equation (7.63) in GERG
          endif
        enddo
        do m=1,nc
        do k=1,nc
          if (x(k).gt.0.d0 .and. x(m).gt.0.d0) then
            dadxij(k,m)=phimxk(k,m)/x(m)/x(k)
            sdadxx(m)=sdadxx(m)+x(k)*dadxij(k,m)
          endif
        enddo
        enddo
        do m=1,nc
        do k=1,nc
          if (x(k).gt.0.d0 .and. x(m).gt.0.d0) then
          d2adxn(k,m)=daddx(m)*(1.d0+dvdn(k))                  !Eq. 7.52
     &             -ar01/rho0*(ddrdxn(k,m)-rho0**2*dr01(m)*dvdn(k))
     &             +dadtx(m)*dtdn(k)
     &             +ar10/t0*(dtrdxn(k,m)-tr01(m)*dtdn(k))
     &             +dadxij(k,m)-dadxi(m)-sdadxx(m)
          endif
        enddo
        enddo

        do k=1,nc
          if (x(k).gt.0.d0) then
            sd2adxn=sd2adxn-x(k)*d2adxn(i,k)
            if (j.eq.k) sd2adxn=sd2adxn+d2adxn(i,k)
          endif
        enddo
        dadnj=ar01*(1.d0+dvdn(j))+ar10*dtdn(j)+dadxi(j)-sdadxi(j)  !Eq. 7.16 in GERG
        d2adnn=dadnj+d2addn*(1.d0+dvdn(j))+d2adtn*dtdn(j)+sd2adxn  !Eq. 7.46
      endif

      RETURN
      end                                            !subroutine PHIDERV
c
c ======================================================================
c
      subroutine FGCTY (t,rho,x,f)
c
c  old routine to compute fugacity for each of the nc components of a mixture
c  by numerical differentiation (using central differences) of the
c  dimensionless residual Helmholtz energy
c
c  based on derivations in E.W. Lemmon, MS Thesis, University of Idaho
c  (1991); section 3.2
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c        f--array (1..nc) of fugacities [kPa]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  12-15-95  MM, original version
c  12-18-95  MM, add pure component fugacity as a special case
c  01-08-96  MM, bug on call to PHIFEQ (wrong arguments)
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                add Zcrit to common /CCON/
c                replace calls to PHIHMX, PHIFEQ with general PHIX, PHIK
c  03-19-19  MM, add dipole moment to /CCON/
c  03-22-96  MM, replace /MODEL/ with /EOSMOD/
c  10-01-97  MM, add compiler switch to allow access by DLL
c  12-16-97  MM, add check for rho = 0; overflow on exponent (set to xerr)
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-22-98 EWL, recalculate R for mixtures based on values for pure fluids
c  05-08-06 EWL, modify how delp and deln are calculated for x>0.9999
c  01-25-07 EWL, change default f(i) from 0 to 1.  Skip calculation if x(i)=0
c  02-26-09 EWL, set f(icomp) equal to a very large number instead of xerr when the variable arg is huge
c  11-20-09 BFT, change deln to delmol in check for x(i).gt.1.0d0-delmol
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: FGCTY
c     dll_export FGCTY
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      dimension x(ncmax),f(ncmax)
      dimension xplus(ncmax),xminus(ncmax)
      character*3 hpheq,heos,hmxeos,hmodcp
      common /NCOMP/ nc,ic
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,xerr,xnotd,xnotc
      delmol=1.0d-4
      call ISPURE (x,icomp)
c
c  fill output fugacity array with zeros (final value for undefined
c  components and insurance against problems for others)
      do i=1,nc
        f(i)=0.0d0
      enddo
      if (t.le.0.d0) return
c
      call RMIX (x)
c  check for zero input density
      if (rho.lt.1.0d-40) RETURN
c
      RTrho=R*t*rho
      if (icomp.ne.0) then
c  pure component
        if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--use BWR-specific routines
          Ar=ABWR(icomp,t,rho)
          p=PBWR(icomp,t,rho)
          f(icomp)=RTrho*exp(Ar/(R*t)+p/RTrho-1.0d0)
        else
c  for other models, use general PHIK routines
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi00=PHIK(icomp,0,0,tau,del)
          phi01=PHIK(icomp,0,1,tau,del)
c  check for potential under- or over-flow (can happen in 2-phase, but
c  the fugacity is meaningless there anyway)
          arg=phi00+phi01
          if (ABS(arg).lt.500.0d0) then
            f(icomp)=RTrho*exp(arg)
          else
            f(icomp)=1.d100
          end if
        end if
      else
c
c  mixture
        do i=1,nc
c  compute positive and negative increments to number of moles
c  general case:  deln < x(i) < 1 - deln
          delp=delmol
          deln=-delmol
          if (x(i).gt.0.d0) then
          if (x(i).lt.delmol) then
c  special case--composition of component i is nearly zero
            deln=-x(i)/2.d0
            delp=-deln
          else if (x(i).gt.1.0d0-delmol) then
c  special case--composition of component i is nearly one (pure fluid)
            delp=(1.0d0-x(i))/2.d0
            deln=-delp
          end if
          delp1=1.0d0/(1.0d0+delp)
          deln1=1.0d0/(1.0d0+deln)
c  since total number of moles is now 1 + (delp or deln), all of the
c  compositions have changed
          do j=1,nc
            xplus(j)=x(j)*delp1
            xminus(j)=x(j)*deln1
          enddo
          xplus(i)=(x(i)+delp)*delp1
          xminus(i)=(x(i)+deln)*deln1
c  derivative is at constant volume, so must adjust density
          Dplus=rho*(1.0d0+delp)
          Dminus=rho*(1.0d0+deln)
c  compute residual Helmholtz at 'plus' and 'minus' density and composition
c  could call subroutine GIBBS here, but more efficient to directly call
c  the core routines (via PHIX)
          call REDX (xplus,t0,rho0)
          tau=t0/t
          del=Dplus/rho0
          Aplus=PHIX(0,0,tau,del,xplus)          !real-gas terms
          call REDX (xminus,t0,rho0)
          tau=t0/t
          del=Dminus/rho0
          Aminus=PHIX(0,0,tau,del,xminus)        !real-gas terms
          dnadn=((1.0d0+delp)*Aplus-(1.0d0+deln)*Aminus)/(delp-deln)
c         write (*,*) ' FGCTY--delp,deln:  ',delp,deln
c         write (*,*) ' FGCTY--i,A+, A-, dAdN: ',i,Aplus,Aminus,dnadn
c  check for potential under- or over-flow (can happen in 2-phase, but
c  the fugacity is meaningless there anyway)
          f(i)=xerr
          if (ABS(dnadn).lt.100.0d0) f(i)=x(i)*RTrho*exp(dnadn) !A is dimensionless
          end if
        enddo
      end if
c
      RETURN
      end                                              !subroutine FGCTY
c
c ======================================================================
c
      subroutine CHEMPOT (t,rho,x,u,ierr,herr)
c
c  compute the chemical potentials for each of the nc components of a
c  mixture
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c        u--array (1..nc) of the chemical potentials [J/mol]
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  03-18-10 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: CHEMPOT
c     dll_export CHEMPOT
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      dimension x(ncmax),u(ncmax)
      common /NCOMP/ nc,ic
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      character*255 herr
c
      ierr=0
      herr=' '
      do i=1,nc
        u(i)=0.d0
      enddo
      if (t.le.0.d0) return
      if (rho.lt.1.0d-40) RETURN
c
      call RMIX (x)
      RT=R*t
      call ISPURE (x,icomp)
      if (icomp.ne.0) then
c  pure component
        call AG (t,rho,x,a,g)
        u(icomp)=g
      else
        do i=1,nc
          if (x(i).gt.0.d0) then
            call PHIDERV (i,0,t,rho,x,dadn,dnadn,ierr,herr)
            da0dn=PHI0K(i,0,0,t,rho)     !ideal-gas terms  (Eq. 7.14 in GERG)
            da0dn=da0dn-href(i)/RT+sref(i)/R+1.D0+log(x(i))
            u(i)=(dnadn+da0dn)*RT
          endif
        enddo
      end if
c
      RETURN
      end                                            !subroutine CHEMPOT
c
c ======================================================================
c
      subroutine FUGCOF (t,rho,x,phi,ierr,herr)
c
c  compute the fugacity coefficient for each of the nc components of a
c  mixture
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c      phi--array (1..nc) of the fugacity coefficients
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  03-21-10 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: FUGCOF
c     dll_export FUGCOF
c
      parameter (ncmax=20)        !max number of components in mixture
      dimension x(ncmax),phi(ncmax),f(ncmax)
      common /NCOMP/ nc,ic
      character*255 herr
c
      ierr=0
      herr=' '
      do i=1,nc
        phi(i)=0.d0
      enddo
      if (t.le.0.d0) return
      if (rho.lt.1.0d-40) RETURN
c
      call ISPURE (x,icomp)
      call FGCTY2 (t,rho,x,f,ierr,herr)
      call PRESS (t,rho,x,p)
      if (p.gt.0.d0) then
        if (icomp.ne.0) then
c  pure component
          phi(icomp)=f(icomp)/p
        else
          do i=1,nc
            if (x(i).gt.0.d0) phi(i)=f(i)/p/x(i)
          enddo
        endif
      endif
c
      RETURN
      end                                             !subroutine FUGCOF
c
c ======================================================================
c
       subroutine ACTVY (t,rho,x,actv,gamma,ierr,herr)
c
c  compute the activity and activity coefficient for each of the nc
c  components of a mixture
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c     actv--array (1..nc) of the activities
c    gamma--array (1..nc) of the activity coefficients
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  03-21-10 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: ACTVY
c     dll_export ACTVY
c
      parameter (ncmax=20)        !max number of components in mixture
      dimension x(ncmax),actv(ncmax),gamma(ncmax),f(ncmax),fp(ncmax)
      dimension xl(ncmax),xv(ncmax)
      common /NCOMP/ nc,ic
      character*255 herr
c
      ierr=0
      herr=' '
      do i=1,nc
        actv(i)=0.d0      !Pure fluid values
        gamma(i)=0.d0
      enddo
      if (t.le.0.d0) return
      if (rho.lt.1.0d-40) RETURN
c
      call ISPURE (x,icomp)
      if (icomp.eq.0) then
        call FGCTY2 (t,rho,x,f,ierr,herr)
        call PRESS (t,rho,x,p)
        ic2=ic
        do i=1,nc
          ic=i
          if (x(i).gt.0.d0) then
            call RMIX (x)
c           call TPRHO(t,p,x,2,0,dp,ierr,herr)
            call TPFLSH (t,p,x,dp,dl,dv,xl,xv,q,e,h,s,cv,cp,w,ierr,herr)
            call FGCTY2 (t,dp,x,fp,ierr,herr)
            if (fp(i).gt.0.d0) actv(i)=f(i)/fp(i)
            gamma(i)=actv(i)/x(i)
          endif
        enddo
        ic=ic2
        call RMIX (x)
      endif
c
      RETURN
      end                                              !subroutine ACTVY
c
c ======================================================================
c
      subroutine VIRB (t,x,b)
c
c  compute second virial coefficient as a function of temperature
c  and composition.
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac]
c  outputs:
c        b--second virial coefficient [L/mol]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  03-27-98 EWL, original version
c  08-30-04 EWL, change rho to 0.00000001
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: VIRB
c     dll_export VIRB
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      b=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      rho=0.00000001d0
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        p=PBWR(icomp,t,rho)
        b=(p/rho/R/t-1.0d0)/rho
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi01=PHIK(icomp,0,1,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi01=PHIX(0,1,tau,del,x)
        end if
        b=phi01/rho
      end if
c
      RETURN
      end                                               !subroutine VIRB
c
c ======================================================================
c
      subroutine DBDT (t,x,dbt)
c
c  compute the 1st derivative of B (B is the second virial coefficient) with
c  respect to T as a function of temperature and composition.
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac]
c  outputs:
c      dbt--1st derivative of B with respect to T [L/mol-K]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  07-30-01 EWL, original version
c  08-30-04 EWL, change rho to 0.00000001
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DBDT
c     dll_export DBDT
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      dbt=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      rho=0.00000001d0
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        p=PBWR(icomp,t,rho)
        dpt=DPTBWR(icomp,t,rho)
        dbt=(dpt - p/t)/rho**2/t/R
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi11=PHIK(icomp,1,1,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi11=PHIX(1,1,tau,del,x)
        end if
        dbt=-phi11/rho/t
      end if
c
      RETURN
      end                                               !subroutine DBDT
c
c ======================================================================
c
      subroutine DBDT2 (t,x,dbt2)
c
c  compute the 2nd derivative of B (B is the second virial coefficient) with
c  respect to T as a function of temperature and composition.
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac]
c  outputs:
c      dbt2--2nd derivative of B with respect to T [L/mol-K^2]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  10-29-07 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DBDT2
c     dll_export DBDT2
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      dbt2=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      rho=0.00000001d0
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        dbt2=0     !not yet implemented
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi21=PHIK(icomp,2,1,tau,del)
          phi11=PHIK(icomp,1,1,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi21=PHIX(2,1,tau,del,x)
          phi11=PHIX(1,1,tau,del,x)
        end if
        dbt2=(phi21+2.d0*phi11)/rho/t**2
      end if
c
      RETURN
      end                                              !subroutine DBDT2
c
c ======================================================================
c
      subroutine VIRC (t,x,c)
c
c  compute the third virial coefficient as a function of temperature
c  and composition.
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac]
c  outputs:
c        c--third virial coefficient [(L/mol)^2]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c   3-27-98 EWL, original version
c  12-02-98 EWL, change rho to 0.0001 to avoid numerical problems in BWR calc.
c  08-30-04 EWL, change rho to 0.00000001
c  05-24-06 EWL, change rho to 0.000001 for the BWR
c  05-07-09 EWL, change rho to 0.0001 for the BWR
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: VIRC
c     dll_export VIRC
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      c=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      rho=0.000001d0
      if (heos.eq.'BWR') then
        rho=0.0001d0
c  pure fluid MBWR equation of state--call BWR-specific routines
        icomp=1
        p=PBWR(icomp,t,rho)
        dpd=DPDBWR(icomp,t,rho)
        c=((dpd-2.0d0*p/rho)/R/t+1.0d0)/rho**2
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi02=PHIK(icomp,0,2,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi02=PHIX(0,2,tau,del,x)
        end if
        c=phi02/rho**2
      end if
c
      RETURN
      end                                               !subroutine VIRC
c
c ======================================================================
c
      subroutine DCDT (t,x,dct)
c
c  compute the 1st derivative of C (C is the third virial coefficient) with
c  respect to T as a function of temperature and composition.
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac]
c  outputs:
c      dct--1st derivative of C with respect to T [(L/mol)^2-K]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  10-29-07 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DCDT
c     dll_export DCDT
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      dct=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      rho=0.00000001d0
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        dct=0.d0      !not yet implemented
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi12=PHIK(icomp,1,2,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi12=PHIX(1,2,tau,del,x)
        end if
        dct=-phi12/rho**2/t
      end if
c
      RETURN
      end                                               !subroutine DCDT
c
c ======================================================================
c
      subroutine DCDT2 (t,x,dct2)
c
c  compute the 2nd derivative of C (C is the third virial coefficient) with
c  respect to T as a function of temperature and composition.
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac]
c  outputs:
c      dct2--2nd derivative of C with respect to T [(L/mol-K)^2]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  10-29-07 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DCDT2
c     dll_export DCDT2
c
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      dct2=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      rho=0.00000001d0
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
        dct2=0     !not yet implemented
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi22=PHIK(icomp,2,2,tau,del)
          phi12=PHIK(icomp,1,2,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi22=PHIX(2,2,tau,del,x)
          phi12=PHIX(1,2,tau,del,x)
        end if
        dct2=(phi22+2.d0*phi12)/rho**2/t**2
      end if
c
      RETURN
      end                                              !subroutine DCDT2
c
c ======================================================================
c
      subroutine VIRD (t,x,d)
c
c  compute the fourth virial coefficient as a function of temperature
c  and composition.
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac]
c  outputs:
c        d--fourth virial coefficient [(L/mol)^3]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  11-26-01 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
      character*3 hpheq,heos,hmxeos,hmodcp
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      d=0.d0
      if (t.le.0.d0) return
      call RMIX (x)
      rho=0.00000001d0
      if (heos.eq.'BWR') then
c  pure fluid MBWR equation of state--call BWR-specific routines
c       icomp=1
c       p=PBWR(icomp,t,rho)
c       dpd=DPDBWR(icomp,t,rho)
c  need to update with correct formula:
c       c=((dpd-2.0d0*p/rho)/R/t+1.0d0)/rho**2
        d=0
c
      else
c  call general PHIK or PHIX routines for all other models
        call ISPURE (x,icomp)
        if (icomp.ne.0) then
c  pure fluid
          tau=tz(icomp)/t
          del=rho/rhoz(icomp)
          phi03=PHIK(icomp,0,3,tau,del)
        else
c  mixture
          call REDX (x,t0,rho0)
          tau=t0/t
          del=rho/rho0
          phi03=PHIX(0,3,tau,del,x)
        end if
        d=phi03/rho**3
      end if
c
      RETURN
      end                                               !subroutine VIRD
c
c ======================================================================
c
      subroutine VIRBA (t,x,ba)
c
c  compute second acoustic virial coefficient as a function of temperature
c  and composition.
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac]
c  outputs:
c        ba--second acoustic virial coefficient [L/mol]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  10-29-07 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: VIRBA
c     dll_export VIRBA
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      ba=0.d0
      if (t.le.0.d0) return
      call VIRB (t,x,b)
      call DBDT (t,x,dbt)
      call DBDT2 (t,x,dbt2)
      cp00=CP0(t,x)
      gpg=cp00/(cp00-R)
c  Trusler and Zarari, J. Chem. Thermodyn., 28:329-335, 1996.
c  Gillis and Moldover, Int. J. Theromphys., 17(6):1305-1324, 1996.
      ba=2.d0*b+2.d0*(gpg-1.d0)*t*dbt+(gpg-1.d0)**2/gpg*t**2*dbt2
c
      RETURN
      end                                              !subroutine VIRBA
c
c ======================================================================
c
      subroutine VIRCA (t,x,ca)
c
c  compute third acoustic virial coefficient as a function of temperature
c  and composition.
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac]
c  outputs:
c        ca--third acoustic virial coefficient [(L/mol)^2]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  10-29-07 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: VIRCA
c     dll_export VIRCA
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      ca=0.d0
      if (t.le.0.d0) return
      call VIRB (t,x,b)
      call VIRBA (t,x,ba)
      call DBDT (t,x,dbt)
      call DBDT2 (t,x,dbt2)
      call VIRC (t,x,c)
      call DCDT (t,x,dct)
      call DCDT2 (t,x,dct2)
      cp00=CP0(t,x)
      gpg=cp00/(cp00-R)
c  Gillis and Moldover, Int. J. Theromphys., 17(6):1305, 1996.
c  Estela-Uribe and Trusler, Int.  J. Theromphys., 21(5):1033, 2000.
      q=b+(2.d0*gpg-1.d0)*t*dbt+(gpg-1.d0)*t**2*dbt2
      ca=(gpg-1.d0)*q**2+(2.d0*gpg+1.d0)*c
      ca=ca+(gpg**2-1.d0)*t*dct+(gpg-1.d0)**2/2.d0*t**2*dct2
      ca=ca/gpg
c  to convert to the pressure expansion form, use this:
c     ca=(el-ba*b)/R/t
c
      RETURN
      end                                              !subroutine VIRCA
c
c ======================================================================
c
      subroutine B12 (t,x,b)
c
c  compute b12 as a function of temperature and composition.
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac]
c  outputs:
c        b--b12 [(L/mol)^2]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  04-19-01 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      if (nc.ne.2) then   !Only calculate b12 for a binary
        b=0
        RETURN
      endif
      call VIRB (t,x,bx)
      ic2=ic
      ic=1
      call VIRB (t,x,b1)
      ic=2
      call VIRB (t,x,b2)
      b=(bx-x(1)**2*b1-x(2)**2*b2)/2.d0/x(1)/x(2)
      ic=ic2
      RETURN
      end                                                !subroutine B12
c
c ======================================================================
c
      subroutine EXCESS (t,p,x,kph,rho,vE,eE,hE,sE,aE,gE,ierr,herr)
c
c  compute excess properties as a function of temperature, pressure,
c  and composition.
c
c  inputs:
c        t--temperature [K]
c        p--pressure [kPa]
c        x--composition [array of mol frac]
c      kph--phase flag:  1 = liquid
c                        2 = vapor
c                        0 = stable phase
c  outputs:
c       rho--molar density [mol/L] (if input less than 0, used as initial guess)
c        vE--excess volume [L/mol]
c        eE--excess energy [J/mol]
c        hE--excess enthalpy [J/mol]
c        sE--excess entropy [J/mol-K]
c        aE--excess Helmholtz energy [J/mol]
c        gE--excess Gibbs energy [J/mol]
c      ierr--error flag:  0 = successful
c                        55 = T,p inputs in different phase for the pure fluids
c      herr--error string (character*255 variable if ierr<>0)
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  04-25-02 EWL, original version
c  11-04-08 EWL, add ierr and herr to argument list
c  11-26-08 EWL, add aE and gE to argument list
c  03-21-10 EWL, add log(x(i)) to sE, gE, and aE
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      common /NCOMP/ nc,ic
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /HCHAR/ htab,hnull
      dimension x(ncmax),xliq(ncmax),xvap(ncmax)
      character*255 herr
c
      ierr=0
      herr=' '
      d=rho
      vE=0.d0
      eE=0.d0
      hE=0.d0
      sE=0.d0
      aE=0.d0
      gE=0.d0
      rho=0.d0
      call ISPURE (x,icomp)
      if (icomp.ne.0) RETURN
c
      kguess=0
      if (d.lt.0) kguess=1
      d=abs(d)
      if (kph.ne.0) then
       call TPRHO (t,p,x,kph,kguess,d,ierr,herr)
       if (ierr.ne.0) then
         d=d*2
         call TPRHO (t,p,x,kph,1,d,ierr,herr)
       endif
       if (ierr.ne.0)
     & call TPFLSH(t,p,x,d,dl,dv,xliq,xvap,q,eE,hE,sE,cv,cp,w,ierr,herr)
       call THERM (t,d,x,pp,eE,hE,sE,cv,cp,w,hjt)
       call AG (t,d,x,aE,gE)
      else
       call TPFLSH(t,p,x,d,dl,dv,xliq,xvap,q,eE,hE,sE,cv,cp,w,ierr,herr)
       call AG (t,d,x,aE,gE)
      endif
      if (ierr.ne.0) RETURN
      rho=d
      if (d.gt.0.d0) vE=1.d0/d
c
      call CRITP (x,tcrit,pcrit,Dcrit,ierr,herr)
      ic2=ic
      do i=1,nc
        ic=i
c...Do not use kph or call TPRHO since the pure fluids could be in either phase
        call RMIX (x)
        RT=R*t
        call INFO (i,wm,ttp,tnbp,tc,pc,dc,zc,acf,dip,rgas)
        call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
        if (d.gt.dc .and. rho.lt.dcrit .and. t.lt.tc) goto 100
        if (d.lt.dc .and. rho.gt.dcrit .and. t.lt.tc) goto 100
        call AG (t,d,x,a,g)
        if (d.gt.0.d0) vE=vE-x(i)/d
        eE=eE-x(i)*e
        hE=hE-x(i)*h
        sE=sE-x(i)*(s- R*log(x(i)))
        aE=aE-x(i)*(a+RT*log(x(i)))
        gE=gE-x(i)*(g+RT*log(x(i)))
      enddo
      ic=ic2
      call RMIX (x)
      RETURN
c
 100  continue
      vE=0.d0
      eE=0.d0
      hE=0.d0
      sE=0.d0
      aE=0.d0
      gE=0.d0
      ierr=55
      herr='[EXCESS error] temperature and pressure inputs are in '//
     &     'different phases for the pure fluids'//hnull
      call ERRMSG (ierr,herr)
      ic=ic2
      call RMIX (x)
      RETURN
      end                                             !subroutine EXCESS
c
c ======================================================================
c
      subroutine FPV (t,rho,p,x,f)
c
c  Compute the supercompressibility factor, Fpv.
c
c  inputs:
c        t--temperature [K]
c        p--pressure [kPa]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c        f--Fpv = sqrt[Z(60 F, 14.73 psia)/Z(T,P)]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  11-07-02 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: FPV
c     dll_export FPV
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
      character*255 herr
c
      tfpv = 288.705555555556d0    !60 F
      pfpv = 101.55977492837d0     !14.73 psia
      call TPRHO (tfpv,pfpv,x,2,0,dfpv,ierr,herr)
      if (p.gt.0.d0) then
        f=SQRT(pfpv/tfpv/dfpv*rho*t/p)
      else
        f=SQRT(pfpv/tfpv/dfpv/R)
      endif

      RETURN
      end                                                !subroutine FPV
Cc
Cc ======================================================================
Cc
C      subroutine SPECGR (t,rho,p,gr)
Cc
Cc  Compute the specific gravity (relative density).
Cc
Cc  inputs:
Cc        t--temperature [K]
Cc        p--pressure [kPa]
Cc      rho--molar density [mol/L]
Cc  outputs:
Cc       gr--specific gravity [dimensionless]
Cc
Cc  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
Cc  11-07-02 EWL, original version
Cc
C      implicit double precision (a-h,o-z)
C      implicit integer (i-k,m,n)
C      parameter (ncmax=20)        !max number of components in mixture
C      parameter (nrefmx=10)       !max number of fluids for transport ECS
C      parameter (n0=-ncmax-nrefmx,nx=ncmax)
C      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
Cc
C      rhoair=1.d0      !Need to add formulation for air here
C      if (rhoair.gt.0.d0) then
C        gr=rho/rhoair
C      else
C        gr=1.d0
C      endif
C
C      RETURN
C      end                                             !subroutine SPECGR
c
c ======================================================================
c
      subroutine RMIX (x)
c
c  inputs:
c        x--composition [array of mol frac]
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  01-19-01 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      call ISPURE (x,icomp)
      if (icomp.ne.0) then
        R=Reos(icomp)
      else
        R=0.0d0
        do i=1,nc
          R=R+x(i)*Reos(i)
        enddo
      endif
      if (R.lt.1.d-10) R=8.314472d0  !Check for bad x(i)
      RETURN
      end                                               !subroutine RMIX
c
c ======================================================================
c
      subroutine RMIX2 (x,Rgas)
c
c  Return the gas "constant" as a combination of the gas constants for
c  the pure fluids
c
c  inputs:
c        x--composition [array of mol frac]
c  outputs:
c     Rgas--gas constant [J/mol-K]
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  03-21-10 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /NCOMP/ nc,ic
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      dimension x(ncmax)
c
      call RMIX (x)
      Rgas=R
      RETURN
      end                                              !subroutine RMIX2
c
c ======================================================================
c
      subroutine THERM3 (t,rho,x,
     &           xkappa,beta,xisenk,xkt,betas,bs,xkkt,thrott,pi,spht)
c
c  Compute miscellaneous thermodynamic properties
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c
c  outputs:
c   xkappa--Isothermal compressibility [1/kPa]
c     beta--Volume expansivity [1/K]
c   xisenk--Isentropic expansion coefficient [-]
c      xkt--Isothermal expansion coefficient [-]
c    betas--Adiabatic compressibility [1/kPa]
c       bs--Adiabatic bulk modulus [kPa]
c     xkkt--Isothermal bulk modulus [kPa]
c   thrott--Isothermal throttling coefficient [L/mol]
c       pi--Internal pressure [kPa]
c     spht--Specific heat input [J/mol]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  06-16-06 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: THERM3
c     dll_export THERM3
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      dimension x(ncmax)
c
      call THERM2 (t,rho,x,p,e,h,s,cv,cp,w,Z,hjt,A,G,xkappa,beta,
     &             dPdrho,d2PdD2,dPT,drhodT,drhodP,
     &             d2PT2,d2PdTD,spare3,spare4)
      wm=WMOL(x)
      xisenk=0.d0
      if (p.le.0.d0) then
        xkt=1
        if (t.gt.0.d0) xisenk=w**2/R/T*wm*0.001d0
      else
        xkt=rho/p*dPdrho               !Isothermal expansion coefficient
        xisenk=w**2*rho/p*wm*0.001d0   !Isentropic expansion coefficient
      endif
      betas=xnotc
      if (rho.gt.0.d0 .and. w.gt.0.d0)
     &  betas=1.d0/rho/w**2/wm*1000.d0 !Adiabatic compressibility
      bs=xisenk*p                      !Adiabatic bulk modulus
      xkkt=xkt*p                       !Isothermal bulk modulus
      thrott=-hjt*cp                   !Isothermal throttling coef.
      pi=t*dpt-p                       !Internal pressure
      if (abs(dpt).gt.1.d-20) then
        spht=rho*cp*dpdrho/dpt         !Specific heat input
      else
        spht=cp*t
      end If
      RETURN
      end                                             !subroutine THERM3
c
c ======================================================================
c
      subroutine VIRBCD (icomp,idel,tau,vir)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (mxtrm=72,mxcrt=20)
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
      common /WCFFQ2/ alpha(n0:nx,mxcrt),beta(n0:nx,mxcrt),
     &                gamma(n0:nx,mxcrt),delta(n0:nx,mxcrt),
     &                eta(n0:nx,mxcrt),eid(n0:nx,mxcrt),eit(n0:nx,mxcrt)
      common /CRTSAV2/ delb(n0:nx,mxcrt),taua(n0:nx,mxcrt),
     &                txp(n0:nx,mxcrt),hxp(n0:nx,mxcrt),
     &                ext(n0:nx,mxcrt),extd(n0:nx,mxcrt),
     &                extt(n0:nx,mxcrt),extdt(n0:nx,mxcrt),
     &                extt2(n0:nx,mxcrt),extd2(n0:nx,mxcrt)
      dimension       phisav(n0:nx,mxtrm)
c     dimension taup(n0:nx,mxtrm)
      if (tau.le.0.0d0) RETURN    !  any and all derivatives
      elntau=log(tau)
C     do j=1,ntp(icomp)
C       taup(icomp,j)=tpower(icomp,j)*elntau
C     enddo
C     do k=1,ntermf(icomp)
C       phisav(icomp,k)=a(icomp,k)*EXP(taup(icomp,itp(icomp,k)))
C     enddo
      do k=1,ntermf(icomp)
        phisav(icomp,k)=a(icomp,k)*EXP(ti(icomp,k)*elntau)
      enddo
c
      phisum=0.0d0
      if (idel.eq.1) then
        do k=1,ntermf(icomp)
          if (di(icomp,k).eq.1) then
            phisum=phisum+phisav(icomp,k)
          endif
        enddo
        phisum=phisum/rhoz(icomp)
      elseif (idel.eq.2) then
        do k=1,ntermf(icomp)
          if (di(icomp,k).eq.2) then
            phisum=phisum+2.d0*phisav(icomp,k)
          elseif (di(icomp,k).eq.1 .and. dli(icomp,k).eq.1) then
            phisum=phisum-2.d0*phisav(icomp,k)
          endif
        enddo
        phisum=phisum/rhoz(icomp)**2
      elseif (idel.eq.3) then
        do k=1,ntermf(icomp)
          if (di(icomp,k).eq.3) then
            phisum=phisum+6.d0*phisav(icomp,k)
          elseif (di(icomp,k).eq.1 .and. dli(icomp,k).eq.2) then
            phisum=phisum-6.d0*phisav(icomp,k)
          elseif (di(icomp,k).eq.2 .and. dli(icomp,k).eq.1) then
            phisum=phisum-6.d0*phisav(icomp,k)
          elseif (di(icomp,k).eq.1 .and. dli(icomp,k).eq.1) then
            phisum=phisum+3.d0*phisav(icomp,k)
          endif
        enddo
        phisum=phisum/rhoz(icomp)**3
      endif
      vir=phisum
      RETURN
      end                                             !subroutine VIRBCD
c
c ======================================================================
c
      subroutine HEAT (t,rho,x,hg,hn,ierr,herr)
c
c  Compute the ideal gas gross and net heating values.
c
c  inputs:
c        t--temperature [K]
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c
c  outputs:
c       hg--gross (or superior) heating value [J/mol]
c       hn--net (or inferior) heating value [J/mol]
c     ierr--error flag:  0 = successful
c                        1 = error in chemical formula
c                        2 = not all heating values available
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  01-09-08 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: HEAT
c     dll_export HEAT
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,hnam80,hsyn1,hsyn2,hcf
      dimension x(ncmax),v(6)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /CNAM80/ hnam80(n0:nx),hsyn1(n0:nx),hsyn2(n0:nx)
      common /HTCM/ hcmbst(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c
      call ISPURE (x,icomp)
      hg=0.d0
      hn=0.d0
      do k=1,nc
        acrb=0.d0   !Number of carbon atoms
        ahyd=0.d0
        aoxy=0.d0
        anit=0.d0
        asul=0.d0
        hgk=0.d0
        hnk=0.d0
c
c  extract chemical formula from hsyn1, use the part in {} if available
        hcf=hsyn1(k)
        i=index(hcf,'{')
        if (i.eq.0) then
          i=index(hcf,'!')
          if (i.ne.0) hcf=hcf(1:i-1)
        else
          hcf=hcf(i+1:255)
          i=index(hcf,'}')
          if (i.ne.0) hcf=hcf(1:i-1)
        endif
c
c  extract the number of carbon, hydrogen, oxygen, ..., atoms from
c  the chemical formula
 10     continue
        if (hcf(2:2).lt.'a' .or. hcf(2:2).gt.'z')
     &      hcf=hcf(1:1)//' '//hcf(2:255)       !Add space in second slot
        do i=3,5
          if (hcf(i:i).lt.'0' .or. hcf(i:i).gt.'9')
     &        hcf=hcf(1:i-1)//' '//hcf(i:255)       !Add space in ith slot
        enddo
c
        read (hcf(3:5),'(i3)') j
        if (j.eq.0) j=1
        if (hcf(1:2).eq.'C ') then           !Carbon
          if (acrb.gt.0) goto 999
          acrb=real(j)
        elseif (hcf(1:2).eq.'H ') then       !Hydrogen
          if (ahyd.gt.0) goto 999
          ahyd=real(j)
        elseif (hcf(1:2).eq.'O ') then       !Oxygen
          if (aoxy.gt.0) goto 999
          aoxy=real(j)
        elseif (hcf(1:2).eq.'N ') then       !Nitrogen
          if (anit.gt.0) goto 999
          anit=real(j)
        elseif (hcf(1:2).eq.'S ') then       !Sulfur
          if (asul.gt.0) goto 999
          asul=real(j)
        endif
        hcf=hcf(6:255)
        if (hcf.ne.' ') goto 10
c
c       if (acrb+ahyd+aoxy+anit+asul.eq.0) goto 999   !Unknown substance
c
        if (abs(hcmbst(k)+1.d0).lt.1.d-12) goto 998
        hgk=hcmbst(k)*1000.d0
        v(1)=-ahyd/2.d0                        !Water produced
        v(2)=acrb+ahyd/4.d0-aoxy/2.d0+asul     !Oxygen needed
        v(3)=-acrb                             !CO2 produced
        v(4)=-anit/2.d0                        !Nitrogen
        v(5)=-asul                             !SO2
c
        t25=298.15d0
        rho0=0.d0
        call ENTHAL (t25,rho0,x,h25)
        call ENTHAL (t,rho0,x,h)
        hgk=hgk-(h25-h)
c
        call ENTHHC (0,t25,t,h)                !Liquid water
        hgk=hgk-v(1)*h
        do i=2,5
          call ENTHHC (i,t25,t,h)
          hgk=hgk-v(i)*h
        enddo
c
        call ENTHHC (0,t,-1.d0,h1)                !Liquid water
        call ENTHHC (1,t,-1.d0,h2)                !Ideal gas water
        hnk=hgk+v(1)*(h2-h1)
c
        if (k.eq.icomp) then
          hg=hgk
          hn=hnk
        else
          hg=hg+x(k)*hgk
          hn=hn+x(k)*hnk
        endif
      enddo
c
      ierr=0
      herr=' '
      RETURN
c
 998  ierr=2
      herr='[HEAT error 2] Heating values are not available for all '//
     &     'species in the mixture'//hnull
      hg=0.d0
      hn=0.d0
      RETURN
c
 999  ierr=1
      herr='[HEAT error 1] Error in chemical formula'//hnull
      hg=0.d0
      hn=0.d0
      RETURN
      end                                               !subroutine HEAT
c
c ======================================================================
c
      subroutine ENTHHC (icmb,t1,t2,h)
c
c  Compute the ideal gas enthalpy difference between temperatures t1 and
c  t2 for several combustion gases.  This is used in conjunction with the
c  heat of combustion subroutine so that the extra fluids do not have
c  loaded in memory.  If the EOS ever changes, then the coefficients
c  here must be updated.
c
c  inputs:
c       t1--temperature [K]
c       t2--temperature [K]
c           if t2 is less than zero, then absolute enthalpy calculated at t1
c     icmb--fluid identifier:
c           0-saturated liquid water
c           1-water
c           2-oxygen
c           3-CO2
c           4-nitrogen
c           5-SO2
c
c  outputs:
c        h--ideal gas enthalpy difference [J/mol]
c           (for i=0, then saturated liquid enthalpy difference of water)
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  01-09-08 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      dimension xkhc(10),cpchc(10)
c
      h=0
      if (icmb.eq.0) then
        ntc=6
        cpchc(1)= 0.8943272d+5     !Enthalpy of saturated liquid water.
        cpchc(2)=-0.6138239d+5     !Fit by EWL on 1/9/2007
        cpchc(3)=-0.4415242d+5     !Matches Wagner & Pruss equation to within
        cpchc(4)=-0.1092092d+5     !0.02% over full saturation range.
        cpchc(5)= 0.6703065d+5     !Extrapolation to 200 K is smooth.
        cpchc(6)=-0.1763063d+6
        xkhc(1) = 0.032d0
        xkhc(2) = 0.078d0
        xkhc(3) = 0.825d0
        xkhc(4) = 4.d0
        xkhc(5) = 9.d0
        xkhc(6) =12.d0
        do i=1,ntc
                          h=h+cpchc(i)*(1.d0-t1/647.096)**xkhc(i)
          if (t2.ge.0.d0) h=h-cpchc(i)*(1.d0-t2/647.096)**xkhc(i)
        enddo
        RETURN
      elseif (icmb.eq.1) then
        ntc=1
        nte=5
        cpchc(1)=0.400632d+1             !Water
        cpchc(2)=0.124360d-1
        cpchc(3)=0.973150d+0
        cpchc(4)=0.127950d+1
        cpchc(5)=0.969560d+0
        cpchc(6)=0.248730d+0
        xkhc(1) =    0.d0
        xkhc(2) =  833.d0
        xkhc(3) = 2289.d0
        xkhc(4) = 5009.d0
        xkhc(5) = 5982.d0
        xkhc(6) =17800.d0
        R=8.314371357587d0
c                           value from above    value from refprop
c                           at 300 K            at 300 K
        if (t2.lt.0.d0) h=(-88471.35339670102d0+45964.71449803960d0)/r
      elseif (icmb.eq.2) then
        ntc=1
        nte=5
        cpchc(1)=3.51808732d0            !Oxygen
        cpchc(2)=0.102323928D+01
        cpchc(3)=0.784357918D+00
        cpchc(4)=0.337183363D-02
        cpchc(5)=-.170864084D-01
        cpchc(6)=0.463751562D-01
        xkhc(1) =0.d0
        xkhc(2) =0.224632440D+04
        xkhc(3) =0.112599763D+05
        xkhc(4) =0.120126209D+04
        xkhc(5) =0.690089445D+02
        xkhc(6) =0.532805445D+04
        R=8.31434d0
        if (t2.lt.0.d0) h=(-56058.60590328777d0+8734.35384436554d0)/r
      elseif (icmb.eq.3) then
        ntc=1
        nte=5
        cpchc(1)=0.35d+01                !CO2
        cpchc(2)=1.99427042d0
        cpchc(3)=0.621052475d0
        cpchc(4)=0.411952928d0
        cpchc(5)=1.04028922d0
        cpchc(6)=0.0832767753d0
        xkhc(1) =   0.d0
        xkhc(2) = 958.49956d0
        xkhc(3) =1858.80115d0
        xkhc(4) =2061.10114d0
        xkhc(5) =3443.89908d0
        xkhc(6) =8238.20035d0
        R=8.31451d0
        if (t2.lt.0.d0) h=(-43458.10456571817d0+22372.0720622156d0)/r
      elseif (icmb.eq.4) then
        ntc=4
        nte=1
        cpchc(1)= 3.5d0                  !Nitrogen
        cpchc(2)= 3.066469d-6
        cpchc(3)= 4.70124d-9
        cpchc(4)=-3.987984d-13
        cpchc(5)= 0.1012941d1
        xkhc(1) = 0.d0
        xkhc(2) = 1.d0
        xkhc(3) = 2.d0
        xkhc(4) = 3.d0
        xkhc(5) = 3364.011d0
        R=8.31451d0
        if (t2.lt.0.d0) h=(-22898.14229042497d0+8723.88255738888d0)/r
      elseif (icmb.eq.5) then
        ntc=2
        nte=2
        cpchc(1)= 4.0d0                  !SO2
        cpchc(2)= 0.72453d-4
        cpchc(3)= 1.0620d0
        cpchc(4)= 1.9401d0
        xkhc(1) =    0.d0
        xkhc(2) =    1.d0
        xkhc(3) =  775.d0
        xkhc(4) = 1851.d0
        R=8.314472d0
        if (t2.lt.0.d0) h=(-28976.84624561172d0+26659.0278471527d0)/r
      endif
c
      do i=1,ntc
        xkhci=xkhc(i)
        xkhc1=xkhci+1.0d0
                        h=h+cpchc(i)*t1**xkhc1/xkhc1
        if (t2.ge.0.d0) h=h-cpchc(i)*t2**xkhc1/xkhc1
      enddo
      do i=1,nte
        j=i+ntc
        expui1=EXP(xkhc(j)/t1)
        expui2=EXP(xkhc(j)/t2)
        if (t2.ge.0.d0)
     &  h=h-cpchc(j)*(-0.5d0*xkhc(j))*(1.0d0+expui2)/(1.0d0-expui2)
        h=h+cpchc(j)*(-0.5d0*xkhc(j))*(1.0d0+expui1)/(1.0d0-expui1)
      enddo
      h=h*R
      end                                             !subroutine ENTHHC
c
c ======================================================================
c
      subroutine ISPURE (x,icomp)
c
c  Determine if the user has requested the properties of a pure fluid.
c  This happens if 1) nc=1, 2) PUREFLD has been called, or 3) one of the
c  compositions in the x array is one.
c
c  inputs:
c        x--composition [array of mol frac]
c
c  outputs:
c    icomp--index set to the pure fluid, which is 1 for nc=1 or 0 for
c           a mixture.  If PUREFLD has been called, icomp=ic.
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  05-01-08 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)        !max number of components in mixture
      common /NCOMP/ nc,ic
      dimension x(ncmax)
c
      icomp=0
      if (nc.eq.1) then
        icomp=1
      elseif (ic.gt.0) then
        icomp=ic
      else
        do i=1,nc
          if (ABS(x(i)-1.d0).lt.1.d-12) then
            icomp=i
            RETURN
          endif
        enddo
      endif
      RETURN
      end                                             !subroutine ISPURE
c
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c ======================================================================
c                                                    end file prop_sub.f
c ======================================================================
