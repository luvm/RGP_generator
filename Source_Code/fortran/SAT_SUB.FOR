c  begin file sat_sub.f
c
c  This file contains routines for saturation properties
c
c  contained here are:
c     subroutine SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr)
c     subroutine SATP (p,x,kph,t,rhol,rhov,xliq,xvap,ierr,herr)
c     subroutine SATD (rho,x,kph,kr,t,p,rhol,rhov,xliq,xvap,ierr,herr)
c     subroutine SATH (h,x,kph,nroot,k1,t1,p1,d1,k2,t2,p2,d2,ierr,herr)
c     subroutine SATE (e,x,kph,nroot,k1,t1,p1,d1,k2,t2,p2,d2,ierr,herr)
c     subroutine SATS (s,x,kph,nroot,k1,t1,p1,d1,k2,t2,p2,d2,
c    &                 k3,t3,p3,d3,ierr,herr)
c     subroutine CSATK (icomp,t,kph,p,rho,csat,ierr,herr)
c     subroutine CV2PK (icomp,t,rho,cv2p,csat,ierr,herr)
c     subroutine DPTSATK (icomp,t,kph,p,rho,csat,dpt,ierr,herr)
c     subroutine TPRHOB (t,p,rho1,rho2,x,rho,ierr,herr)
c     subroutine DLDV (t,p,rhol,rhov,xl,xv,ierr,herr)
c     subroutine SPNDL (t,x,rhol,rhov,ierr,herr)
c     subroutine SATTP (t,p,x,iFlash,iGuess,d,Dl,Dv,xliq,xvap,q,ierr,herr)
c     subroutine SATTEST (t,x,kph,p,x2,ierr,herr)
c     subroutine SATPEST (p,x,kph,t,x2,ierr,herr)
c
c  these routines use the following common blocks from other files
c     common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c     common /NCOMP/ nc,ic
c     common /HCHAR/ htab,hnull
c     common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
c                   ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
c                   tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
c                   wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
c  various arrays are dimensioned with parameter statements
c     parameter (ncmax=20)        !max number of components in mixture
c     parameter (nrefmx=10)       !max number of fluids for transport ECS
c     parameter (n0=-ncmax-nrefmx,nx=ncmax)
c
c ======================================================================
c ======================================================================
c
      subroutine SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr)
c
c  iterate for saturated liquid and vapor states given temperature
c  and the composition of one phase
c
c  inputs:
c        t--temperature [K]
c                       if t is negative, then use other variables as
c                       initial guesses at abs(t)
c        x--composition [array of mol frac] (phase specified by kph)
c      kph--phase flag: 1 = input x is liquid composition (bubble point)
c                       2 = input x is vapor composition (dew point)
c                       3 = input x is liquid composition (freezing point)
c                       4 = input x is vapor composition (sublimation point)
c  outputs:
c        p--pressure [kPa]
c     rhol--molar density [mol/L] of saturated liquid
c     rhov--molar density [mol/L] of saturated vapor
c           For a pseudo pure fluid, the density of the equilibrium phase
c           is not returned.  Call SATT twice, once with kph=1 to get
c           pliq and rhol, and once with kph=2 to get pvap and rhov.
c     xliq--liquid phase composition [array of mol frac]
c     xvap--vapor phase composition [array of mol frac]
c     ierr--error flag:   0 = successful
c                         1 = T < Tmin
c                         8 = x out of range
c                         9 = T and x out of range
c                       120 = CRITP did not converge
c                       121 = T > Tcrit
c                       122 = TPRHO-liquid did not converge (pure fluid)
c                       123 = TPRHO-vapor did not converge (pure fluid)
c                       124 = pure fluid iteration did not converge
c           following 3 error codes are advisory--iteration will either
c           converge on later guess or error out (ierr = 128)
c                      -125 = TPRHO did not converge for parent ph (mix)
c                      -126 = TPRHO did not converge for incipient (mix)
c                      -127 = composition iteration did not converge
c                       128 = mixture iteration did not converge
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  07-11-95  MM, original version
c  09-11-95  MM, add error string to argument list
c  09-25-95  MM, rearrange argument list (outputs in order p, rho, x)
c  10-06-95  MM, use stored acentric factor for pure fluids
c  10-11-95  MM, RETURN if any error detected
c  11-26-95  MM, Raoult's law as first guess for mixture
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  12-18-95  MM, fill xliq, xvap with zeros for undefined components
c  12-19-20  MM, add full mixture iteration using fugacity
c  12-27-95  MM, pratio for new pressure if no converge for TPRHO for mix
c  01-09-96  MM, move check for supercritical outside nc = 1 block
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                add Zcrit to common /CCON/
c  03-06-96  MM, reset warning from TPRHO if loop eventually converges
c  03-19-96  MM, add dipole moment to /CCON/
c  04-05-96  MM, test for supercritical '.ge. tc' rather than '.gt. tc'
c  05-30-96  MM, check input temperature against limits
c  06-03-96  MM, add 'EOS' to calling list for LIMITX
c  06-05-96  MM, refine error numbers and messages;
c                also ensure that all outputs are set on error condition
c  11-14-96  MM, adjust pratio, etc to get closer to critical
c   2-07-96 EWL, add pressure increment/decrement when TPRHO does not converge
c   2-13-96 EWL, add initial guess for densities near the critical point
c                return critical point values if within delta of Tc
c   6-06-96 EWL, return critical point values for failure to converge if
c                Tc - T < 10 mK
c  10-01-97  MM, add compiler switch to allow access by DLL
c  11-13-97 EWL, initialize fpit(j); potential bomb if no value when writing error message
c  11-14-97 EWL, add line following 500 to improve critical region convergence
c  02-09-98  MM, limit delp step, change delp if TPRHO does not converge
c  02-11-98  MM, check that new guess for mix pressure is < p_crit
c  03-24-98 EWL, use critical region initial guesses when converge fails
c                due to liquid and vapor roots being equal
c  04-06-98 EWL, call subroutine AG to get Gibbs energy
c  07-30-98  MM, check for f1=0 in phase 2 iteration; separate kguess for each phase
c  08-03-98 EWL, check for fpit(2)-fpit(1)<>0 before calculation of pit(3)
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-22-98 EWL, zero out all elements of xliq and xvap
c  12-22-98 EWL, use pure fluid algorithm when called with a mixture but x(i)=1
c  02-04-99 EWL, increase tolerance if number of iterations hits 8 or 12
c  02-04-99 EWL, if fugacity is zero, don't restart, just change x2
c  02-04-99 EWL, add check for equal densities
c  02-04-99 EWL, allow pressure to be greater than critical pressure
c  02-04-99 EWL, do not allow outer loop to converge if inner loop has not
c  02-04-99 EWL, if TPRHO fails, set x2 to original values before restarting
c  08-23-99 EWL, change khpsav to kphsav
c  11-16-99 EWL, do not modify x2new(i) after the check for fugacity=0
c                when x2(i)=x1(i)=0 (composition was set to zero by user)
c  02-23-00 EWL, add kph=3,4 as inputs and call melting or sublimation lines
c  03-14-00 MLH, change 'do 100 i=1,ncmax' from ncmax to nc
c  12-18-00 EWL, use different pratio for dew and bubble sides
c  02-27-01 EWL, don't allow too large of jumps on first mixture iteration for p
c  07-16-01 EWL, add calls to ancillary equations to get better estimates
c  07-16-01 EWL, call THERM one last time on the liquid side once converged
c  11-08-01 EWL, add alternative method that converges near the critical point
c  11-20-01 EWL, check for ierr=-16 (t<ttrp)
c  02-25-02 EWL, allow initial guesses to be passed in by negating the pressure
c  05-28-02 EWL, check for bad root in two phase
c  09-19-02 EWL, add ancillary check for liquid pressure
c  09-19-02 EWL, exit after calculating ancillaries for pseudo-pure fluids
c  07-28-03 EWL, add check to keep x2 from bouncing around
c  11-30-04 EWL, add checks to remove crashes (p>1.d5, rhol-rhov<.1)
c  11-16-05 EWL, add check for negative temperature for use in getting initial values
c  07-21-08 EWL, add missing parenthesis in check for delp-delp2
c  02-18-08 EWL, add check for icomp<>old icomp (icsav)
c  04-12-10 EWL, add check for dp/dT negative
c  04-22-10 EWL, rename variables in TSTSAV
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SATT
c     dll_export SATT
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herr3
      character*3 hps,hpsk,hpl,hplk,hdl,hdlk,hdv,hdvk
      character*3 hpheq,heos,hmxeos,hmodcp
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),x2org(ncmax),xs(ncmax)
      dimension x2(ncmax),f1(ncmax),f2(ncmax),x2new(ncmax)
      dimension pit(3),fpit(2)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /prnterr/ iprnterr
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /FSHSAV/ lsatt,lsatp
      common /TSTSAV/ tsavt,psavt,dlsavt,dvsavt,
     &                xsavt(ncmax),xlsavt(ncmax),xvsavt(ncmax),
     &                kphsvt,icsavt
      common /PSMOD/ hps,hpsk(n0:nx)
      common /PLMOD/ hpl,hplk(n0:nx)
      common /DLMOD/ hdl,hdlk(n0:nx)
      common /DVMOD/ hdv,hdvk(n0:nx)
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      iaga=0
      call ISPURE (x,icomp)
c
      herr2=hnull
      if (lsatt) then
        if (ABS(t-tsavt).lt.1.0d-9 .and. kph.eq.kphsvt) then
          lsame=.true.
          if (icomp.ne.icsavt) lsame=.false.
          if (icomp.eq.0) then
            do i=1,nc
              if (ABS(x(i)-xsavt(i)).gt.1.0d-9) lsame=.false.
            enddo
          endif
          if (lsame) then
            p=psavt
            rhol=dlsavt
            rhov=dvsavt
            do i=1,nc
              xliq(i)=xlsavt(i)
              xvap(i)=xvsavt(i)
            enddo
            ierr=0
            herr=' '
            RETURN
          endif
        endif
      endif
c  set tolerance and maximum number of iterations
      tolr=1.0d-6
      itmax=25
c     write (*,*) ' SATT--entering with t,kph = ',t,kph
      ierr=0
      iflag=0        !flag for equal densities close to critical
      herr=' '
      delp2=1.d6
      fpit(1)=0.0d0
      fpit(2)=0.0d0
      initflg=0
c  make the temperature negative to use the other parameters in the call
c  statement as initial guesses
      if (t.lt.0.d0) then
        initflg=1
        t=abs(t)
      else
        p=0.0d0
        rhol=0.0d0
        rhov=0.0d0
        if (icomp.eq.0) then
          do i=1,nc
            xliq(i)=x(i)
            xvap(i)=x(i)
          enddo
        else
          do i=1,nc
            xliq(i)=0
            xvap(i)=0
          enddo
          xliq(icomp)=1.d0
          xvap(icomp)=1.d0
        endif
      endif
c
c  check if melting or sublimation line requested and call appropriate routines.
      if (kph.eq.3) then      !liquid/solid
        call MELTT (t,x,p,ierr2,herr2)
        rhol=dtp(1)
        if (p.gt.1.d-15) call TPRHO (t,p,x,1,1,rhol,ierr,herr)
        if (ierr.eq.0) then
          ierr=ierr2
          herr=herr2
        endif
        goto 900
      else if (kph.eq.4) then      !vapor/solid
        call SUBLT (t,x,p,ierr2,herr2)
        if (t.gt.0.d0) rhov=p/(R*t)
        if (rhov.lt.1.d10) call TPRHO (t,p,x,2,1,rhov,ierr,herr)
        if (ierr.eq.0) then
          ierr=ierr2
          herr=herr2
        endif
        goto 900
      end if
c
c  check that input conditions (in this case t and x) are within limits
c
      Ddum=0.0d0
      pdum=0.0d0
      call LIMITX ('EOS',t,Ddum,pdum,x,tmin,tmax,rhomax,pmax,ierr,herr2)
      if (ierr.gt.0 .or. ierr.eq.-16) then
        if (ierr.eq.1 .and. t.lt.tmin .and. icomp.ne.0) then      !vapor/solid
          ierr2=ierr
          call SUBLT (t,x,p,ierr,herr)
          if (t.gt.0.d0 .and. p.gt.0.d0) then
            rhov=p/(R*t)
            if (rhov.lt.1.d10) call TPRHO (t,p,x,2,1,rhov,ierr3,herr3)
            goto 900
          endif
          ierr=ierr2
        endif
        ierr=abs(ierr)
c  T and/or x are out of bounds, set error flag and return
        write (herr,1000) ierr,herr2(1:238),hnull
 1000   format ('[SATT error',i3,'] ',a238,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      call CRITP (x,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
c  error condition--set outputs, issue warning, and return
        ierr=120
        p=0.0d0
        rhol=0.0d0
        rhov=0.0d0
        do i=1,nc
          xliq(i)=x(i)
          xvap(i)=x(i)
        enddo
        write (herr,1120) herr2(1:237),hnull
 1120   format ('[SATT error 120] ',a237,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      acf=0.d0
      if (icomp.eq.0) then
        do i=1,nc
          acf=acf+x(i)*accen(i)
        enddo
      endif
c
      if (icomp.ne.0) acf=accen(icomp)
      if (t.gt.tc-1.0d-8 .and. icomp.ne.0) then
c  input temperature is equal to or greater than critical,
c  return critical parameters for output pressure and densities
        p=pc
        rhol=rhoc
        rhov=rhoc
        do i=1,nc
          xliq(i)=x(i)
          xvap(i)=x(i)
        enddo
        if (t.gt.tc+1.0d-8) then
c  supercritical temperature as input, set error flag and return
c  critical parameters for output pressure and densities
          ierr=121
          write (herr,1121) t,tc,hnull
          call ERRMSG (ierr,herr)
 1121     format ('[SATT error 121] ',
     &            'temperature input to saturation routine is ',
     &            'greater than critical temperature; T =',g11.5,
     &            ' K, Tcrit =',g11.5,' K.',a1)
        end if
        RETURN
      end if
c
c  generate initial guess for pressure using acentric factor
c
      if (initflg.eq.0) then
        p=pc*10.d0**(-2.333333d0*(1.d0+acf)*(tc/t-1.d0))
      endif
c  for near-critical states, generate initial guesses for density;
c  using correlation developed by E.W. Lemmon, NIST
      theta=(ABS(1.0d0-t/tc))**(1.0d0/3.0d0)*(1.5d0+acf)
      rholi=rhoc*(1.0d0+1.113614d0*theta+0.080400d0*theta**2)
      rhovi=rhoc*(1.0d0-1.078683d0*theta+5.014057d-2*theta**2)
      if (rhovi.lt.0.d0) rhovi=1.d-6
c
      if (heos.eq.'AGA') then
        heos='HMX'
        iaga=1
      endif
c
c     if (icomp.eq.0) then
c       call SATTP(t,p,x,kph,0,d,rhol,rhov,xliq,xvap,q,ierr,herr)
c       if (ierr.eq.0) goto 900
c     endif
c
c  pure fluid iteration
      if (icomp.ne.0) then
        xliq(icomp)=1.d0
        xvap(icomp)=1.d0
c
c  Iterate for saturated liquid and vapor states given temperature using
c  a simple successive substitution method.  The independent variable
c  in the iteration is the vapor pressure.  The convergence criteria
c  is equality of Gibbs free energy in both phases.
c
        if ((t.gt.0.99d0*tc .or. iflag.eq.1) .and. ianc(icomp).eq.0)then
          goto 300
        else
c  assume nothing about densities on initial calls to TPRHO
          kguess=0
        endif
        if (hpsk(icomp).ne.' ' .and. hpsk(icomp).ne.'NBS') then
          call PSATK (icomp,t,p,ierr,herr)
        endif
        if (hplk(icomp).ne.' ' .and. hplk(icomp).ne.'NBS'
     &                         .and. kph.eq.1) then
          call PLSATK (icomp,t,p,ierr,herr)!Get liquid pressure for mixtures
        endif
        if (hdlk(icomp).ne.' ' .and. hdlk(icomp).ne.'NBS') then
          kguess=1
          call DLSATK (icomp,t,rhol,ierr,herr)
        endif
        if (hdvk(icomp).ne.' ' .and. hdvk(icomp).ne.'NBS') then
          kguess=1
          call DVSATK (icomp,t,rhov,ierr,herr)
        endif
        if (ianc(icomp).eq.1 .and. kph.eq.2) then
          call TPRHO (t,p,x,2,kguess,rhov,ierr,herr) !find vapor density
          rhol=rhov   !Don't return liquid density for equilibrium phase
          goto 900    !For pseudo-pure fluid, stop here
        elseif (ianc(icomp).eq.1) then
          call TPRHO (t,p,x,1,kguess,rhol,ierr,herr)!find liquid density
          rhov=rhol   !Don't return vapor density for equilibrium phase
          goto 900    !For pseudo-pure fluid, stop here
        endif
        do 200 it=1,itmax
c       write (*,*) 'SATT--t,p input to TPRHO: ',t,p
        call TPRHO (t,p,x,1,kguess,rhol,ierr,herr2) !find liquid density
        if (ierr.ne.0 .or. rhol.le.0.d0) then
          p=p*1.005d0
          goto 200
        end if
        call TPRHO (t,p,x,2,kguess,rhov,ierr,herr2) !find vapor density
        if (ierr.ne.0 .or. rhov.le.0.d0) then
          p=p*0.95d0
          goto 200
        end if
c  use previous densities as initial guesses for calls to TPRHO after
c  first iteration
        kguess=1
c
        call AG (t,rhol,x,Aliq,Gliq)
        call AG (t,rhov,x,Avap,Gvap)
        ZG=Gliq-Gvap
c       write (*,1014) it,t,rhol,rhov,p,Gliq,Gvap,ZG
c1014   format (1x,'SATT:  it,t,rhol,rhov,p,Gliq,Gvap,ZG: ',
c    &              i4,f8.3,2f12.8,e14.6,2f12.4,e14.6)
c
c  check convergence
c
c  check that liquid and vapor densities are different
        if (ABS(rhol-rhov).lt.1.0d-8 .or.
     &     (ABS(1.d0/rhol-1.d0/rhov).lt.0.1d0 .and. t.lt.tc-1.d0)) then
          ierr=124
          herr='[SATT error 124] density roots equal'
          call ERRMSG (ierr,herr)
          goto 300
        end if
        delp=ZG/(1.d0/rhol-1.d0/rhov)
c  the delp-delp2 check is only important for very low pressures on the
c  liquid surface (propane or R124).  See comments in TPRHO.
        if (abs(delp/p).lt.tolr .or.
     &       abs((delp-delp2)/delp2).lt.1.d-11) then
          p=p-delp
          call TPRHO (t,p,x,2,kguess,rhov,ierr,herr2) !find vap density
          if (ierr.ne.0 .or. p.gt.pc .or. rhov.gt.rhoc) then
            ierr=123
            write (herr,1123) it,herr2(1:147),hnull
 1123       format ('[SATT error 123] vapor density iteration in ',
     &            'saturation routine did not converge for pressure ',
     &            'iteration',i3,'; ',a147,a1)
            call ERRMSG (ierr,herr)
c  return critical parameters if not converged and very close to Tc
            goto 300
          end if
          call TPRHO (t,p,x,1,kguess,rhol,ierr,herr2) !find liq density
c  call THERM again to get exact p at T and rhol.  This is only important
c  for very low pressures.
c...(12-13-06 EWL)  Call THERM with RHOV if p < 1 kPa.  This was put in place
c                   because some of the pressures were coming back erratic
c                   from the PR model (for example, butane at T<160 K)
          if (p.gt.1.d0) then
            call THERM (t,rhol,x,p,e,h,s,cv,cp,w,hjt)
          else
            call THERM (t,rhov,x,p,e,h,s,cv,cp,w,hjt)
          endif
          if (ierr.ne.0 .or. p.gt.pc .or. rhol.lt.rhoc) then
            ierr=122
            write (herr,1122) it,herr2(1:146),hnull
 1122       format ('[SATT error 122] liquid density iteration in ',
     &            'saturation routine did not converge for pressure ',
     &            'iteration',i3,'; ',a146,a1)
c  return critical parameters if not converged and very close to Tc
            goto 300
          end if
c  !debug--next six lines for debug only
c         call GIBBS (t,rhol,x,Aliq,Gliq)
c         call GIBBS (t,rhov,x,Avap,Gvap)
c         ZG=Gliq-Gvap
c         write (*,1015) it,t,rhol,rhov,p,Gliq,Gvap,ZG
c1015   format (1x,'SATT:  it,t,rhol,rhov,p,Gliq,Gvap,ZG: ',
c    &              i4,f8.3,2f12.8,e14.6,2f12.4,e14.6)
          ierr=0
          herr=' '
          goto 900         !normal termination for pure fluid
        end if
        delp2=delp
c
c  continue iteration, define next guess (check that delp
c  will not result in negative [or very small] pressure)
c  02-09-98 MM:  limit step size
        if (ABS(delp).gt.0.4d0*p) then
          do jj=1,10
c           write (*,1198) p,delp
c1198       format(1x,'% SATT advisory; delp gives p<0; p,delp:',2e12.4)
            delp=0.5*delp
            if (ABS(delp).lt.0.4d0*p) goto 110
          enddo
        end if
 110    continue
        p=p-delp
        if (p.gt.1.d5) goto 210   !prevent overflow
c
 200    continue
c  iteration has not converged
 210    continue
        ierr=124
        write (herr,1124) t,hnull
 1124   format ('[SATT error 124] ',
     &          'iteration for saturation state did not converge; ',
     &          'T =',g11.5,' K.',a1)
        call ERRMSG (ierr,herr)
c  return critical parameters if not converged and very close to Tc
        if (t.gt.0.999975*tc) then
          p=pc
          rhol=rhoc
          rhov=rhoc
        end if
c
c  Alternative method for finding saturation boundaries.  The routines works
c  best near critical by finding the spinodal points on the vapor and liquid
c  sides, and uses these points to bound the iteration.
 300    continue
        if (rholi.gt.0.d0 .and. rhovi.gt.0.d0) then
          call SPNDL(t,x,rholi,rhovi,ierr,herr)
          if (ierr.gt.0) goto 390
c
c  Get the pressures at the spinodals and find the densities in the opposite
c  phase at the spinodal pressures.  The liquid pressure could be negative.
c  These new densities will then bound the iteration, i.e., the liquid density
c  will be between rholi and rholj, and likewise for the vapor density.
          call PRESS (t,rholi,x,pl)
          call PRESS (t,rhovi,x,pv)
          CALL TPRHOB (t,pv,rholi,rholi*2.d0,x,rholj,ierr,herr)
          if (pl.le.0.d0) then
            rhovj=0.d0
          else
            CALL TPRHOB (t,pl,0.d0,rhovi,x,rhovj,ierr,herr)
          endif

c  Find the saturation condition.  Use the midpoint of the spinodal pressures
c  for the initial guess.  Call TPRHOB to get both densities, and call AG to
c  get Gibbs energy.  The difference in G is used to get the next pressure.
          it=0
          p=(pl+pv)/2.d0
          if (p.lt.0.d0) p=pv/2.d0
 370      continue
          CALL TPRHOB (t,p,rhovi,rhovj,x,rhov,ierr,herr)
          if (ierr.ne.0) then
            p=p*1.001d0
            goto 380
          endif
          CALL TPRHOB (t,p,rholi,rholj,x,rhol,ierr,herr)
          if (ierr.ne.0) then
            p=p*0.999d0
            goto 380
          endif
          call AG (t,rhol,x,Aliq,Gliq)
          call AG (t,rhov,x,Avap,Gvap)
          ZG=Gliq-Gvap
          delp=ZG/(1.d0/rhol-1.d0/rhov)
          if (p-delp.lt.0.d0) then
            p=p/2.d0
          else
            p=p-delp
          endif
          if (abs(delp/p)*100.lt.tolr) goto 900   !Exit when solved
 380      continue
          it=it+1
          if (it.lt.20) goto 370
 390      continue
          ierr=124
          write (herr,1124) t,hnull
          call ERRMSG (ierr,herr)
          p=pc
          rhol=rhoc
          rhov=rhoc
        endif
c
c  end of pure fluid iteration
c
      else
c
c  begin mixture iteration
c
c  Iterate for the pressure and the composition of the incipient phase
c  (vapor phase for a bubble point calculation, liquid for dew point)
c  given temperature and the composition of the parent phase.  Iteration
c  is generally based on the algorithm given by Smith & Van Ness (Intro
c  to Chem Engr Thermo, McGraw-Hill, 1975); convergence criteria is the
c  equality of fugacity for each component in both phases.
c
        x2sum=0.d0
        tolr=1.d-8
        call SATTEST (t,x,kph,psum,x2,ierr,herr)
        if (psum.gt.pc) psum=0.99d0*pc  !helps critical region conv.
c
c  variable kph2 specifies the state of the incipient phase (x2):
c     1 = liq,  2 = vap
c  it is used in calls to TPRHO
c  pratio is multiplier for pressure to use when TPRHO does not converge
c  different values for liquid and vapor phases, such that new guess for
c  pressure is further into corresponding single-phase region
        if (kph.eq.1) then
          kph2=2
          prtio1=1.05d0          !Choose different pratio for the dew
          prtio2=1.02d0          !and bubble sides to avoid loops
          if (initflg.ne.0) then
            do i=1,nc
              x2(i)=xvap(i)
            enddo
            psum=p
          endif
        else
          kph2=1
          prtio1=0.95d0
          prtio2=0.98d0
          if (initflg.ne.0) then
            do i=1,nc
              x2(i)=xliq(i)
            enddo
            psum=p
          endif
        end if
c
c  begin main mixture iteration--outer loop for pressure,
c  using Raoult's Law result from above as first guess
c
        p=psum
        do i=1,nc
          if (kph.eq.1) then
            xvap(i)=x2(i)
          else
            xliq(i)=x2(i)
          endif
        enddo
c
        kgues1=0              !for first calls to TPRHO
        kgues2=0              !kguess flags for parent & incipient phase
        lppos=.false.         !flags for reguli-falsi iteration
        lpneg=.false.
        pneg=0.0d0
        ppos=0.0d0
        fpneg=0.0d0
        fppos=0.0d0
        if (initflg.eq.0) then
          rho1=0.0d0
          rho2=0.0d0
        else
          kgues1=1
          kgues2=1
          if (kph.eq.1) then
            rho1=rhol
            rho2=rhov
          else
            rho1=rhov
            rho2=rhol
          endif
        endif
        tbad2=0.0d0
        ibad1=0
        ibad2=0
        j=1
        do ii=1,nc
          x2org(ii)=x2(ii)
        enddo
        pit(1)=psum           !first guess for pressure = sum (x1*Pi)
c
c
c
c       write (*,*) ' SATT--begin outer iteration loop for pressure'
        do 400 itp=1,itmax
c  increase tolerance to account for errors in numerical derivatives in FGCTY2
        if (itp.eq.8) tolr=tolr*10
        if (itp.eq.12) tolr=tolr*10
        if (itp.eq.15) tolr=tolr*10
        p=pit(j)
c       write (*,*) ' SATT--pressure iteration',itp,' w/ p =',pit(j)
        lx2con=.false.        !flag for convergence of inner loop
c  compute density and fugacities for parent phase
        i=iprnterr
        iprnterr=0
        call TPRHO (t,p,x,-kph,kgues1,rho1,ierr,herr2)   !parent phase
        iprnterr=i
        kgues1=1
        if (kph.eq.1) then
          if (rho1.lt.rhoc/2.d0) ierr=1
          if (rho1.lt.rhoc*1.2d0 .and.itp.eq.1 .and.initflg.eq.0) ierr=1
        endif
        if (ierr.gt.0) then
          ibad1=ibad1+1
          if (ibad1.lt.6) then
            if (kph.eq.1) then
              call DLDV (t,p,rho1,rho2,x,x2,ierr,herr2)
            else
              call DLDV (t,p,rho2,rho1,x2,x,ierr,herr2)
            endif
            pit(j)=p
          endif
          if (ierr.ne.0) then
            ierr=-125
            write (herr,1125) itp,herr2(1:149),hnull
 1125       format ('[SATT advisory -125] density iteration in ',
     &              'saturation routine did not converge for pressure ',
     &              'iteration',i3,'; ',a149,a1)
            kgues1=0              !do not reuse faulty density as guess
            pit(j)=pit(j)*prtio1 !try another pressure and use up one
            goto 400              !iteration (to prevent infinite loop)
          else
            kgues2=1
          endif
        end if
        call FGCTY2 (t,rho1,x,f1,ierr,herr)
c       write (*,1082) itp,t,p,rho1,(x(i),i=1,2),(f1(i),i=1,2)
c1082   format (1x,' SATT--phase 1:  ',i3,f8.2,2e14.6,2e16.8,2e18.10)
c
c  begin inner iteration loop for composition of phase 2
c
c       write (*,*) 'SATT--begin inner loop for composition of phase 2'
c
        sumdl2=2.0d0
        do itx=1,itmax*2
        do i=1,nc
        xs(i)=x2(i)
        enddo
c  compute density and fugacities for phase 2
        i=iprnterr
        iprnterr=0
        call TPRHO (t,p,x2,-kph2,kgues2,rho2,ierr,herr2)!incipient phase
        if (itx.ne.1) then
          if (ABS(rho2-rhoc).lt.1 .and. ABS(rho1-rhoc).gt.10) then
            rho2o=rho2
c  jump away from critical point when a bad root was found
            rho2=rho2*.5d0
            ierr2=ierr
            call TPRHO (t,p,x2,-kph2,kgues2,rho2,ierr,herr3)
            if (ierr.ne.0 .or. ABS(rho1-rho2).lt.0.1d0) then
              herr2=herr3
              ierr=ierr2
              rho2=rho2o
            endif
          endif
        endif
        iprnterr=i
        kgues2=1
        itx1=0
 420    continue
        if (ierr.gt.0 .or. ABS(rho1-rho2).lt.0.1d0) then
          ibad2=ibad2+1
          if (ibad2.lt.15) then
            if (ABS(rho1-rho2).lt.0.1d0 .and.
     &      ABS(tcrit(1)-tcrit(2)).gt.100.d0) then
              do i=1,nc
                x2(i)=(tbad2*x(i)+x2org(i))/(1.d0+tbad2)
              enddo
              tbad2=tbad2+0.2d0
            endif
            if (kph.eq.1) then
              call DLDV (t,p,rho1,rho2,x,x2,ierr,herr2)
            else
              call DLDV (t,p,rho2,rho1,x2,x,ierr,herr2)
            endif
            call FGCTY2 (t,rho1,x,f1,ierr,herr)
            pit(j)=p
          endif
          if (ABS(f1(1)).gt.1.d10) goto 810
          if (ierr.ne.0 .or. ABS(rho1-rho2).lt.0.1d0) then
            ierr=-126
            write (herr,1126) itx,herr2(1:146),hnull
 1126       format ('[SATT advisory -126] density iteration in ',
     &            'saturation routine did not converge for composition',
     &            ' iteration',i3,'; ',a146,a1)
            do ii=1,nc
              x2(ii)=x2org(ii)
            enddo
            kgues2=0              !do not reuse faulty density as guess
            pit(j)=pit(j)/prtio2 !try another pressure and use up one
            goto 400              !iteration (to prevent infinite loop)
          endif
        endif
        call FGCTY2 (t,rho2,x2,f2,ierr,herr)
c       write (*,1086) itx,rho2,(x2(i),i=1,2),(f2(i),i=1,2)
c1086   format (1x,'       phase 2:  ',i3,22x,e14.6,2e16.8,2e18.10)
c  calculate new x2's by ratio of fugacities; inner loop has converged
c  when x2's change by less than a convergence tolerance
        x2sum=0.0d0
        do i=1,nc
          if (f2(i).gt.1.0d-20 .and. f1(i).gt.0.0d0) then
            x2new(i)=x2(i)*f1(i)/f2(i)
          else
c  in case fugacity is zero; e.g., if x(i)=0, then slightly modify x2(i)
            if (x(i).gt.0.d0) then
              rho2=rho1
              itx1=itx1+1
              if (itx1.lt.10) goto 420
            endif
            x2new(i)=0
c           if (x2(i).gt.0) x2new(i)=x2(i)+0.01d0
c           if (x2new(i).gt.1.0d0) x2new(i)=x2(i)-0.01d0
          end if
          x2sum=x2sum+x2new(i)
        enddo
        if (x2sum.le.0.d0 .or. x2sum.gt.1.d6) goto 810
c  normalize the x2 compositions; this yields next guess for x2 and
c  ensures that the x2 always sum to one
        sumdel=0.0d0
        do i=1,nc
          x2new(i)=x2new(i)/x2sum
          sumdel=sumdel+abs(x2(i)-x2new(i))  !change in compositions
          x2(i)=x2new(i)
        enddo
c       write (*,1560) (x2(i),i=1,nc)
c1560   format (1x,'SATT:  new compositions in inner loop:  ',5f10.5)
        if (sumdel.lt.tolr .or. abs(sumdel-sumdl2).lt.tolr*100.d0) then
c  inner iteration loop has converged
          lx2con=.true.
          ierr=0
          goto 700
        end if
        sumdl2=sumdel
c  if not, continue inner iteration loop
        if (int(itx/10)*10.eq.itx) then
c  occasionally average out the compositions on two successive iterations.
c  sometimes x2 bounces back and forth between two values.
          do i=1,nc
            x2(i)=(x2(i)+xs(i))/2.d0
          enddo
        endif
        enddo
c
c  inner iteration loop has not converged
        ierr=-127
        write (herr,1127) t,sumdel,hnull
 1127   format ('[SATT advisory -127] ',
     &          'iteration for composition in saturation routine ',
     &          'did not converge; T =',g11.5,
     &          ' K; deltaX =',g11.5,' mol frac.',a1)
        call ERRMSG (ierr,herr)
c
c  end of inner (x2) iteration loop
c
 700    continue
        fpit(j)=1.0d0-x2sum
c  outer (pressure) loop has converged when the x2's sum to one, i.e.,
c  when the fugacities of each component in each phase are equal
c       write (*,*) ' SATT--check conv, p, fp: ',pit(j),fpit(j)
        call DPDT (t,rho2,x2,dpt2)
        if (dpt2.le.0.d0) goto 810  !dp/dT should never be negative
        if (ABS(fpit(j)).lt.tolr .and. lx2con) then
          ierr=0
          herr=' '
          goto 850
        else
c  provided that the inner loop has converged, update positive and
c  negative bounds on pressure for possible use in reguli-falsi iteration
          if (lx2con) then
            if (fpit(j).lt.0.0d0) then
              lpneg=.true.
              pneg=pit(j)
              fpneg=fpit(j)
            else
              lppos=.true.
              ppos=pit(j)
              fppos=fpit(j)
            end if
          end if
        end if
c
c  compute new guess for saturation pressure
c
        if (j.eq.1) then
c  for first iteration, new pressure is ratio of old
          if (x2sum.gt.5) x2sum=5
          j=2
          if (kph.eq.1) then
c  bubble point
            pit(2)=pit(1)*x2sum
          else
c  dew point
            pit(2)=pit(1)/x2sum
          end if
        else
c  subsequent iterations--use secant method
          if (ABS(fpit(2)-fpit(1)).gt.1.0d-12)
     &    pit(3)=pit(2)-fpit(2)*(pit(2)-pit(1))/(fpit(2)-fpit(1))
c  check that new pressure is not outside bounds, if so use reguli-falsi
          if (lpneg .and. lppos .and. (pit(3).gt.1.001d0*MAX(ppos,pneg)
     &        .or. pit(3).lt.0.990d0*MIN(ppos,pneg))) then
            pit(3)=ppos-fppos*(ppos-pneg)/(fppos-fpneg)
          end if
c  discard oldest iteration
          pit(1)=pit(2)
          pit(2)=pit(3)
          fpit(1)=fpit(2)
        end if
c       write (*,1799) itp,j,pit(1),pit(2),fpit(2)
c1799   format (1x,' SATT--itp,j,p1,p2,fp2:  ',2i4,3e14.6)
 400    continue
c  outer iteration loop has not converged
 810    continue
        ierr=128
        write (herr,1128) t,hnull
 1128   format ('[SATT error 128] ',
     &          'iteration for saturation state did not converge; T =',
     &          g11.5,' K.',a1)
        call ERRMSG (ierr,herr)
        p=pc
        rho1=rhoc
        rho2=rhoc
c
c  end of outer (pressure) iteration loop
c
 850    continue
c
c  assign final compositions and densities for parent and incipient
c  phases (x and x2, rho1 and rho2, respectively) to outputs
c
        if (kph.eq.1) then
c  bubble point
          rhol=rho1
          rhov=rho2
          do i=1,nc
            xliq(i)=x(i)
            xvap(i)=x2(i)
          enddo
        else
c  dew point
          rhol=rho2
          rhov=rho1
          do i=1,nc
            xliq(i)=x2(i)
            xvap(i)=x(i)
          enddo
        end if
c  call new routine SATTP if ierr>0 and attempt to get convergence
        if (ierr.gt.0 .and. icomp.eq.0) then
          call SATTP(t,p,x,kph,0,d,rhol,rhov,xliq,xvap,q,ierr,herr)
        endif
        if (rhov.gt.rhol*0.99d0) then
          if (ABS(xliq(1)-xvap(1)).lt.0.1d0) then
            p=pc
            rhov=rhoc
            rhol=rhoc
            if (t.lt.tc) then
              ierr=128
              write (herr,1128) t,hnull
              call ERRMSG (ierr,herr)
            else
              ierr=121
              write (herr,1121) t,tc,hnull
              call ERRMSG (ierr,herr)
            endif
          endif
        endif
c  end of mixture iteration
      end if
c
c  save results
 900  continue
      if (ierr.eq.0) then
        if (p.lt.0.00001d0 .and. icomp.ne.0) p=rhov*r*t
        tsavt=t
        psavt=p
        icsavt=icomp
        kphsvt=kph
        dlsavt=rhol
        dvsavt=rhov
        if (icomp.eq.0) then
          do i=1,nc
            xsavt(i)=x(i)
            xlsavt(i)=xliq(i)
            xvsavt(i)=xvap(i)
          enddo
        else
          xsavt(icomp)=1.d0
          xlsavt(icomp)=1.d0
          xvsavt(icomp)=1.d0
        endif
        lsatt=.true.
      endif
      if (iaga.eq.1) then
        heos='AGA'
        if (ierr.eq.0) call TPRHO (t,p,xvap,2,1,rhov,ierr,herr)
      endif
      RETURN
c
      end                                               !subroutine SATT
c
c ======================================================================
c
      subroutine SATTEST (t,x,kph,p,x2,ierr,herr)
c
c  estimate initial values for saturation states given temperature
c  and the composition of one phase
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac] (phase specified by kph)
c      kph--phase flag: 1 = input x is liquid composition (bubble point)
c                       2 = input x is vapor composition (dew point)
c  outputs:
c        p--estimated pressure [kPa]
c       x2--estimated composition of unknown phase [array of mol frac]
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  08-06-09 EWL, original version, taken from code in SATT
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SATTEST
c     dll_export SATTEST
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      dimension x(ncmax),x2(ncmax),pcomp(ncmax)
      character*255 herr
      common /NCOMP/ nc,ic
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
c  generate initial guesses for mixture pressure and compositions
c  using Raoult's law
c
c  estimates for vapor pressures of pure components
      ierr=0
      herr=' '
      do i=1,nc
        pcomp(i)=pcrit(i)
     &        *10.d0**(-2.333333d0*(1.d0+accen(i))*(tcrit(i)/t-1.d0))
      enddo
c
      p=0.d0
      if (kph.eq.1) then
c  bubble point
        do i=1,nc
          p=p+x(i)*pcomp(i)
        enddo
        do i=1,nc
          x2(i)=x(i)*pcomp(i)/p
        enddo
      else
c  dew point
        ysum=0.d0
        do i=1,nc
          ysum=ysum+x(i)/pcomp(i)
        enddo
        do i=1,nc
          x2(i)=x(i)/pcomp(i)/ysum
          p=p+x2(i)*pcomp(i)
        enddo
      end if
      RETURN
c
      end                                            !subroutine SATTEST
c
c ======================================================================
c
      subroutine SATP (p,x,kph,t,rhol,rhov,xliq,xvap,ierr,herr)
c
c  iterate for saturated liquid and vapor states given pressure
c  and the composition of one phase
c
c  inputs:
c        p--pressure [kPa]
c        x--composition [array of mol frac] (phase specified by kph)
c      kph--phase flag:  1 = input x is liquid composition
c                        2 = input x is vapor composition
c                        3 = input x is liquid composition (freezing point)
c                        4 = input x is vapor composition (sublimation point)
c
c  outputs:
c        t--temperature [K]
c     rhol--molar density [mol/L] of saturated liquid
c     rhov--molar density [mol/L] of saturated vapor
c           For a pseudo pure fluid, the density of the equilibrium phase
c           is not returned.  Call SATP twice, once with kph=1 to get
c           tliq and rhol, and once with kph=2 to get tvap and rhov.
c     xliq--liquid phase composition [array of mol frac]
c     xvap--vapor phase composition [array of mol frac]
c     ierr--error flag:  0 = successful
c                        2 = P < Ptp
c                        4 = P < 0
c                        8 = x out of range
c                       12 = P and x out of range
c                      140 = CRITP did not converge
c                      141 = P > Pcrit
c                      142 = TPRHO-liquid did not converge (pure fluid)
c                      143 = TPRHO-vapor did not converge (pure fluid)
c                      144 = pure fluid iteration did not converge
c           following 3 error codes are advisory--iteration will either
c           converge on later guess or error out (ierr = 148)
c                     -144 = Raoult's law (mixture initial guess) did
c                            not converge
c                     -145 = TPRHO did not converge for parent ph (mix)
c                     -146 = TPRHO did not converge for incipient (mix)
c                     -147 = composition iteration did not converge
c                      148 = mixture iteration did not converge
c     herr--error string if ierr<>0 (character*255)
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  07-13-95  MM, original version
c  09-11-95  MM, add error string to argument list
c  09-25-95  MM, rearrange argument list (outputs in order t, rho, x)
c  10-11-95  MM, RETURN if any error detected
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c  12-11-95  MM, Raoult's law as first guess for mixture
c  12-18-95  MM, fill xliq, xvap with zeros for undefined components
c  12-27-95  MM, add full mixture iteration using fugacity, based on SATT
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                add Zcrit to common /CCON/
c  03-19-96  MM, add dipole moment to /CCON/
c  04-05-96  MM, test for supercritical '.ge. pc' rather than '.gt. pc'
c  11-14-96  MM, adjust initial guesses, tratio, etc to get closer to critical
c  02-12-96 EWL, special initial guess for temperature near the critical point
c  07-15-97  MM, add errors/warnings to parallel SATT
c  10-01-97  MM, add compiler switch to allow access by DLL
c  11-13-97 EWL, initialize ft(j); potential bomb if no value when writing error message
c  12-05-97  MM, check that TPRHO gives density within bounds
c                if Raoult's law iteration D.N.C., revert to initial guess
c  02-10-98  MM, add reguli-falsi, quadratic interpolation and bisection
c                to Raoult's law iteration
c  04-06-98 EWL, call subroutine AG to get Gibbs energy
c  07-30-98  MM, check for f1=0 in phase 2 iteration; separate kguess for each phase
c  11-30-98 EWL, change kguess to kgues2 in call to TPRHO in vapor search.
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  12-22-98 EWL, zero out all elements of xliq and xvap
c  12-22-98 EWL, use pure fluid algorithm when called with a mixture but x(i)=1
c  12-23-98 EWL, modify step size for T when TPRHO fails
c  12-23-98 EWL, replace rhov with rhoc/2 when TPRHO finds liquid root
c  08-03-99 EWL, if TPRHO fails, set x2 to original values before restarting
c  08-03-99 EWL, if dpdrho>1d6, then return an error message
c  01-11-00 EWL, remove ierr from lines where 'i3' was not in format statement
c  02-23-00 EWL, add kph=3,4 as inputs and call melting or sublimation lines
c  12-17-01 EWL, call SATT if iteration fails
c  05-28-02 EWL, check for bad root in two phase
c  08-14-02 EWL, add logic if ammonia/water in use
c  09-19-02 EWL, add checks for ancillary routines and exit for pseudo-pures
c  03-24-05 EWL, add check for small rho in liquid phase. (for mixtures, inner loop)
c  09-21-06 EWL, increase itmax, check for ierr=-147, add small value to keep x2sum<>1,
c                add check for nh3+h2o, split tratio into trtio1 and trtio2
c  03-07-07 EWL, add check for very low pcomp(i)
c  02-18-08 EWL, add check for icomp<>old icomp (icsav)
c  04-22-10 EWL, rename variables in PSTSAV
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SATP
c     dll_export SATP
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herr3
      character*3 hps,hpsk,hpl,hplk,hdl,hdlk,hdv,hdvk
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),x2org(ncmax)
      dimension x2(ncmax),f1(ncmax),f2(ncmax),x2new(ncmax)
      dimension tk(4),ft(3)                  !used for T iteration
      character*3 hpheq,heos,hmxeos,hmodcp
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /FSHSAV/ lsatt,lsatp
      common /PSTSAV/ tsavp,psavp,dlsavp,dvsavp,
     &                xsavp(ncmax),xlsavp(ncmax),xvsavp(ncmax),
     &                kphsvp,icsavp
      common /PSMOD/ hps,hpsk(n0:nx)
      common /PLMOD/ hpl,hplk(n0:nx)
      common /DLMOD/ hdl,hdlk(n0:nx)
      common /DVMOD/ hdv,hdvk(n0:nx)
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      iaga=0
      call ISPURE (x,icomp)
      ierr=0
      herr=' '
      herr2=hnull
      if (lsatp) then
        if (ABS(p-psavp).lt.1.d-9 .and. kph.eq.kphsvp
     &          .and. p.gt.1.d-7) then
          lsame=.true.
          if (icomp.ne.icsavp) lsame=.false.
          if (icomp.eq.0) then
            do i=1,nc
              if (ABS(x(i)-xsavp(i)).gt.1.0d-9) lsame=.false.
            enddo
          endif
          if (lsame) then
            t=tsavp
            rhol=dlsavp
            rhov=dvsavp
            do i=1,nc
              xliq(i)=xlsavp(i)
              xvap(i)=xvsavp(i)
            enddo
            RETURN
          endif
        endif
      endif
c  set tolerance and maximum number of iterations
      tolr=1.0d-6
      itmax=100
      ft(1)=0.0d0  !initialize to avoid potential problem
      ft(2)=0.0d0  !when writing error message to GUI
c
c  initialize outputs in event of failure of routines
      t=300.0d0
      rhol=0.0d0
      rhov=0.0d0
      if (icomp.eq.0) then
        do i=1,nc
          xliq(i)=x(i)
          xvap(i)=x(i)
        enddo
      else
        do i=1,nc
          xliq(i)=0
          xvap(i)=0
        enddo
        xliq(icomp)=1.d0
        xvap(icomp)=1.d0
      endif
c
c  check if melting or sublimation line requested and call appropriate routines.
      if (kph.eq.3) then      !liquid/solid
        call MELTP (p,x,t,ierr2,herr2)
        rhol=dtp(1)
        if (t.gt.0.d0) call TPRHO (t,p,x,1,1,rhol,ierr,herr)
        if (ierr.eq.0) then
          ierr=ierr2
          herr=herr2
        endif
        goto 900
      else if (kph.eq.4) then      !vapor/solid
        call SUBLP (p,x,t,ierr2,herr2)
        if (t.gt.0.d0) then
          rhov=p/(R*t)
          if (rhov.lt.1.d10) call TPRHO (t,p,x,2,1,rhov,ierr,herr)
        endif
        if (ierr.eq.0) then
          ierr=ierr2
          herr=herr2
        endif
        goto 900
      end if
c
      call CRITP (x,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
c  error condition--set outputs, issue warning, and return
        ierr=140
        write (herr,1140) herr2(1:237),hnull
 1140   format ('[SATP error 140] ',a237,a1)
        call ERRMSG (ierr,herr2)
        RETURN
      end if
c
c  check that input conditions (in this case p and x) are within limits
c
      Ddum=0.0d0
      tdum=0.8d0*tc
      call LIMITX ('EOS',tdum,Ddum,p,x,tmin,tmax,rhomax,pmax,ierr,herr2)
c     write (*,*) ' SATP--density limit:  ',rhomax
      if (ierr.gt.1) then     !ignore ierr = 1 (t out of range)
c  p and/or x are out of bounds, set error flag and return
        t=0.8d0*tc
        write (herr,1000) ierr,herr2(1:238),hnull
 1000   format ('[SATP error',i3,'] ',a238,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0) then
        if (p.lt.ptp(icomp)) then
          call SUBLP (p,x,t,ierr,herr)        !vapor/solid
          if (t.gt.0.d0) then
            rhov=p/(R*t)
            if (rhov.lt.1.d10) call TPRHO (t,p,x,2,1,rhov,ierr3,herr3)
            goto 900
          endif
          ierr=2
          t=ttp(icomp)
          write (herr,1149) p/1000.0d0,ptp(icomp)/1000.0d0,hnull
          call ERRMSG (ierr,herr)
 1149     format ('[SATP error 2] ',
     &            'pressure less than triple point pressure; P =',
     &          g11.5,' MPa, Ptp =',g11.5,' MPa.',a1)
          return
        endif
      endif
c
      if (p/pc.ge.0.99999999d0 .and. icomp.ne.0) then
c  input pressure is equal to or greater than critical point value,
c  return critical parameters for output temperature and densities
c  output compositions initialized above
        t=tc
        rhol=rhoc
        rhov=rhoc
        if (p/pc.gt.1.000001d0) then
c  supercritical pressure as input, set error flag and return
c  critical parameters for output pressure and densities
          ierr=141
          write (herr,1141) p/1000.0d0,pc/1000.0d0,hnull
          call ERRMSG (ierr,herr)
 1141     format ('[SATP error 141] ',
     &            'pressure input to saturation routine is ',
     &            'greater than critical pressure; P =',g11.5,
     &            ' MPa, Pcrit =',g11.5,' MPa.',a1)
        end if
        RETURN
      end if
c
      if (heos.eq.'AGA') then
        heos='HMX'
        iaga=1
      endif
C
c     if (icomp.eq.0) then
c       call SATTP(t,p,x,kph+2,0,d,rhol,rhov,xliq,xvap,q,ierr,herr)
c       if (ierr.eq.0) goto 900
c     endif
c
      if (icomp.ne.0) then
c
c  pure fluid iteration
        xliq(icomp)=1.d0
        xvap(icomp)=1.d0
c
c  Iterate for saturated liquid and vapor states given pressure using
c  a simple successive substitution method.  The independent variable
c  in the iteration is the temperature.  The convergence criteria is
c  equality of Gibbs free energy in both phases.
c
c  generate initial guess using acentric factor
c
        t=tc/(1.0-0.428571*log10(p/pc)/(1.0+accen(icomp)))
c       write (*,*) 'SATP:  T(0): ',t
c
c  assume nothing about densities on initial calls to TPRHO
        kguess=0
        if (p.gt.0.98*pc) then
c  for near-critical states, generate initial guesses for density;
c  using correlation developed by E.W. Lemmon, NIST
          theta=(1-p/pc)*(1.6d0-accen(icomp))
          t=tc*(1.0d0-0.103947d0*theta-4.108265d-2*theta**2)
          theta=(1-p/pc)**(1.0d0/3.0d0)*(3.0d0+accen(icomp))
          rhov=rhoc*(1.0d0-0.290039d0*theta-7.120197d-3*theta**2)
          rhol=rhoc*(1.0d0+0.298544d0*theta+1.870808d-2*theta**2)
          kguess = 1
        endif
        if (hpsk(icomp).ne.' ' .and. hpsk(icomp).ne.'NBS') then
          call TSATP (p,x,tt,ierr,herr)
          if (tt.gt.0.) t=tt
        endif
        if (hplk(icomp).ne.' ' .and. hplk(icomp).ne.'NBS'
     &                         .and. kph.eq.1) then
          call TSATPL (p,x,tt,ierr,herr)!Get liquid pressure for mixtures
          if (tt.gt.0.) t=tt
        endif
        if (hdlk(icomp).ne.' ' .and. hdlk(icomp).ne.'NBS') then
          kguess=1
          call DLSATK (icomp,t,rhol,ierr,herr)
        endif
        if (hdvk(icomp).ne.' ' .and. hdvk(icomp).ne.'NBS') then
          kguess=1
          call DVSATK (icomp,t,rhov,ierr,herr)
        endif
        stp=1.0002d0
        if (ianc(icomp).eq.1 .and. kph.eq.2) then
          call TPRHO (t,p,x,2,kguess,rhov,ierr,herr) !find vapor density
          rhol=rhov   !Don't return liquid density for equilibrium phase
          goto 900    !For pseudo-pure fluid, stop here
        elseif (ianc(icomp).eq.1) then
          call TPRHO (t,p,x,1,kguess,rhol,ierr,herr)!find liquid density
          rhov=rhol   !Don't return vapor density for equilibrium phase
          goto 900    !For pseudo-pure fluid, stop here
        endif
        do 200 it=1,itmax
        call TPRHO (t,p,x,1,kguess,rhol,ierr,herr2) !find liquid density
        if (ierr.ne.0 .or. rhol.lt.rhoc) then
          t=t/stp
          stp=1.0d0+(stp-1.0d0)/1.2d0
          goto 200
c         herr=' ERROR from SATP:  '//herr2
c         call ERRMSG (ierr,herr)
c         RETURN
        end if
        call TPRHO (t,p,x,2,kguess,rhov,ierr,herr2)  !find vapor density
        if (ierr.ne.0 .or. rhov.gt.rhoc) then
          t=t*1.001d0
          if (t.gt.tc) t=t/1.001d0*1.00005d0
          if (rhov.gt.rhoc) rhov=rhoc/2.0d0
          goto 200
c         herr=' ERROR from SATP:  '//herr2
c         call ERRMSG (ierr,herr)
c         RETURN
        end if
c  use previous densities as initial guesses for calls to TPRHO after
c  first iteration
        kguess=1
c       call GIBBS (t,rhol,x,Aliq,Gliq)
c       call GIBBS (t,rhov,x,Avap,Gvap)
        call AG (t,rhol,x,Aliq,Gliq)
        call AG (t,rhov,x,Avap,Gvap)
        call ENTRO (t,rhol,x,sliq)
        call ENTRO (t,rhov,x,svap)
        ZG=Gliq-Gvap
c       write (*,1014) it,p,rhol,rhov,t,Gliq,Gvap,ZG
c1014   format (1x,'SATP:',i4,e14.6,2f12.8,f10.5,2f12.4,e14.6)
c
c  check convergence
c
        delt=ZG/(sliq-svap)
        if (abs(delt).lt.tolr) then
c  pure component iteration is done (make use of current delt)
          t=t+delt
          call TPRHO (t,p,x,1,kguess,rhol,ierr,herr2) !find liq density
          if (ierr.ne.0 .or. t.gt.tc .or. rhol.lt.rhoc) then
            ierr=142
            write (herr,1142) it,herr2(1:142),hnull
 1142       format ('[SATP error 142] liquid density iteration in ',
     &            'saturation routine did not converge for temperature',
     &            ' iteration',i3,'; ',a142,a1)
            call ERRMSG (ierr,herr)
c  return critical parameters if not converged and very close to Pc
            if (p.gt.0.9999*pc) then
              t=tc
              rhol=rhoc
              rhov=rhoc
            end if
            goto 300
c           RETURN
          end if
          call TPRHO (t,p,x,2,kguess,rhov,ierr,herr2) !find vap density
          if (ierr.ne.0 .or. p.gt.pc .or. rhov.gt.rhoc) then
            ierr=143
            write (herr,1143) it,herr2(1:143),hnull
 1143       format ('[SATP error 143] vapor density iteration in ',
     &            'saturation routine did not converge for temperature',
     &            ' iteration',i3,'; ',a143,a1)
            call ERRMSG (ierr,herr)
c  return critical parameters if not converged and very close to Pc
            if (p.gt.0.9999*pc) then
              t=tc
              rhol=rhoc
              rhov=rhoc
            end if
            goto 300
c           RETURN
          end if
c  !debug--next four lines for debug only
c         call GIBBS (t,rhol,x,Aliq,Gliq)
c         call GIBBS (t,rhov,x,Avap,Gvap)
c         ZG=Gliq-Gvap
c         write (*,1015) it,p,rhol,rhov,t,Gliq,Gvap,ZG
c1015   format (1x,'SATP:',i4,e14.6,2f12.8,f10.5,2f12.4,e14.6)
          ierr=0
          herr=' '
          goto 900         !normal termination for pure fluid
        end if
c
c  continue iteration, define next guess (check that delt
c  will not result in too large a change in temperature)
        if (delt.gt.0.5*t) then
          do j=1,100
c           write (*,1198) t,delt
c1198       format (1x,'% SATP advisory; delt > 0.5*t; t,delt:',2f12.6)
            delt=0.25*delt
            if (delt.lt.0.5*t) goto 110
          enddo
        end if
 110    continue
        t=t+delt
c
 200    continue
c
c  In case of failure, try calling SATT iteratively to find the saturated
c  temperature.  This takes advantage of the alternate method in SATT used
c  at temperatures very close to the critical point.
 300    continue
        i=0
        t1=tc*.999d0
        call SATT (t1,x,kph,p1,rhol,rhov,xliq,xvap,ierr,herr)
        t2=tc*.9995d0
 310    continue
        call SATT (t2,x,kph,p2,rhol,rhov,xliq,xvap,ierr,herr)
        t=t2
        if (ABS(p2-p).lt.tolr .and. ierr.eq.0) goto 900  !Convergence
        i=i+1
        if (ABS(p2-p1).lt.1.d-12) goto 320
        if (i.gt.20) goto 320
        t=t1-(p1-p)/(p2-p1)*(t2-t1)
        if (t.gt.tc .and. t2.gt.t1) t=(t2+tc)/2.d0
        if (t.gt.tc .and. t1.gt.t2) t=(t1+tc)/2.d0
        t1=t2
        p1=p2
        t2=t
        goto 310
c
c  iteration has not converged
 320    continue
        ierr=144
        write (herr,1144) p/1000.0d0,hnull
 1144   format ('[SATP error 144] ',
     &          'iteration for saturation state did not converge; ',
     &          'P =',g11.5,' MPa.',a1)
        call ERRMSG (ierr,herr)
c  return critical parameters if not converged and very close to Tc
        if (p.gt.0.9999*pc) then
          t=tc
          rhol=rhoc
          rhov=rhoc
        end if
        goto 900
c
c  end of pure fluid iteration
c
      else
c
c  begin mixture iteration
c
c  Iterate for the temperature and the composition of the incipient
c  phase (vapor phase for a bubble point calculation, liquid for dew
c  point) given pressure and the composition of the parent phase.
c  Iteration is generally based on the algorithm given by Van Ness
c  & Abbott (Classical Thermodynamics of Nonelectrolyte Solutions with
c  Applications to Phase Equilibria, McGraw-Hill, 1982); convergence
c  criteria is the equality of fugacity for each component in both
c  phases.
c
        call SATPEST (p,x,kph,t,x2,ierr,herr)
        j=2
        tk(j)=t
        if (kph.eq.1) then
c  bubble point
c  variable kph2 specifies the state of the incipient phase (x2): 1=liq, 2=vap
c  trtio is temperature multiplier to use when TPRHO does not converge;
c  different values for liquid and vapor phases, such that new guess
c  is further into corresponding single-phase region
          kph2=2
          trtio1=0.995d0
          trtio2=0.998d0
        else
c  dew point
          kph2=1
          trtio1=1.005d0
          trtio2=1.002d0
        end if
c
c  initial temperature (satisfying Raoult's law) has been found,
c  check that this temperature is not above critical
        tmax=0.998d0*tc
        if (tk(j).gt.tmax .or. tk(j).le.0) tk(j)=tmax
c
c  the following line adjusts the initial guess for ammonia/water mixtures
c  allowing the saturation routines to work substantially better, especially in
c  the critical region.  It may work for other mixtures as well.
        if (iamwat.ne.0) tk(j)=tk(j)*1.05d0
c
c  generate initial guesses for densities & incipient phase composition;
c  the do loop allows for the possibility that TPRHO does not converge,
c  it should normally exit with just one pass
c
        do it=1,itmax
          t=tk(j)
c  first guess for densities (separate flags for each phase)
          kgues1=0
          kgues2=0
          call TPRHO (t,p,x,kph,kgues1,rho1,ierr,herr2)     !parent phase
          kgues1=1
          if (ierr.gt.0) then
            kgues1=0             !do not reuse faulty density as guess
            ierr=-145
            write (herr,1145) herr2(1:146),hnull
 1145       format ('[SATP advisory -145] density iteration in ',
     &            'saturation routine did not converge for the parent',
     &            ' phase; ',a146,a1)
            call ERRMSG (ierr,herr)
            if (t.gt.0.8*tc .and. kph.eq.1) then
c   non-convergence probably because too close to critical
              tnew=0.995*t
            else
              tnew=1.005*t
            end if
            tk(j)=tnew
          else
            call TPRHO (t,p,x2,kph2,kgues2,rho2,ierr,herr2) !incipient ph
            kgues2=1
            if (ierr.gt.0) then
              kgues2=0             !do not reuse faulty density as guess
              ierr=-144
              write (herr,1146) herr2(1:146),hnull
 1146         format ('[SATP advisory -146] density iteration in ',
     &              'saturation routine did not converge for the ',
     &              'incipient phase; ',a146,a1)
              call ERRMSG (ierr,herr)
              if (t.gt.0.8*tc .and. kph2.eq.1) then
c   non-convergence probably because too close to critical
                tnew=0.995*t
              else
                tnew=1.005*t
              end if
              tk(j)=tnew
            else
              goto 550   !both parent and incipient phases have converged
            end if
          end if
        enddo

c
c  main outer iteration loop for mixtures
c  loop for temperature, using Raoult's Law result (above) as first guess
c
 550    continue
c
        t=tk(j)
        do i=1,nc
          if (kph.eq.1) then
            xvap(i)=x2(i)
          else
            xliq(i)=x2(i)
          endif
        enddo
c
        x2sum=0.0d0
        ltpos=.false.         !flags for reguli-falsi iteration
        ltneg=.false.
        tneg=0.0d0
        tpos=0.0d0
        ftneg=0.0d0
        ftpos=0.0d0
        do ii=1,nc
          x2org(ii)=x2(ii)
        enddo
c       write (*,*) ' SATP--start main iteration; j,tk(j):  ',j,tk(j)
        tk(1)=tk(j)           !first guess for temperature from above
        j=1                   !reset iteration flag
        kguess=0              !for first calls to TPRHO
        do 400 itt=1,itmax
        if (int(itt/10)*10.eq.itt .and. itt.gt.30) then
          trtio1=trtio1**2
          trtio2=trtio2**2
        endif
        t=tk(j)
        lx2con=.false.        !flag for convergence of inner loop
c  compute density and fugacities for parent phase
c       write (*,*) ' SATP call TPRHO (parent) for it,T = ',itt,t
        call TPRHO (t,p,x,kph,kguess,rho1,ierr,herr2)   !parent phase
        if (itt.ne.1) then
          if (ABS(rho1-rhoc).lt.1 .and. ABS(rho2-rhoc).gt.4) then
c  jump away from critical point when a bad root was found
            rho1=rho1*1.5d0
            call TPRHO (t,p,x,kph,kguess,rho1,ierr,herr2)
          endif
        endif
        kgues1=1
        if (ierr.gt.0) then
          kgues1=0             !do not reuse faulty density as guess
          ierr=-145
          write (herr,1245) itt,herr2(1:146),hnull
 1245     format ('[SATP advisory -145] density iteration in ',
     &            'saturation routine did not converge for temperature',
     &            ' iteration',i3,'; ',a146,a1)
          tk(j)=tk(j)*trtio1  !try another temperature and use up one
          goto 400             !iteration (to prevent infinite loop)
        else if (rho1.gt.1.1d0*rhomax) then
c  density from TPRHO is out of range, reset
          rho1=rhomax
          ierr=-145
          write (herr,1245) itt,herr2(1:146),hnull
        end if
        call FGCTY2 (t,rho1,x,f1,ierr,herr)
c       write (*,1082) itt,t,p,rho1,(x(i),i=1,2),(f1(i),i=1,2)
c1082   format (1x,' SATP--phase 1:  ',i3,f11.5,2e14.6,2e16.8,2e18.10)
c
c  begin inner iteration loop for composition of phase 2
c
        do itx=1,itmax*2
c  compute density and fugacities for phase 2
          call TPRHO (t,p,x2,kph2,kguess,rho2,ierr,herr2) !incipient phase
          kgues2=1
          if (ierr.gt.0 .or. ABS(rho2-rho1).lt.0.01d0 .or.
     &       (rho2.lt..01d0 .and. kph2.eq.1)) then
            kgues2=0             !do not reuse faulty density as guess
            ierr=-146
            write (herr,1246) itx,herr2(1:146),hnull
 1246       format ('[SATP advisory -146] density iteration in ',
     &            'saturation routine did not converge for composition',
     &            ' iteration',i3,'; ',a146,a1)
            do ii=1,nc
              x2(ii)=x2org(ii)
            enddo
            kguess=0
            tk(j)=tk(j)/trtio2   !try another pressure and use up one
            goto 400             !iteration (to prevent infinite loop)
          else if (rho2.gt.1.2*rhomax) then
c  density from TPRHO is out of range, reset
            rho2=rhomax
            ierr=-146
            write (herr,1246) itt,herr2(1:146),hnull
          end if
          call FGCTY2 (t,rho2,x2,f2,ierr,herr)
c       write (*,1086) itx,rho2,(x2(i),i=1,2),(f2(i),i=1,2)
c1086   format (1x,'       phase 2:  ',i3,25x,e14.6,2e16.8,2e18.10)
c  calculate new x2's by ratio of fugacities; inner loop has converged
c  when x2's change by less than a convergence tolerance
          x2sum=0.0d0
          do i=1,nc
            if (x2(i).gt.0.d0) then
            if (f2(i).gt.1.0d-20 .and. f1(i).gt.0.0d0) then
              x2new(i)=x2(i)*f1(i)/f2(i)
            else
c  in case fugacity is zero; e.g., if x(i)=0 x2(i) does not change
              x2new(i)=x2(i)+0.0001d0 !add a little to keep x2sum<>1
            end if
            x2sum=x2sum+x2new(i)
            endif
          enddo
c  normalize the x2 compositions; this yields next guess for x2 and
c  ensures that the x2 always sum to one
          sumdel=0.0d0
          do i=1,nc
            if (x2(i).gt.0.d0) then
            x2new(i)=x2new(i)/x2sum
            sumdel=sumdel+abs(x2(i)-x2new(i))  !change in compositions
            x2(i)=x2new(i)
            endif
          enddo
          if (sumdel.lt.tolr) then
c  inner iteration loop has converged
            lx2con=.true.
            goto 700
          end if
c  if not, continue inner iteration loop
        enddo
c  inner iteration loop has not converged
        ierr=-147
        write (herr,1147) p/1000.0d0,sumdel,hnull
 1147   format ('[SATP advisory -147] ',
     &          'iteration for composition in saturation routine ',
     &          'did not converge; P =',g11.5,
     &          ' MPa; deltaX =',g11.5,' mol frac.',a1)
        call ERRMSG (ierr,herr)
c
c  end of inner (x2) iteration loop
c
 700    continue
        ft(j)=1.0d0-x2sum
c  outer (temperature) loop has converged when the x2's sum to one, i.e.,
c  when the fugacities of each component in each phase are equal
        if (abs(ft(j)).lt.tolr .and. ierr.ne.-147) then
          goto 850
        else
c  provided that the inner loop has converged, update positive and
c  negative bounds on pressure for possible use in reguli-falsi iteration
          if (lx2con) then
            if (ft(j).lt.0.0d0) then
              ltneg=.true.
              tneg=tk(j)
              ftneg=ft(j)
            else
              ltpos=.true.
              tpos=tk(j)
              ftpos=ft(j)
            end if
          end if
        end if
c
c  compute new guess for saturation temperature
c
        if (j.eq.1) then
c  for first iteration, new temperature is ratio of old
          j=2
          kguess=1    !use previous density as initial guess to TPRHO
c  ratio for next guess of temperature; the 0.10d0 is adjustable
          tnew=1.0d0+0.10d0*(x2sum-1.0d0)
c         write (*,*) ' PSAT--x2sum,tnew:  ',x2sum,tnew
          if (kph.eq.1) then
c  bubble point
            tk(2)=tk(1)/tnew
          else
c  dew point
            tk(2)=tk(1)*tnew
          end if
        else
c  subsequent iterations--use secant method, check for divide by zero
          if (ABS(ft(2)-ft(1)).lt.1.0d-10) then
            tk(3)=0.5*(tk(1)+tk(2))
          else
            tk(3)=tk(2)-ft(2)*(tk(2)-tk(1))/(ft(2)-ft(1))
          end if
c  check that new temperature is not outside bounds, if so use reguli-falsi
          if (ltneg .and. ltpos .and. (tk(3).gt.MAX(tpos,tneg)
     &        .or. tk(3).lt.MIN(tpos,tneg))) then
            tk(3)=tpos-ftpos*(tpos-tneg)/(ftpos-ftneg)
          end if
          if (tk(3).lt.0) tk(3)=tk(2)*.95d0
          if (tk(3).gt.tc*1.5d0) tk(3)=tk(2)*1.05d0
c         write (*,1137) tk(1),tk(2),tk(3),yyyyy,ft(2)
c1137     format (1x,' SATP--tguess_1,2,3; ft_1,2:  ',5e14.6)
c  discard oldest iteration
          tk(1)=tk(2)
          tk(2)=tk(3)
          ft(1)=ft(2)
        end if
 400    continue
c  outer iteration loop has not converged
        ierr=148
        write (herr,1148) p/1000.0d0,hnull
 1148   format ('[SATP error 148] ',
     &          'iteration for saturation state did not converge; P =',
     &          g11.5,' MPa.',a1)
        call ERRMSG (ierr,herr)
c
c  end of outer (temperature) iteration loop
c
 850    continue
c
c  assign final compositions and densities for parent and incipient
c  phases (x and x2, rho1 and rho2, respectively) to outputs
c
        t=tk(j)
        if (kph.eq.1) then
c  bubble point
          rhol=rho1
          rhov=rho2
          do i=1,nc
            xliq(i)=x(i)
            xvap(i)=x2(i)
          enddo
        else
c  dew point
          rhol=rho2
          rhov=rho1
          do i=1,nc
            xliq(i)=x2(i)
            xvap(i)=x(i)
          enddo
        end if
        call DPDD (t,rhol,xliq,dpdrh1)
        call DPDD (t,rhov,xvap,dpdrh2)
        if (dpdrh1.gt.1.0d6 .or. dpdrh2.gt.1.0d6) then
          ierr=148
          write (herr,1148) p/1000.d0,hnull
          call ERRMSG (ierr,herr)
        endif
        if ((rhol.gt.rhoc .and. rhov.gt.rhoc) .or.
     &      (rhol.lt.rhoc .and. rhov.lt.rhoc)) then
c  for some mixtures, both rhol and rhov can be greater than rhoc, so
c  check for p close to pc first before returning error
          if (p.lt.0.95d0*pc .or. icomp.ne.0 .or. kph.eq.2) then
            ierr=148
            write (herr,1148) p/1000.d0,hnull
            call ERRMSG (ierr,herr)
          endif
        endif
c  call new routine SATTP if ierr>0 and attempt to get convergence
        if (ierr.gt.0 .and. icomp.eq.0) then
          call SATTP(t,p,x,kph+2,0,d,rhol,rhov,xliq,xvap,q,ierr,herr)
        endif
c
c  end of mixture iteration
      end if
c
c  save results
 900  continue
      if (ierr.eq.0) then
        tsavp=t
        psavp=p
        icsavp=icomp
        kphsvp=kph
        dlsavp=rhol
        dvsavp=rhov
        if (icomp.eq.0) then
          do i=1,nc
            xsavp(i)=x(i)
            xlsavp(i)=xliq(i)
            xvsavp(i)=xvap(i)
          enddo
        else
          xsavp(icomp)=1.d0
          xlsavp(icomp)=1.d0
          xvsavp(icomp)=1.d0
        endif
        lsatp=.true.
      endif
      if (iaga.eq.1) then
        heos='AGA'
        if (ierr.eq.0) call TPRHO (t,p,xvap,2,1,rhov,ierr,herr)
      endif
      RETURN
c
      end                                               !subroutine SATP
c
c ======================================================================
c
      subroutine SATPEST (p,x,kph,t,x2,ierr,herr)
c
c  estimate initial values for saturation states given pressure
c  and the composition of one phase
c
c  inputs:
c        p--pressure [kPa]
c        x--composition [array of mol frac] (phase specified by kph)
c      kph--phase flag: 1 = input x is liquid composition (bubble point)
c                       2 = input x is vapor composition (dew point)
c  outputs:
c        t--estimated temperature [K]
c       x2--estimated composition of unknown phase [array of mol frac]
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  08-06-09 EWL, original version, taken from code in SATP
c  08-06-09 EWL, remove the check on pcomp(i)<10.  The impact of this may not be good for all situations.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SATPEST
c     dll_export SATPEST
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      dimension x(ncmax),x2(ncmax),pcomp(ncmax)
      dimension tk(4),ft(3)                  !used for T iteration
      character*255 herr
      character*1 htab,hnull
      common /HCHAR/ htab,hnull
      common /NCOMP/ nc,ic
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
      ltpos=.false. !flags for reguli-falsi iteration
      ltneg=.false.
      tneg=0.d0
      tpos=0.d0
      ftneg=0.d0
      ftpos=0.d0
      tolr=1.0d-6
      itmax=100
      j=1                     !flag for first iteration
c
c  generate initial guess using Raoult's law,
c  first generate initial guess for temperature by simple ratio of
c  pure component saturation temperatures
      tk(1)=0.d0
      acf=0.d0              !average acentric factor
      do i=1,nc
        ti=tcrit(i)/(1.0-0.428571*log10(p/pcrit(i))/(1.0+accen(i)))
c         write (*,*) ' SATP:  tsat-guess for component ',i,ti
        tk(1)=tk(1)+x(i)*ti
        acf=acf+x(i)*accen(i)
      enddo
      tzero=tk(1)  !save in case Raoult iteration does not converge
c
c  secant method iteration to find t which satisfies Raoult's law
c
      do it=1,itmax
c  approximate pure component vapor pressures with acentric factor
      do i=1,nc
        pxp=-2.333333*(1.0+accen(i))*(tcrit(i)/tk(j)-1.0)
        if (ABS(pxp).lt.50) pcomp(i)=pcrit(i)*10.0**pxp
c       if (pcomp(i).lt.10.) pcomp(i)=10.d0
      enddo
c
      if (kph.eq.1) then
c  bubble point
        psum=0.d0
        do i=1,nc
          psum=psum+x(i)*pcomp(i)
        enddo
        ft(j)=1.d0-p/psum
        if (abs(ft(j)).lt.tolr) goto 500  !iteration has converged
c  update + & - bounds on temperature for possible use in reguli-falsi
        if (ft(j).lt.0.d0) then
          ltneg=.true.
          tneg=tk(j)
          ftneg=ft(j)
        else
          ltpos=.true.
          tpos=tk(j)
          ftpos=ft(j)
        end if
c  generate next guess
        if (it.eq.1) then
          tratio=1.d0-0.42857d0*log10(psum/p)/(1.d0+acf)
          tk(2)=tk(1)*tratio
          j=2
        else if (it.ge.itmax/3 .and. ltpos .and. ltneg) then
c  if iteration has not converged after many iterations, use bisection
c  (provided that guesses bounding the root are available)
          tk(2)=0.5d0*(tpos+tneg)
          tk(1)=tk(2)
          ft(1)=ft(j)
          j=2
        else
c  use secant method
          if (ABS(ft(2)-ft(1)).gt.1.0d-10) then
            tk(3)=tk(2)-ft(2)*(tk(2)-tk(1))/(ft(2)-ft(1))
          else
            tk(3)=0.5d0*(tk(1)+tk(2))
          end if
c  check that new temperature is not outside bounds, if so use reguli-falsi
          if (ltneg .and. ltpos .and. (tk(j+1).gt.MAX(tpos,tneg)
     &        .or. tk(j+1).lt.MIN(tpos,tneg))) then
            tk(j+1)=tpos-ftpos*(tpos-tneg)/(ftpos-ftneg)
          end if
          tk(1)=tk(2)
          tk(2)=tk(3)
          ft(1)=ft(2)
        end if
      else
c  dew point
        xdamp=1.d0     !damping ratio for secant method
        ypsum=0.d0
        do i=1,nc
          ypsum=ypsum+x(i)*p/pcomp(i)
          ft(j)=1.0-ypsum
        enddo
        if (ABS(ft(j)).lt.1.0d3*tolr) goto 500  !iteration has converged
c  update + & - bounds on pressure for possible use in reguli-falsi
        if (ft(j).lt.0.d0) then
          ltneg=.true.
          tneg=tk(j)
          ftneg=ft(j)
        else
          ltpos=.true.
          tpos=tk(j)
          ftpos=ft(j)
        end if
        if (it.eq.1) then
          tratio=1.d0/(1.0-0.42857*log10(ypsum)/(1.d0+acf))
          tk(2)=tk(1)*tratio
          j=2
        else if (it.eq.2) then
c  secant method for 2nd guess
            if (ABS(ft(2)-ft(1)).gt.1.0d-10) then
              tk(3)=tk(2)-xdamp*ft(2)*(tk(2)-tk(1))/(ft(2)-ft(1))
            else
              tk(3)=0.5d0*(tk(1)+tk(2))
            end if
c  check that new temperature is not outside bounds, if so use reguli-falsi
            if (ltneg .and. ltpos .and. (tk(j+1).gt.MAX(tpos,tneg)
     &          .or. tk(j+1).lt.MIN(tpos,tneg))) then
              tk(j+1)=tpos-ftpos*(tpos-tneg)/(ftpos-ftneg)
            end if
            j=3
        else if (it.le.itmax/3 .or. .not.(ltpos.and.ltneg)) then
c  2nd order secant (inverse quadratic interpolation) for subsequent guesses
c  see Numerical Recipes, p 252
          rr=ft(3)/ft(2)  !these are the R,S,T used in Num. Rec.
          rs=ft(3)/ft(1)
          rt=ft(1)/ft(2)
          rrst=(rt-1.d0)*(rr-1.d0)*(rs-1.d0)
          if (abs(rrst).gt.1.d-20) then
            tk(4)=tk(3)
     &       +rs*(rt*(rr-rt)*(tk(2)-tk(3))-(1.d0-rr)*(tk(3)-tk(1)))/rrst
          else
            tk(j+1)=-1.d0    !Fix so that reguli-falsi will kick in
          endif
c  check that new temperature is not outside bounds, if so use reguli-falsi
          if (ltneg .and. ltpos .and. (tk(j+1).gt.MAX(tpos,tneg)
     &        .or. tk(j+1).lt.MIN(tpos,tneg))) then
            tk(j+1)=tpos-ftpos*(tpos-tneg)/(ftpos-ftneg)
          end if
c  discard oldest iteration
          tk(1)=tk(2)
          tk(2)=tk(3)
          tk(3)=tk(4)
          ft(1)=ft(2)
          ft(2)=ft(3)
        else
c  if iteration has not converged by now, use bisection
          tk(2)=0.5d0*(tpos+tneg)
          tk(1)=tk(2)
          ft(1)=ft(j)
          j=2
        end if
      end if
c       write (*,1006) it,j,tk(j-1),tk(j),ft(j-1)
c1006   format (1x,' SATP Raoult''s: it,j,t1,t2,ft:',2i4,2f12.4,e16.7)
      enddo                                     !next trial for t
c
c  iteration has not converged, issue warning and proceed
      ierr=-144
      herr='[SATP advisory -144] Raoult''s law iteration (to '//
     &     'generate mixture initial guess) has not converged.'//
     &     hnull
      call ERRMSG (ierr,herr)
      if (ABS(ft(j)).lt.1.0d4*tolr .or. ABS(tk(j)-tk(j-1)).lt.1.d0) then
      else
c  if current guess is not even close go back to initial guess
        tk(j)=tzero
      end if
c
 500  continue

      psum=0.d0
      if (kph.eq.1) then
c  bubble point
        do i=1,nc
          psum=psum+x(i)*pcomp(i)
        enddo
        do i=1,nc
          x2(i)=x(i)*pcomp(i)/psum
        enddo
      else
c  dew point
        ysum=0.d0
        do i=1,nc
          ysum=ysum+x(i)/pcomp(i)
        enddo
        do i=1,nc
          x2(i)=x(i)/pcomp(i)/ysum
          psum=psum+x2(i)*pcomp(i)
        enddo
      end if
      t=tk(j)
      RETURN
c
      end                                            !subroutine SATPEST
c
c ======================================================================
c
      subroutine SATD (rho,x,kph,kr,t,p,rhol,rhov,xliq,xvap,ierr,herr)
c
c  iterate for temperature and pressure given a density along the
c  saturation boundary and the composition
c
c  inputs:
c      rho--molar density [mol/L]
c        x--composition [array of mol frac]
c      kph--flag specifying desired root for multi-valued inputs
c           has meaning only for water at temperatures close to its triple point
c          -1 = return middle root (between 0 and 4 C)
c           1 = return highest temperature root (above 4 C)
c           3 = return lowest temperature root (along freezing line)
c  outputs:
c        t--temperature [K]
c        p--pressure [kPa]
c     rhol--molar density [mol/L] of saturated liquid
c     rhov--molar density [mol/L] of saturated vapor
c     xliq--liquid phase composition [array of mol frac]
c     xvap--vapor phase composition [array of mol frac]
c       kr--phase flag: 1 = input state is liquid
c                       2 = input state is vapor in equilibrium with liq
c                       3 = input state is liquid in equilibrium with solid
c                       4 = input state is vapor in equilibrium with solid
c     ierr--error flag:   0 = successful
c                         2 = D > Dmax
c                         8 = x out of range
c                        10 = D and x out of range
c                       160 = CRITP did not converge
c                       161 = SATD did not converge
c     herr--error string (character*255 variable if ierr<>0)
c
c  N.B. kr = 3,4 presently working only for pure components
c
c  either (rhol,xliq) or (rhov,xvap) will correspond to the input state
c  with the other pair corresponding to the other phase in equilibrium
c  with the input state
c
c  written by M. McLinden, NIST Phys & Chem Properties Div, Boulder, CO
c  11-22-99  MM, original version
c  02-08-00 EWL, add pure fluid version using Maxwell criteria
c  02-28-00 EWL, add checks for water
c  08-13-02 EWL, add check for ierr<>0 in mixture routine
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SATD
c     dll_export SATD
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1
      character*3 hsubl,hsublk
      character*3 hmelt,hmeltk
      character*12 hcasn
      character*3 hpheq,heos,hmxeos,hmodcp
      dimension x(ncmax),xliq(ncmax),xvap(ncmax)
      dimension tt(3),ft(2)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /FXPNT/ d72l,d72v,thmax,hmax,htpl,htpv,
     &               stpl,stpv,tsmax,smax,tsmin,smin,tsminm,sminm
      common /MELTMOD/ hmelt,hmeltk(n0:nx)
      common /SUBLMOD/ hsubl,hsublk(n0:nx)
      common /CCAS/ hcasn(n0:nx)
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /prnterr/ iprnterr
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      iaga=0
      call ISPURE (x,icomp)
      ierr=0
      herr=' '
      t0=300.0d0
      p0=100.0d0
      call LIMITX ('EOS',t0,rho,p0,x,tmin,tmax,Dmax,pmax,ierr,herr1)
c  [don't care if t0 or p0 are out of bounds, only rho and x]
      if (ierr.ge.8 .or. rho.gt.Dmax .or. rho.le.0.0d0) then
        write (herr,1001) ierr,herr1(1:235),hnull
 1001   format ('[SATD error',i3,'] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
      call CRITP (x,tc,pc,rhoc,ierr,herr1)
      t=tc
      p=pc
      rhol=rhoc
      rhov=rhoc
      if (ierr.gt.0) then
c  error condition--set outputs, issue warning, and return
        ierr=160
        write (herr,1160) herr1(1:237),hnull
 1160   format ('[SATD error 160] ',a237,a1)
        call ERRMSG (ierr,herr1)
        RETURN
      end if
      if (ABS(rhoc-rho).lt.0.001d0 .and. icomp.ne.0) then
        kr=1
        RETURN
      endif
c  calculate density at triple point if not set
      if (icomp.ne.0 .and. dtpv(icomp).lt.1.d-15) THEN
        call SATT (ttp(icomp),x,2,p1,rhol,rhov,xliq,xvap,ierr,herr)
        dtpv(icomp)=rhov
      endif
c  check if pure water
      iw=0
      rhow=55.504316185d0  !maximum liquid density from Pruss Eq.
      if (hcasn(icomp).eq.'7732-18-5' .and. icomp.ne.0) then
        iw=1
        if (rho.lt.dtp(icomp) .and. (kph.eq.-1 .or. kph.eq.3)) goto 210
        if (rho.gt.rhow .and. (kph.eq.-1 .or. kph.eq.1)) goto 210
      endif
c  determine region
      if (rho.lt.dtpv(icomp) .and. icomp.ne.0) then
        kr=4                         !vapor/solid
      elseif (rho.lt.rhoc) then
        kr=2                         !vapor
      elseif (iw.eq.1 .and. rho.lt.rhow .and. kph.ne.3) then   !water
        kr=1
      elseif (rho.gt.dtp(icomp) .and. abs(dtp(icomp)).gt.1.d-20
     &       .and. icomp.ne.0) then
        kr=3                         !liquid/solid
      else
        kr=1                         !liquid
      end if
c     write (*,*) 'SATD--initial phase:  ',kr
c
      if (heos.eq.'AGA') then
        heos='HMX'
        iaga=1
      endif
c
c  pure fluid VLE iteration
      if (icomp.ne.0 .and. kr.le.2 .and. kph.ne.-1
     &               .and. ianc(icomp).eq.0) then
        tol=1.0d-6
        itmax=20
c  estimate liquid density for any fluid
        if (d72l.lt.1.d-8) then
          t=0.72d0*tc
          call SATT (t,x,1,p,rhol,rhov,xliq,xvap,ierr,herr)
          d72l=rhol
          if (ianc(icomp).eq.1)
     &      call SATT (t,x,2,p,rhol,rhov,xliq,xvap,ierr,herr)
          d72v=rhov
        endif
        if (kr.eq.1) then
          d=rho/rhoc+(1.2d0-d72l/rhoc)*(rho-rhoc)/(d72l-rhoc)-1.d0
          t=9.272d0*LOG(d+1.0d0)-5.195d0*d
          t=tc*(1.0d0-t**(1.d0/0.337d0))
        elseif (kr.eq.2) then
          t=tc
          if (d72v.gt.0) then
            dx=LOG(d72v/rhoc)
            del=LOG(rho/rhoc)
            d=(3.5d0+dx)*(del/dx)**1.3d0-del
            t=0.5d0*LOG(d+1.0d0)-0.026d0*d**0.5D0
          endif
          if (t.lt.0) t=0
          t=tc/(1.d0+t**(1.d0/0.37d0))
        endif
c
        d1=0.001d0
        t2=t+0.00001d0
        stp1=1.00001d0
        stp2=1.00001d0
        do it=1,itmax
c  calculate pressure, if less than zero, increase t and try again
          call PRESS (t,rho,x,p)
          if (p.lt.0) then
            t=t*stp1
            stp1=1.0d0+(stp1-1.0d0)*1.5d0
            goto 140
          endif
          iii=iprnterr
          iprnterr=0
          call TPRHO (t,p,x,3-kr,0,rho2,ierr,herr)
          iprnterr=iii
c  calculate pressure, if error, decrease t and try again
          if (ierr.ne.0) then
            t=t/stp2
            stp2=1.0d0+(stp2-1.0d0)*2d0
            goto 140
          endif
          stp1=1.00001d0
          stp2=1.00001d0
c  use Maxwell criterion to generate next guess for t
          if (ABS(rho-rho2).gt.1d-11) then
            if (kr.eq.1) then
              rhol=rho
              rhov=rho2
            else
              rhol=rho2
              rhov=rho
            endif
            call AG (t,rho,x,a,g1)
            call AG (t,rho2,x,a,g2)
            call ENTRO (t,rho,x,s1)
            call ENTRO (t,rho2,x,s2)
            d2=d1
            if (ABS(s1-s2).gt.1d-11) d1=(g1-g2)/(s1-s2)
            if (ABS(d1-d2).gt.1d-11) f=-(t2-t)/(d2-d1)
            t2=t
            if (ABS(d1).le.tol) goto 900
            if (t+f*d1.le.0) then
              t=t*0.95d0
            else
              t=t+f*d1
            endif
          else
            t=t/stp1
          endif
        enddo
 140    continue
c  method failed (generally at very low temperatures in the liquid), try
c  alternate method:
c       write (*,*) 'Maxwell solution failed'
      endif
c
c  iterate for temperature using a combination of Newton's method
c  and reguli-falsi
c
      if (icomp.ne.0) then
        if (hmeltk(icomp).eq.'NBS' .and. rho.gt.dtp(icomp)) then !no melt line
          ierr=2
          write (herr,1002) ierr,rho,dtp(icomp),hnull
 1002     format('[SATD error',i3,'] density above triple-point density'
     &          ,'; D =',g11.5,' mol/L; Dtp =',g11.5,' mol/L',a1)
          call ERRMSG (ierr,herr)
          goto 900
        endif
        if (hsublk(icomp).eq.'NBS' .and. rho.lt.dtpv(icomp)) then !no subl line
          ierr=2
          write (herr,1003) ierr,rho,dtpv(icomp),hnull
 1003     format('[SATD error',i3,'] density below triple-point density'
     &          ,'; D =',g11.5,' mol/L; Dtp =',g11.5,' mol/L',a1)
          call ERRMSG (ierr,herr)
          goto 900
        endif
      end if
c
      tol=1.0d-7
      itmax=20
      rhowm=rhol
      if (kr.eq.3) then
        tt(1)=ttp(icomp)
        rhol=dtp(icomp)
        if (iw.eq.1) then !check for water
c  calculate p and d at slightly higher than the lowest possible temperature
          tt(1)=251.1650000001d0
          call MLTH2O(tt(1),p,p2)
        else
c  Call melting routine in case dtp(1) is not exactly eq. to rho(ttrp)
          call MELTT (tt(1),x,p,ierr,herr)
        endif
        call TPRHO (tt(1),p,x,1,1,rhol,ierr,herr)
        ft(1)=log(rho/rhol)
      elseif (kr.eq.4) then
        tt(1)=ttp(icomp)
        ft(1)=log(rho/dtpv(icomp))
      else
        tt(1)=tc
        ft(1)=log(rho/rhoc)
        if (iw.eq.1) then
          if (rho.ge.dtp(icomp) .and. kph.eq.-1) then
            tt(1)=ttp(icomp)
            ft(1)=log(dtp(icomp)/rhoc)
          endif
        endif
      endif
c
c  generate second guess for temperature
c
      if (kr.eq.1) then
        tt(2)=0.85*tc
        if (iw.eq.1 .and. kph.eq.-1 .and. rho.ge.dtp(icomp)) tt(2)=275
      else if (kr.eq.2) then
        tt(2)=0.75*tc
      else if (kr.eq.3) then
        tt(2)=0.85*tc
        if (iw.eq.1) tt(2)=273.16d0
      else if (kr.eq.4) then
        tt(2)=0.95*ttp(icomp)
      end if
c
c  initialize iteration flags
c
      lneg=.false.
      lpos=.false.
      ltp=.false.      !flag indicating if SATT called at triple point
      tneg=0.0d0
      tpos=0.0d0
c  store variables for reguli-falsi
      if (ft(1).lt.0.0) then
        lneg=.true.
        tneg=tt(1)
c       ftneg=yyyyy
      else
        lpos=.true.
        tpos=tt(1)
c       ftpos=yyyyy
      end if
      jt=2
c
      do 200 it=1,itmax
      ierr=0
      if (kr.eq.1) then           !liquid
        call SATT (tt(jt),x,1,p,rhol,rhov,xliq,xvap,ierr,herr)
        if (rhol.gt.0.0d0) ft(jt)=log(rho/rhol)
      else if (kr.eq.2 .or. kr.eq.4) then      !vapor
        call SATT (tt(jt),x,kr,p,rhol,rhov,xliq,xvap,ierr,herr)
        if (rhov.gt.0) ft(jt)=log(rho/rhov)
      else if (kr.eq.3) then      !liquid/solid
        if (iw.eq.1) then !check for water
          call MLTH2O(tt(jt),p,p2)
          if (rho.lt.rhowm) p=p2
        else
          call MELTT (tt(jt),x,p,ierr,herr)
        endif
        kguess=1
        rhol=dtp(icomp)
        call TPRHO (tt(jt),p,x,1,kguess,rhol,ierr,herr)
        if (rhol.gt.0) ft(jt)=log(rho/rhol)
      end if
      if (ierr.ne.0) then
        tt(jt)=0.999d0*tt(jt)
        goto 200
      endif
c     write (*,1999) it,jt,tt(jt),ft(jt)
c1999 format (1x,'% it,jt,t,ft: ',2i4,f9.4,f16.10)
c  check for convergence
      t=tt(jt)
      if (abs(ft(jt)).lt.tol) goto 900
c  store variables for reguli-falsi
      if (ft(jt).lt.0.0) then
        lneg=.true.
        tneg=tt(jt)
      else
        lpos=.true.
        tpos=tt(jt)
      end if
c  compute next guess for temperature
      if (jt.le.1) then
        jt=2
      else
        if (abs(ft(2)-ft(1)).gt.1.d-20)
     &  tt(3)=tt(2)-ft(2)*(tt(2)-tt(1))/(ft(2)-ft(1))
c  check if new guess is outside of bounds of previous guesses,
c  if so (and if upper and lower bounds are available) use
c  reguli-falsi
        if (lpos .and. lneg .and.
     &     (tt(3).gt.max(tpos,tneg) .or. tt(3).lt.min(tpos,tneg)))then
           tt(3)=0.5*(tpos+tneg)
c  check against triple point
        elseif (tt(3).lt.ttp(icomp) .and. icomp.ne.0 .and. iw.eq.0) then
          if (.not.ltp) then
            ltp=.true.
            tt(3)=ttp(icomp)
          else if (tt(3).le.0.0d0) then
            tt(3)=0.5d0*tt(2)
          end if
        end if
        tt(1)=tt(2)
        tt(2)=tt(3)
        ft(1)=ft(2)
      end if
 200  continue
 210  continue
      p=pc
      t=tc
      rhol=rhoc
      rhov=rhoc
      ierr=161
      write (herr,1161) rho,hnull
 1161 format ('[SATD error 161] ',
     &        'iteration for saturation state given density did not ',
     &        'converge for D =',g11.5,' mol/L',a1)
      call ERRMSG (ierr,herr)
 900  continue
      if (iaga.eq.1) then
        heos='AGA'
        if (ierr.eq.0) call TPRHO (t,p,xvap,2,1,rhov,ierr,herr)
      endif
      RETURN
c
      end                                               !subroutine SATD
c
c ======================================================================
c
      subroutine SATH (h,x,kph,nroot,k1,t1,p1,d1,k2,t2,p2,d2,ierr,herr)
c
c  iterate for temperature, pressure, and density given enthalpy along
c  the saturation boundary and the composition
c
c  inputs:
c        h--molar enthalpy [J/mol]
c        x--composition [array of mol frac]
c      kph--flag specifying desired root
c           0 = return all roots along the liquid-vapor line
c           1 = return only liquid VLE root
c           2 = return only vapor VLE roots
c           3 = return liquid SLE root (melting line)
c           4 = return vapor SVE root (sublimation line)
c  outputs:
c    nroot--number of roots.  Value is set to one for kph=1,3,4 if ierr=0
c       k1--phase of first root (1-liquid, 2-vapor, 3-melt, 4-subl)
c       t1--temperature of first root [K]
c       p1--pressure of first root [kPa]
c       d1--molar density of first root [mol/L]
c       k2--phase of second root (1-liquid, 2-vapor, 3-melt, 4-subl)
c       t2--temperature of second root [K]
c       p2--pressure of second root [kPa]
c       d2--molar density of second root [mol/L]
c     ierr--error flag:   0 = successful
c                         2 = h < hmin
c                         4 = h > hmax
c                         8 = h > htrp (for subl input)
c                       160 = CRITP did not converge
c                       161 = SATH did not converge for one root
c                       162 = SATH did not converge for both roots
c     herr--error string (character*255 variable if ierr<>0)
c
c  The second root is always set as the root in the vapor at temperatures
c  below the maximum enthalpy on the vapor saturation line.  If kph is
c  set to 2, and only one root is found in the vapor (this occurs when h<hcrit)
c  the state point will be placed in k2,t2,p2,d2.  If kph=0 and this situation
c  occurred, the first root (k1,t1,p1,d1) would be in the liquid (k1=1, k2=2).
c
c  N.B. kph = 3,4 presently working only for pure components
c
c  written by E.W. Lemmon, NIST Phys & Chem Properties Div, Boulder, CO
c  02-22-00 EWL, original version
c  05-08-08 EWL, add iteration count to avoid endless loop
c  09-04-08 EWL, add check for water and adjust tmin to 273.16 K
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SATH
c     dll_export SATH
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1
      character*12 hcasn
      dimension x(ncmax),xliq(ncmax),xvap(ncmax)
      dimension tt(3),ft(2)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /FXPNT/ d72l,d72v,thmax,hmax,htpl,htpv,
     &               stpl,stpv,tsmax,smax,tsmin,smin,tsminm,sminm
      common /FXPNTU/ temax,emax,etpl,etpv
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /ESTATE/ ieflg
      common /CCAS/ hcasn(n0:nx)
      call ISPURE (x,icomp)
      ierr=0
      ierr1=0
      iflag=0
      herr=' '
      k1=0
      t1=0
      p1=0
      d1=0
      k2=0
      t2=0
      p2=0
      d2=0
      nroot=0
      hold=0
      call CRITP (x,tc,pc,rhoc,ierr,herr1)
      if (ierr.gt.0) then
c  error condition--set outputs, issue warning, and return
        ierr=160
        write (herr,1160) herr1(1:237),hnull
 1160   format ('[SATH error 160] ',a237,a1)
        call ERRMSG (ierr,herr1)
        RETURN
      end if
c
c  calculate enthalpy at triple point if not set
      call LIMITX ('EOS',300.0d0,.1d0,.1d0,x,tmin,t,d,p,ierr,herr)
      if (ABS(htpl).lt.1.d-15) THEN
c  check for water
        if (hcasn(icomp).eq.'7732-18-5' .and. icomp.ne.0) tmin=273.16d0
        call SATT (tmin,x,1,p,rhol,rhov,xliq,xvap,ierr,herr)
        call ENTHAL (tmin,rhol,x,htpl)
        call ENERGY (tmin,rhol,x,etpl)
        if (icomp.eq.0 .or. ianc(icomp).eq.1)
     &      call SATT(tmin,x,2,p,rhol,rhov,xliq,xvap,ierr,herr)
        call ENTHAL (tmin,rhov,x,htpv)
        call ENERGY (tmin,rhov,x,etpv)
      endif
      xtpl=htpl
      xtpv=htpv
      call ENTHAL(tc,rhoc,x,xc)
      if (ieflg.eq.1) then
        xtpl=etpl
        xtpv=etpv
        call ENERGY(tc,rhoc,x,xc)
      endif
c
c  find maximum enthalpy along the saturated vapor dome
      xmax=hmax
      if (ieflg.eq.1) xmax=emax
      if (ABS(xmax).lt.1.d-15) THEN
c  set up initial three points
        tt1=tc-tc/10
        call SATT (tt1,x,2,p,rhol,rhov,xliq,xvap,ierr,herr)
        call ENTHAL (tt1,rhov,x,h1)
        if (ieflg.eq.1) call ENERGY (tt1,rhov,x,h1)
        tt2=tt1-tc/10
        call SATT (tt2,x,2,p,rhol,rhov,xliq,xvap,ierr,herr)
        call ENTHAL (tt2,rhov,x,h2)
        if (ieflg.eq.1) call ENERGY (tt2,rhov,x,h2)
        tt3=tt2-tc/10
        call SATT (tt3,x,2,p,rhol,rhov,xliq,xvap,ierr,herr)
        call ENTHAL (tt3,rhov,x,h3)
        if (ieflg.eq.1) call ENERGY (tt3,rhov,x,h3)
        it=0
c  use quadratic solution to find next guess for Tmax
 140    continue
        it=it+1
        b1=(h2-h1)/(tt2-tt1)
        b2=((h3-h2)/(tt3-tt2)-(h2-h1)/(tt2-tt1))/(tt3-tt1)
        b3=b1-b2*tt1-b2*tt2
        txmax=-b3/2.0d0/b2
        if (txmax.gt.tc) txmax=(tc+tt1)/2.0d0
        if (txmax.lt.tmin) txmax=(tmin+tt1)/2.0d0
        call SATT (txmax,x,2,p,rhol,rhov,xliq,xvap,ierr,herr)
        if (ieflg.ne.1) then
          call ENTHAL (txmax,rhov,x,xmax)
          thmax=txmax
          hmax=xmax
        else
          call ENERGY (txmax,rhov,x,xmax)
          temax=txmax
          emax=xmax
        endif
c  discard a point and load in the new one
        if (txmax.lt.tt3) then
          tt1=tt2
          tt2=tt3
          tt3=txmax
          h1=h2
          h2=h3
          h3=xmax
        elseif (txmax.gt.tt1) then
          tt3=tt2
          tt2=tt1
          tt1=txmax
          h3=h2
          h2=h1
          h1=xmax
        elseif (txmax.gt.tt2) then
          tt1=txmax
          h1=xmax
        else
          tt3=txmax
          h3=xmax
        endif
        hdiff=hold-xmax
        hold=xmax
        if (ABS(hdiff).gt.1.d-9 .and. it.lt.50) goto 140  !check for convergence
      endif
c
c  determine region
      if (kph.eq.3 .and. icomp.ne.0) then
        kphs=3                         !liquid/solid
      else if (kph.eq.4 .and. icomp.ne.0) then
        kphs=4                         !vapor/solid
      else
        if (h.lt.xc) then
          kphs=1                         !liquid
        else
          kphs=2                         !vapor
        end if
      endif
c     write (*,*) 'SATH--initial phase:  ',kphs
c
      if (h.gt.xtpv .and. kphs.eq.4) then
        ierr=8
        write (herr,1008) ierr,h,xtpv,hnull
 1008   format('[SATH error',i3,'] enthalpy greater than triple point'
     &     ,' value; h =',g11.5,' J/mol, max =',g11.5,' J/mol',a1)
        call ERRMSG (ierr,herr)
        RETURN
      endif
      if (h.lt.xtpl .and. kphs.ne.4) then
        ierr=2
        write (herr,1002) ierr,h,xtpl,hnull
 1002   format('[SATH error',i3,'] enthalpy below triple-point '
     &      ,'value; h =',g11.5,' J/mol, min =',g11.5,' J/mol',a1)
        call ERRMSG (ierr,herr)
        RETURN
      endif
      if (h.gt.xmax .and. kphs.ne.3) then
        ierr=4
        write (herr,1004) ierr,h,xmax,hnull
 1004   format('[SATH error',i3,'] enthalpy greater than maximum '
     &      ,'value; h =',g11.5,' J/mol, max =',g11.5,' J/mol',a1)
        call ERRMSG (ierr,herr)
        RETURN
      endif
c
c  iterate for temperature using a combination of Newton's method
c  and reguli-falsi
c
      tol=1.0d-6
      itmax=20
      if (kphs.eq.3) then
        tt(1)=tmin
        ft(1)=xtpl-h
      elseif (kphs.eq.4) then
        tt(1)=tmin
        ft(1)=xtpv-h
      else
        tt(1)=tc
        ft(1)=xc-h
      endif
c
c  initialize iteration flags
c
 130  continue
      lneg=.false.
      lpos=.false.
      tpos=0
      tneg=0
c  store variables for reguli-falsi
      if (ft(1).lt.0.0) then
        lneg=.true.
        tneg=tt(1)
        if (kphs.eq.2) then
          lpos=.true.
          tpos=txmax
        endif
      else
        lpos=.true.
        tpos=tt(1)
        if (kphs.eq.2) then
          lneg=.true.
          tneg=txmax
        endif
      end if
c
c  generate second guess for temperature
c
      if (kphs.eq.1) then
        tt(2)=0.95*tc
      else if (kphs.eq.2) then
        tt(2)=(tneg+tpos)/2.0d0
      else if (kphs.eq.3) then
        tt(2)=0.85*tc
      else if (kphs.eq.4) then
        tt(2)=0.95*tmin
      end if
      jt=2
c
      do it=1,itmax
      call SATT (tt(jt),x,kphs,p,rho,rhov,xliq,xvap,ierr,herr)
      if (kphs.eq.2 .or. kphs.eq.4) rho=rhov
      call ENTHAL (tt(jt),rho,x,h1)
      if (ieflg.eq.1) call ENERGY (tt(jt),rho,x,h1)
      ft(jt)=h1-h
c     write (*,1999) it,jt,tt(jt),ft(jt)
c1999 format (1x,'% it,jt,t,ft: ',2i4,f9.4,f16.10)
c  check for convergence
      if (abs(ft(jt)).lt.tol) goto 160
c  store variables for reguli-falsi
      if (ft(jt).lt.0.0) then
        lneg=.true.
        tneg=tt(jt)
      else
        lpos=.true.
        tpos=tt(jt)
      end if
c  compute next guess for temperature
      if (abs(ft(2)-ft(1)).gt.1.d-20)
     &  tt(3)=tt(2)-ft(2)*(tt(2)-tt(1))/(ft(2)-ft(1))
c  check if new guess is outside of bounds of previous guesses,
c  if so (and if upper and lower bounds are available) use
c  reguli-falsi
      if (lpos .and. lneg .and.
     &  (tt(3).gt.max(tpos,tneg) .or. tt(3).lt.min(tpos,tneg))) then
        tt(3)=0.5*(tpos+tneg)
      elseif (tt(3).lt.tmin .and. kphs.ne.4) then
        tt(3)=tmin
      else if (tt(3).le.tmin/10.0d0 .and. icomp.ne.0) then !subl check
        goto 150
      endif
      tt(1)=tt(2)
      tt(2)=tt(3)
      ft(1)=ft(2)
      enddo
 150  continue
      if (ierr1.ne.161) then
        ierr=161
        ierr1=161
        write (herr,1161) h,nroot+1,hnull
 1161   format ('[SATH error 161] ',
     &         'iteration for saturation state did not ',
     &         'converge for h =',g11.5,' J/mol; root = ',i1,a1)
        herr1=herr
      else
        ierr=162
        write (herr,1162) h,hnull
 1162   format ('[SATH error 162] ',
     &         'iteration for saturation state did not ',
     &         'converge for both roots; h =',g11.5,' J/mol',a1)
      endif
      call ERRMSG (ierr,herr)
      tt(jt)=0
      p=0
      rho=0
 160  continue
      t=tt(jt)
      if (kphs.eq.kph .or. kph.eq.0) then
        nroot=nroot+1
        if (iflag.eq.0) then
c  first root found, save values
          t1=t
          p1=p
          d1=rho
          k1=kphs
        else
c  second root found, save values
          t2=t
          p2=p
          d2=rho
          k2=kphs
          if (ierr.eq.0 .and. ierr1.ne.0) then
            ierr=ierr1
            herr=herr1
          endif
        endif
      endif
c  check if second root might exist
      if (h.gt.xtpv .and. h.lt.xmax .and. iflag.eq.0) then
c  only calculate second root if requested
        if (kph.eq.0 .or. kph.eq.2) then
          tt(1)=tmin
          ft(1)=xtpv-h
          kphs=2
          iflag=1
          goto 130
        endif
      endif
      RETURN
c
      end                                               !subroutine SATH
c
c ======================================================================
c
      subroutine SATE (e,x,kph,nroot,k1,t1,p1,d1,k2,t2,p2,d2,ierr,herr)
c
c  iterate for temperature, pressure, and density given energy along
c  the saturation boundary and the composition
c
c  inputs:
c        e--molar energy [J/mol]
c        x--composition [array of mol frac]
c      kph--flag specifying desired root
c           0 = return all roots along the liquid-vapor line
c           1 = return only liquid VLE root
c           2 = return only vapor VLE roots
c           3 = return liquid SLE root (melting line)
c           4 = return vapor SVE root (sublimation line)
c  outputs:
c    see SATH for description of outputs
c
c  written by E.W. Lemmon, NIST Phys & Chem Properties Div, Boulder, CO
c  05-16-05 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SATE
c     dll_export SATE
c
      parameter (ncmax=20)        !max number of components in mixture
      character*255 herr
      dimension x(ncmax)
      common /ESTATE/ ieflg
      ieflg=1
      call SATH (e,x,kph,nroot,k1,t1,p1,d1,k2,t2,p2,d2,ierr,herr)
      ieflg=0
      if (ierr.ne.0) then
        i=index(herr,'SATH')
        if (i.gt.0) herr=herr(1:i+2)//'E'//herr(i+4:255)
        i=index(herr,'h =')
        if (i.gt.0) herr=herr(1:i-1)//'e'//herr(i+1:255)
        i=index(herr,'enthalpy')
        if (i.gt.0) herr=herr(1:i-1)//'energy'//herr(i+8:255)
      endif
      RETURN
c
      end                                               !subroutine SATE
c
c ======================================================================
c
      subroutine SATS (s,x,kph,nroot,k1,t1,p1,d1,k2,t2,p2,d2,
     &                 k3,t3,p3,d3,ierr,herr)
c
c  iterate for temperature, pressure, and density given an entropy along
c  the saturation boundary and the composition
c
c  inputs:
c        s--molar entropy [J/mol-K]
c        x--composition [array of mol frac]
c      kph--flag specifying desired root
c           0 = return all roots along the liquid-vapor line
c           1 = return only liquid VLE root
c           2 = return only vapor VLE roots
c           3 = return liquid SLE root (melting line)
c           4 = return vapor SVE root (sublimation line)
c  outputs:
c    nroot--number of roots.  Set to one for kph=1,3,4 if ierr=0
c       k1--phase of first root (1-liquid, 2-vapor, 3-melt, 4-subl)
c       t1--temperature of first root [K]
c       p1--pressure of first root [kPa]
c       dl--molar density of first root [mol/L]
c       k2--phase of second root (1-liquid, 2-vapor, 3-melt, 4-subl)
c       t2--temperature of second root [K]
c       p2--pressure of second root [kPa]
c       d2--molar density of second root [mol/L]
c       k3--phase of third root (1-liquid, 2-vapor, 3-melt, 4-subl)
c       t3--temperature of third root [K]
c       p3--pressure of third root [kPa]
c       d3--molar density of third root [mol/L]
c     ierr--error flag:   0 = successful
c                         1 = no roots found for specified input phase
c                         2 = s < smin
c                         4 = s > smax
c                         8 = s > strp (for subl input)
c                       160 = CRITP did not converge
c                       161 = SATS did not converge for one root
c                       162 = SATS did not converge for two roots
c                       163 = SATS did not converge for all roots
c     herr--error string (character*255 variable if ierr<>0)
c
c  The second root is always set as the root in the vapor at temperatures
c  below the maximum entropy on the vapor saturation line.  If kph is
c  set to 2, and only one root is found in the vapor (this occurs when s<scrit)
c  the state point will be placed in k2,t2,p2,d2.  If kph=0 and this situation
c  occurred, the first root (k1,t1,p1,d1) would be in the liquid (k1=1, k2=2).
c
c  The third root is the root with the lowest temperature.  For fluids
c  with multiple roots:  When only one root is found in the vapor phase
c  (this happens only at very low temperatures past the region where three
c  roots are located), the value of the root is still placed in
c  k3,t3,p3,d3.  For fluids that never have more than one root (when there
c  is no maximum entropy along the saturated vapor line), the value of the
c  root is always placed in k1,t1,p1,d1.
c
c  N.B. kph = 3,4 presently working only for pure components
c
c  written by E.W. Lemmon, NIST Phys & Chem Properties Div, Boulder, CO
c  02-24-00 EWL, original version
c  04-18-07 EWL, add check for SATT failure on the vapor side when finding stpv
c  03-04-10 EWL, add error message 1001 for root not found in specified phase
c  03-04-10 EWL, check for cases where melting line is not available
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SATS
c     dll_export SATS
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1
      dimension x(ncmax),xliq(ncmax),xvap(ncmax)
      dimension tt(3),ft(2)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /FXPNT/ d72l,d72v,thmax,hmax,htpl,htpv,
     &               stpl,stpv,tsmax,smax,tsmin,smin,tsminm,sminm
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      call ISPURE (x,icomp)
      herr=' '
      ierr=0
      ierr1=0
      iflag=0
      nroot=0
      k1=0
      t1=0.d0
      p1=0.d0
      d1=0.d0
      k2=0
      t2=0.d0
      p2=0.d0
      d2=0.d0
      k3=0
      t3=0.d0
      p3=0.d0
      d3=0.d0
      sold=0.d0
      call CRITP (x,tc,pc,rhoc,ierr,herr1)
      if (ierr.gt.0) then
c  error condition--set outputs, issue warning, and return
        ierr=160
        write (herr,1160) herr1(1:237),hnull
 1160   format ('[SATS error 160] ',a237,a1)
        call ERRMSG (ierr,herr1)
        RETURN
      end if
c
c  calculate entropy at triple point if not set
      call LIMITX ('EOS',300.0d0,.1d0,.1d0,x,tmin,t,d,p,ierr,herr)
      if (ABS(stpl).lt.1.d-15) THEN
        if (icomp.ne.0) then
          call INFO (icomp,wmm,ttrp,tnbpt,tq,pq,Dq,Zc,ac1,dip,Rg1)
          if (tmin.lt.ttrp) tmin=ttrp  !Check for water
        endif
        call SATT (tmin,x,1,p,rhol,rhov,xliq,xvap,ierr,herr)
        call ENTRO (tmin,rhol,x,stpl)
        rhov2=rhov
        if (icomp.eq.0 .or. ianc(icomp).eq.1)
     &        call SATT(tmin,x,2,p,rhol,rhov,xliq,xvap,ierr,herr)
        if (ierr.gt.0) rhov=rhov2
        call ENTRO (tmin,rhov,x,stpv)
      endif
      call ENTRO(tc,rhoc,x,sc)      !entropy at critical point
c
c  find maximum and minimum entropy along the saturated vapor dome, and
c  the minimum entropy along the melting line
      if (ABS(tsmax).lt.1.d-15) THEN
        ii=1
        tsminm=-1
        sminm=stpl
        tsmin=-1
        smin=sc
        tsmax=-1
        smax=stpv
c  check first to see if melting line is double valued
        p=0.d0
        call SATT (tmin,x,3,p,rhol,rhov,xliq,xvap,ierr,herr)
        if (ABS(p).gt.1.d-15) then
          call ENTRO (tmin,rhol,x,s1)
          call SATT (tmin+1.0d0,x,3,p,rhol,rhov,xliq,xvap,ierr,herr)
          call ENTRO (tmin+1.0d0,rhol,x,s2)
          if (s2.lt.s1) ii=0      !melting line is doubled valued if s2<s1
        endif
c  set up initial three points
        j=0
        do i=ii,2
          if (i.eq.1) then
            tt1=tc-tc/20
            tt2=tt1-tc/20
            tt3=tt2-tc/20
          else
            tt3=tmin+tmin/20
            tt2=tt3+tmin/20
            tt1=tt2+tmin/20
          endif
          k=2
          if (i.eq.0) k=3
          call SATT (tt1,x,k,p,rhol,rho,xliq,xvap,ierr,herr)
          if (i.eq.0) rho=rhol
          call ENTRO (tt1,rho,x,s1)
          call SATT (tt2,x,k,p,rhol,rho,xliq,xvap,ierr,herr)
          if (i.eq.0) rho=rhol
          call ENTRO (tt2,rho,x,s2)
          call SATT (tt3,x,k,p,rhol,rho,xliq,xvap,ierr,herr)
          if (i.eq.0) rho=rhol
          call ENTRO (tt3,rho,x,s3)
c  use quadratic solution to find next guess for Tsm
 140      continue
          b1=(s2-s1)/(tt2-tt1)
          b2=((s3-s2)/(tt3-tt2)-(s2-s1)/(tt2-tt1))/(tt3-tt1)
          b3=b1-b2*tt1-b2*tt2
          tsm=-b3/2.0d0/b2
          j=j+1
          if (tsm.lt.tmin .or. j.gt.100) then
            if (i.eq.0) goto 110
            goto 120         !no solution (VLE line is not double valued)
          endif
          if (tsm.gt.tc) tsm=(tc+tt1)/2.0d0
          call SATT (tsm,x,k,p,rhol,rho,xliq,xvap,ierr,herr)
          if (i.eq.0) rho=rhol
          call ENTRO (tsm,rho,x,sm)
c  discard a point and load in the new one
          if (tsm.lt.tt3) then
            tt1=tt2
            tt2=tt3
            tt3=tsm
            s1=s2
            s2=s3
            s3=sm
          elseif (tsm.gt.tt1) then
            tt3=tt2
            tt2=tt1
            tt1=tsm
            s3=s2
            s2=s1
            s1=sm
          elseif (tsm.gt.tt2) then
            tt1=tsm
            s1=sm
          else
            tt3=tsm
            s3=sm
          endif
          sdiff=sold-sm
          sold=sm
          if (ABS(sdiff).gt.1.d-9) goto 140  !check for convergence
c  load appropriate fixed points with the max/min value.
          if (i.eq.1) then
            smax=sm
            tsmax=tsm
          elseif (i.eq.2) then
            smin=sm
            tsmin=tsm
          elseif (i.eq.0) then
            sminm=sm
            tsminm=tsm
          endif
 110      continue
        enddo
      endif
c
c  determine region
 120  continue
      if (kph.eq.3 .and. icomp.ne.0) then
        kphs=3                         !liquid/solid
      else if (kph.eq.4 .and. icomp.ne.0) then
        kphs=4                         !vapor/solid
      else
        if (s.lt.sc) then
          kphs=1                         !liquid
        else
          kphs=2                         !vapor
        end if
      endif
c     write (*,*) 'SATS--initial phase:  ',kphs
c
      if (s.lt.stpv .and. kphs.eq.4) then
        ierr=8
        write (herr,1008) ierr,s,stpv,hnull
 1008   format('[SATS error',i3,'] entropy less than triple point'
     &    ,' entropy; s =',g11.5,' J/mol-K, smin =',g11.5,' J/mol-K',a1)
        call ERRMSG (ierr,herr)
        RETURN
      endif
      if (s.lt.stpl .and. (s.lt.sminm .or. kphs.ne.3)) then
        ierr=2
        smx=stpl
        if (sminm.lt.stpl .and. kphs.eq.3) smx=sminm
        write (herr,1002) ierr,s,smx,hnull
 1002   format('[SATS error',i3,'] entropy below minimum '
     &     ,'entropy; s =',g11.5,' J/mol-K, smin =',g11.5,' J/mol-K',a1)
        call ERRMSG (ierr,herr)
        RETURN
      endif
      if (s.gt.smax .and. s.gt.stpv .and. kphs.le.2) then
        ierr=4
        smx=smax
        if (stpv.gt.smax) smx=stpv
        write (herr,1004) ierr,s,smx,hnull
 1004   format('[SATS error',i3,'] entropy greater than maximum '
     &     ,'entropy; s =',g11.5,' J/mol-K, smax =',g11.5,' J/mol-K',a1)
        call ERRMSG (ierr,herr)
        RETURN
      endif
c
c  iterate for temperature using a combination of Newton's method
c  and reguli-falsi
c
      tol=1.0d-6
      itmax=20
      if (kphs.eq.3) then
        tt(1)=tmin
        ft(1)=stpl-s
        if (tsminm.gt.0) then
          tt(1)=tsminm
          ft(1)=sminm-s
        endif
      elseif (kphs.eq.4) then
        tt(1)=tmin
        ft(1)=stpv-s
      else
        tt(1)=tc
        ft(1)=sc-s
c  Check for situations where three roots exists in the vapor, but the current
c  input is beyond the three root region, and only one root exists.
        if (s.gt.smax .and. tsmax.gt.0) then
          tt(1)=tmin
          ft(1)=stpv-s
          iflag=2
          kphs=2
        endif
      endif
c
c  initialize iteration flags
c
 130  continue
      lneg=.false.
      lpos=.false.
      tpos=0
      tneg=0
c  store variables for reguli-falsi
      if (ft(1).lt.0.0) then
        lneg=.true.
        tneg=tt(1)
        if (kphs.eq.2 .and. tsmax.gt.0) then
          lpos=.true.
          tpos=tsmax
          if (iflag.eq.2) tpos=tsmin
        endif
        if (kphs.eq.3 .and. iflag.eq.1) then
          lpos=.true.
          tpos=tsminm
        endif
      else
        lpos=.true.
        tpos=tt(1)
        if (kphs.eq.2 .and. tsmax.gt.0) then
          lneg=.true.
          tneg=tsmax
          if (iflag.eq.2) tneg=tsmin
        endif
        if (kphs.eq.3 .and. iflag.eq.1) then
          lneg=.true.
          tneg=tsminm
        endif
      end if
c
c  generate second guess for temperature
c
      if (kphs.eq.1) then
        tt(2)=0.95*tc
      else if (kphs.eq.2 .and. tneg.gt.0 .and. tpos.gt.0) then
        tt(2)=(tneg+tpos)/2.0d0
      else if (kphs.eq.2) then
        tt(2)=0.85*tc
      else if (kphs.eq.3 .and. tneg.gt.0 .and. tpos.gt.0) then
        tt(2)=(tneg+tpos)/2.0d0
      else if (kphs.eq.3 .and. tsminm.gt.0) then
        tt(2)=1.05*tsminm
      else if (kphs.eq.3) then
        tt(2)=1.05*tmin
      else if (kphs.eq.4) then
        tt(2)=0.95*tmin
      end if
      jt=2
c
      do it=1,itmax
      call SATT (tt(jt),x,kphs,p,rho,rhov,xliq,xvap,ierr,herr)
      if (p.le.0.d0 .and. ierr.eq.0) goto 200
      if (kphs.eq.2 .or. kphs.eq.4) rho=rhov
      call ENTRO (tt(jt),rho,x,s1)
      ft(jt)=s1-s
c     write (*,1999) it,jt,tt(jt),ft(jt)
c1999 format (1x,'% it,jt,t,ft: ',2i4,f9.4,f16.10)
c  check for convergence
      if (abs(ft(jt)).lt.tol) goto 160
c  store variables for reguli-falsi
      if (ft(jt).lt.0.0) then
        lneg=.true.
        tneg=tt(jt)
      else
        lpos=.true.
        tpos=tt(jt)
      end if
c  compute next guess for temperature
      if (abs(ft(2)-ft(1)).gt.1.d-20)
     &  tt(3)=tt(2)-ft(2)*(tt(2)-tt(1))/(ft(2)-ft(1))
c  check if new guess is outside of bounds of previous guesses,
c  if so (and if upper and lower bounds are available) use
c  reguli-falsi
      if (lpos .and. lneg .and.
     &  (tt(3).gt.max(tpos,tneg) .or. tt(3).lt.min(tpos,tneg))) then
        tt(3)=0.5*(tpos+tneg)
      elseif (tt(3).lt.tmin .and. kphs.ne.4) then
        tt(3)=tmin
      else if (tt(3).le.tmin/10.0d0 .and. icomp.ne.0) then !subl check
        goto 150
      endif
      tt(1)=tt(2)
      tt(2)=tt(3)
      ft(1)=ft(2)
      enddo
 150  continue
      if (ierr1.ne.161 .and. ierr1.ne.162) then
        ierr=161
        ierr1=161
        write (herr,1161) s,nroot+1,hnull
 1161   format ('[SATS error 161] ',
     &         'iteration for saturation state given entropy did not ',
     &         'converge for s =',g11.5,' J/mol-K; root = ',i1,a1)
        herr1=herr
      elseif (ierr1.ne.162) then
        ierr=162
        write (herr,1162) s,hnull
 1162   format ('[SATS error 162] ',
     &         'iteration for saturation state given entropy did not ',
     &         'converge for two of the roots; s =',g11.5,' J/mol-K',a1)
      else
        ierr=163
        write (herr,1163) s,hnull
 1163   format ('[SATS error 163] ',
     &         'iteration for saturation state given entropy did not ',
     &         'converge for all roots; s =',g11.5,' J/mol-K',a1)
      endif
      call ERRMSG (ierr,herr)
      tt(jt)=0
      p=0
      rho=0
 160  continue
      t=tt(jt)
      if (kphs.eq.kph .or. kph.eq.0) then
        nroot=nroot+1
        if (iflag.eq.0) then
c  first root found, save values
          t1=t
          p1=p
          d1=rho
          k1=kphs
        elseif (iflag.eq.1) then
c  second root found, save values
          t2=t
          p2=p
          d2=rho
          k2=kphs
        elseif (iflag.eq.2) then
c  third root found, save values
          t3=t
          p3=p
          d3=rho
          k3=kphs
        endif
        if (ierr.eq.0 .and. ierr1.ne.0) then
          ierr=ierr1
          herr=herr1
        endif
      endif
c  check if second root might exist
c  only calculate second root if requested
      if ((kph.eq.0 .or. kph.eq.2) .and. tsmax.gt.0) then
        if (s.gt.smin .and. s.lt.smax .and. iflag.eq.0) then
          kphs=2
          iflag=iflag+1
          tt(1)=tsmin
          ft(1)=smin-s
          goto 130
        elseif (s.gt.smin .and. s.lt.smax .and. s.lt.stpv
     &         .and. iflag.eq.1) then
          kphs=2
          iflag=iflag+1
          tt(1)=tmin
          ft(1)=stpv-s
          goto 130
        endif
      endif
      if (kph.eq.3 .and. tsminm.gt.0) then
        if (s.gt.sminm .and. s.lt.stpl .and. iflag.eq.0) then
          iflag=iflag+1
          tt(1)=tmin
          ft(1)=smin-s
          goto 130
        endif
      endif
 200  if (k1.eq.0 .and. k2.eq.0 .and. k3.eq.0 .and. ierr.eq.0) then
        ierr=1
        write (herr,1001) ierr,s,hnull
 1001   format('[SATS error',i3,'] no roots found for specified phase;'
     &    ,' s =',g11.5,' J/mol-K', a1)
        call ERRMSG (ierr,herr)
      endif
      RETURN
c
      end                                               !subroutine SATS
c
c ======================================================================
c
      subroutine CSATK (icomp,t,kph,p,rho,csat,ierr,herr)
c
c  compute the heat capacity along the saturation line as a function of
c  temperature for a given component
c
c  csat can be calculated two different ways:
c     Csat = Cp - T(DvDT)(DPDTsat)
c     Csat = Cp - beta/rho*hvap/(vliq - vvap)
c     where beta is the volume expansivity
c
c  inputs:
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c        t--temperature [K]
c      kph--phase flag: 1 = liquid calculation
c                       2 = vapor calculation
c  outputs:
c        p--saturation pressure [kPa]
c      rho--saturation molar density [mol/L]
c     csat--saturation heat capacity [J/mol-K]
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  09-30-98 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: CSATK
c     dll_export CSATK
c
      parameter (ncmax=20)        !max number of components in mixture
      character*255 herr
      dimension x(ncmax),xliq(ncmax),xvap(ncmax)
      common /NCOMP/ nc,ic
c
      ic2=ic
      ic=icomp
      csat=0.d0
      call SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr)
      if (ierr.gt.0 .or. rhol.le.0.d0) then
        ic=ic2
        if (ierr.eq.0) ierr=1
        p=0.d0
        rho=0.d0
        RETURN
      end if
      rho=rhol
      if (kph.eq.2) rho=rhov
      call DPDT (t,rho,x,dpt)
      call DPDD (t,rho,x,dpdrho)
      call CVCPK (icomp,t,rho,cv,cp)
      call ENTHAL (t,rhol,x,hl)
      call ENTHAL (t,rhov,x,hv)
      beta=dpt/dpdrho/rho
      csat=cp-beta/rho*(hl-hv)/(1.d0/rhol-1.d0/rhov)
      ic=ic2
c
      RETURN
      end                                              !subroutine CSATK
c
c ======================================================================
c
      subroutine DPTSATK (icomp,t,kph,p,rho,csat,dpt,ierr,herr)
c
c  compute the heat capacity and dP/dT along the saturation line as a
c  function of temperature for a given component.  See also subroutine CSATK.
c
c  inputs:
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c        t--temperature [K]
c      kph--phase flag: 1 = liquid calculation
c                       2 = vapor calculation
c  outputs:
c        p--saturation pressure [kPa]
c      rho--saturation molar density [mol/L]
c     csat--saturation heat capacity [J/mol-K] (same as that called from CSATK)
c      dpt--dP/dT along the saturation line [kPa/K]
c           (this is not dP/dT "at" the saturation line for the single phase
c            state, but the change in saturated vapor pressure as the
c            saturation temperature changes.)
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  09-25-06 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: DPTSATK
c     dll_export DPTSATK
c
      parameter (ncmax=20)        !max number of components in mixture
      character*255 herr
      dimension x(ncmax),xliq(ncmax),xvap(ncmax)
      common /NCOMP/ nc,ic
c
      ic2=ic
      ic=icomp
      csat=0.d0
      dpt=0.d0
      call SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr)
      if (ierr.gt.0 .or. rhol.le.0.d0) then
        ic=ic2
        if (ierr.eq.0) ierr=1
        p=0.d0
        rho=0.d0
        RETURN
      end if
      rho=rhol
      if (kph.eq.2) rho=rhov
      call DPDT (t,rho,x,dpt)
      call DPDD (t,rho,x,dpdrho)
      call CVCPK (icomp,t,rho,cv,cp)
      call ENTHAL (t,rhol,x,hl)
      call ENTHAL (t,rhov,x,hv)
      beta=dpt/dpdrho/rho
      csat=cp-beta/rho*(hl-hv)/(1.d0/rhol-1.d0/rhov)
      dpt=             (hl-hv)/(1.d0/rhol-1.d0/rhov)/T
      ic=ic2
c
      RETURN
      end                                            !subroutine DPTSATK
c
c ======================================================================
c
      subroutine CV2PK (icomp,t,rho,cv2p,csat,ierr,herr)
c
c  compute the isochoric heat capacity in the two phase (liquid+vapor)
c  region
c
c  inputs:
c    icomp--component number in mixture (1..nc); 1 for pure fluid
c        t--temperature [K]
c      rho--density [mol/l] if known
c           If rho=0, then a saturated liquid state is assumed.
c
c  outputs:
c     cv2p--isochoric two-phase heat capacity [J/mol-K]
c     csat--saturation heat capacity [J/mol-K]
c           (Although there is already a csat routine in REFPROP,
c            it is also returned here.  However, the calculation
c            speed is slower than csat.)
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  03-30-05 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT::CV2PK
c     dll_export CV2PK
c
      parameter (ncmax=20)        !max number of components in mixture
      character*255 herr,herr1,herr2
      dimension x(ncmax),xliq(ncmax),xvap(ncmax)
      common /NCOMP/ nc,ic
c
      ierr=0
      herr=' '
      ic2=ic
      ic=icomp
      cv2p=0.d0
      dt=0.01d0
      i=1
      t1=t+dt
      t2=t-dt
      call SATT (t ,x,i,p ,dl ,dv ,xliq,xvap,ierr ,herr)
      call SATT (t1,x,i,p1,dl1,dv1,xliq,xvap,ierr1,herr1)
      call SATT (t2,x,i,p2,dl2,dv2,xliq,xvap,ierr2,herr2)
      if (rho.le.0.d0) rho=dl
      if (ierr.gt.0 .or. ierr1.gt.0 .or. ierr2.gt.0 .or. dl.le.0.d0)then
        if (ierr.eq.0 .and. ierr1.ne.0) then
          ierr=ierr1
          herr=herr1
        endif
        if (ierr.eq.0 .and. ierr2.ne.0) then
          ierr=ierr2
          herr=herr2
        endif
        if (ierr.eq.0 .and. dl.le.0.d0) ierr=1
        ic=ic2
        RETURN
      endif
c
      call DPDT (t,dl,x,dpt)
      call DPDD (t,dl,x,dpdrho)
      call CVCPK (icomp,t,dl,cv,cp)
      call ENTHAL (t,dl,x,hl)
      call ENTHAL (t,dv,x,hv)
      beta=dpt/dpdrho/dl
      dpdtsat=(hv-hl)/t/(1.d0/dv-1.d0/dl)     ! d(p)/d(T) at sat. liq.
      dddtsat=(dl1-dl2)/2.d0/dt               ! d(rho)/d(T) at sat. liq.
      d2pdtsat=(p2+p1-2.d0*p)/dt**2           ! d^2(p)/d(T)^2 at sat. liq.
      csat=cp-beta/dl*dpdtsat*t
      cv2p=csat+t/dl**2*dddtsat*dpdtsat+t*(1.d0/rho-1.d0/dl)*d2pdtsat
      ic=ic2
c
      RETURN
      end                                              !subroutine CV2PK
c
c ======================================================================
c
      subroutine TPRHOB (t,p,rho1,rho2,x,rho,ierr,herr)
c
c  iterate for density given temperature, pressure and an upper and lower
c  bound for the density.  This routine is only meant to replace TPRHO
c  in special cases near the critical point.  (And is only used by SATT.)
c
c  inputs:
c        t--temperature [K]
c        p--pressure [kPa]
c     rho1--first bound on density [mol/L]
c     rho2--second bound on density [mol/L]
c        x--composition [array of mol frac]
c  outputs:
c      rho--molar density [mol/L]
c     ierr--error flag:   0 = successful
c                         1 = failed to converge
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  11-08-01 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)        !max number of components in mixture
      character*255 herr
      dimension x(ncmax)

      ierr=0
      herr=' '
      tolr=1.d-9
      d1=rho1
      d2=rho2
      if (d1.gt.d2) then
        d1=rho2
        d2=rho1
      endif
      it=1
      rho=(d1+d2)/2.d0
 140  continue
      call PRESS (t,rho,x,p1)
      call DPDD (t,rho,x,dpd)
      rho0=rho
c  use false position to get next root
      if (abs(dpd).gt.1.d-20) rho=rho+(p-p1)/dpd
c  keep within bounds
      if (rho.lt.d1) rho=(rho0+d1)/2.d0
      if (rho.gt.d2) rho=(rho0+d2)/2.d0
      if (ABS(p1-p).lt.tolr) RETURN
      it=it+1
      if (it.eq.10) tolr=tolr*100
      if (it.eq.15) tolr=tolr*100
      if (it.lt.25) goto 140

      ierr=1
      RETURN
      end                                             !subroutine TPRHOB
c
c ======================================================================
c
      subroutine DLDV (t,p,rhol,rhov,xl,xv,ierr,herr)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)        !max number of components in mixture
      character*255 herr
      dimension xl(ncmax),xv(ncmax)

      tolr=1.0d-6
      deld=1.d0
      call CRITP (xl,tc,pc,rhoc,ierr,herr)
      rhol=rhoc*1.1d0
      rhov=rhoc
      it=0
 140  continue
      it=it+1

      it2=0
 110  continue
      call DPDD2 (t,rhol,xl,dpdl2)
      call DPDD (t,rhol,xl,dpdl)
      call PRESS (t,rhol,xl,pl)
      if (abs(rhol).gt.1.d6) then
        ierr=1
        herr=' '
        return
      endif
      if (dpdl2.lt.0.d0 .or. dpdl.lt.0.d0 .or. pl.lt.0.d0) then
        rhol=rhol*1.05d0
        it2=it2+1
        if (it2.lt.100) goto 110
      endif
      if (dpdl.gt.1.d4) then
        if (ABS(pl).gt.pc*1.5d0 .and. rhol-rhoc.lt.5) then
          rhol=rhol*1.5d0
          it2=it2+1
          if (it2.lt.100) goto 110
        endif
      endif


      it2=0
 120  continue
      call DPDD2 (t,rhov,xv,dpdv2)
      call DPDD (t,rhov,xv,dpdv)
      call PRESS (t,rhov,xv,pv)
      if (ABS(dpdv).gt.1.d5) then
        if (ABS(pv).gt.pc*1.5d0) then
          rhov=rhov*0.5d0
          it2=it2+1
          if (it2.lt.100) goto 120
        endif
      endif
      if (dpdv2.gt.0.d0 .or. dpdv.lt.0.d0) then
        rhov=rhov*.95d0
        it2=it2+1
        if (it2.lt.100) goto 120
      endif
      if (it.eq.1) then
        rhov=rhov*.95d0
        call DPDD2 (t,rhov,xv,dpdv2)
        call DPDD (t,rhov,xv,dpdv)
        call PRESS (t,rhov,xv,pv)
      endif



      p=(pl+pv)/2.d0

      if (abs(pv-p).gt.1.d-10 .and. abs(dpdv).gt.1.d-10) then
        deld=1.d0/(-dpdv/(pv-p)+dpdv2/2.d0/dpdv)
      else
        deld=deld/2.d0
      endif
 860  continue
      if (ABS(deld/rhov).gt.0.5d0) then
        deld=deld/10.d0
        goto 860
      endif
      rhov=rhov+deld

      if (abs(pl-p).gt.1.d-10 .and. abs(dpdl).gt.1.d-10) then
        deld=1.d0/(-dpdl/(pl-p)+dpdl2/2.d0/dpdl)
      else
        deld=deld/2.d0
      endif
 870  continue
      if (ABS(deld/rhol).gt.1.d0) then
        deld=deld/10.d0
        goto 870
      endif
      rhol=rhol+deld
      if (it.lt.50 .and. ABS(pl-pv).gt.tolr) goto 140

      ierr=0
      herr=' '
      RETURN
      end                                               !subroutine DLDV
c
c ======================================================================
c
      subroutine SPNDL (t,x,rhol,rhov,ierr,herr)
c
c  Find the spinodal densities for a given temperature.
c  Estimates for the liquid and vapor density are needed (rhol and rhov).
c
c  inputs:
c        t--temperature [K]
c        x--composition [array of mol frac]
c  outputs:
c     rhol--liquid spinodal [mol/L] (initial guess required)
c     rhov--vapor spinodal [mol/L] (initial guess required)
c     ierr--error flag:   0 = successful
c                       124 = failed to converge
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  11-28-06 EWL, original version taken from code in SATT
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)        !max number of components in mixture
      character*255 herr
      dimension x(ncmax)
      character*1 htab,hnull
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull

      ierr=0
      herr=' '
      call CRITP (x,tc,pc,rhoc,ierr,herr)
      call ISPURE (x,icomp)
      if (t.gt.tc+1.0d-8 .and. icomp.ne.0) then
        ierr=121
        write (herr,1121) t,tc,hnull
        call ERRMSG (ierr,herr)
 1121   format ('[SPNDL error 121] ',
     &          'temperature input to spinodal routine is ',
     &          'greater than critical temperature; T =',g11.5,
     &          ' K, Tcrit =',g11.5,' K.',a1)
        RETURN
      end if
      if (rhov.lt.0.d0) rhov=0.1d0
      it=0
      tolr=1.d-8
c  find the liquid spinodal (point where dpdrho=0)
 310  continue
      call DPDD (t,rhol,x,dpd)
      call DPDD2 (t,rhol,x,dpd2)
      d1=rhol
      dp1=dpd
      if (abs(dpd2).gt.1.d-20) rhol=rhol-dpd/dpd2
      it=it+1
      if (it.eq.10) tolr=tolr*100.d0
      if (it.eq.15) tolr=tolr*100.d0
      if (it.gt.30 .or. rhol.lt.rhoc) then
c  in case of failure, use the false position method:
c  (when the nonanalytical terms are used, this may be caused due to the
c   uncalculated dpdd2 part)
        d2=rhol
        call DPDD (t,d2,x,dp2)
c  check for good bounds; if bad, find new ones:
        if (rhol.lt.rhoc .or. dp1*dp2.ge.0.d0) then
          it=0
          d1=rhoc
          d2=rhoc
          call DPDD (t,d1,x,dp1)
 320      continue
          d2=d2*1.02d0
          call DPDD (t,d2,x,dp2)
          it=it+1
          if (it.gt.100) goto 390
          if (dp1*dp2.gt.0) goto 320
        endif
        it=0
        tolr=1.d-8
 330    continue
        if (abs(dp2-dp1).gt.1.d-20) rhol=d1-dp1*(d2-d1)/(dp2-dp1)
        call DPDD (t,rhol,x,dpd)
        if (dpd*dp1.lt.0.d0) then
          d2=rhol
          dp2=dpd
        else
          d1=rhol
          dp1=dpd
        endif
        it=it+1
        if (it.gt.100) goto 390
        if (abs(dpd).gt.tolr) goto 330
      endif
      if (abs(dpd).gt.tolr) goto 310
c
c  find the vapor spinodal
      rhov=rhov*.8d0
      it=0
      tolr=1.d-8
 340  continue
      call DPDD (t,rhov,x,dpd)
      call DPDD2 (t,rhov,x,dpd2)
      d1=rhov
      dp1=dpd
      if (abs(dpd2).gt.1.d-20) rhov=rhov-dpd/dpd2
      it=it+1
      if (it.eq.10) tolr=tolr*100.d0
      if (it.eq.15) tolr=tolr*100.d0
      if (it.gt.30 .or. rhov.gt.rhoc) then
c  false position method:
        d2=rhov
        call DPDD (t,d2,x,dp2)
c  check for good bounds; if bad, find new ones:
        if (d1.gt.rhoc .or. d2.gt.rhoc .or. dp1*dp2.ge.0.d0) then
          it=0
          d1=rhoc
          d2=rhoc
          call DPDD (t,d1,x,dp1)
 350      continue
          d2=d2/1.02d0
          call DPDD (t,d2,x,dp2)
          it=it+1
          if (it.gt.100) goto 390
          if (dp1*dp2.gt.0) goto 350
        endif
        it=0
        tolr=1.d-8
 360    continue
        if (abs(dp2-dp1).gt.1.d-20) rhov=d1-dp1*(d2-d1)/(dp2-dp1)
        call DPDD (t,rhov,x,dpd)
        if (dpd*dp1.lt.0.d0) then
          d2=rhov
          dp2=dpd
        else
          d1=rhov
          dp1=dpd
        endif
        it=it+1
        if (it.gt.100) goto 390
        if (abs(dpd).gt.tolr) goto 360
      endif
      if (abs(dpd).gt.tolr) goto 340
      RETURN
C
 390  continue
      ierr=124
      write (herr,1124) t,hnull
      call ERRMSG (ierr,herr)
      rhol=rhoc
      rhov=rhoc
 1124 format ('[SPNDL error 124] ',
     &        'iteration for spinodals did not converge; ',
     &        'T =',g11.5,' K.',a1)
      RETURN
      end                                              !subroutine SPNDL
c
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c ======================================================================
c                                                     end file sat_sub.f
c ======================================================================
      SUBROUTINE SATTP(t,p,x,iFlash,iGuess,d,Dl,Dv,xliq,xvap,q,ierr,
     & herr)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)   !max number of components in mixture
      parameter (MaxPr=100)
      common /NCOMP/ nc,ic
      common /vletmp1/ NofPrm
      common /vletmp2/ SSQBest,SSQCur,Offset,tcRed,pcRed,dcRed
      common /vletmp3/ xMnBest(MaxPr),xMnCur(MaxPr)
      common /vletmp4/ funs_eval(ncmax*2),StpX(MaxPr)
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),xmn(ncmax)
      character herr*255
      dimension xMnSave(1000,MaxPr)
      dimension SSQSave(1000),sumx(MaxPr)
      dimension xValue(3),yValue(3)
c     dimension i,j,k,Iter,Iter2,iSkp,iStp,ijk,iCurX,icn
c     dimension Tol,var1,Stp,fflp,fm
      dimension Pnts(ncmax*2,ncmax*2),zPnts(ncmax*2)
c     dimension temp1,t1,p1
c     dimension sum,SSQ1,SSQ2
      dimension fl(ncmax),fv(ncmax)

      iCurX=0
      iStp=1
      ierr=0
      ierr=0
      Iter=0
      Stp=0.001
      Tol=0.00000000001
      iSkp=0
      ijk=0
      SSQCur=0
      Offset=0
      fflp=0.02
      NofPrm=2*(nc-1)
      If (iFlash.ne.0) NofPrm=nc
      do i=1,NofPrm+1
        StpX(i)=0.d0
        xMnBest(i)=0.d0
      enddo

C     xmn array is the array being minimized
      If (iGuess.eq.0) Then
        If (iFlash.ne.0) Call Satest(iFlash,t,p,x,xliq,xvap,ierr,herr)
        If (iFlash.eq.0) Call Sat0est(t,p,x,xliq,xvap,ierr,herr)
      End If
      Call CRITP(x,tcRed,pcRed,dcRed,ierr,herr)
      d=dcRed
      Call xmnUpdate(t,p,1,iFlash,x,xliq,xvap,xmn,ierr,herr)
      Call xmnUpdate(t,p,5,iFlash,x,xliq,xvap,xmn,ierr,herr)
      If (t.le.0 .or. p.le.0) Then
        ierr=1
         GoTo 30
       endif

 10   continue
      Call GibbsMin(t,p,x,xmn,iFlash,1,0,Dl,Dv,xliq,xvap,fl,fv,ierr,
     & herr)
      If (ierr.ne.0) Then
        do i=1,NofPrm
          xmn(i)=xMnBest(i)*(1+fflp)
          If (xmn(i).le.0) xmn(i)=xMnBest(i)*(1+Abs(fflp))
        enddo
        Iter=Iter+1
        fflp=-1.2*fflp
        If (ierr.ne.0 .and. Iter.lt.52) GoTo 10
      End If

      SSQ1=SSQCur
C     SSQ2 = SSQCur
      fm=0.0001
      icn=0
      do i=1,30
        t1=t
        p1=p
        Call Newton(t,p,x,xmn,iFlash,Dl,Dv,xliq,xvap,q,ierr,herr)
 60   continue
        If (ierr.gt.0) Then
          Call xmnUpdate(t,p,6,iFlash,x,xliq,xvap,xMnBest,ierr,herr)
          fm=-fm*2
          If (iFlash.eq.0 .or. iFlash.eq.1 .or. iFlash.eq.3) Then
            sum=0
            do j=1,nc
              xvap(j)=Abs(xvap(j)*(1+fm*(-1)**j))
              sum=sum+xvap(j)
            enddo
            do j=1,nc
              xvap(j)=xvap(j)/sum
            enddo
          End If
          If (iFlash.eq.0 .or. iFlash.eq.2 .or. iFlash.eq.4) Then
            sum=0
            do j=1,nc
              xliq(j)=Abs(xliq(j)*(1-fm*(-1)**j))
              sum=sum+xliq(j)
            enddo
            do j=1,nc
              xliq(j)=xliq(j)/sum
            enddo
          End If
          If (iFlash.eq.1 .or. iFlash.eq.2) p=Abs(p*(1-fm/10000))
          If (iFlash.eq.3 .or. iFlash.eq.4) t=Abs(t*(1-fm/10000))
        Else
          If (SSQCur.lt.Tol) GoTo 30
          If (SSQCur.lt.SSQ1/1.05) Then
            SSQ1=SSQCur
          Else
            icn=icn+1
            If (icn.ge.5) Then
              ierr=1
              icn=0
              GoTo 60
            endif
          End If
        End If
      enddo
 50   continue
      q=999
      ierr=1
      herr='[SATTP error 1] iteration failed to converge'
      RETURN

      StpX(1)=Stp

 40   continue
      If (ierr.eq.0) Call xmnUpdate(t,p,1,iFlash,x,xliq,xvap,xmn,ierr,
     &    herr)
      do iCurX=1,50
        iStp=iStp+1
        If (iStp.gt.70) Then
          ierr=1
          herr='SATTP failed to converge'
          GoTo 30
        End If

        Iter=Iter+1
        If (Iter.gt.500) Then
          ierr=1
          GoTo 30
        endif
        do j=1,NofPrm
          xMnSave(iCurX,j)=xmn(j)
          xMnCur(j)=xmn(j)
        enddo
        SSQSave(iCurX)=SSQBest

        xValue(2)=0
        yValue(2)=SSQBest
        Offset=1
        If (iSkp.eq.0) Then
          Call GibbsMin(t,p,x,xmn,iFlash,1,0,Dl,Dv,xliq,xvap,fl,fv,ierr,
     & herr)
          xValue(3)=Offset
          yValue(3)=SSQCur
        End If
C     Use at least from 1 to 2, higher for better convergence along a single line
        do i=1,6
          temp1=0
          j=2
          If (i.eq.1) j=1
          Call Root(j,xValue,yValue,Offset,ierr,herr)
          If (ierr.ne.0) GoTo 30
          If (i.eq.1 .and. iSkp.eq.0) Offset=Offset/10D0
          Call GibbsMin(t,p,x,xmn,iFlash,1,0,Dl,Dv,xliq,xvap,fl,fv,ierr,
     & herr)
C     ierr = 0

          xValue(1)=xValue(2)
          yValue(1)=yValue(2)
          xValue(2)=xValue(3)
          yValue(2)=yValue(3)
          xValue(3)=Offset
          yValue(3)=SSQCur

          If (i.gt.1 .and. yValue(2).ne.0) Then
            If (Abs((yValue(3)-yValue(2))/yValue(2)).lt.0.01) GoTo 15
          End If
          If (i.gt.2 .and. xValue(2).ne.0) Then
            If (Abs((xValue(3)-xValue(2))/xValue(2)).lt.0.01) GoTo 15
          End If
        enddo

C     Use a Newton numerical method for saturation points
        temp1=1
        If (SSQCur.gt.Tol) Then
          If (p.lt.0.93*pcRed) Then
            Call Newton(t,p,x,xmn,iFlash,Dl,Dv,xliq,xvap,q,ierr,herr)
          End If
          temp1=1
        End If

 15   continue
C     Use a Newton numerical method for saturation points
        If (temp1.eq.0) Then
          If (SSQCur.gt.Tol) Then
            If (p.lt.0.93*pcRed) Then
              Call Newton(t,p,x,xmn,iFlash,Dl,Dv,xliq,xvap,q,ierr,herr)
            End If
            temp1=0
          End If
        End If
        Iter2=0

 20   continue
        If (SSQCur.gt.SSQSave(iCurX) .and. Iter2.lt.50 .and. Iter.gt.1)
     &      Then
          Iter2=Iter2+1
          Offset=-Offset/2
          Call GibbsMin(t,p,x,xmn,iFlash,1,0,Dl,Dv,xliq,xvap,fl,fv,ierr,
     & herr)
          If (SSQCur.gt.Tol) Then
            If (p.lt.0.93*pcRed) Then
              Call Newton(t,p,x,xmn,iFlash,Dl,Dv,xliq,xvap,q,ierr,herr)
            End If
            temp1=1
          End If
          If (Iter2.lt.4) GoTo 20
        End If

        If (iCurX.eq.10) Tol=Tol*10
        If (iCurX.eq.20) Tol=Tol*10
        If (iCurX.eq.30) Tol=Tol*10
        If (iCurX.eq.40) Tol=Tol*10
        If (iCurX.eq.50) Tol=Tol*10
        If (iCurX.eq.60) Tol=Tol*10
        If (SSQCur.lt.Tol .and. (SSQBest-SSQCur).lt.0.01) GoTo 30

        do j=1,NofPrm
          If (iCurX.lt.NofPrm) Then
            StpX(j)=0
            StpX(iCurX+1)=Stp
          Else
            i=iCurX-2
            If (i.lt.NofPrm) i=NofPrm
            If (iCurX.eq.NofPrm) i=1
            StpX(j)=(xMnBest(j)-xMnSave(i,j))/xMnBest(j)
            xValue(3)=-1
            yValue(3)=SSQSave(i)
            iSkp=1
            If (StpX(j).eq.0) Then
              StpX(j)=Stp
              iSkp=0
            endif
          End If
          xmn(j)=xMnBest(j)
        enddo

        If (Int(iCurX/3)*3.eq.iCurX .and. iCurX.gt.NofPrm) Then

C VLE.TxtScrn.SelStart = Len(VLE.TxtScrn)
          do i=1,NofPrm
            do j=1,NofPrm
              If (i.eq.1) Then
                Pnts(j,i)=xMnBest(j)
              Else
                Pnts(j,i)=xMnSave(iCurX+2-i,j)
              End If
            enddo
            zPnts(i)=1
          enddo
          Call LUdecomp(NofPrm,Pnts,zPnts,ierr,herr)
          If (ierr.ne.0) GoTo 30
C     Scale them to about the right size
          var1=Abs(StpX(1))*10
          do i=1,NofPrm
C     Normalize by the first one
            StpX(i)=zPnts(i)*var1*Offset/xMnBest(i)
          enddo
          iSkp=0
          If (iCurX.gt.10) Then
            do j=1,NofPrm
              sumx(j)=0
              k=0
              do i=iCurX-3,iCurX
                k=k+1
                sumx(j)=sumx(j)+Abs(xMnSave(i,j)-xMnSave(i-1,j))
              enddo
              sumx(j)=sumx(j)/k
            enddo
            ijk=ijk+1
            If (ijk.gt.NofPrm) ijk=1
            do j=1,NofPrm
              If (sumx(j).ne.0 .and. sumx(ijk).ne.0) Then
                If (Int(iCurX/6)*6.eq.iCurX) Then
                  StpX(j)=StpX(j)*(sumx(ijk)/sumx(j))
                Else
                  StpX(j)=StpX(j)*(sumx(j)/sumx(ijk))
                End If
              End If
            enddo
          End If
        End If
        j=0
        If (iCurX.gt.8) Then
          j=1
          do i=iCurX-8,iCurX
            If (((SSQSave(i)-SSQBest)/SSQBest*100).gt.10) j=2
          enddo
        End If
      enddo

 30   continue

      If (Abs(xvap(1)-xliq(1)).lt.0.00001) Then
        ierr=1
        herr='SATTP error:  Found bad root,compositions are equal'
      endif
      If (Abs(Dl-Dv).lt.0.01) Then
        ierr=1
        herr='SATTP error:  Density roots equal'
      endif
      If (SSQBest.gt.Tol*2) ierr=1
C     VLE.TxtScrn.SelStart = Len(VLE.TxtScrn)

      If (iFlash.eq.1 .or. iFlash.eq.3) Then
        q=0
        d=Dl
      endif
      If (iFlash.eq.2 .or. iFlash.eq.4) Then
        q=1
        d=Dv
      endif
      If (iFlash.eq.0 .and. ierr.eq.0) q=(x(1)-xliq(1))/(xvap(1)-
     &    xliq(1))
      END

      SUBROUTINE GibbsMin(t,p,x,xmn,iFlash,iflag,iknown,Dl,Dv,xliq,xvap,
     & fl,fv,ierr,herr)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)   !max number of components in mixture
      parameter (MaxPr=100)
      common /NCOMP/ nc,ic
      common /vletmp1/ NofPrm
      common /vletmp2/ SSQBest,SSQCur,Offset,tcRed,pcRed,dcRed
      common /vletmp3/ xMnBest(MaxPr),xMnCur(MaxPr)
      common /vletmp4/ funs_eval(ncmax*2),StpX(MaxPr)
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),xmn(ncmax)
      dimension fl(ncmax),fv(ncmax)
      character herr*255
c     dimension i,j,n,Iter3,SumF1
      dimension dgdxl(ncmax),dgdxv(ncmax)

      Iter3=0
      ierr=0
 10   continue

      If (ierr.ne.0 .and. iflag.eq.0) Then
        ierr=1
        GoTo 999
      endif
      do j=1,NofPrm
        xMnCur(j)=xmn(j)
      enddo

      If (iflag.eq.1) Then
        If (ierr.ne.0) Then
          If (Abs(Offset).gt.1E-20) Then
            Offset=Offset/2
          Else
            GoTo 999
          End If
          ierr=0
        End If
        Iter3=Iter3+1
        If (Iter3.gt.100) Then
          ierr=1
          herr='GibbsMin failed to converge'
          GoTo 999
        endif
        do j=1,NofPrm
          xMnCur(j)=xmn(j)*(1+Offset*StpX(j))
          If (xMnCur(j).le.0) Then
            ierr=1
            GoTo 10
          endif
          If (xMnCur(j).ge.5) Then
            ierr=1
            GoTo 10
          endif
          If (xMnCur(j).ge.1) Then
            If (j.lt.NofPrm .or. iFlash.eq.0) Then
              ierr=1
               GoTo 10
             endif
          End If
        enddo
        If (iFlash.eq.1 .or. iFlash.eq.2) p=xMnCur(NofPrm)*pcRed
        If (iFlash.eq.3 .or. iFlash.eq.4) t=xMnCur(NofPrm)*tcRed
        If (t.le.0 .or. p.le.0) Then
          ierr=1
          GoTo 999
        endif
        Call xmnUpdate(t,p,2,iFlash,x,xliq,xvap,xmn,ierr,herr)
        If (ierr.gt.0) GoTo 10
      End If

      If (iknown.ne.1) Then
        Call TPRHO(t,p,xliq,1,0,Dl,ierr,herr)
        If (Dl.gt.dcRed*7 .or. ierr.ne.0 .or. Dl.lt.dcRed/1.1) Then
          Dl=dcRed*3
          Call TPRHO(t,p,xliq,1,1,Dl,ierr,herr)
        End If
        If (Dl.lt.dcRed/1.1 .and. ierr.eq.0) ierr=1
        If (ierr.ne.0) GoTo 10
        Call FGCTY2(t,Dl,xliq,fl,ierr,herr)
      End If

      If (iknown.ne.2) Then
        Call TPRHO(t,p,xvap,2,0,Dv,ierr,herr)
        If (Dv.gt.dcRed*1.5) Then
          Dv=Dv*0.7
          Call TPRHO(t,p,xvap,2,1,Dv,ierr,herr)
        endif
        If (Dv.gt.dcRed*1.5 .and. Abs(Dl-Dv).lt.5 .and. ierr.eq.0)
     &      ierr=1
        If (ierr.ne.0) GoTo 10
        Call FGCTY2(t,Dv,xvap,fv,ierr,herr)
      End If

      If (fv(nc).le.0 .or. fl(nc).le.0) Then
        ierr=1
        GoTo 10
      endif
      do i=1,nc
        If (fv(i).le.0 .or. fl(i).le.0) Then
          ierr=1
          GoTo 10
        endif
        dgdxl(i)=Log(fl(i)/fl(nc))
        dgdxv(i)=Log(fv(i)/fv(nc))
      enddo

      SumF1=0
      funs_eval(nc)=0
      do i=1,nc
        If (i.ne.nc) funs_eval(i)=dgdxv(i)-dgdxl(i)
C     Is this the same as the funs_eval(n) below?
        funs_eval(nc)=funs_eval(nc)+xliq(i)*Log(fl(i))-
     & xvap(i)*Log(fv(i))-(xliq(i)-xvap(i))*dgdxv(i)
        SumF1=SumF1+funs_eval(i)**2
      enddo

      If (iFlash.eq.0 .and. nc.gt.2) Then
        n=nc
        do i=1,nc-2
          n=n+1
          funs_eval(n)=0
          do j=1,nc-1
            If (j.eq.1) Then
              funs_eval(n)=funs_eval(n)+dgdxv(j)-dgdxl(j)
            Else
              If (Abs(xvap(i)-x(i)).lt.1.d-16 .or.
     &            Abs(xliq(i)-x(i)).lt.1.d-16) Then
                ierr=1
                GoTo 10
              End If
              funs_eval(n)=funs_eval(n)+dgdxv(j)*(xvap(j)-
     & x(j))/(xvap(i)-x(i))-dgdxl(j)*(xliq(j)-x(j))/(xliq(i)-x(i))
            End If
          enddo
          SumF1=SumF1+funs_eval(n)**2
        enddo
      End If

      SSQCur=SumF1
      If (ierr.ne.0) GoTo 10
      If ((SSQCur.lt.SSQBest .and. SSQCur.ne.0) .or. SSQBest.eq.0) Call
     &    xmnUpdate(t,p,5,iFlash,x,xliq,xvap,xmn,ierr,herr)
      RETURN

 999  continue
C     If ierr <> 0 And herr = " " Then herr = " "  '"Error in subroutine GibbsMin"
      END

      SUBROUTINE Newton(t,p,x,xmn,iFlash,Dl,Dv,xliq,xvap,q,ierr,herr)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)   !max number of components in mixture
      parameter (MaxPr=100)
      common /NCOMP/ nc,ic
      common /vletmp1/ NofPrm
      common /vletmp2/ SSQBest,SSQCur,Offset,tcRed,pcRed,dcRed
      common /vletmp3/ xMnBest(MaxPr),xMnCur(MaxPr)
      common /vletmp4/ funs_eval(ncmax*2),StpX(MaxPr)
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),xmn(ncmax)
      character herr*255
c     dimension i,j,iLoop
c     dimension cont,w,d,sum,wcRed
      dimension xvar(ncmax*2),xwrk(ncmax),aJac(ncmax*2,ncmax*2),
     & deltaxvar(ncmax*2),anewxvar(ncmax*2)
      dimension fl(ncmax),fv(ncmax)

      iLoop=0
 10   continue
      Call Jacobian(t,p,x,xmn,iFlash,aJac,xvar,Dl,Dv,xliq,xvap,q,ierr,
     & herr)
            If (ierr.gt.0) GoTo 998
C     Call this?
      Call GibbsMin(t,p,x,xmn,iFlash,0,0,Dl,Dv,xliq,xvap,fl,fv,ierr,
     & herr)
                     If (ierr.gt.0) GoTo 998
      do i=1,NofPrm
        deltaxvar(i)=funs_eval(i)
      enddo
      Call LUdecomp(NofPrm,aJac,deltaxvar,ierr,herr)
                                                   If (ierr.gt.0) GoTo
     &                                                 998
      do i=1,NofPrm
        anewxvar(i)=xvar(i)-deltaxvar(i)
      enddo

      If (iFlash.ne.0) Then
        If (iFlash.eq.1 .or. iFlash.eq.2) wcRed=pcRed
        If (iFlash.eq.3 .or. iFlash.eq.4) wcRed=tcRed
        cont=1
 20   continue
        sum=0
        do i=1,nc-1
          xwrk(i)=anewxvar(i)
          sum=sum+xwrk(i)
        enddo
        xwrk(nc)=1-sum
        w=anewxvar(nc)*wcRed
 25   continue
        sum=1
        do i=1,nc
          If (xwrk(i).lt.0 .or. xwrk(i).gt.1) sum=0
        enddo
 30   continue
        If (sum.eq.0 .or. w.lt.0 .or. w/wcRed.gt.5) Then
          cont=cont+1
          If (cont.gt.50) GoTo 998

          If (2*Int(cont/2).eq.cont) Then
            Call xmnUpdate(t,p,4,iFlash,x,xliq,xvap,xwrk,ierr,herr)
            w=xwrk(NofPrm)*wcRed
            GoTo 25
          Else
            do i=1,nc
              anewxvar(i)=xvar(i)-deltaxvar(i)/cont**2
            enddo
            GoTo 20
          End If
        End If
        Call xmnUpdate(t,p,2,iFlash,x,xliq,xvap,xwrk,ierr,herr)
        If (ierr.gt.0) GoTo 999
        If (iFlash.eq.1 .or. iFlash.eq.2) p=w
        If (iFlash.eq.3 .or. iFlash.eq.4) t=w
        Call GibbsMin(t,p,x,xmn,iFlash,0,0,Dl,Dv,xliq,xvap,fl,fv,ierr,
     & herr)
        If (SSQCur.gt.SSQBest*1.2 .or. ierr.ne.0) Then
          sum=0
          GoTo 30
        endif

      Else
        cont=1
 40   continue
        j=0
        sum=0
        do i=1,nc-1
          j=j+1
          xvap(i)=anewxvar(j)
          sum=sum+xvap(i)
        enddo
        xvap(nc)=1-sum
        sum=0
        do i=1,nc-1
          j=j+1
          xliq(i)=anewxvar(j)
          sum=sum+xliq(i)
        enddo
        xliq(nc)=1-sum

        sum=1
        do i=1,nc
          If (xliq(i).lt.0 .or. xliq(i).gt.1 .or. xvap(i).lt.0 .or.
     &        xvap(i).gt.1) sum=0
        enddo

 50   continue
        If (sum.eq.0) Then
          cont=cont+1
          If (cont.gt.10) GoTo 998
          do i=1,2*(nc-1)
            anewxvar(i)=xvar(i)-deltaxvar(i)/cont**5
          enddo
          GoTo 40
        End If
        Call GibbsMin(t,p,x,xmn,iFlash,0,0,Dl,Dv,xliq,xvap,fl,fv,ierr,
     & herr)
        If (SSQCur.gt.SSQBest*1.5 .or. ierr.ne.0) Then
          sum=0
          GoTo 50
        endif
      End If

 998  continue
      If (ierr.ne.0) Then
        iLoop=iLoop+1
        If (iLoop.gt.3) Then
          ierr=2
          GoTo 999
        endif
        Call xmnUpdate(t,p,3,iFlash,x,xliq,xvap,xmn,ierr,herr)
        GoTo 10
      End If
      Call xmnUpdate(t,p,5,iFlash,x,xliq,xvap,xmn,ierr,herr)

 999  continue
C     If ierr <> 0 And herr = " " Then herr = " " '"Error in subroutine Newton"
      END

      SUBROUTINE Jacobian(t,p,x,xmn,iFlash,aJac,xvar,Dl,Dv,xliq,xvap,q,
     & ierr,herr)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)   !max number of components in mixture
      parameter (MaxPr=100)
      common /NCOMP/ nc,ic
      common /vletmp1/ NofPrm
      common /vletmp2/ SSQBest,SSQCur,Offset,tcRed,pcRed,dcRed
      common /vletmp3/ xMnBest(MaxPr),xMnCur(MaxPr)
      common /vletmp4/ funs_eval(ncmax*2),StpX(MaxPr)
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),xmn(ncmax)
      dimension aJac(ncmax*2,ncmax*2),xvar(ncmax*2)
      character herr*255
c     dimension i,j,k,m,iknown
      dimension xvapstart(ncmax),xliqstart(ncmax),xstart(ncmax),
     & xwrk(ncmax),funs_evalaux(ncmax*2,5)
c     dimension d,w,wstart,wcRed,deltax,delta
      dimension fl(ncmax),fv(ncmax)

      ierr=0
C     deltax = 0.00000001
       deltax=0.00000001

      iknown=0
      If (iFlash.ne.0) Then
        If (iFlash.eq.1 .or. iFlash.eq.3) Then
          do i=1,nc
          xstart(i)=xvap(i)
          enddo
        ElseIf (iFlash.eq.2 .or. iFlash.eq.4) Then
          do i=1,nc
          xstart(i)=xliq(i)
          enddo
        End If
        If (iFlash.eq.1 .or. iFlash.eq.2) Then
          wstart=p
          wcRed=pcRed
        ElseIf (iFlash.eq.3 .or. iFlash.eq.4) Then
          wstart=t
          wcRed=tcRed
        End If
        Call xmnUpdate(t,p,1,iFlash,x,xliq,xvap,xvar,ierr,herr)
        If (ierr.gt.0) GoTo 999
        do i=1,nc
          delta=xvar(i)*deltax
          m=0
C     For k = -2 To 2
          do k=-1,1
            If (k.ne.0) Then
              m=m+1
              If (i.eq.nc) Then
                 w=wstart+k*delta*wcRed
                 do j=1,nc
                  xwrk(j)=xstart(j)
                enddo
                iknown=0
              Else
                w=wstart
                do j=1,nc
                  xwrk(j)=xstart(j)
                enddo
                xwrk(i)=xstart(i)+k*delta
                xwrk(nc)=xstart(nc)-k*delta
              End If
              If (iFlash.eq.1 .or. iFlash.eq.2) p=w
              If (iFlash.eq.3 .or. iFlash.eq.4) t=w
              Call xmnUpdate(t,p,2,iFlash,x,xliq,xvap,xwrk,ierr,herr)
              If (ierr.gt.0) GoTo 999
              Call GibbsMin(t,p,x,xmn,iFlash,0,iknown,Dl,Dv,xliq,xvap,
     & fl,fv,ierr,herr)
              If (ierr.gt.0) GoTo 999
              If (iFlash.eq.1 .or. iFlash.eq.3) iknown=1
              If (iFlash.eq.2 .or. iFlash.eq.4) iknown=2
              do j=1,nc
                funs_evalaux(j,m)=funs_eval(j)
              enddo
            End If
          enddo
          do j=1,nc
            aJac(i,j)=(funs_evalaux(j,2)-funs_evalaux(j,1))/(2*delta)
          enddo
          Call xmnUpdate(t,p,2,iFlash,x,xliq,xvap,xstart,ierr,herr)
          If (ierr.gt.0) GoTo 999
          If (iFlash.eq.1 .or. iFlash.eq.2) p=wstart
          If (iFlash.eq.3 .or. iFlash.eq.4) t=wstart
        enddo
      Else
        do i=1,nc
          xvapstart(i)=xvap(i)
          xliqstart(i)=xliq(i)
        enddo
        do i=1,nc-1
          xvar(i)=xvap(i)
          xvar(i+nc-1)=xliq(i)
        enddo
        do i=1,2*(nc-1)
          delta=xvar(i)*deltax
          m=0
          do k=-1,1
            If (k.ne.0) Then
              m=m+1
              do j=1,nc
                xvap(j)=xvapstart(j)
                xliq(j)=xliqstart(j)
              enddo
              If (i.lt.nc) Then
                xvap(i)=xvapstart(i)+k*delta
                xvap(nc)=xvapstart(nc)-k*delta
              Else
                j=i-nc+1
                xliq(j)=xliqstart(j)+k*delta
                xliq(nc)=xliqstart(nc)-k*delta
              End If
              Call GibbsMin(t,p,x,xmn,iFlash,0,0,Dl,Dv,xliq,xvap,fl,fv,
     & ierr,herr)
              If (ierr.gt.0) GoTo 999
              do j=1,2*(nc-1)
                funs_evalaux(j,m)=funs_eval(j)
              enddo
            End If
          enddo
          do j=1,2*(nc-1)
            aJac(i,j)=(funs_evalaux(j,2)-funs_evalaux(j,1))/(2*delta)
          enddo
        enddo
        do i=1,nc
          xvap(i)=xvapstart(i)
          xliq(i)=xliqstart(i)
        enddo
      End If

 999  continue
C     If ierr <> 0 And herr = " " Then herr = " " '"Error in subroutine Jacobian"
      END

      SUBROUTINE Root(n,xValue,yValue,xNewRoot,ierr,herr)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      dimension xValue(3),yValue(3)
      character herr*255
c     dimension xw,xx,xy
c     dimension nn
      ierr=0
      nn=n
 10   continue
      If (nn.eq.1) Then
        If (yValue(2).eq.yValue(3)) Then
          xNewRoot=xNewRoot/2
           RETURN
         endif
        xNewRoot=xValue(2)-yValue(2)/(yValue(2)-yValue(3))*(xValue(2)-
     & xValue(3))
      ElseIf (nn.eq.2) Then
        xw=xValue(2)-xValue(3)
        xx=xValue(3)-xValue(1)
        xy=xValue(1)-xValue(2)
        If (xw*xx*xy.eq.0) GoTo 999
        xw=-(yValue(1)*xw+yValue(2)*xx+yValue(3)*xy)/(xw*xx*xy)
        If (xy.eq.0 .or. xw.eq.0) GoTo 999
        xx=(yValue(1)-yValue(2))/xy-xw*(xValue(1)+xValue(2))
        xNewRoot=-xx/(2*xw)
C     Curvature is negative, so get root using linear solution
        If (xw.lt.0) GoTo 999
      Else
        ierr=2
      End If
      If (xNewRoot.gt.100000) xNewRoot=100000
      If (xNewRoot.lt.-100000) xNewRoot=-100000
      RETURN
 999  continue
      nn=1
      GoTo 10
      END

      SUBROUTINE LUdecomp(n,aMatrix,cMatrix,ierr,herr)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)   !max number of components in mixture
      dimension amatrix(ncmax*2, ncmax*2), cmatrix(ncmax*2)
      dimension iord(ncmax*2), ctemp(ncmax*2), sdecomp(ncmax*2)
      character herr*255
c     dimension i,j,k,sum

      ierr=0
      do i=1,n
        iord(i)=i
        sdecomp(i)=Abs(aMatrix(1,i))
        do j=2,n
          If (Abs(aMatrix(j,i)).gt.sdecomp(i)) sdecomp(i)=Abs(aMatrix(j,
     &        i))
        enddo
C     Singular matrix
        If (sdecomp(i).eq.0) Then
          ierr=1
          herr='Singular matrix'
          RETURN
        endif
      enddo

      j=1
      Call Pivot(n,j,iord,aMatrix,sdecomp)
      If (aMatrix(1,iord(1)).eq.0) Then
        ierr=1
         RETURN
       endif
      do j=2,n
        aMatrix(j,iord(1))=aMatrix(j,iord(1))/aMatrix(1,iord(1))
      enddo
      do j=2,n-1
        do i=j,n
          sum=0
          do k=1,j-1
            sum=sum+aMatrix(k,iord(i))*aMatrix(j,iord(k))
          enddo
          aMatrix(j,iord(i))=aMatrix(j,iord(i))-sum
        enddo
        Call Pivot(n,j,iord,aMatrix,sdecomp)
        do k=j+1,n
          sum=0
          do i=1,j-1
            sum=sum+aMatrix(i,iord(j))*aMatrix(k,iord(i))
          enddo
          If (aMatrix(j,iord(j)).eq.0) aMatrix(j,iord(j))=1E+20
          aMatrix(k,iord(j))=(aMatrix(k,iord(j))-sum)/aMatrix(j,iord(j))
        enddo
      enddo
      sum=0
      do k=1,n-1
        sum=sum+aMatrix(k,iord(n))*aMatrix(n,iord(k))
      enddo
      aMatrix(n,iord(n))=aMatrix(n,iord(n))-sum
      If (aMatrix(n,iord(n)).eq.0) aMatrix(n,iord(n))=1E+20

      cMatrix(iord(1))=cMatrix(iord(1))/aMatrix(1,iord(1))
      do i=2,n
        sum=0
        do j=1,i-1
          sum=sum+aMatrix(j,iord(i))*cMatrix(iord(j))
        enddo
        cMatrix(iord(i))=(cMatrix(iord(i))-sum)/aMatrix(i,iord(i))
      enddo

      do i=n-1,1,-1
        sum=0
        do j=i+1,n
          sum=sum+aMatrix(j,iord(i))*cMatrix(iord(j))
        enddo
        cMatrix(iord(i))=cMatrix(iord(i))-sum
      enddo
      do i=1,n
        ctemp(i)=cMatrix(iord(i))
      enddo
      do i=1,n
        cMatrix(i)=ctemp(i)
      enddo
      END

      SUBROUTINE Pivot(n,j,iord,aMatrix,sdecomp)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)   !max number of components in mixture
      dimension iord(ncmax*2), amatrix(ncmax*2, ncmax*2),
     & sdecomp(ncmax*2)
c     dimension i,idummy,ipivt,big,dummy
      ipivt=j
      big=Abs(aMatrix(j,iord(j))/sdecomp(iord(j)))
      do i=j+1,n
        dummy=Abs(aMatrix(j,iord(i))/sdecomp(iord(i)))
        If (dummy.gt.big) Then
          big=dummy
          ipivt=i
        endif
      enddo
      idummy=iord(ipivt)
      iord(ipivt)=iord(j)
      iord(j)=idummy
      END

      SUBROUTINE Satest(iFlash,t,p,x,xliq,xvap,ierr,herr)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)   !max number of components in mixture
      common /NCOMP/ nc,ic
      dimension x(ncmax),xliq(ncmax),xvap(ncmax)
      character herr*255
c     dimension tini,ti,Tol,told,fnold
c     dimension i,j,k,dummy
c     dimension amag,fn,fnprime,tmp2
c     dimension tmp3,tmp4,taux1,taux2,taux3
c     dimension fnaux1,fnaux2,fnaux3,fnini
c     dimension tx,alpha
      dimension pc2(ncmax),tc2(ncmax),acf2(ncmax),wmm2(ncmax)
c     dimension wmm,ttrp,tnbpt,tc,pc,dc,zc,acf,dip,Rgas

      ierr=0
      do i=1,nc
        Call INFO(i,wmm2(i),ttrp,tnbpt,tc2(i),pc2(i),dc,zc,acf2(i),dip,
     & Rgas)
        xliq(i)=x(i)
        xvap(i)=x(i)
      enddo
      alpha=1

      If (t.eq.0 .or. iFlash.eq.3 .or. iFlash.eq.4) Then
        If (p.le.0) Then
          ierr=1
          RETURN
        endif
        tini=0
        do i=1,nc
          tini=tini+x(i)*tc2(i)/(1-
     & 0.428571*(Log(p/pc2(i))/Log(10D0))/(1+acf2(i)))
        enddo
      Else
        tini=t
      End If
      Tol=0.00000001

C     Do not remove parenthesis around p on the calls to PTest.  It needs to send the value, not the memory location.
      If (iFlash.eq.1) Then
        Call PTest(1,t,(p),alpha,x,tc2,pc2,acf2,p,xvap,ierr,herr)
        Call PTest(5,t,(p),alpha,x,tc2,pc2,acf2,p,xvap,ierr,herr)

      ElseIf (iFlash.eq.2) Then
        Call PTest(2,t,(p),alpha,x,tc2,pc2,acf2,p,xliq,ierr,herr)
        Call PTest(6,t,(p),alpha,x,tc2,pc2,acf2,p,xliq,ierr,herr)

      ElseIf (iFlash.eq.3 .or. iFlash.eq.4) Then
        do j=1,7
          dummy=0.7-(j-1)/10
          Call PTest(iFlash,tini,p,alpha,x,tc2,pc2,acf2,fnini,xvap,ierr,
     & herr)
          taux1=dummy*tini
          Call PTest(iFlash,taux1,p,alpha,x,tc2,pc2,acf2,fnaux1,xvap,
     & ierr,herr)
          taux2=tini/dummy
          Call PTest(iFlash,taux2,p,alpha,x,tc2,pc2,acf2,fnaux2,xvap,
     & ierr,herr)
          If (fnini*fnaux1.lt.0 .or. fnini*fnaux2.lt.0) GoTo 10
        enddo
 10   continue
        If (fnini*fnaux1.lt.0) Then
          taux2=(tini+taux1)/2
          Call PTest(iFlash,taux2,p,alpha,x,tc2,pc2,acf2,fnaux2,xvap,
     & ierr,herr)
          If (fnaux2*fnini.lt.0) Then
            taux1=(taux2+tini)/2
            Call PTest(iFlash,taux1,p,alpha,x,tc2,pc2,acf2,fnaux1,xvap,
     & ierr,herr)
            t=(taux1+taux2)/2
            If (fnaux1*fnini.lt.0) t=(taux1+tini)/2
          Else
            taux3=(taux2+taux1)/2
            Call PTest(iFlash,taux3,p,alpha,x,tc2,pc2,acf2,fnaux3,xvap,
     & ierr,herr)
            t=(taux3+taux1)/2
            If (fnaux3*fnaux2.lt.0) t=(taux3+taux2)/2
          End If
        ElseIf (fnini*fnaux2.lt.0) Then
          taux1=(tini+taux2)/2
          Call PTest(iFlash,taux1,p,alpha,x,tc2,pc2,acf2,fnaux1,xvap,
     & ierr,herr)
          If (iFlash.eq.3) Then
            t=(taux1+taux2)/2
            If (fnaux1*fnini.lt.0) t=(taux1+tini)/2
          ElseIf (iFlash.eq.4) Then
            If (fnaux1*fnini.lt.0) Then
              taux2=(taux1+tini)/2
              Call PTest(iFlash,taux2,p,alpha,x,tc2,pc2,acf2,fnaux2,
     & xvap,ierr,herr)
              t=(taux2+taux1)/2
              If (fnaux2*fnini.lt.0) t=(taux2+tini)/2
            Else
              taux3=(taux1+taux2)/2
              Call PTest(iFlash,taux3,p,alpha,x,tc2,pc2,acf2,fnaux3,
     & xvap,ierr,herr)
              t=(taux3+taux2)/2
              If (fnaux3*fnaux1.lt.0) t=(taux3+taux1)/2
            End If
          End If
        End If
        amag=2*Tol
        tmp3=100
        tmp4=100
        j=0
 30   continue
        If (amag.gt.Tol .or. tmp3.gt.1 .or. tmp4.gt.1) Then
          j=j+1
          If (j.gt.100) Then
            ierr=1
            RETURN
          endif
          Call PTest(iFlash,t,p,alpha,x,tc2,pc2,acf2,fn,xvap,ierr,herr)
          tmp2=0
          do i=1,nc
            tx=5.373*tc2(i)
            If (iFlash.eq.3) tmp2=tmp2+(Exp((-tx*acf2(i)-tx)/t+
     &          5.373*acf2(i)))*(tx*pc2(i)*acf2(i)*x(i)+tx*pc2(i)*x(i))
            If (iFlash.eq.4) tmp2=tmp2+(Exp((tx*acf2(i)+tx)/t-
     &          5.373*acf2(i)))*((-tx*acf2(i)*x(i))/pc2(i)-
     &          (tx*x(i))/pc2(i))
          enddo
          fnprime=215.508424158*tmp2*((fn-1)/t)**2/p
          If (iFlash.eq.4) fnprime=-0.004640189839*p*tmp2/(t**2)
          told=t
          If (fnprime.ne.0) t=told-(fn/fnprime)
          fnold=fn
          Call PTest(iFlash,t,p,alpha,x,tc2,pc2,acf2,fn,xvap,ierr,herr)
          amag=Abs(fn)
          tmp3=(Abs(told-t)/told)*100
          tmp4=(Abs(fnold-fn)/fnold)*100
          GoTo 30
        End If
        If (iFlash.eq.3) Call PTest(5,t,p,alpha,x,tc2,pc2,acf2,fn,xvap,
     &      ierr,herr)
        If (iFlash.eq.4) Call PTest(6,t,p,alpha,x,tc2,pc2,acf2,fn,xliq,
     &      ierr,herr)
      End If
      END

      SUBROUTINE Sat0est(t,p,x,xliq,xvap,ierr,herr)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)   !max number of components in mixture
      common /NCOMP/ nc,ic
      dimension x(ncmax),xliq(ncmax),xvap(ncmax)
      character herr*255
c     dimension alphaini,ti,Tol,alphaold,fnold
c     dimension i,j,k,dummy,Pb,Tb,Pd,Td
c     dimension temp1,amag,fn,fnprime,temp2
c     dimension temp3,temp4,alpha1,alpha2,alpha3,alpha
c     dimension fnaux1,fnaux2,fnaux3,fnini
      dimension pc2(ncmax),tc2(ncmax),acf2(ncmax),wmm2(ncmax)
c     dimension wmm,ttrp,tnbpt,tc,pc,dc,zc,acf,dip,Rgas

      ierr=0
      do i=1,nc
        Call INFO(i,wmm2(i),ttrp,tnbpt,tc2(i),pc2(i),dc,zc,acf2(i),dip,
     & Rgas)
      enddo

C     Check if (T,P) is in the two-phase region
      Tb=0.d0
      Td=0.d0
      Call Satest(1,t,Pb,x,xliq,xvap,ierr,herr)
      Call Satest(2,t,Pd,x,xliq,xvap,ierr,herr)
      Call Satest(3,Tb,p,x,xliq,xvap,ierr,herr)
      Call Satest(4,Td,p,x,xliq,xvap,ierr,herr)
      If ((t-Td)*(t-Tb).gt.0 .or. (p-Pb)*(p-Pd).gt.0) Then
        do i=1,nc
        xvap(i)=x(i)/2
        xliq(i)=x(i)
        enddo
        RETURN
      End If

      alphaini=0.5
      alpha1=0.000000001
      alpha2=0.999999999
      Call PTest(4,t,p,alphaini,x,tc2,pc2,acf2,fnini,xvap,ierr,herr)
      Call PTest(4,t,p,alpha1,x,tc2,pc2,acf2,fnaux1,xvap,ierr,herr)
      Call PTest(4,t,p,alpha2,x,tc2,pc2,acf2,fnaux2,xvap,ierr,herr)

      If (fnini*fnaux1.lt.0) Then
        alpha2=(alphaini+alpha1)/2
        Call PTest(4,t,p,alpha2,x,tc2,pc2,acf2,fnaux2,xvap,ierr,herr)
        If (fnaux2*fnini.lt.0) Then
          alpha1=(alpha2+alphaini)/2
          Call PTest(4,t,p,alpha1,x,tc2,pc2,acf2,fnaux1,xvap,ierr,herr)
          alpha=(alpha1+alpha2)/2
          If (fnaux1*fnini.lt.0) alpha=(alpha1+alphaini)/2
        Else
          alpha3=(alpha2+alpha1)/2
          Call PTest(4,t,p,alpha3,x,tc2,pc2,acf2,fnaux3,xvap,ierr,herr)
          alpha=(alpha3+alpha1)/2
          If (fnaux3*fnaux2.lt.0) alpha=(alpha3+alpha2)/2
        End If
      ElseIf (fnini*fnaux2.lt.0) Then
        alpha1=(alphaini+alpha2)/2
        Call PTest(4,t,p,alpha1,x,tc2,pc2,acf2,fnaux1,xvap,ierr,herr)
        If (fnaux1*fnini.lt.0) Then
          alpha2=(alpha1+alphaini)/2
          Call PTest(4,t,p,alpha2,x,tc2,pc2,acf2,fnaux2,xvap,ierr,herr)
          alpha=(alpha2+alpha1)/2
          If (fnaux2*fnini.lt.0) alpha=(alpha2+alphaini)/2
        Else
          alpha3=(alpha1+alpha2)/2
          Call PTest(4,t,p,alpha3,x,tc2,pc2,acf2,fnaux3,xvap,ierr,herr)
          alpha=(alpha3+alpha2)/2
          If (fnaux3*fnaux1.lt.0) alpha=(alpha3+alpha1)/2
        End If
      End If

      Tol=0.00001
      amag=2*Tol
      temp3=100
      temp4=100
      j=0
 30   continue
      If (amag.gt.Tol .or. temp3.gt.1 .or. temp4.gt.1) Then
        j=j+1
        Call PTest(4,t,p,alpha,x,tc2,pc2,acf2,fn,xvap,ierr,herr)
        temp1=0
        do i=1,nc
          temp2=Exp((5.373*tc2(i)*acf2(i)+5.373*tc2(i))/t)
          temp3=pc2(i)*215.508424158**(1+acf2(i))
          temp1=temp1-(x(i)*0.004640189839*temp2*(p*temp2-
     & temp3))/((alpha*(p*temp2-temp3)-p*temp2)**2)
        enddo
        fnprime=p*temp1*215.508424158
        alphaold=alpha
        alpha=alphaold-(fn/fnprime)
        If (alpha.gt.1E+20 .or. j.gt.100) Then
          ierr=1
          RETURN
        endif
        fnold=fn
        Call PTest(4,t,p,alpha,x,tc2,pc2,acf2,fn,xvap,ierr,herr)
        amag=Abs(fn)
        temp3=(Abs(alphaold-alpha)/alphaold)*100
        temp4=(Abs(fnold-fn)/fnold)*100
        GoTo 30
      End If

      Call PTest(6,t,p,alpha,x,tc2,pc2,acf2,fn,xliq,ierr,herr)
      do i=1,nc
        xvap(i)=(x(i)-xliq(i)*(1-alpha))/alpha
      enddo
      END

      SUBROUTINE PTest(inp,t,p,alpha,x,tc2,pc2,acf2,calc,xout,ierr,herr)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)   !max number of components in mixture
      common /NCOMP/ nc,ic
      dimension x(ncmax),xout(ncmax),pc2(ncmax),tc2(ncmax),acf2(ncmax)
      character herr*255

c     dimension i
c     dimension sum
      dimension ptemp(ncmax)

      ierr=0

      If (t.le.0) RETURN
      do i=1,nc
        ptemp(i)=pc2(i)*Exp(5.373*(1+acf2(i))*(1-tc2(i)/t))
        If (ptemp(i).gt.pc2(i)) ptemp(i)=pc2(i)+1*(ptemp(i)-pc2(i))
      enddo
      sum=0
      If (inp.eq.1) Then
        do i=1,nc
          sum=sum+x(i)*ptemp(i)
        enddo
        calc=sum
      ElseIf (inp.eq.2) Then
        do i=1,nc
          sum=sum+x(i)/ptemp(i)
        enddo
        calc=1/sum
      ElseIf (inp.eq.3) Then
        do i=1,nc
          sum=sum+x(i)*ptemp(i)
        enddo
        calc=1-p/sum
      ElseIf (inp.eq.4) Then
        do i=1,nc
          sum=sum+x(i)/(1-alpha+alpha*ptemp(i)/p)
        enddo
        calc=1-sum
      ElseIf (inp.eq.5) Then
        do i=1,nc
          xout(i)=x(i)/p*ptemp(i)
          If (xout(i).lt.0.000001) xout(i)=0.000001
          If (xout(i).gt.0.999999) xout(i)=0.999999
          sum=sum+xout(i)
        enddo
        do i=1,nc
          xout(i)=xout(i)/sum
        enddo
      ElseIf (inp.eq.6) Then
        do i=1,nc
          xout(i)=x(i)/(1-alpha+alpha*ptemp(i)/p)
          If (xout(i).lt.0.000001) xout(i)=0.000001
          If (xout(i).gt.0.999999) xout(i)=0.999999
          sum=sum+xout(i)
        enddo
        do i=1,nc
          xout(i)=xout(i)/sum
        enddo
      Else
      End If
      END

      SUBROUTINE xmnUpdate(t,p,ixmn,iFlash,x,xliq,xvap,xmn,ierr,herr)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      parameter (ncmax=20)   !max number of components in mixture
      parameter (MaxPr = 100)
      common /NCOMP/ nc,ic
      common /vletmp1/ NofPrm
      common /vletmp2/ SSQBest,SSQCur,Offset,tcRed,pcRed,dcRed
      common /vletmp3/ xMnBest(MaxPr),xMnCur(MaxPr)
      common /vletmp4/ funs_eval(ncmax*2),StpX(MaxPr)
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),xmn(ncmax)
      character herr*255
c     dimension i,j,sum
      j=0
      ierr=0

C     Put info into the minimizing array (xmn), compositions from 1 to nc-1 and p or t in the final position
      If (ixmn.eq.1) Then
        If (iFlash.eq.0 .or. iFlash.eq.1 .or. iFlash.eq.3) Then
          do i=1,nc-1
            j=j+1
            xmn(j)=xvap(i)
          enddo
        End If
        If (iFlash.eq.0 .or. iFlash.eq.2 .or. iFlash.eq.4) Then
          do i=1,nc-1
            j=j+1
            xmn(j)=xliq(i)
          enddo
        End If
        j=j+1
        If (iFlash.eq.1 .or. iFlash.eq.2) xmn(j)=p/pcRed
        If (iFlash.eq.3 .or. iFlash.eq.4) xmn(j)=t/tcRed

C     Get liq and vap compositions from the minimizing array, but not p and t.  x(nc) set from 1-sum[x(1..nc-1)]
      ElseIf (ixmn.eq.2) Then
        If (iFlash.eq.2 .or. iFlash.eq.4) Then
          do i=1,nc
            xvap(i)=x(i)
          enddo
        Else
          sum=0
          do i=1,nc-1
            j=j+1
            xvap(i)=xmn(j)
            sum=sum+xvap(i)
          enddo
          If (sum.gt.1) Then
            ierr=1
            RETURN
          endif
          xvap(nc)=1-sum
        End If
        If (iFlash.eq.1 .or. iFlash.eq.3) Then
          do i=1,nc
            xliq(i)=x(i)
          enddo
        Else
          sum=0
          do i=1,nc-1
            j=j+1
            xliq(i)=xmn(j)
            sum=sum+xliq(i)
          enddo
          If (sum.gt.1) Then
            ierr=1
            RETURN
          endif
          xliq(nc)=1-sum
        End If

C     Replace values with half way between current value and best value
C     Used when an iteration fails and the program backs up to a good solution
      ElseIf (ixmn.eq.3) Then
        If (iFlash.eq.0 .or. iFlash.eq.1 .or. iFlash.eq.3) Then
          sum=0
          do i=1,nc-1
            j=j+1
            xvap(i)=(xvap(i)+xMnBest(j))/2
            sum=sum+xvap(i)
          enddo
          If (sum.gt.1) Then
            ierr=1
            RETURN
          endif
          xvap(nc)=1-sum
        End If
        If (iFlash.eq.0 .or. iFlash.eq.2 .or. iFlash.eq.4) Then
          sum=0
          do i=1,nc-1
            j=j+1
            xliq(i)=(xliq(i)+xMnBest(j))/2
            sum=sum+xliq(i)
          enddo
          If (sum.gt.1) Then
            ierr=1
            RETURN
          endif
          xliq(nc)=1-sum
        End If
        j=j+1
        If (iFlash.eq.1 .or. iFlash.eq.2) p=(p+xMnBest(j)*pcRed)/2
        If (iFlash.eq.3 .or. iFlash.eq.4) t=(t+xMnBest(j)*tcRed)/2

C     Replace values with half way between current value and best value
C     Used when an iteration fails and the program backs up to a good solution
      ElseIf (ixmn.eq.4) Then
        do i=1,NofPrm
          xmn(i)=(xmn(i)+xMnBest(i))/2
        enddo

C     Save info in xMnBest array for possible use later
      ElseIf (ixmn.eq.5) Then
        If (iFlash.eq.0 .or. iFlash.eq.1 .or. iFlash.eq.3) Then
          do i=1,nc-1
            j=j+1
            xMnBest(j)=xvap(i)
            xMnCur(j)=xvap(i)
          enddo
        End If
        If (iFlash.eq.0 .or. iFlash.eq.2 .or. iFlash.eq.4) Then
          do i=1,nc-1
            j=j+1
            xMnBest(j)=xliq(i)
            xMnCur(j)=xliq(i)
          enddo
        End If
        j=j+1
        If (iFlash.eq.1 .or. iFlash.eq.2) xMnBest(j)=p/pcRed
        If (iFlash.eq.3 .or. iFlash.eq.4) xMnBest(j)=t/tcRed
        xMnCur(j)=xMnBest(j)
        SSQBest=SSQCur

C     Get info from xmn array
      ElseIf (ixmn.eq.6) Then
        If (iFlash.eq.0 .or. iFlash.eq.1 .or. iFlash.eq.3) Then
          sum=0
          do i=1,nc-1
            j=j+1
            xvap(i)=xmn(j)
            sum=sum+xvap(i)
          enddo
          xvap(nc)=1-sum
        End If
        If (iFlash.eq.0 .or. iFlash.eq.2 .or. iFlash.eq.4) Then
          sum=0
          do i=1,nc-1
            j=j+1
            xliq(i)=xmn(j)
            sum=sum+xliq(i)
          enddo
          xliq(nc)=1-sum
        End If
        j=j+1
        If (iFlash.eq.1 .or. iFlash.eq.2) p=xmn(j)*pcRed
        If (iFlash.eq.3 .or. iFlash.eq.4) t=xmn(j)*tcRed
      Else
      End If
      END
