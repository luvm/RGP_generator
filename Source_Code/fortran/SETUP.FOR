c  begin file setup.for
c
c  This file contains the routines which initialize the models, fluid-
c  specific parameters and coefficients, etc.  The subroutine SETUP must
c  be called before any of the other routines (except SETMOD) are called.
c  Call(s) to SETMOD and SETREF are optional and may be used to specify
c  non-standard models and reference states.
c
c  Subroutine SETUP0 calls subroutine SETUP using the same techniques as
c  those used by the graphical interface and the Excel spreadsheet.
c
c  Subroutine SETMIX reads the *.MIX files and makes the appropriate call
c  to subroutine SETUP.
c
c  Subroutine SETPATH sets the path where the *.fld files can be found
c
c  Subroutine PUREFLD allows the user to calculate the properties of a pure
c  fluid when a mixture has been loaded and the fluid is one of the
c  constituents in the mixture.
c
c  contained here are:
c     subroutine SETUP (nc,hfiles,hfmix,hrf,ierr,herr)
c     subroutine SETUP0 (i,hfld,hfm,hrf,ierr,herr)
c     subroutine SETFLD (icomp,hfile,ierr,herr)
c     subroutine SETMOD (nc,htype,hmix,hcomp,ierr,herr)
c     subroutine SETREF (hrf,ixflag,x0,h0,s0,t0,p0,ierr,herr)
c     subroutine SETMIX (hmxnme,hfmix,hrf,ncc,hfiles,x,ierr,herr)
c     subroutine SETPATH (hpth)
c     subroutine SETNC (ncomp)
c     subroutine RFFILE (hfilei,hflref)
c     subroutine OPENFL (nread,hfl,idir,ierr,herr)
c     subroutine PUREFLD (icomp)
c     subroutine GERG04 (nc,iflag,ierr,herr)
c     subroutine RESETA
c     block data BDSET
c
c  these routines set the values in the following common blocks
c     common /NCOMP/ ncomp,ic
c     common /CCAS/ hcas(n0:nx)
c     common /CNAM/ hname(n0:nx)
c     common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
c     common /TRNMOD/ heta,hetak(nrf0:nx),htcx,htcxk(nrf0:nx)
c     common /STNMOD/ hsten,hstenk(n0:nx)
c     common /DEMOD/ hdiel,hdielk(n0:nx)
c     common /MELTMOD/ hmelt,hmeltk(n0:nx)
c     common /SUBLMOD/ hsubl,hsublk(n0:nx)
c     common /PSMOD/ hps,hpsk(n0:nx)
c     common /PLMOD/ hpl,hplk(n0:nx)
c     common /DLMOD/ hdl,hdlk(n0:nx)
c     common /DVMOD/ hdv,hdvk(n0:nx)
c     common /REFST/ hsvrfs,hsvrfd(n0:nx)
c     common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
c     common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
c                   ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
c                   tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
c                   wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c
c  various arrays are dimensioned with parameter statements
c     parameter (ncmax=20)        !max number of components in mixture
c     parameter (nrefmx=10)       !max number of fluids for transport ECS
c     parameter (n0=-ncmax-nrefmx,nx=ncmax)
c     parameter (nrf0=n0)    !lower limit for transport ref fluid arrays
c
c ======================================================================
c ======================================================================
c
      subroutine SETUP (nc,hfiles,hfmix,hrf,ierr,herr)
c
c  define models and initialize arrays
c
c  A call to this routine is required.
c
c  inputs:
c       nc--number of components (1 for pure fluid) [integer]
c           if called with nc=-1, the version number*100 will be returned in ierr
c   hfiles--array of file names specifying fluid/mixture components
c           [character*255 variable] for each of the nc components;
c           e.g., :fluids:r134a.fld (Mac) or fluids\r134a.fld (DOS) or
c           [full_path]/fluids/r134a.fld (UNIX)
c    hfmix--mixture coefficients [character*255]
c           file name containing coefficients for mixture model,
c           if applicable
c           e.g.,  fluids\hmx.bnc
c      hrf--reference state for thermodynamic calculations [character*3]
c           'DEF':  default reference state as specified in fluid file
c                   is applied to each pure component
c           'NBP':  h,s = 0 at pure component normal boiling point(s)
c           'ASH':  h,s = 0 for sat liquid at -40 C (ASHRAE convention)
c           'IIR':  h = 200, s = 1.0 for sat liq at 0 C (IIR convention)
c           other choices are possible, but these require a separate
c           call to SETREF
c  outputs:
c     ierr--error flag:  0 = successful
c                      101 = error in opening file
c                      102 = error in file or premature end of file
c                     -103 = unknown model encountered in file
c                      104 = error in setup of model
c                      105 = specified model not found
c                      111 = error in opening mixture file
c                      112 = mixture file of wrong type
c                      114 = nc<>nc from setmod
c                     -117 = binary pair not found, all parameters will be estimated
c                      117 = No mixture data are available, mixture is outside the
c                            range of the model and calculations will not be made
c     herr--error string (character*255 variable if ierr<>0)
c     [fluid parameters, etc. returned via various common blocks]
c
c  explanation of parameters (used to dimension arrays)
c     ncmax:    maximum number of mixture components
c        n0:    lower bound on component arrays (-ncmax to accommodate
c               multiple ECS-thermo reference fluids)
c        nx:    same as ncmax
c      nrf0:    lower bound on arrays associated with transport props
c               (element 0 stores reference fluid information)
c
c  explanation of commons
c    /NCOMP/    nc--number of components
c
c    /EOSMOD/   equation of state (thermodynamic) models
c        hpheq--model for phase equilibria calcs (not currently used)
c         heos--equation of state model currently loaded (equal to hmxeos(1) if nc=1)
c       hmxeos--models for mixture components
c       hmodcp--ideal heat capacity model to use for component i
c               (set by call to SETeos routine)
c
c    /TRNMOD/   transport property models
c         heta--viscosity model for mixture
c        hetak--viscosity models for mixture components
c         htcx--thermal conductivity model for mixture
c        htcxk--conductivity models for mixture components
c
c    /STNMOD/   surface tension models
c        hsten--surface tension model for mixture
c       hstenk--surface tension models for mixture components
c
c    /DEMOD/    dielectric constant models
c    /MELTMOD/  melting line models
c    /SUBLMOD/  sublimation line models
c    /PSMOD/    vapor pressure equations
c    /PLMOD/    liquid pressure equations
c    /DLMOD/    saturated liquid density equations
c    /DVMOD/    saturated vapor density equations
c
c    /REFST/
c       hsvrfs--reference state (as specified by input argument hrf)
c       hsvrfd--default reference state for each component
c
c    /CCON/     constants for each of the NCMAX components; these values
c               are taken from the coefficients/array of selected model
c      tc(i):     critical temperature (K)
c      pc(i):     critical pressure (kPa)
c      rhoc(i):   critical density (mol/L)
c      zcrit(i):  critical compressibility factor (Pc/(R*Tc*rhoc))
c      ttp(i):    triple point temperature (K)
c      ptp(i):    triple point pressure (kPa)
c      dtp(i):    triple point density (mol/L)
c      tnbp(i):   normal boiling point temperature (K)
c      wm(i):     molecular mass (g/mol)
c      accen(i):  acentric factor
c      dipole(i): dipole moment (debye) at normal boiling pt
c
c    /CREF/     reference state for each of the NCMAX components
c      tref(i):   reference temperature for enthalpy and entropy
c      rhoref(i): reference density for enthalpy and entropy
c      href(i):   enthalpy at tref(i), rhoref(i) for component i
c      sref(i):   entropy at tref(i), rhoref(i) for component i
c
c    /CCAS/     CAS numbers for the n0..ncmax components
c      hcas(i):   "i" is mixture component number (n0..ncmax)
c
c    /CNAM/     short-hand names for the nc components
c      hname(i):  "i" is mixture component number
c
c    /HCHAR/    characters used to delimit/terminate output
c      htab:      tab (or other character) to delimit output tables
c      hnull:     ASCII null character to terminate error strings
c                 (for compatibility with mixed-language DLLs)
c
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  11-19-94  MM, original version
c  07-21-95  MM, restructure to add subordinate (model-specific)
c                set-up routines
c  07-24-95  MM, fluid const (e.g. Tc) set in SETeos rather than here
c  09-10-95  MM, add common /HCHAR/ and define contents
c  09-13-95  MM, add ierr, herr to argument list
c  09-25-95  MM, new argument list to SATT (outputs in order p, rho, x)
c  10-03-95  MM, read coefficients from files, change argument list
c  11-02-95  MM, expand to include mixture Helmholtz (HMX) model
c                add common block CNAM with fluid names
c  11-08-95  MM, initialize hmxeos array, even if heos not HMX
c  11-29-95  MM, variable lower limit on coefficient/constant arrays
c                to accommodate ECS reference fluid
c                kludge to load R134a as ECS reference fluid
c  12-08-95  MM, move file read to separate subroutine SETFLD
c  12-12-95  MM, remove R134a reference fluid kludge
c  01-09-96  MM, add call to SETHMX (read HMX.bnc file)
c  01-10-96  MM, add check for inputs same as previous call
c  01-11-96  MM, move reference state calculation to SETREF
c  01-19-96  MM, reset flag ksetrf to zero on call to SETUP
c  01-31-96  MM, implement 'NBS' option for EOS (and transport model)
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                add Zcrit to common /CCON/
c  02-29-96  MM, cover case of pure-fluid model specified for mixture
c  03-13-96  MM, changes to merge with Klein's transport code,
c                add common /TRNMOD/
c  03-19-96  MM, add dipole moment to /CCON/
c  03-20-96  MM, reduce argument list, move model specification to SETMOD
c                add mix models to /TRNMOD/
c  03-21-96  MM, replace /MODEL/ with /EOSMOD/, /STNMOD/, /REFST/
c  03-27-96  MM, fix bug when nc = 1 but overall model is a mix model
c  05-08-96  MM, add /MXINFO/ load hfmix into hmfile(0)
c  05-14-96  MM, add call to SETPH0 model (Helmholtz form)
c                disable check if inputs same as previous (!temporary)
c  11-19-96  MM, add commons related to critical lines, print parameters
c  11-25-96  MM, move debug print of mixture pars to SETHMX
c  02-20-97  MM, add hsvrfd to /REFST/, new common /CREFDF/ (for default ref st)
c  02-24-97  MM, add /CREMOD/ (pointer to transport critical enhancement models)
c  03-25-97  MM, bug fix:  x0 used but not dimensioned
c  03-28-97  MM, bug:  models initialized to 'NBS' only for first call to SETUP,
c                subsequent calls find past settings or 'NUL' if ncomp has increased;
c                add flag in /RESETM/ to fix
c  05-27-97  MM, reset gas constant when lreset=.true.
c  07-03-97  MM, use 'DEF' reference state if 'OTH' is specified
c  10-01-97  MM, add compiler switch to allow access by DLL
c  03-30-98  MM, initialize ierr1
c  04-10-98  MM, do not set component models to 'NBS' if mix model is 'NBS'
c                (this had effect of undoing calls to SETMOD)
c  12-01-98 EWL, add routine to convert hfiles*80 to hfiles*255 if necessary
c  12-01-98 EWL, set hfmix2 to hfmix to avoid problems with *80 to *255 change
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  03-11-99 EWL, reset lsatt and lsatp to .false. when called
c  06-03-99  MM, change default mixture surface tension model to STH
c  11-01-99 EWL, return from SETFEQ when ierr>0
c  11-02-99 EWL, make sure that selected equation of state has been loaded
c  11-02-99 EWL, set lreset to true before RETURN's
c  01-26-00  MM, add initialization of models & citations in /CITE/
c  03-16-00 EWL, initialize variables in FXPNT
c  12-18-00 MLH, set mixture reference fluid for transport to nitrogen
c  01-23-02 EWL, split common block CITE into 2 pieces for some compilers
c  01-23-02 EWL, split common block MXINFO into 2 pieces for some compilers
c  05-27-02 EWL, check for ammonia/water mixture, adjust R, Href, Sref
c  07-29-02 EWL, changed third argument of SETREF to x0 instead of 0.d0
c  08-20-02 EWL, zero out the iamwat flag
c  09-19-02 EWL, add the ianc flag
c  09-19-02 EWL, add liquid pressure ancillary setup
c  10-31-02 EWL, add AGA8 equation of state setup
c  11-26-02 EWL, add OT0 reference state
c  10-12-04 MLH, add PR setup
c  11-19-04 EWL, read UN number
c  01-03-05 EWL, remove AGA8 equation of state setup, just call setaga after call to setup
c  07-22-05 MLH, add multiple reference fluids for transport properties model
c  09-28-05  DT, change occurrences of a backslash to 'char(92)' for compatibility
c                with some compilers using F77
c  10-07-05 EWL, only modify ammonia/water reference state if user has not redefined it
c  06-13-06 EWL, modified most of the fortran files to remove the line numbers
c                in the do...continue loops.
c  06-13-06 EWL, modified every fortran file throughout all of Refprop to enable
c                the user to load multiple fluids (as if it were a mixture) and
c                then set which pure fluid is being used (though subroutine PureFld)
c                in the calculations.  Common block NCOMP was modified everywhere
c                to include the variable ic, which sets the pure fluid currently
c                in use.  When ic is zero, then mixture calculations are made.
c  07-18-06 EWL, add linit to set strings to 'NUL' on first call
c  07-24-06 EWL, rework lreset logic
c  07-26-06 EWL, remove excess code, rework setup
c  08-10-06 EWL, add action='READ' to all open statements
c  09-25-06 EWL, changed all occurrences in every routine from herr(1:1)=hnull to herr=' '
c  03-30-07 EWL, add dll_export switches to all top level REFPROP subroutines
c  09-13-09 EWL, change loop to call SETPRCO from icomp=0 to icomp=1
c  02-04-10 EWL, change the ianc variable into an array
c  03-03-10 EWL, return version number when called with nc=-1
c  10-28-10 EWL, initialize more variables relating to transport properties
c
c  Compiler switches to allow access to the routines in the DLL.
c  The DEC$ switch is for use with Digital Visual Fortran; this should be
c    treated as a comment by all other compilers.
c  The dll_export switch is for use with Lahey Fortran and needs to be
c    uncommented throughout all the fortran files in order for it to work.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SETUP
c     dll_export SETUP
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (nrf0=n0)     !lower limit for transport ref fluid arrays
      parameter (nrefluids=4) ! number of ref fluids available for transport properties model
      character*1 htab,hnull
      character*3 hrf,hsvrfs,hrefdf,hpheq,heos,hmxeos,hmodcp
      character*3 heta,hetak,htcx,htcxk,hetacr,htcxcr,hsten,hstenk
      character*3 hdiel,hdielk,hmelt,hmeltk,hsubl,hsublk
      character*3 hps,hpsk,hpl,hplk,hdl,hdlk,hdv,hdvk
      character*3 hsvph,hsveqn,hsveqk,hsvvis,hsvvik,hsvcnd,hsvcnk
      character*3 hsvdil,hsvdik,hsvsrf,hsvsrk
      character*12 hcas,hname
      character*255 hfiles(ncmax),hfile(n0:nx),hfmix,hfmix2,rflnam
      character*255 hsvfld,hsvmix
      character*255 herr,herr1,ucase
      character*255 hieos,hicp0,hieta,hietac,hitcx,hitcxc,histn,
     &              hidiel,himelt,hisubl,hips,hipl,hidl,hidv
      logical lreset,lrst,linit
      dimension x0(ncmax)
c
      common /NCOMP/ ncomp,ic
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /TRNMOD/ heta,hetak(nrf0:nx),htcx,htcxk(nrf0:nx)
      common /CREMOD/ hetacr(nrf0:ncmax),htcxcr(nrf0:ncmax)
      common /STNMOD/ hsten,hstenk(n0:nx)
      common /DEMOD/  hdiel,hdielk(n0:nx)
      common /MELTMOD/ hmelt,hmeltk(n0:nx)
      common /SUBLMOD/ hsubl,hsublk(n0:nx)
      common /PSMOD/ hps,hpsk(n0:nx)
      common /PLMOD/ hpl,hplk(n0:nx)
      common /DLMOD/ hdl,hdlk(n0:nx)
      common /DVMOD/ hdv,hdvk(n0:nx)
      common /REFST/ hsvrfs,hrefdf(n0:nx)
      common /CCAS/ hcas(n0:nx)
      common /CNAM/ hname(n0:nx)
      common /HCHAR/ htab,hnull
      common /SETSAV/ hsvfld(n0:nx),hsvmix
      common /IRFSAV/ ixfsav,ksetrf
      common /MODSAV/ hsvph,hsveqn,hsveqk(nx),hsvvis,hsvvik(nx),
     &            hsvcnd,hsvcnk(nx),hsvsrf,hsvsrk(nx),hsvdil,hsvdik(nx)
      common /RESETM/ lreset,lrst,linit,ncset !flag indicating need to reset models
      common /CITE/  hieos(n0:nx),hicp0(n0:nx),
     &               hieta(nrf0:nx),hietac(nrf0:nx),
     &               hitcx(nrf0:nx),hitcxc(nrf0:nx),histn(n0:nx)
      common /CITE2/ hidiel(n0:nx),himelt(n0:nx),hisubl(n0:nx),
     &               hips(n0:nx),hipl(n0:nx),hidl(n0:nx),hidv(n0:nx)
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /WLMSTN/ tmstn(n0:nx),txstn(n0:nx)
      common /WLMTCX/ tmtcx(nrf0:nx),txtcx(nrf0:nx),ptcx(nrf0:nx),
     &                Dtcx(nrf0:nx)
      common /WLMETA/ tmeta(nrf0:nx),txeta(nrf0:nx),peta(nrf0:nx),
     &                Deta(nrf0:nx)
      common /WLMTRN/ tmecs(nrf0:nx),txecs(nrf0:nx),pecs(nrf0:nx),
     &                Decs(nrf0:nx)
c
c  initialize variables upon first call to setup
      if (nc.eq.-1) then
        ierr=900   !Send back version number*100.
        RETURN
      endif

      if (linit) then
        hsveqn='NUL'
        hsvvis='NUL'
        hsvcnd='NUL'
        hsvsrf='NUL'
        hsvdil='NUL'
        hsvmix='NUL'
        hsvrfs='NUL'
        hsvph ='NUL'
        do i=1, nx
          hsveqk(i)='NUL'
          hsvvik(i)='NUL'
          hsvcnk(i)='NUL'
          hsvsrk(i)='NUL'
          hsvdik(i)='NUL'
        enddo
        do i=n0, nx
          hsvfld(i)='NUL'
          hrefdf(i)='NUL'
        enddo
        linit=.FALSE.
      endif
c
c  reset fixed points
c
      htab=CHAR(9)         !tab character for output tables
      hnull=CHAR(0)        !null character to terminate error strings
      ierr=0
      ierr1=0
      herr=' '
      hfmix2=hfmix
c  used to check for version 6.0 inputs:

c     if(hfmix(79:80).eq.'  '.and.hfmix(81:82).ne.'  ')hfmix2=hfmix(1:80)
c
c  compare inputs to previous values (stored in /xxxMOD/ and /SETSAV/)
c  if nothing has changed SETUP can be bypassed
c
      lsame=.false.
      if (nc.eq.ncomp .and. hsvrfs.eq.hrf) then
c       write (*,*) ' SETUP:  ncomp and ref state same as previous'
        if((hsveqn.eq.heos  .or. hsveqn.eq.'NBS') .and. !EOS model
     &     (hsvvis.eq.heta  .or. hsvvis.eq.'NBS') .and. !viscosity model
     &     (hsvcnd.eq.htcx  .or. hsvcnd.eq.'NBS') .and. !therm cond model
     &     (hsvsrf.eq.hsten .or. hsvsrf.eq.'NBS') .and. !surface ten model
C    &     (hsvdil.eq.hdiel .or. hsvdil.eq.'NBS') .and. !dielectric model
     &      hsvph.eq.hpheq) then      !phase equil method, not implemented
          lsame=.true.
c         write (*,*) ' SETUP:  after check of models lsame =  ',lsame
          do i=1,nc
c  check that .fld file and component models are the same as last call
            if(hfiles(i).ne.hsvfld(i) .or.
     &        (hsveqk(i).ne.hmxeos(i)  .and. hsveqk(i).ne.'NBS').or.
     &        (hsvsrk(i).ne.hstenk(i)  .and. hsvsrk(i).ne.'NBS').or.
C    &        (hsvdik(i).ne.hdielk(i)  .and. hsvdik(i).ne.'NBS').or.
     &        (hsvvik(i).ne.hetak(i)   .and. hsvvik(i).ne.'NBS').or.
     &        (hsvcnk(i).ne.htcxk(i)   .and. hsvcnk(i).ne.'NBS')) then
              lsame=.false.
            end if
c         write (*,*) ' SETUP:  after component models lsame = ',lsame
          enddo
c  check for previous call to setref with ixflag=2
          if (ixfsav.eq.2) lsame=.false.
          if (hsvmix.ne.hfmix2) lsame=.false.   !mixture coeff file
c         write (*,*) ' SETUP:  after mixture file lsame =     ',lsame
        end if
c  if lrst has been set equal to true, then both SETMOd and SETUP were called
c  simultaneously earlier, and this call to SETUP must reset the earlier call
c  to SETMOD.
        if (lrst) lsame=.false.
        if (lsame) RETURN
      end if
c
c
c-------start of main code--------
c
c
      lrst=.false.
      if (.not.lreset) lrst=.true.  !Flag to indicate that SETMOD was called outside of SETUP
      hsvmix='x'  !Delete saved input in case setup called with nc=0 to reset system
c  check state of flag lreset; a value of .true. indicates need to reset
c  all models to 'NBS'; this occurs when SETUP is called a second time
c  in a given application
      if (lreset) call SETMOD (nc,'NBS',heos,hmxeos,ierr,herr)
      lreset=.true.
      iamwat=0
      iwat=0
      do i=0,ncmax
        ianc(i)=0
      enddo
      call PUREFLD(0)
      if (ierr.gt.0) RETURN
      if (nc.ne.ncset) then
        ierr=114
        herr='[SETUP error 114] number of components does not match '//
     &  'value sent to SETMOD'//hnull
        call ERRMSG (ierr,herr)
        RETURN
      endif
      R=8.314472d0   !gas constant, Mohr and Taylor, JPCRD, 28(6):1713-1852, 1999.
c     R=8.314510d0   !gas constant, CODATA recommended value
c
c     write (*,*)    'Input values'
c     write (*,1002) heos, (hmxeos(j),j=0,12),
c    &               heta, (hetak(j), j=0,12),
c    &               htcx, (htcxk(j), j=0,12),
c    &               hsten,(hstenk(j),j=0,12)
c     write (*,*)    'Saved values'
c     write (*,1002) hsveqn,(hsveqk(j),j=0,12),
c    &               hsvvis,(hsvvik(j),j=0,12),
c    &               hsvcnd,(hsvcnk(j),j=0,12),
c    &               hsvsrf,(hsvsrk(j),j=0,12)
c1002 format ('    heos, hmxeos: ',2(a3,3x),12(a3,1x)/
c    &        '     visc models: ',2(a3,3x),12(a3,1x)/
c    &        '      tcx models: ',2(a3,3x),12(a3,1x)/
c    &        ' surf ten models: ',2(a3,3x),12(a3,1x))
c
c  set pointers to models
c
c  special case if equation of state is NIST recommendation
c  only mixture choice at present is the HMX model
      if (heos.eq.'NBS' .or. heos.eq.'nbs') then
        heos='HMX'
        if (nc.eq.1) heos=hmxeos(1)
      else if (heos.eq.'HMX') then
c  mixture model is HMX, use component models as set in SETMOD
      else
        if (nc.gt.1) then
c  this condition should not be accessed; a mix model other than HMX
c  has been specified, but only choice at present is the HMX model
          ierr=-105
          herr='[SETUP warning -105] unknown mixture model specified,'
     &       //' the HMX model will be used.'//hnull
          call ERRMSG (ierr,herr)
          heos='HMX'
          do i=1,nc
            hmxeos(i)='NBS'
          enddo
        end if
      end if
c
c  special case if viscosity model is NIST recommendation
c  only mixture choice at present is the ECS model
      if (heta.eq.'NBS' .or. heta.eq.'nbs') then
        heta='ECS'
        if (nc.eq.1) heta=hetak(1)
      else if (heta.eq.'ECS') then
c  mixture model is ECS, use component models as set in SETMOD
      else
        if (nc.gt.1) then
c  this condition should not be accessed; a mix model other than ECS
c  has been specified, but only choice at present is the ECS model
c         write (*,*) ' SETUP--unknown mix viscosity model, ECS used'
          heta='ECS'
          do i=1,nc
            hetak(i)='NBS'
          enddo
        end if
      end if
c
c  special case if thermal conductivity model is NIST recommendation
c  only mixture choice at present is the ECS model
      if (htcx.eq.'NBS' .or. htcx.eq.'nbs') then
        htcx='ECS'
        if (nc.eq.1) htcx=htcxk(1)
      else if (htcx.eq.'ECS') then
c  mixture model is ECS, use component models as set in SETMOD
      else
        if (nc.gt.1) then
c  this condition should not be accessed; a mix model other than ECS
c  has been specified, but only choice at present is the ECS model
c         write (*,*) ' SETUP--unknown mix conductivity model, ECS used'
          htcx='ECS'
          do i=1,nc
            htcxk(i)='NBS'
          enddo
        end if
      end if
c
c  special case if surface tension model is NIST recommendation
      if (hsten.eq.'NBS' .or. hsten.eq.'nbs') then
        hsten='STH'
        if (nc.eq.1) hsten=hstenk(1)
c     else if (hsten.eq.'STM' .or. hsten.eq.'STX') then
c  mixture model is STM or STX, both are contained in routine STN;
c  use component models as specified in SETMOD
      end if
c
c  special case if dielectric constant model is NIST recommendation
      if (hdiel.eq.'NBS' .or. hdiel.eq.'nbs') then
        hdiel='DEM'
        if (nc.eq.1) hdiel=hdielk(1)
c     else if (hdiel.eq.'DEM' .or. hdiel.eq.'DEX') then
c  mixture model is DEM or DEX, both are contained in routine DE;
c  use component models as specified in SETMOD
      end if
c
c  old code for compatibility with version 6.0
c     if (hfiles(1)(79:80).eq.'  ' .and. hfiles(1)(81:82).ne.'  ') then
c       hfile(1)=hfiles(1)(1:80)
c       if (nc.ge.2) hfile(2)=hfiles(1)(81:160)
c       if (nc.ge.3) hfile(3)=hfiles(1)(161:240)
c       if (nc.ge.4) hfile(4)=hfiles(1)(241:255)//hfiles(2)(1:65)
c       if (nc.ge.5) hfile(5)=hfiles(2)(66:145)
c     endif
c
      ksetrf=0
      hsvph='NBS'
      hsvmix=hfmix2
      ncomp=nc
c
c  fill up info arrays for unused/undefined components/reference fluid
c
      do i=n0,nc
        hfile(i)='NUL'
        hcas(i)='not_defined'
        hname(i)='not_defined'
        hieos(i)='NUL'//hnull
        hicp0(i)='NUL'//hnull
        histn(i)='NUL'//hnull
        hidiel(i)='NUL'//hnull
        himelt(i)='NUL'//hnull
        hisubl(i)='NUL'//hnull
        hips(i)='NUL'//hnull
        hipl(i)='NUL'//hnull
        hidl(i)='NUL'//hnull
        hidv(i)='NUL'//hnull
        hitcx(i)='NUL'//hnull
        hieta(i)='NUL'//hnull
        hitcxc(i)='NUL'//hnull
        hietac(i)='NUL'//hnull
        tmstn(i)=300.d0
        txstn(i)=300.d0
        tmtcx(i)=300.d0
        txtcx(i)=300.d0
        tmeta(i)=300.d0
        txeta(i)=300.d0
        tmecs(i)=300.d0
        txecs(i)=300.d0
        ptcx(i)=0.d0
        Dtcx(i)=0.d0
        peta(i)=0.d0
        Deta(i)=0.d0
        pecs(i)=0.d0
        Decs(i)=0.d0
      enddo
      do i=1,nc
        hfile(i)=hfiles(i)
        hsvfld(i)=hfiles(i)
      enddo
      hetak(0)='NUL'
      htcxk(0)='NUL'
      hmxeos(0)='FEQ'     !Set default in case .fld file has no transport model
c
c  store fluid constants for each of the nc components
c
      do icomp=1,nc
        if (INDEX(UCASE(hfile(icomp),255),'PPF').ne.0) ianc(icomp)=1
        call SETFLD (icomp,hfile,ierr,herr)
        if (ierr.gt.0) RETURN  !error in opening file--further processing pointless
      enddo
c
c  store fluid constants for ECS reference fluids (if any)
c
      do icomp=-ncmax,0    !ecs ref fluids from -ncmax to 0 (not n0 to 0)
        if (hfile(icomp).ne.'NUL') then
          if (icomp.eq.0 .and. nc.gt.1) then
c           no need to call SETFLD here, it gets called below for i=0
          else
            call SETFLD (icomp,hfile,ierr,herr)
            if (ierr.gt.0) RETURN  !error in opening file--further processing pointless
          endif
        end if
      enddo
c
c
c  for a mixture, set reference fluids for mixture calculations; currently
c  hfile(0) corresponds to ref fluid of last pure fluid loaded; reset to nitrogen.
c  presently there are 4 possible reference fluids for mixtures
c  slot      0: nitrogen
c     -ncmax-1: propane
c     -ncmax-2: dodecane
c     -ncmax-3: R134a
c
c     write (*,*) 'Reference fluid:  ',hfile(0)(1:60)
      if (nc.gt.1) then
        ipos=0
        do ii=1,255
c  EWL change from hfile(0) to hfile(1)
          if (hfile(1)(ii:ii).eq.'/' .or. hfile(1)(ii:ii).eq.':'
     &        .or. hfile(1)(ii:ii).eq.char(92)) then
            ipos=ii
          end if
        enddo
        do ij=0,nrefluids-1
          il=-ncmax-ij
          if (ij.eq.0) then
            il=0
            rflnam='nitrogen.fld'
          elseif (ij.eq.1) then
            rflnam='propane.fld'
          elseif (ij.eq.2) then
            rflnam='c12.fld'
          elseif (ij.eq.3) then
            rflnam='r134a.fld'
          endif
          if (ipos.eq.0) then
            hfile(il)=rflnam
          else
            hfile(il)=hfile(1)(1:ipos)//rflnam !append ref fluid name
          endif
          call SETFLD (il,hfile,ierr,herr)
          if (ierr.gt.0) RETURN  !error in opening file--further processing pointless
        enddo
      endif
c
c
c  set up mixture model, if applicable
c
      if (heos.eq.'HMX') then
        call SETHMX (hfmix2,ierr1,herr1)
c       write (*,*) 'ierr from SETHMX:  ',ierr1,':  ',herr1
      end if
c
c  set reference state
c
      ixflag=1    !only 'pure fluid' reference state can be called here
      if (hrf(1:2).eq.'OT') then
        ierr=-105
        herr='[SETUP warning -105] must use routine SETREF for (OTH) '
     &      //'reference state choice; will use default (DEF) choice'
     &      //hnull
        call ERRMSG (ierr,herr)
        hrf='DEF'        !use the default reference state
        if (ierr1.eq.0) then
          ierr1=ierr
          herr1=herr
        end if
      end if
      call SETREF (hrf,ixflag,x0,h0,s0,t0,p0,ierr,herr)
c     write (*,*) ' SETUP--ierr returned from SETREF:  ',ierr
c  any error from SETREF takes precedence over error from SETHMX
      if (ierr.eq.0 .and. ierr1.ne.0) then
        ierr=ierr1
        herr=herr1
      end if
c
c  save current models into /MODSAV/--check on subsequent calls to
c  SETUP, if nothing has changed, can bypass call
c
      if (ierr.eq.0) then
        hsvph=hpheq
        hsveqn=heos
        hsvvis=heta
        hsvcnd=htcx
        hsvsrf=hsten
        hsvdil=hdiel
        do k=1,nc
          hsveqk(k)=hmxeos(k)
          hsvvik(k)=hetak(k)
          hsvcnk(k)=htcxk(k)
          hsvsrk(k)=hstenk(k)
          hsvdik(k)=hdielk(k)
        enddo
      endif
c
c     write (*,1442) heos,(hmxeos(j),j=-5,5),
c    &                    (hmodcp(j),j=-5,5),
c    &               heta,(hetak(j),j=0,5),
c    &                    (hetacr(j),j=0,ncomp)
c     write (*,1443) htcx,(htcxk(j),j=0,5),
c    &                    (htcxcr(j),j=0,ncomp)
c     write (*,1444) hsten,(hstenk(j),j=-5,5)
c1442 format ('   SETUP--heos, hmxeos: ',a3,3x,5(a3,1x),a5,3x,5(a3,1x)/
c    &        '   (exit)     Cp0 mods: ',3x,3x,5(a3,1x),a5,3x,5(a3,1x)/
c    &        '           visc models: ',a3,25x,a3,3x,5(a3,1x)/
c    &        '        visc crit mods: ',3x,25x,a3,3x,5(a3,1x))
c1443 format ('            tcx models: ',a3,25x,a3,3x,5(a3,1x)/
c    &        '        t.c. crit mods: ',3x,25x,a3,3x,5(a3,1x))
c1444 format ('       surf ten models: ',a3,3x,5(a3,1x),a5,3x,5(a3,1x))
c     write (*,*) 'ierr at SETUP exit: ',ierr,':  ',herr(1:40)
c
      call AMH2OR  !Setup ammonia/water reference state
      do icomp=1,nc
        call SETPRCO(icomp)
      enddo
c
      RETURN
c
      end                                              !subroutine SETUP
c
c ======================================================================
c
      subroutine SETUP0 (i,hfld,hfm,hrf,ierr,herr)
c
c     call the SETUP routine with the same inputs as SETUP except for
c     the hfld variable.  This subroutine is generally used in calls
c     to the REFPROP DLL since the call cannot handle a string array.
c     The hfld variable is a string of length 10000.  For a pure fluid,
c     it simply contains the name of the fluid file (with a path if needed).
c     For a mixture, it contains the names of the constituents in the
c     mixture separated by a |.  For the air mixture, this would be
c     something like (depends on the need for paths):
c     hfld='fluids\nitrogen.fld|fluids\argon.fld|fluids\oxygen.fld|'
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  04-07-05 EWL, original version
c  01-16-06 EWL, check for missing '|' at end of line
C
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SETUP0
c     dll_export SETUP0
c
      parameter (ncmax=20)
      character hfld*10000,hfm*255,hrf*3,herr*255,hfmix*255
      character*255 hf(ncmax)
      do j=1,ncmax
        hf(j)=' '
      enddo
      if (i.eq.1) then
        hf(1)=hfld
        j=index(hf(1),'|')
        if (j.ne.0) hf(1)(j:)=' '
      else
        m=1
        do k=1,i
          j=index(hfld(m:10000),'|')
          if (j.ne.0) then
            j=j+m-1
            hf(k)=hfld(m:j-1)
            m=j+1
c  check for inputs without a '|' at the end
          elseif (k.eq.i .and. hfld(m:).ne.' ') then
            hf(k)=hfld(m:)
          else
            ierr=1
            herr='[SETUP0 error 1] not enough fluid names on input line'
            RETURN
          endif
        enddo
      endif
      hfmix=hfm
      call SETUP (i,hf,hfmix,hrf,ierr,herr)
      end                                             !subroutine SETUP0
c
c ======================================================================
c
      subroutine SETFLD (icomp,hfile,ierr,herr)
c
c  open a fluid file and read model coefficients (or get from block data)
c
c  inputs:
c    icomp--pointer specifying component number
c           zero and negative values are used for ECS reference fluid(s)
c    hfile--array of file names specifying fluid/mixture components
c           [character*255 variable] for each of the components;
c           --or--
c           when hf(i) is of the form:
c             BDATA:nn-nn-nn
c           use coefficients stored in block data for fluid with CAS
c           number specified by nn-nn-nn;
c           e.g. to use stored formulation for R134a,
c           hf(i) = 'BDATA:811-97-2'
c  outputs:
c     ierr--error flag:  0 = successful
c                      101 = error in opening file
c                      102 = error in file or premature end of file
c                     -103 = unknown model encountered in file
c                      104 = error in setup of model
c                      105 = specified model not found
c                      106 = cp0 equation not found
c     herr--error string (character*255 variable if ierr<>0)
c     [fluid parameters, etc. returned via various common blocks]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  12-08-95  MM, original version--extracted from old subroutine SETUP
c  12-12-95  MM, add calls to SETECS, SETETA, SETTCX
c  12-14-95  MM, add call to RFFILE (add directory to ref fluid file)
c  01-12-96  MM, do not reset hmxeos(0) if in file 0 (reference fluid)
c  01-23-96  MM, always call SETETA, SETTCX, etc. if found in .fld file
c  01-31-96  MM, implement 'NBS' option for EOS (and transport model)
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                add Zcrit to common /CCON/
c  02-29-96  MM, check only hmxeos(i) [not heos] for call to SETxxx
c  03-13-96  MM, changes to merge with Klein's transport code
c                add common /TRNMOD/; add calls to SETTRN, SETVSi, SETTCi
c  03-15-96  MM, always load ECS-transport model (unless htran = 'NUL')
c  03-19-96  MM, read dipole moment and add to /CCON/
c  03-20-96  MM, add mixture models to /TRNMOD/
c  03-21-96  MM, replace /MODEL/ with /EOSMOD/, /STNMOD/
c  03-27-96  MM, add calls to set up surface tension
c  03-28-96  MM, add flags ltcx,leta to cover case of fluid-specific
c                model available for only one of the transport props
c  06-17-96  MM, check only 'CP' rather than 'CPP' to allow CP1
c  06-18-96  MM, add calls to SETVS4-6 and SETTC2-4 for future use
c  10-16-96  MM, add call to SETCI2 (collision integral AUX function)
c  10-30-96  MM, add calls to SETTC5, SETTC6
c  11-06-96  MM, check that i.ge.nrf0 before calling collision integral
c  01-21-97  MM, add call to SETCI1 (collision integral AUX function)
c  02-06-97  MM, add /CNAM80/ to store full chemical name
c  02-20-97  MM, read in default reference state from fluid file,
c                add hsvrfd to /REFST/, new common /CREFDF/
c  02-24-97  MM, add calls to transport critical models and /CREMOD/
c  03-26-97 EWL, minor changes for Fortran 90 compatibility
c  03-27-97 EWL, if file not found, search in likely directories
c  08-22-97  MM, use double backslash on file open to avoid problem with Unix machines
c  10-28-97  MM, put fluid file version no. into common
c  12-01-97  MM, read in synonyms, add to /CNAM80/
c  03-30-98  MM, initialize ierrvs,ierrtc,ierrst
c  04-08-98 EWL, initialize ierr1
c  10-30-98 EWL, change hmxeos to FEQ when called by SETMOD with FE1, etc.
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  06-22-99 EWL, comment out the "ierr1=104" in the AUX section
c  11-01-99 EWL, add error if cp0 equation not found
c                load CPP as default (in case flag to CPP occurs after coefs. in file)
c  01-14-00 EWL, allow use of 'VIS' as well as 'ETA'
c  01-26-00  MM, read in citation in addition to model specification
c                add common /CITE/ to pass this information to GETMOD
c  02-15-00 EWL, check for case sensitive fluid names
c  03-07-00 EWL, add TC0 and VS0 models
c  05-24-00 EWL, initialize hetemp and httemp
c  08-23-00 EWL, change heos from FE1, etc. to FEQ, and similarly for BWR and ECS
c  01-17-01 EWL, remove double backslash, which caused problems in Windows
c  01-23-02 EWL, split common block CITE into 2 pieces for some compilers
c  07-08-02 EWL, revise checks on location of fld file
c  01-13-05 WLJ, add new checks for open failure (ORNL)
c  07-22-05 MLH, allow for transport reference fluids
c  09-28-05  DT, change occurrences of a backslash to 'char(92)' for compatibility
c                with some compilers using F77
c  04-18-06 EWL, make call to "inquire" to check for multiple use of nread=12
c  07-24-06 EWL, do not call setfld again if @FEQ found after #FE1, etc.
c  08-10-06 EWL, create new subroutine to open file and move appropriate code
c  10-30-06 EWL, check for FEK
c  11-08-06 EWL, after the call to STFLD2, add href to hfiles() in the negative position so that the transport routines of the reference fluid will be loaded.
c  04-18-07 EWL, add checks for UN parameters
c  12-12-07 MLH, allow FEX for generalized model of Xiang and Deiters
c  01-09-08 EWL, read in heat of combustion
c  02-26-09 EWL, check for missing '.fld' on file names
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (nrf0=n0)         !lower limit for transport ref fluid arrays
      character*1 htab,hnull
      character*1 h1,hstar
      character*3 hsveqn,heta,htcx,hetak,htcxk
      character*3 hetemp,httemp
      character*3 hpheq,heos,hmxeos,hmodcp
      character*3 hsten,hstenk
      character*3 hdiel,hdielk
      character*3 hmelt,hmeltk
      character*3 hsubl,hsublk
      character*3 hps,hpsk,hpl,hplk,hdl,hdlk,hdv,hdvk
      character*3 hmaux
      character*3 hflag,htype
      character*3 hsvrfs,hrefdf
      character*3 hetacr,htcxcr
      character*12 hcasn(n0:nx),hcas,hname
      character*255 hfile(n0:nx),href,hnam80,hsyn1,hsyn2
      character*255 family,UNNumb
      character*255 herr,herr1
      character*251 hcite
      character*255 hieos,hicp0,hieta,hietac,hitcx,hitcxc,histn,
     &              hidiel,himelt,hisubl,hips,hipl,hidl,hidv
      character*255 hstr
      integer unflag
      logical lbdata,leta,ltcx,leosfl,lprfl
c
      common /NCOMP/ ncomp,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /TRNMOD/ heta,hetak(nrf0:nx),htcx,htcxk(nrf0:nx)
      common /CREMOD/ hetacr(nrf0:ncmax),htcxcr(nrf0:ncmax)
      common /STNMOD/ hsten,hstenk(n0:nx)
      common /DEMOD/ hdiel,hdielk(n0:nx)
      common /MELTMOD/ hmelt,hmeltk(n0:nx)
      common /SUBLMOD/ hsubl,hsublk(n0:nx)
      common /PSMOD/ hps,hpsk(n0:nx)
      common /PLMOD/ hpl,hplk(n0:nx)
      common /DLMOD/ hdl,hdlk(n0:nx)
      common /DVMOD/ hdv,hdvk(n0:nx)
      common /REFST/ hsvrfs,hrefdf(n0:nx)
      common /CREFDF/ tdef(n0:nx),pdef(n0:nx),hdef(n0:nx),sdef(n0:nx)
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /CCAS/ hcas(n0:nx)
      common /CNAM/ hname(n0:nx)
      common /CNAM80/ hnam80(n0:nx),hsyn1(n0:nx),hsyn2(n0:nx)
      common /HCHAR/ htab,hnull
      common /VERS/ verfl(n0:nx),vermx    !fluid & mix file version nos.
      common /FAML/ family(n0:nx),UNNumb(n0:nx)  !fluid family type
      common /CITE/  hieos(n0:nx),hicp0(n0:nx),
     &               hieta(nrf0:nx),hietac(nrf0:nx),
     &               hitcx(nrf0:nx),hitcxc(nrf0:nx),histn(n0:nx)
      common /CITE2/ hidiel(n0:nx),himelt(n0:nx),hisubl(n0:nx),
     &               hips(n0:nx),hipl(n0:nx),hidl(n0:nx),hidv(n0:nx)
      common /CONPR/ prcoef(0:nx,20),iprflag(0:nx)
      common /CONUN/ tminUN(n0:nx),tmaxUN(n0:nx),pmaxUN(n0:nx),
     &               rhomUN(n0:nx),prmUN(n0:nx,10),ntrmUN(n0:nx),
     &               unflag(n0:nx)
      common /HTCM/ hcmbst(n0:nx)
      common /ALTFEQ/ igenfl(n0:nx)
c
 10   continue
      ierr1=0
      herr=' '
      herr1=hnull
      i=icomp
      icpflg=0
      lprfl=.false.
      ltcx=.false.
      leta=.false.
      leosfl=.false.
      hetemp=' '
      httemp=' '
      igenfl(i)=0
      UNNumb(i)=' '
      family(i)=' '
      hcmbst(i)=-1.d0
c     initialize
      if (i.ge.0) unflag(i)=0
c     do not set pr coef for transport/ecs ref fluids other than zero slot
      if (i.ge.0) then
        do j=1,20
          prcoef(i,j)=0.d0
        enddo
        iprflag(i)=0
      endif
c     write (*,1102) i,hfile(i)
c1102 format (/1x,'input file for fluid ',i3,': (',a40,')')
c     write (*,1003) heos,(hmxeos(j),j=n0,ncmax),
c    &               heta,(hetak(j),j=nrf0,nx),htcx,(htcxk(j),j=nrf0,nx)
c    &              ,hsten,(hstenk(j),j=n0,ncmax)
c1003 format ('  SETFLD--heos, hmxeos: ',a3,3x,5(a3,1x),a5,3x,5(a3,1x)/
c    &        '  (input)  visc models: ',a3,25x,a3,3x,5(a3,1x)/
c    &        '            tcx models: ',a3,25x,a3,3x,5(a3,1x)/
c    &        '       surf ten models: ',a3,3x,5(a3,1x),a5,3x,5(a3,1x))
      if (hfile(i)(1:5).eq.'BDATA' .or. hfile(i)(1:5).eq.'bdata') then
c  get coefficients from common blocks
        lbdata=.true.
        nread=0
        hcasn(i)=hfile(i)(7:18)
        hcas(i)=hcasn(i)
c       write (*,1104) hcasn(i)
c1104   format (1x,'coefficients from block data; CAS no.: (',a12,')')
      else
c  read coefficients from file
c       write (*,1105) hfile(i)
c1105   format (1x,'reading coefficients from file: (',a40,')')
        lbdata=.false.
        nread=12             !logical unit for file reads
        if (index(hfile(i),'.').eq.0) then   !Check for missing '.fld'
          j=index(hfile(i),' ')
          if (j.ne.0 .and. hfile(i)(j:j+3).eq.'    ')
     &        hfile(i)=hfile(i)(1:j-1)//'.fld'
        endif
        call OPENFL (nread,hfile(i),1,ierr,herr)
        if (ierr.ne.0) goto 998
c
        read (nread,2012) hname(i)           !short name
        read (nread,2012) hcasn(i)           !CAS number
        hcas(i)=hcasn(i)
        read (nread,2080) hnam80(i)          !full chemical name
c       write (*,*) ' SETUP--full name:  ',hnam80(i)
        read (nread,2080) hsyn1(i)           !synonym 1
        read (nread,2080) hsyn2(i)           !synonym 2
c  note:  read in generic fluid constants; these may be reset by
c         some SETeos routines to correspond with the values used
c         in the respective models
        read (nread,*) wm(i)           !molecular weight [g/mol]
        read (nread,*) ttp(i)          !triple point temperature [K]
        read (nread,*) tnbp(i)         !normal boiling point [K]
        read (nread,*) tc(i)           !critical temperature [K]
        read (nread,*) pc(i)           !critical pressure [kPa]
        read (nread,*) rhoc(i)         !critical density [mol/L]
        read (nread,*) accen(i)        !acentric factor [-]
        read (nread,*) dipole(i)       !dipole moment [debye]
        read (nread,2003) hrefdf(i)    !default reference state
        if (hrefdf(i)(1:2).eq.'OT') then
c  for 'OTHer' reference state read in reference T,p,h,s
          read (nread,*)  tdef(i),pdef(i),hdef(i),sdef(i)
        end if
        read (nread,*) verfl(i)        !version number
c
 150    continue
        read (nread,2080) hstr
        if (hstr.ne.' ') then
          j=INDEX(hstr,'!')
          if (j.ne.0) then
            if (hstr(j:j+2).eq.'!UN') then
              UNNumb(i)=hstr(1:j-1)
            elseif (hstr(j:j+2).eq.'!fa') then
              family(i)=hstr(1:j-1)
            elseif (hstr(j:j+2).eq.'!he') then
              read (hstr(1:j-1),*) hcmbst(i)
            endif
          endif
          goto 150
        endif
c       write (*,*) ' SETUP--end of constants; default = ',hsvrfd(i)
      end if
c
c  search for key characters in cols 1 and 2-4 (if reading from file)
c  and branch to subsidiary setup routines for specified models
c
 100  continue
      if (.not.lbdata) then           !read from file
        read (nread,2013,end=199,err=199) hstar,hflag
        if (hstar.eq.'!') goto 100
c       write (*,*) ' SETUP--input line:   ',hstar,hflag
c  hstar = '#' or '@' indicates start of model specification
c        # indicates NIST-recommended model
c        @ indicates other model(s)
c  hflag indicates type of model:
c        'EOS' = equation of state
c        'TRN' = transport property model (i.e. ECS model)
c        'ETA' = pure fluid viscosity model
c        'TCX' = pure fluid thermal conductivity model
c        'AUX' = auxiliary model, such as ideal heat capacity
        if (hstar.eq.'#' .or. hstar.eq.'@') then
          read (nread,2083) htype,hcite
c  htype is 3-letter key for particular models + 1-line citation, for example:
c         'NBS' = NIST-recommended model
c         'BWR' = modified Benedict-Webb-Rubin equation of state
c         'FEQ' = fundamental (Helmholtz) equation of state
c         'ECS' = extended corresponding states model
c         'HMX' = mixture Helmholtz model
c         'CPP' = polynomial fit of ideal gas heat capacity
c         'VSi' = pure fluid viscosity model #i
c         'TCi' = pure fluid thermal conductivity model #i
c         write (*,*) ' SETUP--hstar,hflag:  ',hstar,hflag
c         write (*,*) ' SETUP--htype:        ',htype
c
          do k=1,1000             !skip over source comments
            read (nread,2001) h1
            if (h1.ne.'?') goto 120   !comment block terminated by '!'
          enddo
 120      continue
        end if
      end if
c
      if (lbdata) then
        hflag='EOS'
        htype=heos
      end if
c
      if (hflag.eq.'EOS' .or. hflag.eq.'eos') then
c
c  set up equation of state
c
        if (hstar.eq.'#') then
c  special case for NIST-recommended model
          if (ncomp.eq.1 .and. heos.eq.'NBS') heos=htype
          if (i.ge.1 .and. hmxeos(i).eq.'NBS') hmxeos(i)=htype
c  allow for transport reference fluids
          if (i.lt.-ncmax .and. hmxeos(i).eq.'NBS') hmxeos(i)=htype
        end if
c
        if (htype(1:2).eq.'BW') then
c  modified Benedict-Webb-Rubin equation of state found in file
c         if (heos.eq.'BWR' .or. hmxeos(i).eq.'BWR') then
          if (hmxeos(i).eq.htype) then
c           write (*,*) ' SETUP--about to call SETBWR'
            if (.not.leosfl) then
              hmxeos(i)='BWR'
              if (heos(1:2).eq.'BW') heos='BWR'
              call SETBWR (nread,i,hcasn(i),ierr1,herr1)
              hieos(i)=htype//hcite//hnull
              leosfl=.true.
            end if
          end if
        else if (htype(1:2).eq.'FE') then
c  fundamental (Helmholtz) equation of state found in file
c         if (heos.eq.'FEQ' .or. hmxeos(i).eq.'FEQ') then
          if (htype(1:3).eq.'FEX') igenfl(i)=1
          if (hmxeos(i).eq.htype) then
            if (.not.leosfl) then
c           write (*,*) ' SETUP--about to call SETFEQ'
              if (heos(1:2).eq.'FE') heos='FEQ'
              call SETFEQ (nread,i,hcasn(i),ierr1,herr1)
              hieos(i)=htype//hcite//hnull
              leosfl=.true.
            end if
          end if
        else if (htype(1:2).eq.'QU') then
c  fundamental Quintic (Helmholtz) equation of state found in file
c         if (heos.eq.'QUI' .or. hmxeos(i).eq.'QUI') then
          if (hmxeos(i).eq.htype) then
            if (.not.leosfl) then
c           write (*,*) ' SETUP--about to call SETQUI'
              if (heos(1:2).eq.'QU') heos='QUI'
              call SETQUI (nread,i,hcasn(i),ierr1,herr1)
              hieos(i)=htype//hcite//hnull
              leosfl=.true.
            end if
          end if
        else if (htype(1:2).eq.'EC') then
c  extended corresponding states (ECS) thermo model found in file
c         if (heos.eq.'ECS' .or. hmxeos(i).eq.'ECS') then
          if (hmxeos(i).eq.htype) then
c           write (*,*) ' SETUP--about to call SETECS (thermo)'
            if (.not.leosfl) then
              hmxeos(i)='ECS'
              if (heos(1:2).eq.'EC') heos='ECS'
              call SETECS (nread,i,hcasn(i),href,hsveqn,ierr1,herr1)
c  add directory information to reference fluid file
              call RFFILE (hfile(i),href)
              iref=-i
              hfile(iref)=href
              hmxeos(iref)=hsveqn
              hieos(i)=htype//hcite//hnull
              leosfl=.true.
            end if
          end if
        else if (htype(1:2).eq.'PR') then
c  cubic (Peng-Robinson, etc.) equation of state found in file
          if (hmxeos(i)(1:2).eq.htype(1:2)) then
            hmxeos(i)='PR'
            if (heos(1:2).eq.'PR') heos='PR'
            call SETPR (nread,i,hcasn(i),ierr1,herr1)
            hieos(i)=htype//hcite//hnull
            leosfl=.true.
          else
c  get tpr coefficients from fluid file, but nothing else
            call SETPR (-nread,i,hcasn(i),ierr1,herr1)
          end if
          lprfl=.true.
        else
c  unidentified model found in file
          ierr=-103
          write (herr,2103) htype,i,hnull
          call ERRMSG (ierr,herr)
        end if
        if (ierr1.ne.0) then
          ierr=104
          write (herr,2104) hmxeos(i),i,herr1(1:182),hnull
          call ERRMSG (ierr,herr)
        end if
      end if
c
      if (lbdata) then
        hflag='AUX'
        htype=hmodcp(i)
      end if
c
      if (hflag.eq.'AUX') then
c
c  set up auxiliary model(s)
c  (model(s) are specified in call to appropriate SETmod routine)
c
c       write (*,*) ' SETUP--aux model (',htype,') found in fld file '
c       write (*,*) ' SETUP--i,htype,hmodcp(i): ',i,htype,hmodcp(i)
        ierr1=0
        hmaux=hmodcp(i)    !possible use in error message
        if (htype(1:2).eq.'CP' .and. hmodcp(i).eq.htype) then
c  ideal gas heat capacity function
c         write (*,*) ' SETUP--about to call SETCPP'
          call SETCPP (nread,i,hcasn(i),ierr1,herr1)
          hicp0(i)=htype//hcite//hnull
          icpflg=2
        else if (htype(1:2).eq.'PH' .and. hmodcp(i).eq.htype) then
c  Helmholtz for ideal gas heat state
c         write (*,*) ' SETUP--about to call SETPH0'
          call SETPH0 (nread,i,hcasn(i),ierr1,herr1)
          hicp0(i)=htype//hcite//hnull
          icpflg=2
        elseif (htype(1:3).eq.'CPP' .and. (hmodcp(i)(1:2).eq.'CP'
     &         .or. hmodcp(i).eq.'NBS') .and. icpflg.eq.0) then
          call SETCPP (nread,i,hcasn(i),ierr1,herr1)
          hicp0(i)=htype//hcite//hnull
          icpflg=1
        elseif (htype(1:2).eq.'UN') then
          read (nread,*) tminUN(i)             !lower temperature limit
          read (nread,*) tmaxUN(i)             !upper temperature limit
          read (nread,*) pmaxUN(i)             !upper pressure limit
          read (nread,*) rhomUN(i)             !upper density limit
          read (nread,*) unflag(i)             !flag to allow calcs
          read (nread,*) ntrmUN(i)
          jterm=0
          if (ntrmUN(i).ge.1) then
            do j=1,ntrmUN(i)           !read coefficients
              jterm=jterm+1
              read (nread,*) prmUN(i,jterm)
            enddo
          end if
        end if
        if (ierr1.ne.0) then
          ierr=104
          write (herr,2104) hmaux,i,herr1(1:182),hnull
          call ERRMSG (ierr,herr)
        end if
      end if
c
      if (lbdata) hflag='TRN'
c
c  If the transport models, dielectric models, and surface tension eqs.
c  (vis, tcx, st, de) are not required, and only the default values in the
c  fluid files are used, the following files do not have to be compiled and
c  linked with the other source code:  CORE_DE, CORE_STN, SETUP2,
c  TRNS_ECS, TRNS_TCX, TRNS_VIS, and TRNSP.  In addition, the
c  following two lines should be commented out:
      call STFLD2 (nread,i,hcasn,hcite,hsveqn,hetemp,hfile,hflag,
     &              href,hstar,httemp,htype,leta,ltcx,ierr,herr)
      if (href.ne.' ' .and. i.gt.0) then
        iref=-i
        if (hfile(iref).eq.'NUL') then
          hfile(iref)=href
          hmxeos(iref)=hsveqn
        endif
      endif
c
c  do not call the following models if component number less than 0
      if (i.ge.nrf0 .and. hstar.eq.'#') then
c
c  melting line models
        ierr1=0
        if (hflag.eq.'MLT' .or. hflag.eq.'mlt') then
          hmelt=htype
          hmeltk(i)=htype
          call SETMLT (nread,i,ierr1,herr1)
          himelt(i)=htype//hcite//hnull
c  sublimation line models
        elseif (hflag.eq.'SBL' .or. hflag.eq.'sbl') then
          hsubl=htype
          hsublk(i)=htype
          call SETSBL (nread,i,ierr1,herr1)
          hisubl(i)=htype//hcite//hnull
c  vapor pressure equations
        elseif (hflag.eq.'PS ' .or. hflag.eq.'ps ') then
          hps=htype
          hpsk(i)=htype
          call SETPS (nread,i,ierr1,herr1)
          hips(i)=htype//hcite//hnull
c  liquid pressure equations (for pseudo-pure fluids)
        elseif (hflag.eq.'PL ' .or. hflag.eq.'pl ') then
          hpl=htype
          hplk(i)=htype
          call SETPL (nread,i,ierr1,herr1)
          hipl(i)=htype//hcite//hnull
c  saturated liquid density equations
        elseif (hflag.eq.'DL ' .or. hflag.eq.'dl ') then
          hdl=htype
          hdlk(i)=htype
          call SETDL (nread,i,ierr1,herr1)
          hidl(i)=htype//hcite//hnull
c  saturated vapor density equations
        elseif (hflag.eq.'DV ' .or. hflag.eq.'dv ') then
          hdv=htype
          hdvk(i)=htype
          call SETDV (nread,i,ierr1,herr1)
          hidv(i)=htype//hcite//hnull
        end if
        if (ierr1.ne.0 .and. ierr.eq.0) then
          ierr=104
          write (herr,2104) hflag,i,herr1(1:182),hnull
          call ERRMSG (ierr,herr)
        end if
      end if
c
      if (lbdata) RETURN                  !next component
c
      if (hflag.eq.'END') then
c  if fluid-specific 'ETA' and/or 'TCX' models have not been loaded, use
c  the 'TRN' model (e.g. the ECS model)
        if (i.ge.nrf0) then
          if (.not.leta) hetak(i)=hetemp
          if (.not.ltcx) htcxk(i)=httemp
        end if
c  close the file, return to setup and process the next component
        rewind (nread)
        close (unit=nread,err=998)
c       write (*,1190) hfile(i)
c1190   format ('  SETUP--file closed: ',a80)
        if (.not.leosfl) then
        if (i.ge.-ncmax) then
c  if a pure fluid equation for the GERG-2004 model was requested, but does
c  not exist, use the default (FEQ) without returning an error message
          if (hmxeos(i).eq.'FEK') then
            hmxeos(i)='FEQ'
            goto 10
          endif
          ierr=105
          write (herr,2105) i,hfile(i)(1:100),hnull
 2105     format ('[SETUP error 105] selected equation of state',
     &            ' not found for component #',i3,':  ',a100,a1)
          call ERRMSG (ierr,herr)
          RETURN
        endif
        if (icpflg.eq.0 .or. (icpflg.eq.1.and.hmodcp(i).ne.'CPP')) then
          ierr=106
          write (herr,2106) i,hnull
 2106     format ('[SETUP error 106] ideal gas heat capacity equation',
     &            ' not found for component #',i3,'.',a1)
          call ERRMSG (ierr,herr)
        endif
        endif
c  call SETPR to put constants into PR common block if PRT not in FLD file
        if (.not.lprfl) call SETPR(999,i,hcasn(i),ierr1,herr1)
        RETURN
      else
        goto 100
      end if
 199  continue
      ierr=102
      write (herr,2102) i,hnull
 2102 format ('[SETUP error 102] error in file or premature end of',
     &        ' file for component #',i3,'.',a1)
      call ERRMSG (ierr,herr)
c  rewind and close the file
      rewind (nread)
      close (unit=nread,err=998)
c     write (*,1199) hfile(i)
c1199 format (1x,'SETFLD--premature end of file: ',a80)
c  if fluid-specific 'ETA' and/or 'TCX' models have not been loaded, use
c  the 'TRN' model (e.g. the ECS model)
      if (i.ge.nrf0) then
        if (.not.leta) hetak(i)=hetemp
        if (.not.ltcx) htcxk(i)=httemp
      end if
      RETURN
c
 998  continue
      ierr=101
      write (herr,2101) i,hfile(i)(1:80),hnull
 2101 format ('[SETUP error 101] error in opening file for component #',
     &        i3,'; filename = (',a80,').',a1)
      call ERRMSG (ierr,herr)
      RETURN
c
 2001 format (a1)
 2003 format (a3)
 2012 format (a12)
 2013 format (a1,a3)
 2080 format (a255)
 2083 format (a3,a251)
 2103 format ('[SETUP warning -103] unknown model (',a3,
     &        ') encountered in file for component #',i3,'.',a1)
 2104 format ('[SETUP error 104] Error in setup of (',a3,
     &        ') model for component #',i3,':  ',a182,a1)
c
      end                                             !subroutine SETFLD
c
c ======================================================================
c
      subroutine SETMOD (nc,htype,hmix,hcomp,ierr,herr)
c
c  set model(s) other than the NIST-recommended ('NBS') ones
c
c  This subroutine must be called before SETUP; it need not be called
c  at all if the default (NIST-recommended) models are desired.
c
c  inputs:
c       nc--number of components (1 for pure fluid) [integer]
c    htype--flag indicating which models are to be set [character*3]
c           'EOS':  equation of state for thermodynamic properties
c           'ETA':  viscosity
c           'TCX':  thermal conductivity
c           'STN':  surface tension
c           'NBS':  reset all of the above model types and all
c                   subsidiary component models to 'NBS';
c                   values of hmix and hcomp are ignored
c     hmix--mixture model to use for the property specified in htype [character*3];
c           ignored if number of components = 1
c           some allowable choices for hmix:
c             'NBS':  use NIST recommendation for specified fluid/mixture
c             'HMX':  mixture Helmholtz model for thermodynamic properties
c             'ECS':  extended corresponding states for viscosity or therm. cond.
c             'STX':  surface tension mixture model
c    hcomp--component model(s) to use for property specified in htype [array (1..nc) of character*3]
c             'NBS':  NIST recommendation for specified fluid/mixture
c           some allowable choices for an equation of state:
c             'FEQ':  Helmholtz free energy model
c             'BWR':  pure fluid modified Benedict-Webb-Rubin (MBWR)
c             'ECS':  pure fluid thermo extended corresponding states
c           some allowable choices for viscosity:
c             'ECS':  extended corresponding states (all fluids)
c             'VS1':  the 'composite' model for R134a, R152a, NH3, etc.
c             'VS2':  Younglove-Ely model for hydrocarbons
c             'VS4':  Generalized friction theory of Quinones-Cisneros and Deiters
c             'VS5':  Chung et al. (1988) predictive model
c           some allowable choices for thermal conductivity:
c             'ECS':  extended corresponding states (all fluids)
c             'TC1':  the 'composite' model for R134a, R152a, etc.
c             'TC2':  Younglove-Ely model for hydrocarbons
c             'TC5':  Chung et al. (1988) predictive model
c           some allowable choices for surface tension:
c             'ST1':  surface tension as f(tau); tau = 1 - T/Tc
c
c  outputs:
c     ierr--error flag:  0 = successful
c                      113 = nc outside of bounds
c     herr--error string (character*255 variable if ierr<>0)
c     [fluid parameters, etc. returned via various common blocks]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  03-20-96  MM, original version
c  03-21-96  MM, replace /MODEL/ with /EOSMOD/, /STNMOD/
c  03-28-96  MM, add option if htype='NBS' reset all models
c  03-28-97  MM, add initialization of hnull (repeat of SETUP initialization);
c                add /RESETM/ to check for initialization of models
c  03-31-97  MM, reset all models if lreset true, even for htype<>'NBS'
c  04-10-97  MM, add and initialize /CREMOD/ (critical enhancements)
c  10-01-97  MM, add compiler switch to allow access by DLL
c  10-30-98 EWL, change heos to FEQ when called by SETMOD with FE1, etc.
c  03-11-99 EWL, reset lsatt and lsatp to .false. when called
c  01-14-00 EWL, allow use of 'VIS' as well as 'ETA'
c  07-24-06 EWL, remove excess code
c  07-24-06 EWL, only use hmix if nc<>1
c  12-12-07 MLH, add generalized Helmholtz model FEX of Xiang and Deiters (2007)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SETMOD
c     dll_export SETMOD
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (nrf0=n0)    !lower limit for transport ref fluid arrays
      character*1 htab,hnull
      character*3 htype,hmix,hcomp(1:ncmax)
      character*3 heos,hpheq,hmxeos,hsten,hstenk,hmodcp
      character*3 hdiel,hdielk,hmelt,hmeltk,hsubl,hsublk
      character*3 hps,hpsk,hpl,hplk,hdl,hdlk,hdv,hdvk
      character*3 heta,hetak,htcx,htcxk
      character*3 hetacr,htcxcr
      character*3 hmdeta,hmdtcx
      character*3 hsav,hsveos
      character*255 herr
      integer unflag
      logical lreset,lrst,linit
c
      common /NCOMP/ ncomp,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /TRNMOD/ heta,hetak(nrf0:nx),htcx,htcxk(nrf0:nx)
      common /CREMOD/ hetacr(nrf0:ncmax),htcxcr(nrf0:ncmax)
      common /STNMOD/ hsten,hstenk(n0:nx)
      common /DEMOD/ hdiel,hdielk(n0:nx)
      common /MELTMOD/ hmelt,hmeltk(n0:nx)
      common /SUBLMOD/ hsubl,hsublk(n0:nx)
      common /PSMOD/ hps,hpsk(n0:nx)
      common /PLMOD/ hpl,hplk(n0:nx)
      common /DLMOD/ hdl,hdlk(n0:nx)
      common /DVMOD/ hdv,hdvk(n0:nx)
      common /HCHAR/ htab,hnull
      common /RESETM/ lreset,lrst,linit,ncset !flag indicating need to reset models
      common /FSHSAV/ lsatt,lsatp
      common /OMGMOD/ hmdeta(nrf0:nx),hmdtcx(nrf0:nx)
      common /PRHSAV/ hsav,hsveos(n0:nx)
      common /ALTFEQ/ igenfl(n0:nx)
      common /CONUN/ tminUN(n0:nx),tmaxUN(n0:nx),pmaxUN(n0:nx),
     &               rhomUN(n0:nx),prmUN(n0:nx,10),ntrmUN(n0:nx),
     &               unflag(n0:nx)
c
c  repeat initialization of following from SETUP in case SETMOD is called
c  before first call to SETUP
      hnull=CHAR(0)        !null character to terminate error strings
      hsav=' '
      do i=0,nc
        ntrmUN(i)=0
        unflag(i)=0
      enddo
      do i=n0,nx
        hsveos(i)=' '
      enddo
      lsatt=.false.
      lsatp=.false.
      ierr=0
      herr=' '
      if (nc.lt.1 .or. nc.gt.ncmax) then
        ierr=113
        herr='[SETUP error 113] number of components out of bounds'
        call ERRMSG (ierr,herr)
        RETURN
      endif
      ncset=nc
c
c     write (*,1000) nc,htype,hmix,(hcomp(i),i=1,nc)
c1000 format (1x,' SETMOD--input nc,htype,hmix,hcomp:  ',i3,7(1x,a3))
      if (htype.eq.'NBS' .or. htype.eq.'nbs' .or. lreset) then
c  reset all models
c       write (*,*) ' SETMOD--about to reset all models to NBS'
        heos='NBS'
        hpheq='NBS'   !phase equilibria model, not currently implemented
        heta='NBS'
        htcx='NBS'
        hsten='NBS'
        hdiel='NBS'
        hmelt='NBS'
        hsubl='NBS'
        hps='NBS'
        hpl='NBS'
        hdl='NBS'
        hdv='NBS'
        do i=n0,ncmax
          hmxeos(i)='NBS'
          hmodcp(i)='NBS'
          hstenk(i)='NBS'
          hdielk(i)='NBS'
          hmeltk(i)='NBS'
          hsublk(i)='NBS'
          hpsk(i)='NBS'
          hplk(i)='NBS'
          hdlk(i)='NBS'
          hdvk(i)='NBS'
        enddo
        do i=nrf0,ncmax
          hetak(i)='NBS'
          hetacr(i)='NUL'
          htcxk(i)='NBS'
          htcxcr(i)='NUL'
          hmdeta(i)='NUL'
          hmdtcx(i)='NUL'
        enddo
        lreset=.false.
      end if
      if (htype.eq.'EOS' .or. htype.eq.'eos') then
c  equation of state specification
        heos=hmix
        if (heos(1:3).eq.'FEX') igenfl(i)=1
        if (nc.eq.1) heos=hcomp(1)
        if (heos(1:2).eq.'FE') heos='FEQ'
        if (heos(1:2).eq.'QU') heos='QUI'
        if (heos(1:2).eq.'BW') heos='BWR'
        if (heos(1:2).eq.'EC') heos='ECS'
        do i=1,nc
          hmxeos(i)=hcomp(i)
        enddo
      else if (htype.eq.'ETA' .or. htype.eq.'eta' .or.
     &         htype.eq.'VIS' .or. htype.eq.'vis') then
c  viscosity specification
        heta=hmix
        if (nc.eq.1) heta=hcomp(1)
        do i=1,nc
          hetak(i)=hcomp(i)
        enddo
      else if (htype.eq.'TCX' .or. htype.eq.'tcx') then
c  viscosity specification
        htcx=hmix
        if (nc.eq.1) htcx=hcomp(1)
        do i=1,nc
          htcxk(i)=hcomp(i)
        enddo
      else if (htype.eq.'STN' .or. htype.eq.'stn') then
c  surface tension specification
        hsten=hmix
        if (nc.eq.1) hsten=hcomp(1)
        do i=1,nc
          hstenk(i)=hcomp(i)
        enddo
      end if
c     write (*,1003) heos,(hmxeos(j),j=n0,ncmax),
c    &               heta,(hetak(j),j=nrf0,nx),
c    &               htcx,(htcxk(j),j=nrf0,nx),
c    &               hsten,(hstenk(j),j=n0,ncmax)
c1003 format ('  SETMOD--heos, hmxeos: ',a3,3x,5(a3,1x),a5,3x,5(a3,1x)/
c    &        '  (exit)   visc models: ',a3,25x,a3,3x,5(a3,1x)/
c    &        '           t.c. models: ',a3,25x,a3,3x,5(a3,1x)/
c    &        '       surf ten models: ',a3,3x,5(a3,1x),a5,3x,5(a3,1x))
c
      RETURN
c
      end                                             !subroutine SETMOD
c
c ======================================================================
c
      subroutine GERG04 (nc,iflag,ierr,herr)
c
c  set the pure model(s) to those used by the GERG 2004 formulation.
c
c  This subroutine must be called before SETUP; it need not be called
c  at all if the default (NIST-recommended) models are desired.
c  To turn off the GERG settings, call this routine again with iflag=0,
c  and then call the SETUP routine to reset the parameters of the equations
c  of state.
c
c  inputs:
c       nc--number of components (1 for pure fluid)
c    iflag--set to 1 to load the GERG 2004 equations, set to 0 for defaults
c  outputs:
c     ierr--error flag:  0 = successful
c     herr--error string returned from SETMOD
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  01-24-07 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: GERG04
c     dll_export GERG04
c
      parameter (ncmax=20)        !max number of components in mixture
      character*255 herr
      character*3 htype,hmix,hcomp(1:ncmax)
      common /GERG2004/ iGERG04
c
      htype='EOS'
      hmix='HMX'
      do i=1,nc
        hcomp(i)='DEF'
        if (iflag.eq.1) hcomp(i)='FEK'
      enddo
      call SETMOD (nc,htype,hmix,hcomp,ierr,herr)
      iGERG04=iflag
      RETURN
      end                                             !subroutine GERG04
c
c ======================================================================
c
      subroutine SETREF (hrf,ixflag,x0,h0,s0,t0,p0,ierr,herr)
c
c  set reference state enthalpy and entropy
c
c  This subroutine must be called after SETUP; it need not be called at
c  all if the reference state specified in the call to SETUP is to be
c  used.
c
c  inputs:
c      hrf--reference state for thermodynamic calculations [character*3]
c           'NBP':  h,s = 0 at normal boiling point(s)
c           'ASH':  h,s = 0 for sat liquid at -40 C (ASHRAE convention)
c           'IIR':  h = 200, s = 1.0 for sat liq at 0 C (IIR convention)
c           'DEF':  default reference state as specified in fluid file
c                   is applied to each component (ixflag = 1 is used)
c           'OTH':  other, as specified by h0, s0, t0, p0 (real gas state)
c           'OT0':  other, as specified by h0, s0, t0, p0 (ideal gas state)
c           '???':  change hrf to the current reference state and exit.
c   ixflag--composition flag:  1 = ref state applied to pure components
c                              2 = ref state applied to mixture x0
c  following input has meaning only if ixflag = 2
c       x0--composition for which h0, s0 apply; array(1:nc) [mol frac]
c           this is useful for mixtures of a predefined composition,
c           e.g. refrigerant blends such as R410A
c  following inputs have meaning only if hrf = 'OTH'
c       h0--reference state enthalpy at t0,p0 {x0} [J/mol]
c       s0--reference state entropy at t0,p0 {x0} [J/mol-K]
c       t0--reference state temperature [K]
c           t0 = -1 indicates saturated liquid at normal boiling point
c                   (bubble point for a mixture)
c       p0--reference state pressure [kPa]
c           p0 = -1 indicates saturated liquid at t0 {and x0}
c           p0 = -2 indicates saturated vapor at t0 {and x0}
c  outputs:
c     ierr--error flag:  0 = successful
c                       22 = Tmin > Tref for IIR reference state
c                       23 = Tcrit < Tref for IIR reference state
c                       24 = Tmin > Tref for ASHRAE reference state
c                       25 = Tcrit < Tref for ASHRAE reference state
c                       26 = Tmin > Tnbp for NBP reference state
c                       27 = Tref, Pref for OTH ref state outside limits
c                      -28 = can't apply 'DEF' to mixture;
c                            will apply to pure components
c                      -29 = unknown reference state specified;
c                            will use 'DEF'
c     herr--error string (character*255 variable if ierr<>0)
c     [fluid parameters, etc. returned via various common blocks]
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  01-11-96  MM, original version--extracted from old subroutine SETUP
c  01-16-96  MM, add check for inputs same as previous call
c  01-19-96  MM, bug fixes and set flag ksetrf to 1 when SETREF called
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c                add Zcrit to common /CCON/
c  03-19-96  MM, add dipole moment to /CCON/
c  03-21-96  MM, replace /MODEL/ with /REFST/
c  06-18-96  MM, check ref state temperature against LIMITX
c  01-07-97  MM, check limits for case of NBP ref state (re: CO2 bug)
c  01-10-97  MM, pass error flag from limits errors to calling routine
c  02-20-97  MM, add 'DEF' option (default ref state read from fld file)
c                add hsvrfd to /REFST/, new common /CREFDF/
c  03-25-97  MM, xliq,xvap used but not dimensioned; change to xl,xv
c  03-26-97 EWL, changes for Fortran 90 compatibility
c  05-14-97  MM, fix bug in checking for inputs same as last call
c  10-01-97  MM, add compiler switch to allow access by DLL
c  10-16-97  MM, allow ixflag = 2 for 'DEF' when all fluids are same
c  11-24-97  MM, check for errors on calls to SATT, SATP, TPFLSH
c  04-24-98  MM, move check for hrf='DEF' with ixflag=2;
c                add check for 'DEF' within mixture-adjust calc
c  08-13-98  MM, initialize "negative" components only for ECS model
c                (do 200 loop); add common /EOSMOD/
c  08-14-98  MM, pass tref to LIMITX (avoid undefined variable error)
c  08-17-98  MM, initialize Dx,px for ixflag=2 case
c  08-21-98  MM, initialize ierrsv,herrsv
c  12-01-98 EWL, add Reos and triple point pressure and density to /CCON/
c  11-25-02 EWL, reset TSAV so that cp0 is recalculated (because Tref changed)
c  11-26-02 EWL, add OTO input when h0 and s0 are given for the ideal gas
c  08-31-04 EWL, only do ixflg=2 option if ncomp.gt.1
c  08-30-05 EWL, add check for OTH when storing reference entropy and enthalpy in common
c  05-31-06 EWL, add ability to find out what the reference state is when hrf='???'
c  09-06-06 EWL, update logic and fix several small bugs
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SETREF
c     dll_export SETREF
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*3 hrf
      character*3 hsvrfs,hrefdf,hrefi(n0:nx),hrdef
      character*3 hpheq,heos,hmxeos,hmodcp
      character*255 herr,herr2,herrsv
      dimension x(ncmax),x0(ncmax),xl(ncmax),xv(ncmax),prefi(n0:nx)
c
      common /NCOMP/ ncomp,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /REFST/ hsvrfs,hrefdf(n0:nx)
      common /HCHAR/ htab,hnull
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
      common /CREFDF/ tdef(n0:nx),pdef(n0:nx),hdef(n0:nx),sdef(n0:nx)
      common /REFSAV/ x0sav(nx),h0sav,s0sav,t0sav,p0sav
      common /IRFSAV/ ixfsav,ksetrf
      common /CPPSAV/ cp0sav(n0:nx),cpisav(n0:nx),cptsav(n0:nx),
     &                tsav(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
c
c  initialize error codes and strings
      ierr=0
      ierrsv=0
      herr=' '
      herrsv=hnull
      do i=1,ncomp
        x(i)=0.d0
      enddo
c     write (*,1000) hrf,ixflag,ncomp
c1000 format (/1x,' SETREF--entering reference state = ',a3,
c    &           '; with ixflag, ncomp =',2i4)
c     write (*,1001) 'input',hrf,ixflag,h0,s0,t0,p0,(x0(i),i=1,ncomp)
c     write (*,1001) 'saved',hsvrfs,ixfsav,h0sav,s0sav,t0sav,p0sav,
c    &                (x0sav(i),i=1,ncomp)
c1001 format ('  SETREF--',a5,' href,xflag,h0,s0,t0,p0: ',a3,i3,4f10.4/
c    &        '          ',5x,'                     x0: ',6x,5f10.4)
c
c  compare inputs to previous values (stored in /REFST/ and /REFSAV/)
c  if nothing has changed SETREF can be bypassed
c
      x(1)=1.d0
      l1same=.false.    !flag for initialization of pure components
      l2same=.false.    !flag for application of ref state to mixture
      ixflg=ixflag      !input value changed if input href not valid
      if (hrf.eq.'???') then
        hrf=hrefdf(1)      !Return the default reference state for a pure fluid
        RETURN
      endif
      if (hrf.eq.hsvrfs) then
        l1same=.true.
        if (hrf(1:2).eq.'OT') then
c  check h0, s0, t0, p0 only if 'OTHer' reference state is used
          if (abs(s0-s0sav).gt.1.0d-10 .or.
     &        abs(h0-h0sav).gt.1.0d-10 .or.
     &        abs(t0-t0sav).gt.1.0d-10 .or.
     &        abs(p0-p0sav).gt.1.0d-10) then
            l1same=.false.
          end if
        end if
        if (ixflag.eq.ixfsav) then
          l2same=.true.
          if (ixflag.eq.2) then
c  check x0 only if ixflag = 2
            do i=1,ncomp
              if (abs(x0(i)-x0sav(i)).gt.1.0d-10) then
                l2same=.false.
              end if
            enddo
          end if
        end if
        if (l1same .and. l2same .and. ksetrf.eq.1) then
c  return if all inputs are same and also SETUP has not been executed
c  since last call to SETREF (as indicated by ksetrf)
c         write (*,*)
c         write (*,*) ' SETREF--inputs are identical, return'
          RETURN
        end if
      end if
c
      if ((hrf.eq.'DEF' .or. hrf.eq.'def') .and. ixflag.eq.2 .and.
     &    ncomp.gt.1) then
c  check if all fluids have same reference state
        ldef=.true.
c       write (*,2001) (hsvrfd(i),i=1,ncomp)
c2001   format (1x,' SETREF--default ref states:  ',5(a3,1x))
        do i=1,ncomp
          if (hrefdf(1).ne.hrefdf(i) .or. hrefdf(i)(1:2).eq.'OT') then
            ldef=.false.
          end if
        enddo
        if (.not.ldef) then
          ierrsv=-28
          herrsv='[SETREF warning -28] cannot apply reference '
     &         //'state (DEF) at a specified mixture composition '
     &         //'when components have different default reference '
     &         //'states; will apply defaults to each pure component'
     &         //hnull
          call ERRMSG (ierrsv,herrsv)
          ixflg=1
        end if
      end if
c
c     write (*,*)
      if (l1same .and. ksetrf.eq.1 .and. ixflg.ge.ixfsav) then
c  can skip initialization of pure components if:
c  --reference state is the same (as indicated by l1same)
c  --SETUP has not been called since last call to SETREF (indicated by ksetrf)
c  --value of ixflag is >= ixfsav, i.e. this call applies ref state to mixture
c    while last call applied to the pures, or if both applied to pures,
c    but not if last call was applied to mix, this call to pures
c       write (*,*) ' SETREF--skip initialization of pure components'
      else
c  must initialize pure component reference states
        do i=n0,ncomp
c  following is temporary assignment, needed for  SATP and/or THERM
c  (SATP calls ENTRO, and ENTRO needs rhoref, but SATP is used to
c  calculate rhoref; thus need to put in temporary value here)
          rhoref(i)=1.0d0
        enddo
c
c  initialize reference temperature, entropy, and enthalpy
c
c       write (*,*) ' SETREF--240 loop from/to: ',-ncomp,ncomp
        do i=-ncomp,ncomp
c  "negative" components refer to ECS reference fluids, for other models
c  the "negative" components are not defined
c       write (*,*) 'SETREF--start init for i,model: ',i,hmxeos(i)
          if (i.ge.1 .or. hmxeos(i).ne.'NUL') then
c         write (*,*) 'SETREF--start init for i,model: ',i,hmxeos(i)
            if (hrf.eq.'IIR' .or. hrf.eq.'ASH' .or. hrf.eq.'NBP' .or.
     &          hrf(1:2).eq.'OT' .or. hrf.eq.'NA') then
              hrefi(i)=hrf
            else         !DEFault or unknown reference state
              hrefi(i)=hrefdf(i)
              ixflg=1    !must apply ref state to pures only
              if (hrf.ne.'DEF' .and. hrf.ne.'def') then
                ierr=-29
                herr='[SETREF warning -29] reference state not found, '
     &             //'using default reference state'//hnull
                call ERRMSG (ierr,herr)
              end if
            end if
c
c         write (*,*) ' SETREF--i,hrefi:  ',i,' ',hrefi(i)
            href(i)=0.0d0
            sref(i)=0.0d0
            if (hrefi(i).eq.'IIR') then
c  IIR reference state
c           write (*,*) ' SETREF--IIR ref state with ixflag =',ixflg
              tref(i)=273.15d0
              prefi(i)=-1.0d0   !code for reference state of sat liquid
              href(i)=200.0d0*wm(i) !reference enthalpy is 200.0 kJ/kg
              sref(i)=wm(i)         !reference entropy is 1.0 kJ/(kg.K)
            else if (hrefi(i).eq.'ASH') then
c  ASHRAE reference state
c           write (*,*) ' SETREF--ASH ref state with ixflag =',ixflg
              tref(i)=233.15d0
              prefi(i)=-1.0d0   !code for reference state of sat liquid
            else if (hrefi(i).eq.'NBP') then
c  normal boiling point reference state
c           write (*,*) ' SETREF--NBP ref state with ixflag =',ixflg
              tref(i)=tnbp(i)
              prefi(i)=-1.0d0   !code for ref state of saturated liquid
            else if (hrefi(i)(1:2).eq.'OT') then
c  'other' reference state--use input values of h0, s0, t0, p0
c           write (*,*) ' SETREF--OTHer ref state with ixflag =',ixflg
              if (hrf.eq.'DEF' .or. hrf.eq.'def') then
c  if input to SETREF is 'DEF' reference t,p,h,s come from fluid file
                tref(i)=tdef(i)
                prefi(i)=pdef(i)
                href(i)=hdef(i)
                sref(i)=sdef(i)
              else
c  otherwise, reference t,p,h,s are inputs to SETREF
                tref(i)=t0
                prefi(i)=p0
                href(i)=h0
                sref(i)=s0
              end if
            else if (hrefi(i).eq.'NUL') then
c  this case should be encountered only for i<0 when any model other
c  than the ECS model is specified; fill arrays with dummy values
              tref(i)=300.0
              prefi(i)=1.0d-6
            else if (hrefi(i).eq.'NA') then
c  this case is for bypassing the reference state in equation development
              tref(i)=300.0
              prefi(i)=1.0d-6
            else
c  unknown reference state (use normal boiling point)
              ierr=-29
              write (herr,1238) hrefi(i),i,hnull
 1238         format ('[SETREF warning -29] reference state (',a3,') ',
     &           'for component ',i3,' not found; will use ',
     &           'NBP reference state. ',a1)
              call ERRMSG (ierr,herr)
              tref(i)=tnbp(i)
              prefi(i)=-1.0d0   !code for ref state of saturated liquid
            end if
          end if
        enddo
c
c  check limits of equation of state against reference state temperature
c
        do 400 i=1,ncomp
        Dx=0.0d0
        px=0.0d0
        t=tref(i)
        call LIMITK ('EOS',i,t,Dx,px,tmin,tmax,Dmax,pmax,ierr2,herr2)
c
        if (hrefi(i).eq.'IIR') then
c  IIR reference state
          if (tmin.gt.273.15d0) then
c  minimum temperature of EOS is above IIR reference state temperature
c  issue error message
            ierr=22
            write (herr,1022) i,tmin,hrefdf(i),hnull
 1022       format ('[SETUP warning 22] ',
     &        'minimum temperature for component ',i3,' is greater ',
     &        'than IIR reference temperature; Tmin =',g11.5,
     &        ' K, Tref = 273.15 K.; will use default (',a3,'). ',a1)
          else
            if (tc(i).lt.274.15d0) then
c  critical temperature is too low to use IIR reference state
c  issue error message
              ierr=23
              write (herr,1023) i,tc(i),hrefdf(i),hnull
 1023         format ('[SETUP warning 23] ',
     &        'critical temperature for component ',i3,' is less ',
     &        'than IIR reference temperature; Tcrit =',g11.5,
     &        ' K, Tref = 273.15 K; will use default (',a3,'). ',a1)
            end if
          end if
c
        else if (hrefi(i).eq.'ASH') then
c  ASHRAE reference state
          if (tmin.gt.233.15d0) then
c  minimum temperature of EOS is above reference state temperature
            ierr=24
            write (herr,1024) i,tmin,hrefdf(i),hnull
 1024       format ('[SETUP warning 24] ',
     &        'minimum temperature for component ',i3,' is greater ',
     &        'than ASHRAE reference temperature; Tmin =',g11.5,
     &        ' K, Tref = 233.15 K; will use default (',a3,'). ',a1)
          else
            if (tc(i).lt.234.15d0) then
c  critical temperature is too low to use IIR reference state
              ierr=25
              write (herr,1025) i,tc(i),hrefdf(i),hnull
 1025         format ('[SETUP warning 25] ',
     &        'critical temperature for component ',i3,' is less ',
     &        'than ASHRAE reference temperature; Tcrit =',g11.5,
     &        ' K, Tref = 233.15 K; will use default (',a3,'). ',a1)
            end if
          end if
        else if (hrefi(i).eq.'NBP') then
c  Normal Boiling Point reference state
          if (tmin.gt.tnbp(i)) then
c  minimum temperature of EOS is above boiling point temperature
c  issue error message (this is a problem with CO2)
            ierr=26
            write (herr,1026) i,tmin,tnbp(i),hrefdf(i),hnull
 1026       format ('[SETUP warning 26] ',
     &        'minimum temperature for component ',i3,' is greater ',
     &        'than boiling point temperature; Tmin =',g11.5,
     &        ' K, Tnbp =',g11.5,' K; will use default (',a3,'). ',a1)
          end if
        end if
c
        if (ABS(ierr).ge.22) then
c  if specified reference state is outside limits for component i,
c  use default state (read in from fluid file); this code assumes
c  that the specification in the fluid file is valid
          call ERRMSG (ierr,herr)
          href(i)=0.0d0
          sref(i)=0.0d0
          if (hrefdf(i).eq.'IIR') then
c           write (*,1282) i,hrefi(i),hsvrfd(i)
c1282       format (2x,'SETREF--specified reference state for component'
c    &                ,i4,' (',a3,'), is not valid; will apply default '
c    &                ,'reference state (',a3,') to pures. ')
            tref(i)=273.15d0
            prefi(i)=-1.0d0
            href(i)=200.0d0*wm(i) !reference enthalpy is 200.0 kJ/kg
            sref(i)=wm(i)         !reference entropy is 1.0 kJ/(kg.K)
          else if (hrefdf(i).eq.'ASH') then
c           write (*,1282) i,hrefi(i),hsvrfd(i)
            tref(i)=233.15d0
            prefi(i)=-1.0d0
          else if (hrefdf(i).eq.'NBP') then
c           write (*,1282) i,hrefi(i),hsvrfd(i)
            tref(i)=tnbp(i)
            prefi(i)=-1.0d0
          else if (hrefdf(i)(1:2).eq.'OT') then
c           write (*,1282) i,hrefi(i),hsvrfd(i)
            tref(i)=tdef(i)
            prefi(i)=pdef(i)
            href(i)=hdef(i)
            sref(i)=sdef(i)
          else     !unknown reference state; set to reasonable values
c           write (*,*) ' SETREF--default reference state is unknown'
            tref(i)=0.8d0*tc(i)
            prefi(i)=-1.0d0
          end if
          hrefi(i)=hrefdf(i)
          ixflg=1                 !must apply ref state to pures only
        end if
c  save any error flags/messages to pass up the chain; then reset
        if (ierr.ne.0) then
          ierrsv=ierr
          herrsv=herr
        endif
        ierr=0
        herr=' '
 400    continue
c
c  apply reference enthalpy and entropy to each component; this for
c  case of ixflag =1; also used as initialization for ixflag = 2
c
        do 600 i=1,ncomp
c       write (*,1002) i,tref(i),href(i),sref(i)
c1002   format (1x,' SETREF--i,tref,href,sref: ',i3,f7.2,f10.1,f10.3)
        call PUREFLD(i)       !if mixture, calculate pure component i
        if (hrefi(i).eq.'NBP') then
c  find bubble point at one atmosphere for each pure component
          pnbp=101.325d0
          kbub=1
c         write (*,*) ' SETREF--i,x(i): ',i,x(i)
          call SATP (pnbp,x,kbub,tbub,rhol,rhov,xl,xv,ierr,herr)
          rhoref(i)=rhol
          tref(i)=tbub
        elseif (hrefi(i).eq.'NA') then
c  Don't do anything if deactivated
        else if (prefi(i).lt.0.0d0) then
c  find saturated liquid or vapor density (includes ASH and IIR states)
          if (prefi(i).lt.-1.99d0) then
            kbub=2
            call SATT (tref(i),x,kbub,p,rhol,rhov,xl,xv,ierr,herr)
            rhoref(i)=rhov
          else
            kbub=1
            call SATT (tref(i),x,kbub,p,rhol,rhov,xl,xv,ierr,herr)
            rhoref(i)=rhol
          end if
        else if (hrefi(i).eq.'OTH') then
c  find density at specified t,p (either input to SETREF or default value)
          t=tref(i)
          p=prefi(i)
          call TPFLSH (t,p,x,D,Dl,Dv,xl,xv,q,e,h,s,cv,cp,w,ierr,herr)
          rhoref(i)=D
        else if (hrefi(i).eq.'OT0') then
c  find density at specified t,p (either input to SETREF or default value)
          t=tref(i)
          p=prefi(i)
          call RMIX(x0)
          rhoref(i)=p/R/T
        end if
        if (ierr.ge.1) then
c  error encountered in calculating reference state--set reference
c  density to critical density?
          rhoref(i)=rhoc(i)
          ierr=119
          herr2=herr
          write (herr,1119) i,herr2(1:162),hnull
          call ERRMSG (ierr,herr)
 1119     format ('[SETREF error 119] convergence failure in ',
     &            'calculating reference state for component #',i3,
     &            ':  ',a162,a1)
        end if
        if (hrefi(i).ne.'NA') then
          tsav(i)=0.d0  !reset tsav to force recalculation of cp0
          if (hrefi(i).eq.'OT0') then
            call THERM0(tref(i),rhoref(i),x,ptherm,e,h,s,cv,cp,w,a,g)
          else
            call THERM (tref(i),rhoref(i),x,ptherm,e,h,s,cv,cp,w,hjt)
          endif
          if (ixflg.ne.2 .or. ncomp.eq.1 .or.
     &      hrefi(i)(1:2).ne.'OT') then
            sref(i)=s              !store reference entropy in common
            href(i)=h              !ditto for enthalpy
          endif
          tsav(i)=0.d0
        endif
c       write (*,1004) i,tref(i),href(i),sref(i),rhoref(i)
c1004   format (1x,' SETREF--after SATT, etc.: ',i3,f7.2,f10.1,2f10.3)
 600    continue
      end if     !end of initialization for pure components
c
      call PUREFLD(0)
      if (ixflg.eq.2 .and. ncomp.gt.1) then
c
c  reference state applied to mixture of specified composition
c  this option is not available if the 'DEFault' state is chosen
c
c       write (*,1005) hrf,(x0(i),i=1,ncomp)
c1005   format (1x,' SETREF--about to apply ref state "',a3,
c    &             '" to mixture w/ x = ',5f10.5)
        if (hrf(1:2).eq.'OT') then
          tmix=t0
          pref=p0
          smix=s0
          hmix=h0
        else
          hrdef='NUL'   !initialize pointer to default ref state
          if (hrf.eq.'DEF' .or. hrf.eq.'def') then
c  DEFault reference state--should get here only if all components have the
c  same default ref state stored in their respective fluid files
            hrdef=hrefdf(1)
          end if
          hmix=0.0d0
          smix=0.0d0
          if (hrf.eq.'IIR' .or.hrdef.eq.'IIR') then
c  IIR reference state
            tmix=273.15d0
            pref=-1.0d0       !code for ref state of saturated liquid
            smix=1.0d0        !reference entropy is 1.0 kJ/(kg.K)
            hmix=200.0d0      !reference enthalpy is 200. kJ/kg
          else if (hrf.eq.'ASH' .or.hrdef.eq.'ASH') then
c  ASHRAE reference state
            tmix=233.15d0
            pref=-1.0d0       !code for reference state of sat liquid
          else
c  boiling point or unknown ref state
c  find bubble point of specified mixture
            pnbp=101.325d0
            call SATP (pnbp,x0,1,tbub,Dmix,rhov,xl,xv,ierr,herr)
            tmix=tbub
            pref=-1.0d0       !code for reference state of sat liquid
          end if
c  check mix reference temperature against limits, if outside limits
c  take average of component Tref's (would be different only if one or
c  more components were outside limits)
          Dx=0.0d0     !initialize (check only temperature)
          px=0.0d0
          call LIMITX ('EOS',tmix,Dx,px,x0,tmin,tmax,Dmx,pmx,ie2,herr2)
          if (tmix.lt.tmin .or. tmix.gt.tmax) then
            tmix=0.0d0
            do i=1,ncomp
              tmix=tmix+x0(i)*tref(i)
            enddo
          end if
c  predefined reference states are on mass basis
          wmm=WMOL(x0)
          smix=smix*wmm
          hmix=hmix*wmm
        end if
c
        if (pref.lt.0.0d0) then
c  find liquid or vapor density (includes IIR and ASHrae states)
          if (pref.lt.-1.99d0) then
            kbub=2
            call SATT (tmix,x0,kbub,p,rhol,Dmix,xl,xv,ierr,herr)
          else
            kbub=1
            call SATT (tmix,x0,kbub,p,Dmix,rhov,xl,xv,ierr,herr)
          end if
        else
c  find density at specified t0, p0
          call LIMITX ('EOS',t0,Dx,p0,x0,tmin,tmax,Dmx,pmx,ierr,herr2)
          if (ierr.ge.1) then
            ierr=27
            herr='[SETREF error 27] specified reference state for '//
     &           'mixture is outside equation limits: '//herr2(1:170)//
     &           hnull
            call ERRMSG (ierr,herr)
            RETURN
          end if
          call TPFLSH (t0,p0,x0,D,Dl,Dv,xl,xv,q,e,h,s,cv,cp,w,ierr,herr)
          tmix=t0
          Dmix=D
        end if
c
        call THERM (tmix,Dmix,x0,ptherm,e,h,s,cv,cp,w,hjt)
        do i=1,ncomp
c  modify sref, href by amount needed for specified mixture
c  to agree with s0, h0
c         write (*,1006) s,smix,h,hmix
c1006     format (1x,' SETREF  s,smix,h,hmix:  ',2f10.4,2f10.1)
          sref(i)=sref(i)+s-smix
          href(i)=href(i)+h-hmix
c         write (*,1008) i,tref(i),href(i),sref(i),rhoref(i)
c1008 format (1x,' SETREF--after mix adjust: ',i3,f7.2,f10.1,2f10.3)
        enddo
      end if
c
c  variable ksetrf is flag to indicate that SETREF has been called
c  it is reset to 0 each time SETUP is called (indicating that SETREF
c  must be recalled)
      ksetrf=1
c
c  copy inputs to common block
c
      hsvrfs=hrf
      ixfsav=ixflg
      if (hrf(1:2).eq.'OT') then
c  save h0, s0, t0, p0 only if 'OTHer' reference state is used
        h0sav=h0
        s0sav=s0
        t0sav=t0
        p0sav=p0
      end if
      if (ixflg.eq.2) then
c  save x0 only if ixflg = 2
        do i=1,ncomp
          x0sav(i)=x0(i)
        enddo
      end if
c
c  set value of hrf if 'DEF' was input value
c
c     if (hrf.eq.'DEF') then
c       hrf=hrefi(1)
c       do i=1,ncomp
c         if (hrefi(i).ne.hrf) then
c           hrf='MIX'
c         end if
c       enddo
c     end if
c
c  reset error flag/message (the more significant errors are those which
c  arise from out-of-bounds conditions, not any potential errors in calls
c  to SATT, etc.)
      if (ierr.eq.0 .and. ierrsv.ne.0) then
        ierr=ierrsv
        herr=herrsv
      end if
c
      RETURN
c
      end                                             !subroutine SETREF
c
c ======================================================================
c
      subroutine SETMIX (hmxnme,hfmix,hrf,ncc,hfiles,x,ierr,herr)
c
c  open a mixture file (e.g., R410A.mix) and read constituents and
c  mole fractions
c
c  inputs:
c    hmxnme--mixture file name to be read in [character*255]
c    hfmix--mixture coefficients [character*255]
c           file name containing coefficients for mixture model
c    hrf--reference state for thermodynamic calculations [character*3]
c         (see info in subroutine setup for specifics)
c  outputs:
c     ncc--number of fluids in mixture
c     hfiles--array of file names specifying mixture components
c             that were used to call setup. [character*255 variable]
c     x--array of mole fractions for the specified mixture
c     ierr--error flag:  0 = successful
c                      101 = error in opening file
c                     -102 = mixture file contains mixing parameters
c                     -103 = composition not equal to one
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  03-27-06 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SETMIX
c     dll_export SETMIX
c
      parameter (ncmax=20)        !max number of components in mixture
      character*3 hrf
      character*255 hmxnme,hfmix,hfiles(ncmax),herr,herr2,MixName
      dimension x(ncmax)
c
      ierr2=0
      nread=12             !logical unit for file reads
      call OPENFL (nread,hmxnme,2,ierr,herr)
      if (ierr.ne.0) goto 998
      read (nread,2000,err=998,end=998) MixName
      read (nread,*,err=998,end=998) xMixMolarMass,xMixTc,xMixPc,xMixDc
      read (nread,*,err=998,end=998) ncc
      do i=1,ncc
        read (nread,2000,err=998,end=998) hfiles(i)
        call RFFILE(hfmix,hfiles(i))
      enddo
      sum=0
      do i=1,ncc
        read (nread,*,err=998,end=998) x(i)
        sum=sum+x(i)
      enddo
      if (ABS(sum-1.d0).gt.0.0000001.and.sum.gt.0.d0.and.ncc.gt.1) then
        ierr2=-103
        herr2='[SETMIX error -103] '//
     &       'The composition of the mixture does not sum '//
     &       'to one. The compositions have been normalized.'
        call ERRMSG (ierr2,herr2)
        do i=1,ncc
          x(i)=x(i)/sum
        enddo
      endif
      read (nread,*,err=998,end=998) i
      if (i.ne.0) then
        ierr2=-102
        herr2='[SETMIX error -102] '//
     &       'The mixture file contains updated binary interaction '//
     &       'parameters which have been ignored.'
        call ERRMSG (ierr2,herr2)
      endif
      close (nread)
      call SETUP (ncc,hfiles,hfmix,hrf,ierr,herr)
      if (ierr.eq.0 .and. ierr2.ne.0) then
        ierr=ierr2
        herr=herr2
      endif
      RETURN
 998  continue
      ierr=101
      herr='[SETMIX error 101] error in opening mixture file:  '//hmxnme
      call ERRMSG (ierr,herr)
      RETURN
c
 2000 format (a255)
      end                                             !subroutine SETMIX
c
c ======================================================================
c
      subroutine RFFILE (hfilei,hflref)
c
c  provide a full path specification for the '.fld' file for the ECS
c  reference fluid given a specification for a component '.fld' file
c  e.g. 'c:\REFPROP\fluids\r123.fld' and the reference fluid pointer
c  which is stored in that file, e.g. 'r134a.fld'
c
c  N.B.  This routine assumes that the reference fluid '.fld' file
c        is in the same directory as the component file(s)
c
c  inputs:
c   hfilei--file name specifying mixture component (character*255)
c   hflref--file name specifying reference fluid (character*255)
c  output:
c   hflref--file name specifying reference fluid (character*255)
c           the output value of hflref consists of the directory
c           information gleaned from hfilei concatenated with the
c           input value of hflref
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  12-14-95  MM, original version
c  12-18-95  MM, use htemp to concatenate character strings (cannot
c                concatenate string with itself on some compilers)
c  02-22-07 EWL, check for null character in input string
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      character*255 hfilei,hflref,htemp
c
      if (hfilei(1:5).eq.'BDATA' .or. hfilei(1:5).eq.'bdata') then
c  special case if component data is read from block data;
c  assume that reference fluid is in same directory as main program
c  i.e., hflref is simply fluid name + '.fld', e.g. 'r134a.fld'
c  --or--
c  the pointer to the reference fluid in block data is also to
c  block data, e.g. 'BDATA:811-97-2' for R134a
      else
c  usual case of fluid read from '.fld' file
        ncfile=index(hfilei,char(0))
        if (ncfile.eq.0) ncfile=LEN(hfilei)
        jpoint=0
        do i=1,ncfile
          j=1+ncfile-i
c  search for directory delimiters
c  (\,:,/  for DOS, Macintosh, and Unix, respectively)
          jchar=ICHAR(hfilei(j:j))
          if (jchar.eq.92 .or. jchar.eq.58 .or. jchar .eq.47) then
            jpoint=j
            goto 100
          end if
        enddo
 100    continue             !file delimiter character has been found
c       write (*,1004) ncfile,jpoint
c1004   format (1x,' RFFILE--ncfile, jpoint: ',2i4)
        if (jpoint.le.0) then
c  no path name provided in specification of component .fld file
c  assume that reference fluid is in same directory as main program
c  i.e., hflref is simply fluid name + '.fld', e.g. r134a.fld
        else
          htemp=hfilei(1:jpoint)//hflref
          hflref=htemp
        end if
      end if
c     write (*,1005) hflref
c1005 format (1x,' RFFILE--full path for reference fld file: ',a255)
c
      RETURN
      end                                             !subroutine RFFILE
c
c ======================================================================
c
      subroutine OPENFL (nread,hfl,idir,ierr,herr)
c
c  open a file for input
c
c  inputs:
c    nread--unit number
c      hfl--file name [character*255 variable]
c  outputs:
c     ierr--error flag:  0 = successful
c                      101 = error in opening file
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  08-10-06 EWL, original version (with code taken from old SETFLD)
c  11-10-06 EWL, add additional checks on the input file name
c  01-04-07 EWL, add even more checks
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      character*26 lower,upper
      character*255 herr,hfl,hfilei,hfilej,hpath
      character*1 bs
      common /FLDPTH/ hpath
      logical lread
c
      lower = 'abcdefghijklmnopqrstuvwxyz'
      upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      ierr=0
      herr=' '
c
c  check for ^M in filename caused by binary download to UNIX machines
      j=INDEX(hfl,CHAR(13))
      if (j.ne.0) hfl=hfl(1:j-1)//hfl(j+1:255)
c
 10   continue
      INQUIRE(unit=nread,opened=lread)  !Check that unit #nread is not
      if (lread) then                   !already being used.
        nread=nread+1
        goto 10
      endif
c
c  find last occurrence of a character other than a space
      kp=0
      if (hpath.ne.' ') then
        do j=1,255
          if (hpath(j:j).ne.' ') kp=j
        enddo
      endif
c
      do m=1,2
        bs=char(92)
        if (m.eq.2) bs='/'      !Try same options for a Unix system
c  add the file path if SETPATH called
        hfilei=hfl
        if (kp.ne.0) then
          if (hpath(kp:kp).ne.'/' .and. hpath(kp:kp).ne.char(92)) then
            hfilei=hpath(1:kp)//bs//hfl
          else
            hfilei=hpath(1:kp)//hfl
          endif
        endif
c  if your compile reports an error message on the following line,
c  remove the action='READ' piece and try it again.
c  There are two more similar occurrences below.
        open (unit=nread,file=hfilei,status='old',iostat=ioerr,
     &        action='READ')
        if (ioerr.eq.0) RETURN
c
c  try the file path with 'fluids\' added to it if SETPATH called
        if (kp.ne.0) then
          if (idir.ne.2) then
            hfilei=hpath(1:kp)//'fluids'//bs//hfl  !Look in hpath\fluids directory
            if (hpath(kp:kp).ne.'/' .and. hpath(kp:kp).ne.char(92))
     &          hfilei=hpath(1:kp)//bs//'fluids'//bs//hfl
          else
            hfilei=hpath(1:kp)//'mixtures'//bs//hfl
            if (hpath(kp:kp).ne.'/' .and. hpath(kp:kp).ne.char(92))
     &          hfilei=hpath(1:kp)//bs//'mixtures'//bs//hfl
          endif
          open (unit=nread,file=hfilei,status='old',iostat=ioerr,
     &          action='READ')
          if (ioerr.eq.0) RETURN
        endif
      enddo
c
c  Loop over possible paths in DOS/Windows and UNIX machines
c  and try upper/lower case scenarios.
c  01-13-05 WLJ at ORNL
      do m=1,4
        bs=char(92)
        if (m.ge.3) bs='/'      !Try same options for a Unix system
        do k=1,6
          if (idir.ne.2) then
            if (k.eq.1) hfilei=' '                             ! Current directory
            if (k.eq.2) hfilei='fluids'//bs                    ! Fluids subdirectory
            if (k.eq.3) hfilei=bs//'refprop'//bs//'fluids'//bs ! From \Refprop\Fluids\
            if (k.eq.4) hfilei='..'//bs//'fluids'//bs          ! From ..\Fluids\
            if (k.eq.5) hfilei=':fluids:'                      ! Mac
            if (k.eq.6) hfilei='c:'//bs//'Program Files'//bs// ! C:\Program files\Refprop\Fluids
     &                         'Refprop'//bs//'fluids'//bs
          else
            if (k.eq.1) hfilei=' '
            if (k.eq.2) hfilei='mixtures'//bs
            if (k.eq.3) hfilei=bs//'refprop'//bs//'mixtures'//bs
            if (k.eq.4) hfilei='..'//bs//'mixtures'//bs
            if (k.eq.5) hfilei=':mixtures:'
            if (k.eq.6) hfilei='c:'//bs//'Program Files'
     &                         //bs//'Refprop'//bs//'mixtures'//bs
          endif
c  find last occurrence of a character other than a space
          kl=1
          do j=1,255
            if (hfilei(j:j).ne.' ') kl=j+1
          enddo
          hfilei(kl:)=hfl
          hfilej=hfilei
c  find last occurrence of bs so that only the part after bs gets upper/lower cased when k=1
          km=1
          do j=1,255
            if (hfilei(j:j).eq.bs) km=j+1
          enddo
          if (k.eq.1 .and. km.ne.1) kl=km
c  checks for case sensitive machines
          do kk=1,5
c  add the path from SETPATH to the beginning
            if (m.eq.2 .or. m.eq.4) then
              if (kp.ne.0) then
                if(hpath(kp:kp).ne.'/'.and.hpath(kp:kp).ne.char(92)
     &         .and. hfilei(1:1).ne.'/'.and.hfilei(1:1).ne.char(92))then
                  hfilei=hpath(1:kp)//bs//hfilei
                else
                  hfilei=hpath(1:kp)//hfilei
                endif
              endif
            endif
            open (unit=nread,file=hfilei,status='old',iostat=ioerr,
     &            action='READ')
            if (ioerr.eq.0) RETURN
            hfilei=hfilej
c  try all upper case
            if (kk.eq.1) then
              do j=kl,255
                kc=INDEX(lower, hfilei(j:j))
                if (kc.ne.0) hfilei(j:j)=upper(kc:kc)
              enddo
c  try lower case '.fld':
            elseif (kk.eq.2) then
              j=INDEX(hfilei,'.FLD')
              if (j.ne.0) then
                hfilei=hfilej(1:j-1)//'.fld'//hfilej(j+4:255)
              endif
              j=INDEX(hfilei,'.BNC')
              if (j.ne.0) then
                hfilei=hfilej(1:j-1)//'.bnc'//hfilej(j+4:255)
              endif
              j=INDEX(hfilei,'.MIX')
              if (j.ne.0) then
                hfilei=hfilej(1:j-1)//'.mix'//hfilej(j+4:255)
              endif
c  try all lower case
            elseif (kk.eq.3) then
              do j=kl,255
                kc=INDEX(upper, hfilei(j:j))
                if (kc.ne.0) hfilei(j:j)=lower(kc:kc)
              enddo
c  try upper case '.FLD'
            elseif (kk.eq.4) then
              j=INDEX(hfilei,'.fld')
              if (j.ne.0) then
                hfilei=hfilej(1:j-1)//'.FLD'//hfilej(j+4:255)
              endif
              j=INDEX(hfilei,'.bnc')
              if (j.ne.0) then
                hfilei=hfilej(1:j-1)//'.BNC'//hfilej(j+4:255)
              endif
              j=INDEX(hfilei,'.mix')
              if (j.ne.0) then
                hfilei=hfilej(1:j-1)//'.MIX'//hfilej(j+4:255)
              endif
            end if
          enddo
c         write (*,1106) hfl(i),hfl
c1106     format (1x,'initial file specification failed:  ',a80/
c    &            1x,'reading coefficients from file:     ',a80)
        enddo
      enddo
      ierr=101
      RETURN
      end                                             !subroutine OPENFL
c
c ======================================================================
c
      subroutine SETPATH (hpth)
c
c  set the path where the fluid files are located
c
c  inputs:
c     hpth--location of the fluid files [character*255 variable]
c           The path does not need to contain the ending "\" and it can
c           point directly to the location where the DLL is stored if a
c           fluids subdirectory (with the corresponding fluid files) is
c           located there.
c           example:  hpth='C:\Program Files\Refprop'
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  09-22-06 EWL, original version
c  03-07-07 EWL, check for null character
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SETPATH
c     dll_export SETPATH
c
      character hpath*255,hpth*(*)
      common /FLDPTH/ hpath
      hpath=hpth
      i=index(hpath,char(0))
      if (i.ne.0) hpath(i:255)=' '
      RETURN
      end                                            !subroutine SETPATH
c
c ======================================================================
c
      subroutine PUREFLD (icomp)
c
c  Change the standard mixture setup so that the properties of one fluid
c  can be calculated as if SETUP had been called for a pure fluid.
c  Calling this routine will disable all mixture calculations.
c  To reset the mixture setup, call this routine with icomp=0.
c
c  inputs:
c   icomp--fluid number in a mixture to use as a pure fluid
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  06-13-06 EWL, original version
c  01-07-10 EWL, move most variables to subroutine RESETA
c  02-04-10 EWL, allow icomp to be set to values greater than nc, but not
c                greater than ncmax.  This is used for loading a mixture
c                and a ppf file, but changing nc to that of the mixture
c                so that the ppf is located at nc+1 (see subroutine SETNC)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
      parameter (ncmax=20)        !max number of components in mixture
c
cDEC$ ATTRIBUTES DLLEXPORT :: PUREFLD
c     dll_export PUREFLD
c
      common /NCOMP/ nc,ic
c
      call RESETA
      ic=0
      if (icomp.ge.0 .and. icomp.le.ncmax) ic=icomp
c
      RETURN
      end                                            !subroutine PUREFLD
c
c ======================================================================
c
      subroutine SETNC (ncomp)
c
c  Allow the user to modify the value of nc (the number of components in
c  a mixture) so that a subset of the full setup can be used.
c  For example, a 5 component mixture could be set up, but nc could
c  be set to 3 to calculate properties for just the first three components.
c  This also allows a mixture to be setup with one or more pure fluids
c  loaded after the components in the mixture.  For example, R407C could
c  be loaded by calling setup with R32.fld, R125.fld, R134a.fld, and R407C.ppf.
c  The number of components would be set to 3, and PUREFLD would be called
c  to access the properties from component 4, the pseudo-pure fluid equation.
c
c  inputs:
c   ncomp--number of components in the mixture
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  02-04-10 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SETNC
c     dll_export SETNC
c
      common /NCOMP/ nc,ic
c
      nc=ncomp
c
      RETURN
      end                                              !subroutine SETNC
c
c ======================================================================
c
      subroutine RESETA
c
c  Reset all variables that store saved properties so that new
c  new calculations are made regardless of previous inputs.
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  01-07-10 EWL, original version
c  02-04-10 EWL, change the ianc variable into an array
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (mxtrm=72)
c
      common /NCOMP/ nc,ic
      common /FSHSAV/ lsatt,lsatp
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /FXPNT/ d72l,d72v,thmax,hmax,htpl,htpv,
     &               stpl,stpv,tsmax,smax,tsmin,smin,tsminm,sminm
      common /FEQSAV/ phisav(n0:nx,mxtrm),delsav(n0:nx),tausav(n0:nx),
     &                taup(n0:nx,mxtrm),delp(n0:nx,mxtrm),
     &                delli(n0:nx,mxtrm),drvsav(n0:nx,16)
      common /TSTSAV/ tsavt,psavt,dlsavt,dvsavt,
     &                xsavt(ncmax),xlsavt(ncmax),xvsavt(ncmax),
     &                kphsvt,icsavt
      common /PSTSAV/ tsavp,psavp,dlsavp,dvsavp,
     &                xsavp(ncmax),xlsavp(ncmax),xvsavp(ncmax),
     &                kphsvp,icsavp
      common /WLMPS/ pstmin(n0:nx),pstmax(n0:nx)
      common /WLMPL/ pltmin(n0:nx),pltmax(n0:nx)
      common /WLMDL/ dltmin(n0:nx),dltmax(n0:nx)
      common /WLMDV/ dvtmin(n0:nx),dvtmax(n0:nx)
      common /ESTATE/ ieflg
c
c  Do not include variables such as ianc here which should not be reset
c  in general.  They belong in the call to SETUP.
c
      d72l=0.0d0
      d72v=0.0d0
      thmax=0.0d0
      hmax=0.0d0
      htpl=0.0d0
      htpv=0.0d0
      stpl=0.0d0
      stpv=0.0d0
      tsmax=0.0d0
      smax=0.0d0
      tsmin=0.0d0
      smin=0.0d0
      tsminm=0.0d0
      sminm=0.0d0
      lsatt=.false.
      lsatp=.false.
      ic=0
      ieflg=0
      do i=n0,nx
        delsav(i)=0.0d0
        tausav(i)=0.0d0
      enddo
c initialize saturated values
      tsavt=0.d0
      psavt=0.d0
      dlsavt=0.d0
      dvsavt=0.d0
      kphsvt=0
      icsavt=0
      tsavp=0.d0
      psavp=0.d0
      dlsavp=0.d0
      dvsavp=0.d0
      kphsvp=0
      icsavp=0
      do i=1,ncmax
        xsavt (i)=0.d0
        xlsavt(i)=0.d0
        xvsavt(i)=0.d0
        xsavp (i)=0.d0
        xlsavp(i)=0.d0
        xvsavp(i)=0.d0
        pstmin(i)=0.d0
        pltmin(i)=0.d0
        dltmin(i)=0.d0
        dvtmin(i)=0.d0
      enddo
c
      RETURN
      end                                             !subroutine RESETA
c
c ======================================================================
c
      block data BDSET
c
c  This block data initializes the common blocks which save input
c  arguments to SETUP and SETREF.  This is done to check if arguments
c  are unchanged since previous call; if so, routine can be bypassed.
c
c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
c  01-10-96  MM, original version
c  01-16-96  MM, add saves for SETREF
c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
c  03-20-96  MM, changes to accommodate SETMOD:  initialize /TRNMOD/
c                initialize heos, etc to 'NBS'
c  03-21-96  MM, replace /MODEL/ with /EOSMOD/, /STNMOD/, /REFST/
c  03-22-96  MM, initialize heta,htcx,hsten to 'NBS', rather than 'NUL'
c  03-27-96  MM, ditto for component models, hmxeos,heta,htcxk,hstenk
c  02-20-97  MM, initialize default ref states in /REFST/
c  03-28-97  MM, add /RESETM/ and set lreset=.true. so that models are
c                initialized by call to SETMOD rather than here
c  07-15-97  MM, move /FLAGS/ here from ftn_pas.f
c  11-25-97  MM, insert 2 extra 9's into the funny numbers in /FLAGS/
c  12-16-97  MM, set last open variable in /FLAGS/ to 'error'
c  12-01-98 EWL, add default for iprnterr
c  07-10-06 EWL, move any variables set to 'NUL' to beginning of SETUP
c                to comply with ANSI standards
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nx=ncmax)
      logical lreset,lrst,linit
      character*255 hpath
c
      common /NCOMP/ nc,ic
      common /REFSAV/ x0sav(nx),h0sav,s0sav,t0sav,p0sav
      common /IRFSAV/ ixfsav,ksetrf
      common /RESETM/ lreset,lrst,linit,ncset !flag indicating need to reset models
      common /FLDPTH/ hpath
c  flags indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,xerr,xnotd,xnotc
      common /prnterr/ iprnterr
      common /GERG2004/ iGERG04
c
      data nc,ic /2*0/
      data lreset /.true./ !reset all models on 1st call to SETUP/SETMOD
      data lrst /.false./
      data linit /.true./
      data x0sav /nx*0.0d0/
      data h0sav,s0sav,t0sav,p0sav /4*0.0d0/
      data ixfsav,ksetrf /2*0/
      data xnota /-9.99991d6/ !flag indicating 'not applicable'
      data x2ph  /-9.99992d6/ !flag indicating '2-phase'
      data xsubc /-9.99993d6/ !flag indicating 'subcooled'
      data xsuph /-9.99994d6/ !flag indicating 'superheated'
      data xsupc /-9.99995d6/ !flag indicating 'supercritical'
      data xinf  /-9.99996d6/ !flag indicating 'infinite'
      data xerr  /-9.99997d6/ !flag indicating 'error'
      data xnotd /-9.99998d6/ !flag indicating property 'not defined'
      data xnotc /-9.99999d6/ !flag indicating property 'not calculated'
      data iprnterr /0/
      data iGERG04 /0/
c
      data hpath /' '/
      end                                              !block data BDSET
c
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c ======================================================================
c                                                       end file setup.f
c ======================================================================
