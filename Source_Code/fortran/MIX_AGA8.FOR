c  begin file Mix_AGA8.f
c
c  This file contains the routines implementing the AGA8 equation of
c  state for hydrocarbon mixtures
c
c  contained here are:
c     function PHIAGA (itau,idel,tau,del,x)
c     subroutine SETAGA (ierr,herr)
c     subroutine UNSETAGA
c     subroutine SETAG (x)
c     block data AGA8CF
c
c ======================================================================
c ======================================================================
c
      function PHIAGA (itau,idel,tau,del,x)
c
c  compute reduced Helmholtz energy or a derivative as functions
c  of dimensionless temperature and density for the AGA8
c  equation of state
c
c  based on the DETAIL compressibility factor equation of:
c  Starling, K.E. and Savidge, J.L.
c  Compressibility Factors of Natural Gas and Other Related Hydrocarbon
c  Gases,
c  Transmission Measurement Committee Report No. 8, Catalog No. XQ9212,
c  American Gas Association, 1994.
c
c  inputs:
c     itau--flag specifying order of temperature derivative to calc
c     idel--flag specifying order of density derivative to calculate
c           when itau = 0 and idel = 0, compute A/RT
c           when itau = 0 and idel = 1, compute 1st density derivative
c           when itau = 1 and idel = 1, compute cross derivative
c           etc.
c      tau--dimensionless temperature (To/T)
c      del--dimensionless density (D/Do)
c  output (as function value):
c      phi--residual (real-gas) part of the AGA8 equation, or one
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
c  10-31-02 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      parameter (ncmax=20)        !max number of components in mixture
      double precision K
      dimension bs(18),cns(58)
      dimension an(58),un(58)
      dimension a0(64),x(ncmax)
c
      common /VAR/ K,bs,cns
      common /EOSCOEF/ an,un
c
      call SETAG(x)
      call REDX (x,t0,rho0)
      t=t0/tau
      rho=del*rho0
      D=K**3*rho
      d1=D
      d2=D*D
      d3=d2*D
      d4=d3*D
      d5=d4*D
      d6=d5*D
      d7=d6*D
      d8=d7*D
      d9=d8*D
      e0=exp(0.d0)
      e1=exp(-d1)
      e2=exp(-d2)
      e3=exp(-d3)
      e4=exp(-d4)
      t5    =sqrt(t)
      t15   =t5*t
      t10   =t
      t20   =t*t
      t30   =t20*t
      t40   =t30*t
      t50   =t40*t
      t60   =t50*t
      t70   =t60*t
      t80   =t70*t
      t90   =t80*t
      t35   =t30*t5
      t45   =t40*t5
      t75   =t70*t5
      t95   =t90*t5
      t110  =t90*t20
      t120  =t110*t
      t125  =t120*t5
      t210  =t120*t90
      t220  =t210*t
      t230  =t220*t
      tn5   =1.d0/t5
      tn130 =1.d0/t120/t
      tn10  =1.d0/t
      tn60  =1.d0/t60
c
      a0(1) =D* bs(1)
      a0(2) =D* bs(2)          /t5
      a0(3) =D* bs(3)          /t10
      a0(4) =D* bs(4)          /t35
      a0(5) =D* bs(5)          /tn5
      a0(6) =D* bs(6)          /t45
      a0(7) =D* bs(7)          /t5
      a0(8) =D* bs(8)          /t75
      a0(9) =D* bs(9)          /t95
      a0(10)=D* bs(10)         /t60
      a0(11)=D* bs(11)         /t120
      a0(12)=D* bs(12)         /t125
      a0(59)=D*(bs(13)-cns(13))/tn60
      a0(60)=D*(bs(14)-cns(14))/t20
      a0(61)=D*(bs(15)-cns(15))/t30
      a0(62)=D*(bs(16)-cns(16))/t20
      a0(63)=D*(bs(17)-cns(17))/t20
      a0(64)=D*(bs(18)-cns(18))/t110
c
      a0(13)=cns(13)/tn60 *d1*e3
      a0(14)=cns(14)/t20  *d1*e2
      a0(15)=cns(15)/t30  *d1*e2
      a0(16)=cns(16)/t20  *d1*e2
      a0(17)=cns(17)/t20  *d1*e4
      a0(18)=cns(18)/t110 *d1*e4
      a0(19)=cns(19)/tn5  *d2*e0
      a0(20)=cns(20)/t5   *d2*e0
      a0(21)=cns(21)      *d2*e2
      a0(22)=cns(22)/t40  *d2*e2
      a0(23)=cns(23)/t60  *d2*e2
      a0(24)=cns(24)/t210 *d2*e4
      a0(25)=cns(25)/t230 *d2*e4
      a0(26)=cns(26)/t220 *d2*e4
      a0(27)=cns(27)/tn10 *d2*e4
      a0(28)=cns(28)/tn5  *d3*e0
      a0(29)=cns(29)/t70  *d3*e1
      a0(30)=cns(30)/tn10 *d3*e1
      a0(31)=cns(31)/t60  *d3*e2
      a0(32)=cns(32)/t40  *d3*e2
      a0(33)=cns(33)/t10  *d3*e3
      a0(34)=cns(34)/t90  *d3*e3
      a0(35)=cns(35)/tn130*d3*e4
      a0(36)=cns(36)/t210 *d3*e4
      a0(37)=cns(37)/t80  *d3*e4
      a0(38)=cns(38)/tn5  *d4*e0
      a0(39)=cns(39)      *d4*e0
      a0(40)=cns(40)/t20  *d4*e2
      a0(41)=cns(41)/t70  *d4*e2
      a0(42)=cns(42)/t90  *d4*e2
      a0(43)=cns(43)/t220 *d4*e4
      a0(44)=cns(44)/t230 *d4*e4
      a0(45)=cns(45)/t10  *d5*e0
      a0(46)=cns(46)/t90  *d5*e2
      a0(47)=cns(47)/t30  *d5*e2
      a0(48)=cns(48)/t80  *d5*e4
      a0(49)=cns(49)/t230 *d5*e4
      a0(50)=cns(50)/t15  *d6*e0
      a0(51)=cns(51)/t50  *d6*e2
      a0(52)=cns(52)/tn5  *d7*e0
      a0(53)=cns(53)/t40  *d7*e2
      a0(54)=cns(54)/t70  *d8*e1
      a0(55)=cns(55)/t30  *d8*e2
      a0(56)=cns(56)      *d8*e2
      a0(57)=cns(57)/t10  *d9*e2
      a0(58)=cns(58)      *d9*e2
c
      if (idel.eq.1) then
        a0(19)=a0(19)*2d0
        a0(20)=a0(20)*2d0
        a0(28)=a0(28)*3d0
        a0(38)=a0(38)*4d0
        a0(39)=a0(39)*4d0
        a0(45)=a0(45)*5d0
        a0(50)=a0(50)*6d0
        a0(52)=a0(52)*7d0
        a0(13)=a0(13)*(1d0-3d0*d3)
        a0(14)=a0(14)*(1d0-2d0*d2)
        a0(15)=a0(15)*(1d0-2d0*d2)
        a0(16)=a0(16)*(1d0-2d0*d2)
        a0(17)=a0(17)*(1d0-4d0*d4)
        a0(18)=a0(18)*(1d0-4d0*d4)
        a0(21)=a0(21)*(2d0-2d0*d2)
        a0(22)=a0(22)*(2d0-2d0*d2)
        a0(23)=a0(23)*(2d0-2d0*d2)
        a0(24)=a0(24)*(2d0-4d0*d4)
        a0(25)=a0(25)*(2d0-4d0*d4)
        a0(26)=a0(26)*(2d0-4d0*d4)
        a0(27)=a0(27)*(2d0-4d0*d4)
        a0(29)=a0(29)*(3d0-    d1)
        a0(30)=a0(30)*(3d0-    d1)
        a0(31)=a0(31)*(3d0-2d0*d2)
        a0(32)=a0(32)*(3d0-2d0*d2)
        a0(33)=a0(33)*(3d0-3d0*d3)
        a0(34)=a0(34)*(3d0-3d0*d3)
        a0(35)=a0(35)*(3d0-4d0*d4)
        a0(36)=a0(36)*(3d0-4d0*d4)
        a0(37)=a0(37)*(3d0-4d0*d4)
        a0(40)=a0(40)*(4d0-2d0*d2)
        a0(41)=a0(41)*(4d0-2d0*d2)
        a0(42)=a0(42)*(4d0-2d0*d2)
        a0(43)=a0(43)*(4d0-4d0*d4)
        a0(44)=a0(44)*(4d0-4d0*d4)
        a0(46)=a0(46)*(5d0-2d0*d2)
        a0(47)=a0(47)*(5d0-2d0*d2)
        a0(48)=a0(48)*(5d0-4d0*d4)
        a0(49)=a0(49)*(5d0-4d0*d4)
        a0(51)=a0(51)*(6d0-2d0*d2)
        a0(53)=a0(53)*(7d0-2d0*d2)
        a0(54)=a0(54)*(8d0-    d1)
        a0(55)=a0(55)*(8d0-2d0*d2)
        a0(56)=a0(56)*(8d0-2d0*d2)
        a0(57)=a0(57)*(9d0-2d0*d2)
        a0(58)=a0(58)*(9d0-2d0*d2)
c
      elseif (idel.eq.2) then
        do n=1,12
          a0(n)=0
        enddo
        do n=59,64
          a0(n)=0
        enddo
        a0(19)=a0(19)*2d0
        a0(20)=a0(20)*2d0
        a0(28)=a0(28)*6d0
        a0(38)=a0(38)*12d0
        a0(39)=a0(39)*12d0
        a0(45)=a0(45)*20d0
        a0(50)=a0(50)*30d0
        a0(52)=a0(52)*42d0
        a0(13)=a0(13)*((1d0-3d0*d3)*(   -3d0*d3)- 9d0*d3)
        a0(14)=a0(14)*((1d0-2d0*d2)*(   -2d0*d2)- 4d0*d2)
        a0(15)=a0(15)*((1d0-2d0*d2)*(   -2d0*d2)- 4d0*d2)
        a0(16)=a0(16)*((1d0-2d0*d2)*(   -2d0*d2)- 4d0*d2)
        a0(17)=a0(17)*((1d0-4d0*d4)*(   -4d0*d4)-16d0*d4)
        a0(18)=a0(18)*((1d0-4d0*d4)*(   -4d0*d4)-16d0*d4)
        a0(21)=a0(21)*((2d0-2d0*d2)*(1d0-2d0*d2)- 4d0*d2)
        a0(22)=a0(22)*((2d0-2d0*d2)*(1d0-2d0*d2)- 4d0*d2)
        a0(23)=a0(23)*((2d0-2d0*d2)*(1d0-2d0*d2)- 4d0*d2)
        a0(24)=a0(24)*((2d0-4d0*d4)*(1d0-4d0*d4)-16d0*d4)
        a0(25)=a0(25)*((2d0-4d0*d4)*(1d0-4d0*d4)-16d0*d4)
        a0(26)=a0(26)*((2d0-4d0*d4)*(1d0-4d0*d4)-16d0*d4)
        a0(27)=a0(27)*((2d0-4d0*d4)*(1d0-4d0*d4)-16d0*d4)
        a0(29)=a0(29)*((3d0-    d1)*(2d0-    d1)-     d1)
        a0(30)=a0(30)*((3d0-    d1)*(2d0-    d1)-     d1)
        a0(31)=a0(31)*((3d0-2d0*d2)*(2d0-2d0*d2)- 4d0*d2)
        a0(32)=a0(32)*((3d0-2d0*d2)*(2d0-2d0*d2)- 4d0*d2)
        a0(33)=a0(33)*((3d0-3d0*d3)*(2d0-3d0*d3)- 9d0*d3)
        a0(34)=a0(34)*((3d0-3d0*d3)*(2d0-3d0*d3)- 9d0*d3)
        a0(35)=a0(35)*((3d0-4d0*d4)*(2d0-4d0*d4)-16d0*d4)
        a0(36)=a0(36)*((3d0-4d0*d4)*(2d0-4d0*d4)-16d0*d4)
        a0(37)=a0(37)*((3d0-4d0*d4)*(2d0-4d0*d4)-16d0*d4)
        a0(40)=a0(40)*((4d0-2d0*d2)*(3d0-2d0*d2)- 4d0*d2)
        a0(41)=a0(41)*((4d0-2d0*d2)*(3d0-2d0*d2)- 4d0*d2)
        a0(42)=a0(42)*((4d0-2d0*d2)*(3d0-2d0*d2)- 4d0*d2)
        a0(43)=a0(43)*((4d0-4d0*d4)*(3d0-4d0*d4)-16d0*d4)
        a0(44)=a0(44)*((4d0-4d0*d4)*(3d0-4d0*d4)-16d0*d4)
        a0(46)=a0(46)*((5d0-2d0*d2)*(4d0-2d0*d2)- 4d0*d2)
        a0(47)=a0(47)*((5d0-2d0*d2)*(4d0-2d0*d2)- 4d0*d2)
        a0(48)=a0(48)*((5d0-4d0*d4)*(4d0-4d0*d4)-16d0*d4)
        a0(49)=a0(49)*((5d0-4d0*d4)*(4d0-4d0*d4)-16d0*d4)
        a0(51)=a0(51)*((6d0-2d0*d2)*(5d0-2d0*d2)- 4d0*d2)
        a0(53)=a0(53)*((7d0-2d0*d2)*(6d0-2d0*d2)- 4d0*d2)
        a0(54)=a0(54)*((8d0-    d1)*(7d0-    d1)-     d1)
        a0(55)=a0(55)*((8d0-2d0*d2)*(7d0-2d0*d2)- 4d0*d2)
        a0(56)=a0(56)*((8d0-2d0*d2)*(7d0-2d0*d2)- 4d0*d2)
        a0(57)=a0(57)*((9d0-2d0*d2)*(8d0-2d0*d2)- 4d0*d2)
        a0(58)=a0(58)*((9d0-2d0*d2)*(8d0-2d0*d2)- 4d0*d2)
c
      elseif (idel.eq.3) then
        do n=1,12
          a0(n)=0
        enddo
        do n=59,64
          a0(n)=0
        enddo
        a0(19)=0.d0
        a0(20)=0.d0
        a0(28)=a0(28)*6d0
        a0(38)=a0(38)*24d0
        a0(39)=a0(39)*24d0
        a0(45)=a0(45)*60d0
        a0(50)=a0(50)*120d0
        a0(52)=a0(52)*210d0
        d12=d6*d6
        a0(13)=a0(13)*(  0d0 - d3*24d0  + d6*81d0  - 27d0*d9)
        a0(14)=a0(14)*(  0d0 - d2*6d0   + d4*24d0  - 8d0 *d6)
        a0(15)=a0(15)*(  0d0 - d2*6d0   + d4*24d0  - 8d0 *d6)
        a0(16)=a0(16)*(  0d0 - d2*6d0   + d4*24d0  - 8d0 *d6)
        a0(17)=a0(17)*(  0d0 - d4*60d0  + d8*192d0 - 64d0*d12)
        a0(18)=a0(18)*(  0d0 - d4*60d0  + d8*192d0 - 64d0*d12)
        a0(21)=a0(21)*(  0d0 - d2*24d0  + d4*36d0  - 8d0 *d6)
        a0(22)=a0(22)*(  0d0 - d2*24d0  + d4*36d0  - 8d0 *d6)
        a0(23)=a0(23)*(  0d0 - d2*24d0  + d4*36d0  - 8d0 *d6)
        a0(24)=a0(24)*(  0d0 - d4*120d0 + d8*240d0 - 64d0*d12)
        a0(25)=a0(25)*(  0d0 - d4*120d0 + d8*240d0 - 64d0*d12)
        a0(26)=a0(26)*(  0d0 - d4*120d0 + d8*240d0 - 64d0*d12)
        a0(27)=a0(27)*(  0d0 - d4*120d0 + d8*240d0 - 64d0*d12)
        a0(29)=a0(29)*(  6d0 - d1*18d0  + d2*9d0   - 1d0 *d3)
        a0(30)=a0(30)*(  6d0 - d1*18d0  + d2*9d0   - 1d0 *d3)
        a0(31)=a0(31)*(  6d0 - d2*54d0  + d4*48d0  - 8d0 *d6)
        a0(32)=a0(32)*(  6d0 - d2*54d0  + d4*48d0  - 8d0 *d6)
        a0(33)=a0(33)*(  6d0 - d3*114d0 + d6*135d0 - 27d0*d9)
        a0(34)=a0(34)*(  6d0 - d3*114d0 + d6*135d0 - 27d0*d9)
        a0(35)=a0(35)*(  6d0 - d4*204d0 + d8*288d0 - 64d0*d12)
        a0(36)=a0(36)*(  6d0 - d4*204d0 + d8*288d0 - 64d0*d12)
        a0(37)=a0(37)*(  6d0 - d4*204d0 + d8*288d0 - 64d0*d12)
        a0(40)=a0(40)*( 24d0 - d2*96d0  + d4*60d0  - 8d0 *d6)
        a0(41)=a0(41)*( 24d0 - d2*96d0  + d4*60d0  - 8d0 *d6)
        a0(42)=a0(42)*( 24d0 - d2*96d0  + d4*60d0  - 8d0 *d6)
        a0(43)=a0(43)*( 24d0 - d4*312d0 + d8*336d0 - 64d0*d12)
        a0(44)=a0(44)*( 24d0 - d4*312d0 + d8*336d0 - 64d0*d12)
        a0(46)=a0(46)*( 60d0 - d2*150d0 + d4*72d0  - 8d0 *d6)
        a0(47)=a0(47)*( 60d0 - d2*150d0 + d4*72d0  - 8d0 *d6)
        a0(48)=a0(48)*( 60d0 - d4*444d0 + d8*384d0 - 64d0*d12)
        a0(49)=a0(49)*( 60d0 - d4*444d0 + d8*384d0 - 64d0*d12)
        a0(51)=a0(51)*(120d0 - d2*216d0 + d4*84d0  - 8d0 *d6)
        a0(53)=a0(53)*(210d0 - d2*294d0 + d4*96d0  - 8d0 *d6)
        a0(54)=a0(54)*(336d0 - d1*168d0 + d2*24d0  - 1d0 *d3)
        a0(55)=a0(55)*(336d0 - d2*384d0 + d4*108d0 - 8d0 *d6)
        a0(56)=a0(56)*(336d0 - d2*384d0 + d4*108d0 - 8d0 *d6)
        a0(57)=a0(57)*(504d0 - d2*486d0 + d4*120d0 - 8d0 *d6)
        a0(58)=a0(58)*(504d0 - d2*486d0 + d4*120d0 - 8d0 *d6)
      endif
c
      if (itau.eq.1) then
        do n=1,58
          a0(n)=a0(n)*un(n)
        enddo
        do n=59,64
          a0(n)=a0(n)*un(n-46)
        enddo
      elseif (itau.eq.2) then
        do n=1,58
          a0(n)=a0(n)*un(n)*(un(n)-1.d0)
        enddo
        do n=59,64
          a0(n)=a0(n)*un(n-46)*(un(n-46)-1.d0)
        enddo
      endif
c
      ar=0.d0
      do n=1,64
        ar=ar+a0(n)
      enddo
      phiaga=ar
c
      RETURN
      end                                               !function PHIAGA
c
c ======================================================================
c
      subroutine SETAGA (ierr,herr)
c
c  set up working arrays for use with AGA8 equation of state
c
c  input:
c  outputs:
c     ierr--error flag:  0 = successful
c                        1 = error (e.g. fluid not found)
c     herr--error string (character*255 variable if ierr<>0)
c     [fluid parameters, etc. returned via various common blocks]
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  10-31-02 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: SETAGA
c     dll_export SETAGA
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (ncppmx=20)       !max number of Cp0 terms
      character*3 hpheq,heos,hmxeos,hmodcp,hagasv
      character*1 htab,hnull
      character*12 hcas
      character*255 herr
      dimension an(58),un(58)
      dimension qb(21),fb(21),sb(21),wb(21)
      dimension ifp(ncmax)
      double precision eijs(21,21),uij(21,21),kij(21,21),gijs(21,21)
      double precision mrb(21),eb(21),kb(21),gb(21)
      double precision kb2(21),eb2(21),kij2(21,21),uij2(21,21),
     &       gij2(21,21),bs2(18,21,21)
      double precision  acp0(21),bcp0(21),ccp0(21),dcp0(21),ecp0(21),
     &         fcp0(21),gcp0(21),hcp0(21),icp0(21),jcp0(21),kcp0(21)
c
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCAS/ hcas(n0:nx)
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /EOSCOEF/ an,un
      common /FLDCOEF/ mrb,eb,kb,gb,qb,fb,sb,wb,eijs,uij,kij,gijs
      common /INTMCOEF/ kb2,eb2,kij2,uij2,gij2,bs2,ifp
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /CP0COEF/ acp0,bcp0,ccp0,dcp0,ecp0,fcp0,
     &                 gcp0,hcp0,icp0,jcp0,kcp0
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),nCOSH(n0:nx),
     &                nSINH(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WRDCPP/ tred(n0:nx),Cred(n0:nx)
      common /WLMCPP/ tmin(n0:nx),tmax(n0:nx),pmax(n0:nx),rhomax(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
c
      common /AGASV1/ hagasv
      common /AGASV2/ ntermcsv(ncmax),ntermesv(ncmax),nCOSHsv(ncmax),
     &                nSINHsv(ncmax)
      common /AGASV3/ Reossv(ncmax),wmsv(ncmax),tminsv(ncmax),
     &                tmaxsv(ncmax),tredsv(ncmax),Credsv(ncmax),
     &                trefsv(ncmax),rhorefsv(ncmax),hrefsv(ncmax),
     &                srefsv(ncmax),xksv(ncmax,5),xthsv(ncmax,5),
     &                xhsv(ncmax,5),cpcsv(ncmax,5),cphsv(ncmax,5),Rsv
c
      ierr=0
      herr=' '
c
      do i=1,nc
        ifp(i)=0
        if (hcas(i).eq.'74-82-8')    ifp(i)=1   !Methane
        if (hcas(i).eq.'7727-37-9')  ifp(i)=2   !Nitrogen
        if (hcas(i).eq.'124-38-9')   ifp(i)=3   !Carbon Dioxide
        if (hcas(i).eq.'74-84-0')    ifp(i)=4   !Ethane
        if (hcas(i).eq.'74-98-6')    ifp(i)=5   !Propane
        if (hcas(i).eq.'7732-18-5')  ifp(i)=6   !Water
        if (hcas(i).eq.'7783-06-4')  ifp(i)=7   !Hydrogen Sulfide
        if (hcas(i).eq.'1333-74-0')  ifp(i)=8   !Hydrogen
        if (hcas(i).eq.'1333-74-0p') ifp(i)=8   !Hydrogen (para)
        if (hcas(i).eq.'630-08-0')   ifp(i)=9   !Carbon Monoxide
        if (hcas(i).eq.'7782-44-7')  ifp(i)=10  !Oxygen
        if (hcas(i).eq.'75-28-5')    ifp(i)=11  !Isobutane
        if (hcas(i).eq.'106-97-8')   ifp(i)=12  !Butane
        if (hcas(i).eq.'78-78-4')    ifp(i)=13  !Isopentane
        if (hcas(i).eq.'109-66-0')   ifp(i)=14  !Pentane
        if (hcas(i).eq.'110-54-3')   ifp(i)=15  !Hexane
        if (hcas(i).eq.'142-82-5')   ifp(i)=16  !Heptane
        if (hcas(i).eq.'111-65-9')   ifp(i)=17  !Octane
        if (hcas(i).eq.'111-84-2')   ifp(i)=18  !Nonane
        if (hcas(i).eq.'124-18-5')   ifp(i)=19  !Decane
        if (hcas(i).eq.'7440-59-7')  ifp(i)=20  !Helium
        if (hcas(i).eq.'7440-37-1')  ifp(i)=21  !Argon
        if (hcas(i).eq.'463-82-1')   ifp(i)=13  !Set neopentane as isopentane
        if (hcas(i).eq.'108-88-3')   ifp(i)=16  !Set toluene as heptane
        if (hcas(i).eq.'71-43-2')    ifp(i)=15  !Set benzene as hexane
        if (hcas(i).eq.'74-85-1')    ifp(i)=4   !Set ethylene as ethane
        if (hcas(i).eq.'115-07-1')   ifp(i)=5   !Set propylene as propane
        if (hcas(i).eq.'106-98-9')   ifp(i)=12  !Set butene as butane
        if (ifp(i).eq.0) then
          ierr=1
          herr='[SETAGA error 1] Not all requested fluids are '//
     &         'available in AGA8, '//hcas(i)//hnull
          call ERRMSG (ierr,herr)
          RETURN
        endif
      enddo
c
      if (heos.ne.'AGA') then
        hagasv=heos
        Rsv=R
        do i=1,nc
          Reossv(i)    = Reos(i)
          wmsv(i)      = wm(i)
          tminsv(i)    = tmin(i)
          tmaxsv(i)    = tmax(i)
          tredsv(i)    = tred(i)
          Credsv(i)    = Cred(i)
          ntermcsv(i)  = ntermc(i)
          ntermesv(i)  = nterme(i)
          nCOSHsv(i)   = nCOSH(i)
          nSINHsv(i)   = nSINH(i)
          trefsv(i)    = tref(i)
          rhorefsv(i)  = rhoref(i)
          hrefsv(i)    = href(i)
          srefsv(i)    = sref(i)
          do j=1,5
            xksv(i,j)  = xk(i,j)
            xthsv(i,j) = xth(i,j)
            xhsv(i,j)  = xh(i,j)
            cpcsv(i,j) = cpc(i,j)
            cphsv(i,j) = cph(i,j)
          enddo
        enddo
      endif
c
      heos='AGA'
      R=8.31451d0
      do i=1,nc
        Reos(i)=8.31451d0
      enddo
c
      do i=1,21
        kb2(i)=kb(i)**2.5d0
        eb2(i)=eb(i)**2.5d0
      enddo
      do i=1,20
        do j=i+1,21
          kij2(i,j)=2.d0*(kij(i,j)**5-1.d0)*kb2(i)*kb2(j)
          uij2(i,j)=2.d0*(uij(i,j)**5-1.d0)*eb2(i)*eb2(j)
          gij2(i,j)=(gijs(i,j)-1.d0)*(gb(i)+gb(j))
          kij2(j,i)=kij2(i,j)
          uij2(j,i)=uij2(i,j)
          gij2(j,i)=gij2(i,j)
        enddo
      enddo
      do n=1,18
        do i=1,21
        do j=i,21
          eij=eijs(i,j)*sqrt(eb(i)*eb(j))
          bb=1.d0
          if (n.eq.5 .or. n.eq.6)    bb=gijs(i,j)*(gb(i)+gb(j))/2.d0
          if (n.eq.7 .or. n.eq.16)   bb=qb(i)*qb(j)
          if (n.eq.8 .or. n.eq.9)    bb=sb(i)*sb(j)
          if (n.ge.10 .and. n.le.12) bb=wb(i)*wb(j)
          if (n.eq.13)               bb=fb(i)*fb(j)
          if (i.ne.j) bb=2.d0*bb
          bs2(n,i,j)=bb*eij**un(n)*(kb(i)*kb(j))**1.5d0
          bs2(n,j,i)=bs2(n,i,j)
        enddo
        enddo
      enddo
      do i=1,nc
        wm(i)=mrb(ifp(i))
        tmin(i)=100.d0
        tmax(i)=1000.d0
        tred(i)=1.d0
        Cred(i)=4.184d0
        ntermc(i)=1
        nterme(i)=0
        nCOSH(i)=2
        nSINH(i)=2
        cpc(i,1)=bcp0(ifp(i))
        xk(i,1)=0.d0
        do j=2,5
          xk(i,j)=-2.d0
          xth(i,j)=-1.d0
          xh(i,j)=-2.d0
        enddo
        j=ifp(i)
        cpc(i,2)=ecp0(j)*fcp0(j)**2
        cph(i,2)=fcp0(j)
        cpc(i,3)=icp0(j)*jcp0(j)**2
        cph(i,3)=jcp0(j)
        cpc(i,4)=ccp0(j)*dcp0(j)**2
        cph(i,4)=dcp0(j)
        cpc(i,5)=gcp0(j)*hcp0(j)**2
        cph(i,5)=hcp0(j)
        tref(i)=298.15d0
        rhoref(i)=101.325d0/r/tref(i)
        t=tref(i)
        h=acp0(j)+bcp0(j)*t
     &   +ccp0(j)*dcp0(j)/TANH(dcp0(j)/t)
     &   -ecp0(j)*fcp0(j)*TANH(fcp0(j)/t)
     &   +gcp0(j)*hcp0(j)/TANH(hcp0(j)/t)
     &   -icp0(j)*jcp0(j)*TANH(jcp0(j)/t)
        s=kcp0(j)+bcp0(j)*LOG(t)
     &   +ccp0(j)*(dcp0(j)/t/TANH(dcp0(j)/t)-LOG(SINH(dcp0(j)/t)))
     &   -ecp0(j)*(fcp0(j)/t*TANH(fcp0(j)/t)-LOG(COSH(fcp0(j)/t)))
     &   +gcp0(j)*(hcp0(j)/t/TANH(hcp0(j)/t)-LOG(SINH(hcp0(j)/t)))
     &   -icp0(j)*(jcp0(j)/t*TANH(jcp0(j)/t)-LOG(COSH(jcp0(j)/t)))
        href(i)=-h*4.184d0
        sref(i)=-s*4.184d0
      enddo
c
      RETURN
      end                                             !subroutine SETAGA
c
c ======================================================================
c
      subroutine UNSETAGA
c
c  Load original values into arrays changed in the call to SETAGA.  This
c  routine resets the values back to those loaded when SETUP was called.
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  01-07-10 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: UNSETAGA
c     dll_export UNSETAGA
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport ECS
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      parameter (ncppmx=20)       !max number of Cp0 terms
      character*3 hpheq,heos,hmxeos,hmodcp,hagasv
c
      common /NCOMP/ nc,ic
      common /EOSMOD/ hpheq,heos,hmxeos(n0:nx),hmodcp(n0:nx)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /CCON/ tc(n0:nx),pc(n0:nx),rhoc(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /WNTCPP/ ntermc(n0:nx),nterme(n0:nx),nCOSH(n0:nx),
     &                nSINH(n0:nx),nsp1(n0:nx),nsp2(n0:nx),nsp3(n0:nx)
      common /WRDCPP/ tred(n0:nx),Cred(n0:nx)
      common /WLMCPP/ tmin(n0:nx),tmax(n0:nx),pmax(n0:nx),rhomax(n0:nx)
      common /WCPCPP/ cpc(n0:nx,ncppmx),xk(n0:nx,ncppmx),
     &                cph(n0:nx,ncppmx),xth(n0:nx,ncppmx),
     &                                  xh(n0:nx,ncppmx)
      common /CREF/ tref(n0:nx),rhoref(n0:nx),href(n0:nx),sref(n0:nx)
c
      common /AGASV1/ hagasv
      common /AGASV2/ ntermcsv(ncmax),ntermesv(ncmax),nCOSHsv(ncmax),
     &                nSINHsv(ncmax)
      common /AGASV3/ Reossv(ncmax),wmsv(ncmax),tminsv(ncmax),
     &                tmaxsv(ncmax),tredsv(ncmax),Credsv(ncmax),
     &                trefsv(ncmax),rhorefsv(ncmax),hrefsv(ncmax),
     &                srefsv(ncmax),xksv(ncmax,5),xthsv(ncmax,5),
     &                xhsv(ncmax,5),cpcsv(ncmax,5),cphsv(ncmax,5),Rsv
c
      call RESETA
      if (heos.eq.'AGA') then
        heos=hagasv
        R=Rsv
        do i=1,nc
          Reos(i)   = Reossv(i)
          wm(i)     = wmsv(i)
          tmin(i)   = tminsv(i)
          tmax(i)   = tmaxsv(i)
          tred(i)   = tredsv(i)
          Cred(i)   = Credsv(i)
          ntermc(i) = ntermcsv(i)
          nterme(i) = ntermesv(i)
          nCOSH(i)  = nCOSHsv(i)
          nSINH(i)  = nSINHsv(i)
          tref(i)   = trefsv(i)
          rhoref(i) = rhorefsv(i)
          href(i)   = hrefsv(i)
          sref(i)   = srefsv(i)
          do j=1,5
            xk(i,j) = xksv(i,j)
            xth(i,j)= xthsv(i,j)
            xh(i,j) = xhsv(i,j)
            cpc(i,j)= cpcsv(i,j)
            cph(i,j)= cphsv(i,j)
          enddo
        enddo
      endif
c
      RETURN
      end                                           !subroutine UNSETAGA
c
c ======================================================================
c
      subroutine SETAG (x)
c
c  set up working arrays for use with the AGA8 equation of state
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  10-31-02 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      parameter (ncmax=20)        !max number of components in mixture
      dimension an(58),un(58)
      dimension qb(21),fb(21),sb(21),wb(21)
      dimension ifp(ncmax)
      double precision K,K3
      double precision mrb(21),eb(21),kb(21),gb(21)
      double precision eijs(21,21),uij(21,21),kij(21,21),gijs(21,21)
      double precision kb2(21),eb2(21),kij2(21,21),uij2(21,21),
     &       gij2(21,21),bs2(18,21,21)
      dimension bs(18),cns(58)
      dimension x(ncmax),z(ncmax)
c
      common /NCOMP/ nc,ic
      common /EOSCOEF/ an,un
      common /FLDCOEF/ mrb,eb,kb,gb,qb,fb,sb,wb,eijs,uij,kij,gijs
      common /INTMCOEF/ kb2,eb2,kij2,uij2,gij2,bs2,ifp
      common /VAR/ K,bs,cns
c
      sum=0.d0
      do i=1,nc
        sum=sum+x(i)
      enddo
      if (sum.le.0.d0) RETURN
      do i=1,nc
        z(i)=x(i)/sum
      enddo
      K=0.d0
      F=0.d0
      U=0.d0
      Q=0.d0
      G=0.d0
      do i=1,nc
        K=K+z(i)*kb2(ifp(i))
        U=U+z(i)*eb2(ifp(i))
        F=F+z(i)**2*fb(ifp(i))
        Q=Q+z(i)*qb(ifp(i))
        G=G+z(i)*gb(ifp(i))
      enddo
      K=K**2
      U=U**2
      do i=1,nc-1
        do j=i+1,nc
          xij=z(i)*z(j)
          K=K+xij*kij2(ifp(i),ifp(j))
          U=U+xij*uij2(ifp(i),ifp(j))
          G=G+xij*gij2(ifp(i),ifp(j))
        enddo
      enddo
      K=K**0.2d0
      U=U**0.2d0
c
      Q=Q**2
      K3=K**3
      do n=1,18
        bs(n)=0.d0
        do i=1,nc
        do j=i,nc
          bs(n)=bs(n)+bs2(n,ifp(i),ifp(j))*z(i)*z(j)
        enddo
        enddo
        bs(n)=bs(n)*an(n)/K3
      enddo
      do n=13,58
        cns(n)=an(n)*U**un(n)
      enddo
      cns(13)=cns(13)*F
      cns(16)=cns(16)*Q
      cns(25)=cns(25)*G
      cns(26)=cns(26)*Q
      cns(27)=cns(27)*F
      cns(28)=cns(28)*Q
      cns(29)=cns(29)*G
      cns(30)=cns(30)*F
      cns(32)=cns(32)*G
      cns(33)=cns(33)*G
      cns(34)=cns(34)*G
      cns(35)=cns(35)*F
      cns(37)=cns(37)*Q
      cns(42)=cns(42)*Q
      cns(47)=cns(47)*Q
      cns(49)=cns(49)*Q
      cns(51)=cns(51)*G
      cns(52)=cns(52)*Q
      cns(54)=cns(54)*G
      cns(56)=cns(56)*G
      cns(58)=cns(58)*Q
c
      RETURN
      end                                              !subroutine SETAG
c
c ======================================================================
c
      block data AGA8CF
c
c  AGA8 equation of state coefficients taken from the DETAIL
c  compressibility factor equation of:
c
c  Starling, K.E. and Savidge, J.L.
c  Compressibility Factors of Natural Gas and Other Related Hydrocarbon
c  Gases,
c  Transmission Measurement Committee Report No. 8, Catalog No. XQ9212,
c  American Gas Association, 1994.
c
c  written by E.W. Lemmon, NIST Thermophysics Division, Boulder, Colorado
c  10-31-02 EWL, original version
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      dimension an(58),un(58)
      dimension qb(21),fb(21),sb(21),wb(21)
      double precision  acp0(21),bcp0(21),ccp0(21),dcp0(21),ecp0(21),
     &         fcp0(21),gcp0(21),hcp0(21),icp0(21),jcp0(21),kcp0(21)
      double precision mrb(21),eb(21),kb(21),gb(21)
      double precision eijs(21,21),uij(21,21),kij(21,21),gijs(21,21)
c
      common /EOSCOEF/ an,un
      common /FLDCOEF/ mrb,eb,kb,gb,qb,fb,sb,wb,eijs,uij,kij,gijs
      common /CP0COEF/ acp0,bcp0,ccp0,dcp0,ecp0,fcp0,
     &                 gcp0,hcp0,icp0,jcp0,kcp0
c
      data an/0.1538326d0,1.341953d0,-2.998583d0,-0.04831228d0,
     & 0.3757965d0,-1.589575d0,-0.05358847d0,0.88659463d0,-0.71023704d0,
     & -1.471722d0,1.32185035d0,-0.78665925d0,2.29129D-9,0.1576724d0,
     & -0.4363864d0,-0.04408159d0,-0.003433888d0,0.03205905d0,
     & 0.02487355d0,0.07332279d0,-0.001600573d0,0.6424706d0,
     & -0.4162601d0,-0.06689957d0,0.2791795d0,-0.6966051d0,
     & -0.002860589d0,-0.008098836d0,3.150547d0,0.007224479d0,
     & -0.7057529d0,0.5349792d0,-0.07931491d0,-1.418465d0,
     & -5.99905D-17,0.1058402d0,0.03431729d0,-0.007022847d0,
     & 0.02495587d0,0.04296818d0,0.7465453d0,-0.2919613d0,7.294616d0,
     & -9.936757d0,-0.005399808d0,-0.2432567d0,0.04987016d0,
     & 0.003733797d0,1.874951d0,0.002168144d0,-0.6587164d0,
     & 0.000205518d0,0.009776195d0,-0.02048708d0,0.01557322d0,
     & 0.006862415d0,-0.001226752d0,0.002850908d0/
c     data bn/18*1,9*2,10*3,7*4,5*5,2*6,2*7,3*8,2*9/
c     data cn/12*0,6*1,2*0,7*1,0,9*1,2*0,5*1,0,4*1,0,1,0,6*1/
c     data kn/12*0,3,3*2,2*4,2*0,3*2,4*4,0,2*1,2*2,2*3,3*4,2*0,3*2,2*4,
c    & 0,2*2,2*4,0,2,0,2,1,4*2/
      data un/0d0,0.5d0,1d0,3.5d0,-0.5d0,4.5d0,0.5d0,7.5d0,9.5d0,6d0,
     & 12d0,12.5d0,-6d0,2d0,3d0,2d0,2d0,11d0,-0.5d0,0.5d0,0d0,4d0,6d0,
     & 21d0,23d0,22d0,-1d0,-0.5d0,7d0,-1d0,6d0,4d0,1d0,9d0,-13d0,21d0,
     & 8d0,-0.5d0,0d0,2d0,7d0,9d0,22d0,23d0,1d0,9d0,3d0,8d0,23d0,1.5d0,
     & 5d0,-0.5d0,4d0,7d0,3d0,0d0,1d0,0d0/
c     data gn/4*0,2*1,18*0,1,3*0,1,2*0,3*1,16*0,1,2*0,1,0,1,2*0/
c     data qn/6*0,1,8*0,1,9*0,1,0,1,8*0,1,4*0,1,4*0,1,0,1,2*0,1,5*0,1/
c     data fn/12*0,1,13*0,1,2*0,1,4*0,1,23*0/
c     data sn/7*0,2*1,49*0/
c     data wn/9*0,3*1,46*0/
c
      data mrb/16.043d0,28.0135d0,44.01d0,30.07d0,44.097d0,18.0153d0,
     & 34.082d0,2.0159d0,28.01d0,31.9988d0,58.123d0,58.123d0,72.15d0,
     & 72.15d0,86.177d0,100.204d0,114.231d0,128.258d0,142.285d0,
     & 4.0026d0,39.948d0/
      data eb/151.3183d0,99.73778d0,241.9606d0,244.1667d0,298.1183d0,
     & 514.0156d0,296.355d0,26.95794d0,105.5348d0,122.7667d0,324.0689d0,
     & 337.6389d0,365.5999d0,370.6823d0,402.636293d0,427.72263d0,
     & 450.325022d0,470.840891d0,489.558373d0,2.610111d0,119.6299d0/
      data kb/0.4619255d0,0.4479153d0,0.4557489d0,0.5279209d0,
     & 0.583749d0,0.3825868d0,0.4618263d0,0.3514916d0,0.4533894d0,
     & 0.4186954d0,0.6406937d0,0.6341423d0,0.6738577d0,0.6798307d0,
     & 0.7175118d0,0.7525189d0,0.784955d0,0.8152731d0,0.8437826d0,
     & 0.3589888d0,0.4216551d0/
      data gb/0.d0,0.027815d0,0.189065d0,0.0793d0,0.141239d0,0.3325d0,
     & 0.0885d0,0.034369d0,0.038953d0,0.021d0,0.256692d0,0.281835d0,
     & 0.332267d0,0.366911d0,0.289731d0,0.337542d0,0.383381d0,
     & 0.427354d0,0.469659d0,0.d0,0.d0/
      data qb/2*0.d0,0.69d0,2*0.d0,1.06775d0,0.633276d0,14*0.d0/
      data fb/7*0.d0,1.d0,13*0.d0/
      data sb/5*0.d0,1.5822d0,0.39d0,14*0.d0/
      data wb/5*0.d0,1.d0,15*0.d0/
c
c  Binary interaction parameter values
      data (eijs(1,j),j=1,21)/1.d0,0.97164d0,0.960644d0,1.d0,0.994635d0,
     & 0.708218d0,0.931484d0,1.17052d0,0.990126d0,1.d0,1.01953d0,
     & 0.989844d0,1.00235d0,0.999268d0,1.107274d0,0.88088d0,0.880973d0,
     & 0.881067d0,0.881161d0,1.d0,1.d0/
      data (eijs(2,j),j=1,21)/2*1.d0,1.02274d0,0.97012d0,0.945939d0,
     & 0.746954d0,0.902271d0,1.08632d0,1.00571d0,1.021d0,0.946914d0,
     & 0.973384d0,0.95934d0,0.94552d0,7*1.d0/
      data (eijs(3,j),j=1,21)/3*1.d0,0.925053d0,0.960237d0,0.849408d0,
     & 0.955052d0,1.28179d0,1.5d0,1.d0,0.906849d0,0.897362d0,0.726255d0,
     & 0.859764d0,0.855134d0,0.831229d0,0.80831d0,0.786323d0,0.765171d0,
     & 2*1.d0/
      data (eijs(4,j),j=1,21)/4*1.d0,1.02256d0,0.693168d0,0.946871d0,
     & 1.16446d0,3*1.d0,1.01306d0,1.d0,1.00532d0,7*1.d0/
      data (eijs(5,j),j=1,21)/7*1.d0,1.034787d0,3*1.d0,1.0049d0,9*1.d0/
      data (eijs(6,j),j=1,21)/21*1.d0/
      data (eijs(7,j),j=1,21)/14*1.d0,1.008692d0,1.010126d0,1.011501d0,
     & 1.012821d0,1.014089d0,2*1.d0/
      data (eijs(8,j),j=1,21)/8*1.d0,1.1d0,1.d0,1.3d0,1.3d0,9*1.d0/
      data (eijs(9,j),j=1,21)/21*1.d0/
      data (eijs(10,j),j=1,21)/21*1.d0/
      data (eijs(11,j),j=1,21)/21*1.d0/
      data (eijs(12,j),j=1,21)/21*1.d0/
      data (eijs(13,j),j=1,21)/21*1.d0/
      data (eijs(14,j),j=1,21)/21*1.d0/
      data (eijs(15,j),j=1,21)/21*1.d0/
      data (eijs(16,j),j=1,21)/21*1.d0/
      data (eijs(17,j),j=1,21)/21*1.d0/
      data (eijs(18,j),j=1,21)/21*1.d0/
      data (eijs(19,j),j=1,21)/21*1.d0/
      data (eijs(20,j),j=1,21)/21*1.d0/
      data (eijs(21,j),j=1,21)/21*1.d0/
      data (uij(1,j),j=1,21)/1.d0,0.886106d0,0.963827d0,1.d0,0.990877d0,
     & 1.d0,0.736833d0,1.15639d0,3*1.d0,0.992291d0,1.d0,1.00367d0,
     & 1.302576d0,1.191904d0,1.205769d0,1.219634d0,1.233498d0,2*1.d0/
      data (uij(2,j),j=1,21)/2*1.d0,0.835058d0,0.816431d0,0.915502d0,
     & 1.d0,0.993476d0,0.408838d0,3*1.d0,0.993556d0,9*1.d0/
      data (uij(3,j),j=1,21)/3*1.d0,0.96987d0,2*1.d0,1.04529d0,1.d0,
     & 0.9d0,5*1.d0,1.066638d0,1.077634d0,1.088178d0,1.098291d0,
     & 1.108021d0,2*1.d0/
      data (uij(4,j),j=1,21)/4*1.d0,1.065173d0,1.d0,0.971926d0,
     & 1.61666d0,2*1.d0,4*1.25d0,7*1.d0/
      data (uij(5,j),j=1,21)/21*1.d0/
      data (uij(6,j),j=1,21)/21*1.d0/
      data (uij(7,j),j=1,21)/14*1.d0,1.028973d0,1.033754d0,1.038338d0,
     & 1.042735d0,1.046966d0,2*1.d0/
      data (uij(8,j),j=1,21)/21*1.d0/
      data (uij(9,j),j=1,21)/21*1.d0/
      data (uij(10,j),j=1,21)/21*1.d0/
      data (uij(11,j),j=1,21)/21*1.d0/
      data (uij(12,j),j=1,21)/21*1.d0/
      data (uij(13,j),j=1,21)/21*1.d0/
      data (uij(14,j),j=1,21)/21*1.d0/
      data (uij(15,j),j=1,21)/21*1.d0/
      data (uij(16,j),j=1,21)/21*1.d0/
      data (uij(17,j),j=1,21)/21*1.d0/
      data (uij(18,j),j=1,21)/21*1.d0/
      data (uij(19,j),j=1,21)/21*1.d0/
      data (uij(20,j),j=1,21)/21*1.d0/
      data (uij(21,j),j=1,21)/21*1.d0/
c
      data (kij(1,j),j=1,21)/1.d0,1.00363d0,0.995933d0,1.d0,1.007619d0,
     & 1.d0,1.00008d0,1.02326d0,3*1.d0,0.997596d0,1.d0,1.002529d0,
     & 0.982962d0,0.983565d0,0.982707d0,0.981849d0,0.980991d0,2*1.d0/
      data (kij(2,j),j=1,21)/2*1.d0,0.982361d0,1.00796d0,2*1.d0,
     & 0.942596d0,1.03227d0,13*1.d0/
      data (kij(3,j),j=1,21)/3*1.d0,1.00851d0,2*1.d0,1.00779d0,7*1.d0,
     & 0.910183d0,0.895362d0,0.881152d0,0.86752d0,0.854406d0,2*1.d0/
      data (kij(4,j),j=1,21)/4*1.d0,0.986893d0,1.d0,0.999969d0,
     & 1.02034d0,13*1.d0/
      data (kij(5,j),j=1,21)/21*1.d0/
      data (kij(6,j),j=1,21)/21*1.d0/
      data (kij(7,j),j=1,21)/14*1.d0,0.96813d0,0.96287d0,0.957828d0,
     & 0.952441d0,0.948338d0,2*1.d0/
      data (kij(8,j),j=1,21)/21*1.d0/
      data (kij(9,j),j=1,21)/21*1.d0/
      data (kij(10,j),j=1,21)/21*1.d0/
      data (kij(11,j),j=1,21)/21*1.d0/
      data (kij(12,j),j=1,21)/21*1.d0/
      data (kij(13,j),j=1,21)/21*1.d0/
      data (kij(14,j),j=1,21)/21*1.d0/
      data (kij(15,j),j=1,21)/21*1.d0/
      data (kij(16,j),j=1,21)/21*1.d0/
      data (kij(17,j),j=1,21)/21*1.d0/
      data (kij(18,j),j=1,21)/21*1.d0/
      data (kij(19,j),j=1,21)/21*1.d0/
      data (kij(20,j),j=1,21)/21*1.d0/
      data (kij(21,j),j=1,21)/21*1.d0/
c
c
      data (gijs(1,j),j=1,21)/2*1.d0,.807653d0,4*1.d0,1.95731d0,13*1.d0/
      data (gijs(2,j),j=1,21)/2*1.d0,0.982746d0,18*1.d0/
      data (gijs(3,j),j=1,21)/3*1.d0,0.370296d0,1.d0,1.67309d0,15*1.d0/
      data (gijs(4,j),j=1,21)/21*1.d0/
      data (gijs(5,j),j=1,21)/21*1.d0/
      data (gijs(6,j),j=1,21)/21*1.d0/
      data (gijs(7,j),j=1,21)/21*1.d0/
      data (gijs(8,j),j=1,21)/21*1.d0/
      data (gijs(9,j),j=1,21)/21*1.d0/
      data (gijs(10,j),j=1,21)/21*1.d0/
      data (gijs(11,j),j=1,21)/21*1.d0/
      data (gijs(12,j),j=1,21)/21*1.d0/
      data (gijs(13,j),j=1,21)/21*1.d0/
      data (gijs(14,j),j=1,21)/21*1.d0/
      data (gijs(15,j),j=1,21)/21*1.d0/
      data (gijs(16,j),j=1,21)/21*1.d0/
      data (gijs(17,j),j=1,21)/21*1.d0/
      data (gijs(18,j),j=1,21)/21*1.d0/
      data (gijs(19,j),j=1,21)/21*1.d0/
      data (gijs(20,j),j=1,21)/21*1.d0/
      data (gijs(21,j),j=1,21)/21*1.d0/
c
c  Cp0 coefficients given in:
c  McFall, R.L., M.S. Thesis, University of Oklahoma, 1984.
c  Aly, F.A. and Lee, L.L., Fluid Phase Equilib., 6:169, 1981.
c
      data acp0/-29776.4d0,-3495.34d0,20.7307d0,-37524.4d0,-56072.1d0,
     & -13773.1d0,-10085.4d0,-5565.6d0,-2753.49d0,-3497.45d0,-72387.d0,
     & -72674.8d0,-91505.5d0,-83845.2d0,-94982.5d0,-103353.d0,
     & -109674.d0,-122599.d0,-133564.d0,0d0,0d0/
      data bcp0/7.95454d0,6.95587d0,6.96237d0,7.98139d0,8.14319d0,
     & 7.97183d0,7.9468d0,6.66789d0,6.95854d0,6.96302d0,17.8143d0,
     & 18.6383d0,21.3861d0,22.5012d0,26.6225d0,30.4029d0,34.0847d0,
     & 38.5014d0,42.7143d0,4.968d0,4.968d0/
      data ccp0/43.9417d0,0.272892d0,2.68645d0,24.3668d0,37.0629d0,
     & 6.27078d0,-0.0838d0,2.33458d0,2.02441d0,2.40013d0,58.2062d0,
     & 57.4178d0,74.341d0,69.5789d0,80.3819d0,90.6941d0,100.253d0,
     & 111.446d0,122.173d0,0.d0,0.d0/
      data dcp0/1037.09d0,662.738d0,500.371d0,752.32d0,735.402d0,
     & 2572.63d0,433.801d0,2584.98d0,1541.22d0,2522.05d0,1787.39d0,
     & 1792.73d0,1701.58d0,1719.58d0,1718.49d0,1669.32d0,1611.55d0,
     & 1646.48d0,1654.85d0,100.d0,100.d0/
      data ecp0/1.56373d0,-0.291318d0,-2.56429d0,3.5399d0,9.38159d0,
     & 2.0501d0,2.85539d0,0.749019d0,0.096774d0,2.21752d0,40.7621d0,
     & 38.6599d0,47.0587d0,46.2164d0,55.6598d0,63.2028d0,69.7675d0,
     & 80.5015d0,90.2255d0,0.d0,0.d0/
      data fcp0/813.205d0,-680.562d0,-530.443d0,272.846d0,247.19d0,
     & 1156.72d0,843.792d0,559.656d0,3674.81d0,1154.15d0,808.645d0,
     & 814.151d0,775.899d0,802.174d0,802.069d0,786.001d0,768.847d0,
     & 781.588d0,785.564d0,100.d0,100.d0/
      data gcp0/-24.9027d0,1.7898d0,3.91921d0,8.44724d0,13.4556d0,
     & 0.d0,6.31595d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     & 0.d0,0.d0,0.d0,0.d0/
      data hcp0/1019.98d0,1740.06d0,500.198d0,1020.13d0,1454.78d0,
     & 100.d0,1481.43d0,100.d0,100.d0,100.d0,100.d0,100.d0,100.d0,
     & 100.d0,100.d0,100.d0,100.d0,100.d0,100.d0,100.d0,100.d0/
      data icp0/-10.1601d0,0.d0,2.1329d0,-13.2732d0,-11.7342d0,
     & 0.d0,-2.88457d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     & 0.d0,0.d0,0.d0,0.d0,0.d0/
      data jcp0/1070.14d0,100.d0,2197.22d0,869.51d0,984.518d0,
     & 100.d0,1102.23d0,100.d0,100.d0,100.d0,100.d0,100.d0,100.d0,
     & 100.d0,100.d0,100.d0,100.d0,100.d0,100.d0,100.d0,100.d0/
      data kcp0/-20.0615d0,4.49823d0,5.81381d0,-22.401d0,-24.0426d0,
     & -3.24989d0,-0.51551d0,-7.94821d0,6.23387d0,9.19749d0,-44.1341d0,
     & -46.1938d0,-60.2474d0,-62.2197d0,-77.5366d0,-92.0164d0,
     & -106.149d0,-122.444d0,-138.006d0,1.8198d0,8.6776d0/
      end
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c ======================================================================
c                                                    end file MIX_AGA8.f
c ======================================================================
      SUBROUTINE CHARGS(METHOD,HV,GR,X,TH,TD,PD,TGR,PGR,ierr,herr)
c
c  Setup the values required in the SGERG equation of state.
c
c  inputs:
c       hg--gross (or superior) heating value [J/mol]
c        x--composition [array of mol frac]
c
c  outputs:
c     ierr--error flag:  0 = successful
c                        1 = iteration failed
c     herr--error string (character*255 variable if ierr<>0)
c
c  written by E.W. Lemmon, NIST Physical & Chem Properties Div, Boulder, CO
c  01-13-10 EWL, original version
c


C       METHOD - Option number for selecting the method: (Input)
C         METHOD 1 -- Use gross caloric value (HV), relative density (Dr)
C                     and the mole fraction of carbon dioxide.
C         METHOD 2 -- Use relative density (Dr) and mole fractions
C                     of nitrogen and carbon dioxide.
C       HV     - Gross calorific heating value for the gas mixture
C                in kJ/dm^3. (Input for Method 1)
C       GR     - Relative density (specific gravity). (Input)
C       X      - An array of 5 elements containing the mole factions of:
C                X(1) - The equivalent hydrocarbon (Output)
C                X(2) - Nitrogen         (Input for Method 2)
C                X(3) - Carbon Dioxide   (Input for Methods 1 and 2)
C                X(4) - Hydrogen         (Input for Methods 1 and 2)
C                X(5) - Carbon Monoxide  (Input for Methods 1 and 2)
C       TH     - Reference temperature for heating value (K). (Input)
C       TD     - Reference temperature for molar density (K). (Input)
C       TGR    - Reference temperature for relative density (K). (Input)
C       PD     - Reference pressure for molar density (MPa). (Input)
C       PGR    - Reference pressure for relative density (MPa). (Input)
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      character*255 herr
      double precision B0(5,5),B1(5,5),B2(5,5),BB(0:2,0:3),Bm1(5,5)
      double precision Cm1(5,5),EB(5,5),EC(5,5),EBB(0:3),ECC(0:3)
      double precision C0(5,5),C1(5,5),C2(5,5),CC(0:2,0:3),CMT(1:6,0:1)
      COMMON/VIRIAL2/ BB,B0,B1,B2,CMT,CC,C0,C1,C2,TOLD,BBMIX,CCMIX,
     &                G1,G2,Bm1,Cm1,EB,EC,EBB,ECC

      double precision MW(5)
      COMMON/CONSTANTS/ RGAS,MW

      INTEGER FLAG
      double precision X(5),MR

      TOLD   = 0
      ierr = 0
      VIR    = -0.12527D0 + 5.91D-4*TGR - 6.62D-7*TGR**2
      D0AIR  = 28.96256D0/(RGAS*TGR/PGR + VIR)

c     G1     = -2.709328D0
c     G2     = 0.021062199D0
      HTV4   = 285.83D0
      HTV5   = 282.98D0

C-----------------------------------------------------------------------
C     Method 1 - Given the caloric value, specific gravity and
C                the mole fraction of CO2
C-----------------------------------------------------------------------
      IF (METHOD.EQ.1) THEN
        Z0     = 1.D0
        Z0TDPD = 1.D0
 300    HN0   = HV*Z0TDPD*RGAS*TD/PD*(1.D0 + 1.027D-4*(TH - 298.15D0))
        MR    = GR*Z0*RGAS*TGR/PGR*D0AIR
        SUM   = X(3)*(MW(2) - MW(3))
     &        + X(4)*(MW(2) - MW(4) + G2*HTV4)
     &        + X(5)*(MW(2) - MW(5) + G2*HTV5)
        X(1)  = (MR - G2*HN0 - MW(2) + SUM)/(G1 - MW(2))
        X(2)  = 1.D0 - X(1) - X(3) - X(4) - X(5)

        FLAG = 0
        IF (X(2).LT.0) THEN
          FLAG = 1
          X(2) = 0
          X(1) = 1.D0 - X(2) - X(3) - X(4) - X(5)
        ENDIF

        HCH   = (HN0 - X(4)*HTV4 - X(5)*HTV5)/X(1)
        IF (HCH.LT.0) HCH = 0
        MW(1) = G1 + G2*HCH
        BCH   = BB(0,0)+TD*(BB(0,1)+BB(0,2)*TD)+BB(0,3)/TD
     &        +(BB(1,0)+TD*(BB(1,1)+BB(1,2)*TD)+BB(1,3)/TD)*HCH
     &        +(BB(2,0)+TD*(BB(2,1)+BB(2,2)*TD)+BB(2,3)/TD)*HCH**2.D0
        CALL VIRGS(TD,X,BMIX,TEMP,BCH,1,ierr,herr)
        Z0TDPD = 1.D0 + BMIX*PD/RGAS/TD
        BCH   = BB(0,0)+TGR*(BB(0,1)+BB(0,2)*TGR)+BB(0,3)/TGR
     &        +(BB(1,0)+TGR*(BB(1,1)+BB(1,2)*TGR)+BB(1,3)/TGR)*HCH
     &        +(BB(2,0)+TGR*(BB(2,1)+BB(2,2)*TGR)+BB(1,3)/TGR)*HCH**2.D0
        CALL VIRGS(TGR,X,BMIX,TEMP,BCH,1,ierr,herr)
        Z0NEW = 1.D0 + BMIX*PGR/RGAS/TGR
        IF (ABS(Z0/Z0NEW - 1.D0).GT.0.5D-10) THEN
          Z0 = Z0NEW
          GOTO 300
        ENDIF

        X(2) = 1.D0 - X(1) - X(3) - X(4) - X(5)
        IF (X(2).LT.0 .OR. FLAG.EQ.1) THEN
          WRITE (*,*) 'CONFLICTING VALUES OF RELATIVE DENSITY, ',
     &                'HEATING VALUE, AND '
          write (*,*) 'CARBON DIOXIDE CONTENT'
          ierr = 3
          RETURN
        ENDIF

C-----------------------------------------------------------------------
C     Method 2 - Given the specific gravity and the mole fractions of
C                N2 and CO2
C-----------------------------------------------------------------------
      ELSEIF (METHOD.EQ.2) THEN
        Z0    = 1.D0
        X(1)  = 1.D0 - X(2) - X(3) - X(4) - X(5)
        i=1
 100    MR    = GR*Z0*RGAS*TGR/PGR*D0AIR
        MW(1) = (MR - X(2)*MW(2) - X(3)*MW(3)
     &        - X(4)*MW(4) - X(5)*MW(5))/X(1)
        HCH   = (MW(1) - G1)/G2
        BCH   = BB(0,0)+TGR*(BB(0,1)+BB(0,2)*TGR)+BB(0,3)*TGR**EB(1,1)
     &        +(BB(1,0)+TGR*(BB(1,1)+BB(1,2)*TGR)+BB(1,3)*TGR**EB(1,1))
     &         *HCH
     &        +BB(2,0)*HCH**EBB(0)+ TGR*(BB(2,1)*HCH**EBB(1)
     &        +BB(2,2)*HCH**EBB(2)*TGR)+ BB(2,3)*HCH**EBB(3)
     &         *TGR**EB(1,1)
c    &        +(BB(2,0)+TGR*(BB(2,1)+BB(2,2)*TGR)+BB(2,3)*TGR**EB(1,1))
c    &         *HCH**2.D0
        CALL VIRGS(TGR,X,BMIX,TEMP,BCH,1,ierr,herr)
        Z0NEW = 1.D0 + BMIX*PGR/RGAS/TGR
        i=i+1
        if (z0new.gt.100 .or. z0new.lt.0 .or. i.gt.100) then
          herr='bad Z0'
          ierr=1
          return
        endif
        IF (ABS(Z0/Z0NEW - 1.D0).GT.0.5D-10) THEN
          Z0 = Z0NEW
          GOTO 100
        ENDIF
      ENDIF

      IF (ierr.NE.0) RETURN

C.....Find the virial coefficient constants for pure hydrocarbon
C.....using the caloric value.
      B0(1,1) = BB(0,0) + HCH*BB(1,0) + BB(2,0)*HCH**EBB(0)
      B1(1,1) = BB(0,1) + HCH*BB(1,1) + BB(2,1)*HCH**EBB(1)
      B2(1,1) = BB(0,2) + HCH*BB(1,2) + BB(2,2)*HCH**EBB(2)
      Bm1(1,1)= BB(0,3) + HCH*BB(1,3) + BB(2,3)*HCH**EBB(3)
      C0(1,1) = CC(0,0) + HCH*CC(1,0) + CC(2,0)*HCH**ECC(0)
      C1(1,1) = CC(0,1) + HCH*CC(1,1) + CC(2,1)*HCH**ECC(1)
      C2(1,1) = CC(0,2) + HCH*CC(1,2) + CC(2,2)*HCH**ECC(2)
      Cm1(1,1)= CC(0,3) + HCH*CC(1,3) + CC(2,3)*HCH**ECC(3)
c     C0(1,1) = -0.302488D0  + HCH*( 0.646422D-3 - 0.332805D-06*HCH)
c     C1(1,1) =  0.195861D-2 + HCH*(-0.422876D-5 + 0.223160D-08*HCH)
c     C2(1,1) = -0.316302D-5 + HCH*( 0.688157D-8 - 0.367713D-11*HCH)

      END
c
c ======================================================================
c
      SUBROUTINE PARAMGS
C
C     PURPOSE:
C       Sets up constants used by the GERG model.
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      double precision MW(5)
      COMMON/CONSTANTS/ RGAS,MW

      double precision B0(5,5),B1(5,5),B2(5,5),BB(0:2,0:3),Bm1(5,5)
      double precision Cm1(5,5),EB(5,5),EC(5,5),EBB(0:3),ECC(0:3)
      double precision C0(5,5),C1(5,5),C2(5,5),CC(0:2,0:3),CMT(1:6,0:1)
      COMMON/VIRIAL2/ BB,B0,B1,B2,CMT,CC,C0,C1,C2,TOLD,BBMIX,CCMIX,
     &                G1,G2,Bm1,Cm1,EB,EC,EBB,ECC

C.....Store fluid information - Nitrogen is index 2, CO2 is index 3,
C     Hydrogen is index 4, CO is index 5, and the hydrocarbon is index 1.
      RGAS  = 8.314472D-3
      MW(2) = 28.0134D0
      MW(3) = 44.0095D0
      MW(4) = 2.01588D0
      MW(5) = 28.0101D0
C.....Assign virial coefficient constants for N2 and CO2.
      open (unit=11,file='sgerg.cof')
      read (11,*) B0(2,2)
      read (11,*) B1(2,2)
      read (11,*) B2(2,2)
      read (11,*) B0(3,3)
      read (11,*) B1(3,3)
      read (11,*) B2(3,3)
      read (11,*) B0(4,4)
      read (11,*) B1(4,4)
      read (11,*) B2(4,4)
      read (11,*) B0(5,5)
      read (11,*) B1(5,5)
      read (11,*) B2(5,5)
      read (11,*) B0(2,3)
      read (11,*) B1(2,3)
      read (11,*) B2(2,3)
      read (11,*) B0(1,4)
      read (11,*) B1(1,4)
      read (11,*) B2(1,4)
      read (11,*) B0(1,5)
      read (11,*) B1(1,5)
      read (11,*) B2(1,5)
      read (11,*) B0(2,4)
      read (11,*) BB(0,0)
      read (11,*) BB(1,0)
      read (11,*) BB(2,0)
      read (11,*) BB(0,1)
      read (11,*) BB(1,1)
      read (11,*) BB(2,1)
      read (11,*) BB(0,2)
      read (11,*) BB(1,2)
      read (11,*) BB(2,2)
      read (11,*) C0(2,2)
      read (11,*) C1(2,2)
      read (11,*) C2(2,2)
      read (11,*) C0(3,3)
      read (11,*) C1(3,3)
      read (11,*) C2(3,3)
      read (11,*) C0(2,3)
      read (11,*) C1(2,3)
      read (11,*) C2(2,3)
      read (11,*) C0(3,2)
      read (11,*) C1(3,2)
      read (11,*) C2(3,2)
      read (11,*) C0(4,4)
      read (11,*) C1(4,4)
      read (11,*) C2(4,4)
      read (11,*) C0(1,5)
      read (11,*) C1(1,5)
      read (11,*) C2(1,5)
      read (11,*) CC(0,0)
      read (11,*) CC(1,0)
      read (11,*) CC(2,0)
      read (11,*) CC(0,1)
      read (11,*) CC(1,1)
      read (11,*) CC(2,1)
      read (11,*) CC(0,2)
      read (11,*) CC(1,2)
      read (11,*) CC(2,2)
      read (11,*) CMT(1,0)
      read (11,*) CMT(1,1)
      read (11,*) CMT(2,0)
      read (11,*) CMT(2,1)
      read (11,*) CMT(3,0)
      read (11,*) CMT(3,1)
      read (11,*) CMT(4,0)
      read (11,*) CMT(4,1)
      read (11,*) CMT(5,0)
      read (11,*) CMT(5,1)
      read (11,*) CMT(6,0)
      read (11,*) CMT(6,1)
      read (11,*) B0(1,2)
      read (11,*) B1(1,2)
      read (11,*) B2(1,2)
      read (11,*) B0(1,3)
      read (11,*) G1
      read (11,*) G2
      read (11,*) Bm1(2,2)
      read (11,*) Bm1(3,3)
      read (11,*) Bm1(2,3)
      read (11,*) BB(0,3)
      read (11,*) BB(1,3)
      read (11,*) BB(2,3)
      read (11,*) CC(0,3)
      read (11,*) CC(1,3)
      read (11,*) CC(2,3)
      read (11,*) Cm1(2,2)
      read (11,*) Cm1(3,3)
      read (11,*) Cm1(2,3)
      read (11,*) Cm1(3,2)
      read (11,*) EB(1,1)
      read (11,*) EB(2,2)
      read (11,*) EB(3,3)
      read (11,*) EB(2,3)
      read (11,*) EC(1,1)
      read (11,*) EC(2,2)
      read (11,*) EC(3,3)
      read (11,*) EC(2,3)
      read (11,*) EC(3,2)
      read (11,*) EBB(0)
      read (11,*) EBB(1)
      read (11,*) EBB(2)
      read (11,*) EBB(3)
      read (11,*) ECC(0)
      read (11,*) ECC(1)
      read (11,*) ECC(2)
      read (11,*) ECC(3)
      close (11)
c     B0(2,2) = -0.144600D0
c     B1(2,2) =  0.740910D-3
c     B2(2,2) = -0.911950D-6
c     B0(3,3) = -0.868340D0
c     B1(3,3) =  0.403760D-2
c     B2(3,3) = -0.516570D-5
c     B0(4,4) = -0.110596D-2
c     B1(4,4) =  0.813385D-4
c     B2(4,4) = -0.987220D-7
c     B0(5,5) = -0.130820D0
c     B1(5,5) =  0.602540D-3
c     B2(5,5) = -0.644300D-6
c     B0(2,3) = -0.339693D0
c     B1(2,3) =  0.161176D-2
c     B2(2,3) = -0.204429D-5
c     B0(1,4) = -0.521280D-1
c     B1(1,4) =  0.271570D-3
c     B2(1,4) = -0.250000D-6
c     B0(1,5) = -0.687290D-1
c     B1(1,5) = -0.239381D-5
c     B2(1,5) =  0.518195D-6
c     B0(2,4) =  0.012D0
c     BB(0,0) = -0.425468D0
c     BB(1,0) =  0.877118D-3
c     BB(2,0) = -0.824747D-6
c     BB(0,1) =  0.286500D-2
c     BB(1,1) = -0.556281D-5
c     BB(2,1) =  0.431436D-8
c     BB(0,2) = -0.462073D-5
c     BB(1,2) =  0.881510D-8
c     BB(2,2) = -0.608319D-11
c     C0(2,2) =  0.784980D-2
c     C1(2,2) = -0.398950D-4
c     C2(2,2) =  0.611870D-7
c     C0(3,3) =  0.205130D-2
c     C1(3,3) =  0.348880D-4
c     C2(3,3) = -0.837030D-7
c     C0(2,3) =  0.552066D-2
c     C1(2,3) = -0.168609D-4
c     C2(2,3) =  0.157169D-7
c     C0(3,2) =  0.358783D-2
c     C1(3,2) =  0.806674D-5
c     C2(3,2) = -0.325798D-7
c     C0(4,4) =  0.104711D-2
c     C1(4,4) = -0.364887D-5
c     C2(4,4) =  0.467095D-8
c     C0(1,5) =  0.736748D-2
c     C1(1,5) = -0.276578D-4
c     C2(1,5) =  0.343051D-7
      END
c
c ======================================================================
c
      SUBROUTINE VIRGS(T,X,BMIX,CMIX,BCH,iopt,ierr,herr)
C
C     PURPOSE:
C       Calculates the second and third virial coefficients in the GERG
C       model at the given temperature.  The coefficient constants are
C       stored in arrays and are mixed using the combining rules.
C
C     DESCRIPTION OF ARGUMENTS:
C       T      - Temperature in kelvins. (Input)
C       BMIX   - Second virial coefficient of the mixture. (Output)
C       CMIX   - Third virial coefficient of the mixture. (Output)
C       BCH    - Binary CH-CH interaction coefficient. (Input/Output)
C       iopt    - Option number: (Input)
C         iopt = 0  -- Calculate BCH.
C         iopt = 1  -- Use BCH from input.
C       ierr - Error Flag: (Output)
C         ierr = 0       -- No error.
C         ierr = 1 or 2  -- Iteration failed to converge.
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      double precision B0(5,5),B1(5,5),B2(5,5),BB(0:2,0:3),Bm1(5,5)
      double precision Cm1(5,5),EB(5,5),EC(5,5),EBB(0:3),ECC(0:3)
      double precision C0(5,5),C1(5,5),C2(5,5),CC(0:2,0:3),CMT(1:6,0:1)
      COMMON/VIRIAL2/ BB,B0,B1,B2,CMT,CC,C0,C1,C2,TOLD,BBMIX,CCMIX,
     &                G1,G2,Bm1,Cm1,EB,EC,EBB,ECC

      character*255 herr
      double precision X(5)
      ierr = 0

      bmix=0
      cmix=0

C.....Calculate the second virial coefficient
      IF (T.EQ.TOLD) THEN
        BMIX = BBMIX
        CMIX = CCMIX
        RETURN
      ENDIF
      X11 = X(1)*X(1)
      X22 = X(2)*X(2)
      X33 = X(3)*X(3)
      X44 = X(4)*X(4)
      X55 = X(5)*X(5)
      IF (iopt.EQ.0) BCH = B0(1,1) + T*(B1(1,1) + B2(1,1)*T)
     &                   + Bm1(1,1)*T**EB(1,1)
      BN2  = B0(2,2) + T*(B1(2,2) + B2(2,2)*T) + Bm1(2,2)*T**EB(2,2)
      BCO2 = B0(3,3) + T*(B1(3,3) + B2(3,3)*T) + Bm1(3,3)*T**EB(3,3)
      BH2  = B0(4,4) + T*(B1(4,4) + B2(4,4)*T)
      BCO  = B0(5,5) + T*(B1(5,5) + B2(5,5)*T)
      IF (BCO2*BCH.LT.0.D0) THEN
        herr='VIRGS:  SQRT NEGATIVE'
        ierr = 1
        BMIX=1.D9
        RETURN
      ENDIF
      B12  = (B0(1,2)+ T*(B1(1,2) + B2(1,2)*T))*(BN2 + BCH)/2.D0
      B13  = B0(1,3)*DSQRT(BCO2*BCH)
c     B12  = (0.72D0 + 1.875D-5*(320.D0-T)*(320.D0-T))*(BN2 + BCH)/2.D0
c     B13  = -0.865D0*DSQRT(BCO2*BCH)
      B14  = B0(1,4) + T*(B1(1,4) + B2(1,4)*T)
      B15  = B0(1,5) + T*(B1(1,5) + B2(1,5)*T)
      B23  = B0(2,3) + T*(B1(2,3) + B2(2,3)*T) + Bm1(2,3)*T**EB(2,3)
      B24  = B0(2,4)
      BMIX = BCH*X11 + BN2*X22 + BCO2*X33 + BH2*X44 + BCO*X55
     &     + 2.D0*B12*X(1)*X(2) + 2.D0*B13*X(1)*X(3)
     &     + 2.D0*B14*X(1)*X(4) + 2.D0*B15*X(1)*X(5)
     &     + 2.D0*B23*X(2)*X(3) + 2.D0*B24*X(2)*X(4)

C.....Since methods 1 and 2 change X(1) and BCH, TOLD must not
C       be initialized until after the two methods have been set up.
C       During their setup, iopt is equal to one.
      IF (iopt.EQ.1) RETURN
      TOLD = T
      BBMIX = BMIX

C.....Calculate the third virial coefficient
c     E    = 0.92D0 + 0.0013D0*(T - 270.D0)
      F    = 1.D0/3.D0
      C11  = C0(1,1) + T*(C1(1,1) + C2(1,1)*T) + Cm1(1,1)*T**EC(1,1)
      C22  = C0(2,2) + T*(C1(2,2) + C2(2,2)*T) + Cm1(2,2)*T**EC(2,2)
      C33  = C0(3,3) + T*(C1(3,3) + C2(3,3)*T) + Cm1(3,3)*T**EC(3,3)
      C44  = C0(4,4) + T*(C1(4,4) + C2(4,4)*T)
      IF (C11.LT.0 .OR. C22.LT.0 .OR. C33.LT.0) THEN
        herr='INVALID TERM IN VIRGS'
        ierr = 1
        CMIX=1.D9
        RETURN
      ENDIF
      C15  = 3.D0*(C0(1,5) + T*(C1(1,5) + C2(1,5)*T)             )
      C23  = 3.D0*(C0(2,3)+T*(C1(2,3)+C2(2,3)*T)+Cm1(2,3)*T**EC(2,3))
      C32  = 3.D0*(C0(3,2)+T*(C1(3,2)+C2(3,2)*T)+Cm1(3,2)*T**EC(3,2))
      CMIX = C11*X11*X(1) + C22*X22*X(2) + C33*X33*X(3) + C44*X44*X(4)
      CMIX=CMIX+ C23*X22*X(3) + C32*X33*X(2) + C15*X11*X(5)
      CMIX=CMIX+ (CMT(1,0)+CMT(1,1)*(T - 270.D0))*3.D0*X11*X(2)*
     &           (C11*C11*C22)**F
      CMIX=CMIX+ (CMT(2,0)+CMT(2,1)*(T - 270.D0))*3.D0*X22*X(1)*
     &           (C11*C22*C22)**F
      CMIX=CMIX+ (CMT(3,0)+CMT(3,1)*(T - 270.D0))*3.D0*X11*X(3)*
     &           (C11*C11*C33)**F
      CMIX=CMIX+ (CMT(4,0)+CMT(4,1)*(T - 270.D0))*3.D0*X33*X(1)*
     &           (C11*C33*C33)**F
      CMIX=CMIX+ (CMT(5,0)+CMT(5,1)*(T - 270.D0))*3.D0*X11*X(4)*
     &           (C11*C11*C44)**F
      CMIX=CMIX+ (CMT(6,0)+CMT(6,1)*(T - 270.D0))*6.D0*X(1)*X(2)*X(3)*
     &           (C11*C22*C33)**F
c     CMIX = C11*X11*X(1) + C22*X22*X(2) + C33*X33*X(3) + C44*X44*X(4)
c    &     + C23*X22*X(3) + C32*X33*X(2) + C15*X11*X(5)
c    &     +      E*3.D0*X11* X(2)     *(C11*C11*C22)**F
c    &     +      E*3.D0*X22* X(1)     *(C11*C22*C22)**F
c    &     + 0.92D0*3.D0*X11* X(3)     *(C11*C11*C33)**F
c    &     + 0.92D0*3.D0*X33* X(1)     *(C11*C33*C33)**F
c    &     + 1.20D0*3.D0*X11* X(4)     *(C11*C11*C44)**F
c    &     + 1.10D0*6.D0*X(1)*X(2)*X(3)*(C11*C22*C33)**F
      CCMIX = CMIX
      END
c
c ======================================================================
c
      FUNCTION PGROSS (T,D,X)
C
C     PURPOSE:
C       Calculates the pressure from the GERG model as a function of
C       density and temperature.
C
C     DESCRIPTION OF ARGUMENTS:
C       D      - Molar density in mol/dm^3. (Input)
C       T      - Temperature in kelvins. (Input)
C       PGROSS - Pressure in MPa. (Output)
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      double precision MW(5),X(5)
      COMMON/CONSTANTS/ RGAS,MW
      character*255 herr

      Z = 0
      PGROSS = 0
      CALL VIRGS(T,X,BMIX,CMIX,TEMP,0,ierr,herr)
      IF (ierr.NE.0) RETURN
      Z = 1.D0 + BMIX*D + CMIX*D*D
      PGROSS = D*RGAS*T*Z
      END
