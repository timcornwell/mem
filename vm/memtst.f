      PROGRAM MEMTEST
C
      INTEGER IMSZ,IMCNTR,ITER,NITER
      REAL MEM(64), DEF(64), GCH(64), STP(64), MODEL(64), DAT(64),
     $   WT(64), SNORM
      REAL TFLUX, TCHISQ, TOL, ENTRPY, FLUX, CHISQ, NRMGRD, SIGMA
      EXTERNAL GRADCH
C
      COMMON /ODDS/ DAT
C
      WRITE(*,*) ' '
      WRITE(*,*) '   MAXIMUM ENTROPY DECONVOLUTION PROGRAM 5.0'
      WRITE(*,*) '   -----------------------------------------'
      WRITE(*,*) ' '
      WRITE(*,*) '                 WRITTEN BY'
      WRITE(*,*) '             TIM CORNWELL, NRAO'
      WRITE(*,*) ' '
C
C Initial values
C
      Q = 1.0
      TOL = 0.1
      IMSZ = 64
      SIGMA = 5.0
      TCHISQ = IMSZ
      NITER = 100
C
C Make model and convolve
C
      DO 1 IMCNTR=1,64
    1 MODEL(IMCNTR)=0.0
      MODEL(33)=900.0
      DO 2 IMCNTR=34,48
    2 MODEL(IMCNTR)=20.0
      MODEL(49)=200.0
      CALL CONVOLVEARRAY(IMSZ,MODEL,DAT)
      do imcntr=1,64
         write(*,*) imcntr, dat(imcntr)
      end do

C
C Add noise: Gaussian + photon
C
      ISEED = 2000001
      TFLUX = 0.0
      DO 11 IMCNTR = 1,64
C         WT(IMCNTR) = 1.0 / (SIGMA**2 + DAT(IMCNTR))
         WT(IMCNTR) = 1.0
C         DAT(IMCNTR) = DAT(IMCNTR) + SNORM(ISEED) / SQRT(WT(IMCNTR))
         TFLUX = TFLUX + DAT(IMCNTR)
 11   CONTINUE
C
C Set up default image and starting image
C
      DO 10 IMCNTR = 1,64
         DEF(IMCNTR) = ABS(TFLUX)/64.0
         MEM(IMCNTR) = DEF(IMCNTR)
  10  CONTINUE
C
C Now iterate
C
      DO 20 ITER = 1, NITER
         CALL MEMONEIT (ITER.EQ.1, IMSZ, MEM, DEF, DEF(1),
     $      STP, GCH, WT, Q, TFLUX, TCHISQ, TOL, FLUX, CHISQ, ENTRPY,
     $      NRMGRD, GRADCH)
         WRITE(6,*) ITER,' Fit = ', SQRT(CHISQ/TCHISQ),' NRMGRD = ',
     $      NRMGRD,' FLUX = ', FLUX
   20 CONTINUE
C
      END
C----------------------------------------------------------------------
      SUBROUTINE GRADCH(IMSZ, MEM, WT, GCH, CHISQ)
C
C  Calculate CHISQ and its gradient.
C
      INTEGER IMSZ
      REAL MEM(*), GCH(*), WT(*), CHISQ
      INTEGER IMCNTR
      REAL DAT(64)
      COMMON /ODDS/ DAT
C      
      CALL CONVOLVEARRAY(IMSZ,MEM,GCH)
      CHISQ = 0.0
      DO 10 IMCNTR=1,IMSZ
         CHISQ=CHISQ+WT(IMCNTR)*(DAT(IMCNTR)-GCH(IMCNTR))**2
         GCH(IMCNTR) = 2.0 *(GCH(IMCNTR)-DAT(IMCNTR))*WT(IMCNTR)
   10 CONTINUE
      CALL CONVOLVEARRAY(IMSZ,GCH,GCH)
      RETURN
      END
C
C-------------------------------------------------------------------
      SUBROUTINE CONVOLVEARRAY(N,S,T)
C
C  Convolve circularly
C
      REAL T(*),S(*),TEMP(512)
      INTEGER N,J
      TEMP(1)=(S(N-1)+S(N)+S(1)+S(2)+S(3))
      TEMP(2)=(S(N)+S(1)+S(2)+S(3)+S(4))
      DO 1  J=3,N-2
    1 TEMP(J)=(S(J-2)+S(J-1)+S(J)+S(J+1)+S(J+2))
      TEMP(N-1)=(S(N-3)+S(N-2)+S(N-1)+S(N)+S(1))
      TEMP(N)=(S(N-2)+S(N-1)+S(N)+S(1)+S(2))
      DO 2 J=1,N
    2 T(J) = 0.2*TEMP(J)
      RETURN
      END
C
c     algorithm 599, collected algorithms from acm.
c     algorithm appeared in acm-trans. math. software, vol.9, no. 2,
c     jun., 1983, p. 255-257.
c**********************************************************************c
c**********************************************************************c
c**********************************************************************c
c                                                                      c
c                                                                      c
c                                                                      c
c     f o r t r a n  software package for random number generation     c
c                                                                      c
c                                                                      c
c                                                                      c
c**********************************************************************c
c**********************************************************************c
c**********************************************************************c
c
c
c
c     contents:
c
c     1) sunif  -  0,1 -uniform distribution
c
c     2) sexpo  - (standard-) exponential distribution
c
c     3) snorm  - (standard-) normal distribution
c
c     4) sgamma - (standard-) gamma distribution
c
c     5) kpoiss - poisson distribution
c
c
c     this package constitutes a fortran-77 documentation of a set of
c     assembler functions for sampling from the above distributions.
c     all routines make ample use of binary representations of numbers,
c     they are among the most accurate and fast sampling functions
c     known. the fortran programs below yield the same random number
c     sequences as the ones from our assembler package, but they are
c     of course much slower (by factors 5-8 on our siemens 7760
c     computer.)
c     the set of routines will also be acceptable to fortran iv
c     compilers which allow data statements for arrays without
c     implicit do-loops.
c
c
c     remarks:
c
c     -  no care is taken to ensure that the parameter values lie
c        in the allowed range (e.g. a/mu > 0.0 for sgamma/kpoiss).
c
c     -  the parameter 'ir' must be set to some  4*k+1 > 0  before
c        the first call of any of the generators. thereafter ir
c        must not be altered until a new initialization is desired.
c
c     -  the package provides random deviates of 6-7 digits accuracy.
c        on more accurate computers the constants in sexpo, snorm,
c        sgamma and kpoiss ought to be adjusted according to local
c        comments or with the aid of the tables in the literature
c        quoted at the beginning of each function.
c
c
c**********************************************************************c
c**********************************************************************c
c                                                                      c
c                                                                      c
c       0 , 1   - u n i f o r m  distribution                          c
c                                                                      c
c                                                                      c
c**********************************************************************c
c**********************************************************************c
c                                                                      c
c     for details see:                                                 c
c                                                                      c
c               ahrens, j.h., dieter, u. and grube, a.                 c
c               pseudo-random numbers:  a new proposal                 c
c                     for the choice of multiplicators                 c
c               computing, 6 (1970), 121 - 138                         c
c                                                                      c
c**********************************************************************c
c
      real function sunif(ir)
      double precision r,factor,two28
c
c     factor - integer of the form 8*k+5 as close as possible
c              to  2**26 * (sqrt(5)-1)/2     (golden section)
c     two28  = 2**28  (i.e. 28 significant bits for deviates)
c
      data factor /41475557.0d0/, two28 /268435456.0d0/
c
c     returns sample u from the  0,1 -uniform distribution
c     by a multiplicative congruential generator of the form
c        r := r * factor (mod 1) .
c     in the first call r is initialized to
c        r := ir / 2**28 ,
c     where ir must be of the form  ir = 4*k+1.
c     then r assumes all values  0 < (4*k+1)/2**28 < 1 during
c     a full period 2**26 of sunif.
c     the parameter ir is used only in the first call for
c     initialization of sunif. thereafter (when negative)
c     ir becomes a dummy variable.
c
      if (ir .ge. 0) go to 1
c
c     standard case:  sampling
c
      r=dmod(r*factor,1.0d0)
      sunif=sngl(r)
      return
c
c     first call: initialization
c
 1    r=dble(float(ir))/two28
      r=dmod(r*factor,1.0d0)
      sunif=sngl(r)
      ir=-1
      return
      end
c
c**********************************************************************c
c**********************************************************************c
c                                                                      c
c                                                                      c
c     (standard-)  e x p o n e n t i a l   distribution                c
c                                                                      c
c                                                                      c
c**********************************************************************c
c**********************************************************************c
c                                                                      c
c     for details see:                                                 c
c                                                                      c
c               ahrens, j.h. and dieter, u.                            c
c               computer methods for sampling from the                 c
c               exponential and normal distributions.                  c
c               comm. acm, 15,10 (oct. 1972), 873 - 882.               c
c                                                                      c
c     all statement numbers correspond to the steps of algorithm       c
c     'sa' in the above paper (slightly modified implementation)       c
c                                                                      c
c**********************************************************************c
c
      real function sexpo(ir)
      dimension q(8)
      equivalence (q(1),q1)
c
c     q(n) = sum(alog(2.0)**k/k!)    k=1,..,n ,      the highest n
c     (here 8) is determined by q(n)=1.0 within standard precision
c
      data q/.6931472,.9333737,.9888778,.9984959,
     ,.9998293,.9999833,.9999986,.9999999/
c
   1  a=0.0
      u=sunif(ir)
      go to 2
   3  a=a+q1
   2  u=u+u
      if (u.le.1.0) go to 3
   4  u=u-1.0
      if (u.gt.q1) go to 6
   5  sexpo=a+u
      return
   6  i=1
      ustar=sunif(ir)
      umin=ustar
   7  ustar=sunif(ir)
      if (ustar.lt.umin) umin=ustar
   8  i=i+1
      if (u.gt.q(i)) go to 7
   9  sexpo=a+umin*q1
      return
      end
c
c**********************************************************************c
c**********************************************************************c
c                                                                      c
c                                                                      c
c     (standard-)  n o r m a l  distribution                           c
c                                                                      c
c                                                                      c
c**********************************************************************c
c**********************************************************************c
c                                                                      c
c     for details see:                                                 c
c                                                                      c
c               ahrens, j.h. and dieter, u.                            c
c               extensions of forsythe's method for random             c
c               sampling from the normal distribution.                 c
c               math. comput., 27,124 (oct. 1973), 927 - 937.          c
c                                                                      c
c     all statement numbers correspond to the steps of algorithm 'fl'  c
c     (m=5) in the above paper     (slightly modified implementation)  c
c                                                                      c
c**********************************************************************c
c
      real function snorm(ir)
      dimension a(32),d(31),t(31),h(31)
c
c     the definitions of the constants a(k), d(k), t(k) and
c     h(k) are according to the abovementioned article
c
      data a/0.0,.3917609e-1,.7841241e-1,.1177699,.1573107,
     ,.1970991,.2372021,.2776904,.3186394,.3601299,.4022501,
     ,.4450965,.4887764,.5334097,.5791322,.6260990,.6744898,
     ,.7245144,.7764218,.8305109,.8871466,.9467818,1.009990,
     ,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,
     ,1.675940,1.862732,2.153875/
      data d/5*0.0,.2636843,.2425085,.2255674,.2116342,.1999243,
     ,.1899108,.1812252,.1736014,.1668419,.1607967,.1553497,
     ,.1504094,.1459026,.1417700,.1379632,.1344418,.1311722,
     ,.1281260,.1252791,.1226109,.1201036,.1177417,.1155119,
     ,.1134023,.1114027,.1095039/
      data t/.7673828e-3,.2306870e-2,.3860618e-2,.5438454e-2,
     ,.7050699e-2,.8708396e-2,.1042357e-1,.1220953e-1,.1408125e-1,
     ,.1605579e-1,.1815290e-1,.2039573e-1,.2281177e-1,.2543407e-1,
     ,.2830296e-1,.3146822e-1,.3499233e-1,.3895483e-1,.4345878e-1,
     ,.4864035e-1,.5468334e-1,.6184222e-1,.7047983e-1,.8113195e-1,
     ,.9462444e-1,.1123001,.1364980,.1716886,.2276241,.3304980,
     ,.5847031/
      data h/.3920617e-1,.3932705e-1,.3950999e-1,.3975703e-1,
     ,.4007093e-1,.4045533e-1,.4091481e-1,.4145507e-1,.4208311e-1,
     ,.4280748e-1,.4363863e-1,.4458932e-1,.4567523e-1,.4691571e-1,
     ,.4833487e-1,.4996298e-1,.5183859e-1,.5401138e-1,.5654656e-1,
     ,.5953130e-1,.6308489e-1,.6737503e-1,.7264544e-1,.7926471e-1,
     ,.8781922e-1,.9930398e-1,.1155599,.1404344,.1836142,.2790016,
     ,.7010474/
c
   1  u=sunif(ir)
      s=0.0
      if (u.ge.0.5) s=1.0
      u=u+u-s
   2  u=32.0*u
      i=int(u)
      if (i.eq.0) go to 9
c
c                                start center
c
   3  ustar=u-float(i)
      aa=a(i)
   4  if (ustar.le.t(i)) go to 5
      w=(ustar-t(i))*h(i)
c
c                                exit   (both cases)
c
  17  y=aa+w
      snorm=y
      if (s.eq.1.0) snorm=-y
      return
c
c                                center continued
c
   5  u=sunif(ir)
      w=u*(a(i+1)-aa)
      tt=(0.5*w+aa)*w
      go to 6
   8  tt=u
      ustar=sunif(ir)
   6  if (ustar.gt.tt) go to 17
   7  u=sunif(ir)
      if (ustar.ge.u) go to 8
      ustar=sunif(ir)
      go to 4
c
c                                start tail
c
   9  i=6
      aa=a(32)
      go to 10
  11  aa=aa+d(i)
      i=i+1
  10  u=u+u
      if (u.lt.1.0) go to 11
  12  u=u-1.0
  13  w=u*d(i)
      tt=(0.5*w+aa)*w
      go to 14
  16  tt=u
  14  ustar=sunif(ir)
      if (ustar.gt.tt) go to 17
  15  u=sunif(ir)
      if (ustar.ge.u) go to 16
      u=sunif(ir)
      go to 13
      end
c
c**********************************************************************c
c**********************************************************************c
c                                                                      c
c                                                                      c
c     (standard-)  g a m m a  distribution                             c
c                                                                      c
c                                                                      c
c**********************************************************************c
c**********************************************************************c
c                                                                      c
c               parameter  a >= 1.0  !                                 c
c                                                                      c
c**********************************************************************c
c                                                                      c
c     for details see:                                                 c
c                                                                      c
c               ahrens, j.h. and dieter, u.                            c
c               generating gamma variates by a                         c
c               modified rejection technique.                          c
c               comm. acm, 25,1 (jan. 1982), 47 - 54.                  c
c                                                                      c
c     step numbers correspond to algorithm 'gd' in the above paper     c
c                                 (straightforward implementation)     c
c                                                                      c
c**********************************************************************c
c                                                                      c
c               parameter  0.0 < a < 1.0  !                            c
c                                                                      c
c**********************************************************************c
c                                                                      c
c     for details see:                                                 c
c                                                                      c
c               ahrens, j.h. and dieter, u.                            c
c               computer methods for sampling from gamma,              c
c               beta, poisson and binomial distributions.              c
c               computing, 12 (1974), 223 - 246.                       c
c                                                                      c
c     (adapted implementation of algorithm 'gs' in the above paper)    c
c                                                                      c
c**********************************************************************c
c
      real function sgamma(ir,a)
c
c     input:  ir=current state of basic random number generator
c             a =parameter (mean) of the standard gamma distribution
c     output: sgamma = sample from the gamma-(a)-distribution
c
c     coefficients q(k) - for q0 = sum(q(k)*a**(-k))
c     coefficients a(k) - for q = q0+(t*t/2)*sum(a(k)*v**k)
c     coefficients e(k) - for exp(q)-1 = sum(e(k)*q**k)
c
      data q1,q2,q3,q4,q5,q6,q7 /.04166669,.02083148,
     ,.00801191,.00144121,-.00007388,.00024511,.00024240/
      data a1,a2,a3,a4,a5,a6,a7 /.3333333,-.2500030,
     ,.2000062,-.1662921,.1423657,-.1367177,.1233795/
      data e1,e2,e3,e4,e5 /1.,.4999897,.1668290,.0407753,.0102930/
c
c     previous a pre-set to zero - aa is a', aaa is a"
c     sqrt32 is the squareroot of 32 = 5.656854249492380
c
      data aa /0.0/, aaa /0.0/, sqrt32 /5.656854/
c
      if (a .eq. aa) go to 1
      if (a .lt. 1.0) go to 12
c
c     step  1:  recalculations of s2,s,d if a has changed
c
      aa=a
      s2=a-0.5
      s=sqrt(s2)
      d=sqrt32-12.0*s
c
c     step  2:  t=standard normal deviate,
c               x=(s,1/2)-normal deviate.
c               immediate acceptance (i)
c
   1  t=snorm(ir)
      x=s+0.5*t
      sgamma=x*x
      if (t .ge. 0.0) return
c
c     step  3:  u= 0,1 -uniform sample. squeeze acceptance (s)
c
      u=sunif(ir)
      if (d*u .le. t*t*t) return
c
c     step  4:  recalculations of q0,b,si,c if necessary
c
      if (a .eq. aaa) go to 4
      aaa=a
      r=1.0/a
      q0=((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r
c
c               approximation depending on size of parameter a
c               the constants in the expressions for b, si and
c               c were established by numerical experiments
c
      if (a .le. 3.686) go to 3
      if (a .le. 13.022) go to 2
c
c               case 3:  a .gt. 13.022
c
      b=1.77
      si=.75
      c=.1515/s
      go to 4
c
c               case 2:  3.686 .lt. a .le. 13.022
c
   2  b=1.654+.0076*s2
      si=1.68/s+.275
      c=.062/s+.024
      go to 4
c
c               case 1:  a .le. 3.686
c
   3  b=.463+s-.178*s2
      si=1.235
      c=.195/s-.079+.016*s
c
c     step  5:  no quotient test if x not positive
c
   4  if (x .le. 0.0) go to 7
c
c     step  6:  calculation of v and quotient q
c
      v=t/(s+s)
      if (abs(v) .le. 0.25) go to 5
      q=q0-s*t+0.25*t*t+(s2+s2)*alog(1.0+v)
      go to 6
   5  q=q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
c
c     step  7:  quotient acceptance (q)
c
   6  if (alog(1.0-u) .le. q) return
c
c     step  8:  e=standard exponential deviate
c               u= 0,1 -uniform deviate
c               t=(b,si)-double exponential (laplace) sample
c
   7  e=sexpo(ir)
      u=sunif(ir)
      u=u+u-1.0
      t=b+sign(si*e,u)
c
c     step  9:  rejection if t .lt. tau(1) = -.71874483771719
c
      if (t .lt. (-.7187449)) go to 7
c
c     step 10:  calculation of v and quotient q
c
      v=t/(s+s)
      if (abs(v) .le. 0.25) go to 8
      q=q0-s*t+0.25*t*t+(s2+s2)*alog(1.0+v)
      go to 9
   8  q=q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
c
c     step 11:  hat acceptance (h) (if q not positive go to step 8)
c
   9  if (q .le. 0.0) go to 7
      if (q .le. 0.5) go to 10
      w=exp(q)-1.0
      go to 11
  10  w=((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q
c
c               if t is rejected, sample again at step 8
c
  11  if (c*abs(u) .gt. w*exp(e-0.5*t*t)) go to 7
      x=s+0.5*t
      sgamma=x*x
      return
c
c     alternate method for parameters a below 1  (.3678794=exp(-1.))
c
  12  aa=0.0
      b=1.0+.3678794*a
  13  p=b*sunif(ir)
      if (p .ge. 1.0) go to 14
      sgamma=exp(alog(p)/a)
      if (sexpo(ir) .lt. sgamma) go to 13
      return
  14  sgamma=-alog((b-p)/a)
      if (sexpo(ir) .lt. (1.0-a)*alog(sgamma)) go to 13
      return
      end
c
c**********************************************************************c
c**********************************************************************c
c                                                                      c
c                                                                      c
c     p o i s s o n  distribution                                      c
c                                                                      c
c                                                                      c
c**********************************************************************c
c**********************************************************************c
c                                                                      c
c     for details see:                                                 c
c                                                                      c
c               ahrens, j.h. and dieter, u.                            c
c               computer generation of poisson deviates                c
c               from modified normal distributions.                    c
c               acm trans. math. software, 8,2 (june 1982), 163 - 179. c
c                                                                      c
c     (slightly modified version of the program in the above article)  c
c                                                                      c
c**********************************************************************c
c
      integer function kpoiss(ir,mu)
c
c     input:  ir=current state of basic random number generator
c             mu=mean mu of the poisson distribution
c     output: kpoiss=sample from the poisson-(mu)-distribution
c
      real mu, muprev, muold
c
c     muprev=previous mu, muold=mu at last execution of step p or b.
c     tables: coefficients a0-a7 for step f. factorials fact
c     coefficients a(k) - for px = fk*v*v*sum(a(k)*v**k)-del
c
      dimension fact(10), pp(35)
      data muprev,muold /0.,0./
      data a0,a1,a2,a3,a4,a5,a6,a7 /-.5,.3333333,-.2500068,
     ,.2000118,-.1661269,.1421878,-.1384794,.1250060/
      data fact /1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880./
c
c     separation of cases a and b
c
      if (mu .eq. muprev) go to 1
      if (mu .lt. 10.0) go to 12
c
c     c a s e  a. (recalculation of s,d,l if mu has changed)
c
      muprev=mu
      s=sqrt(mu)
      d=6.0*mu*mu
c
c             the poisson probabilities pk exceed the discrete normal
c             probabilities fk whenever k >= m(mu). l=ifix(mu-1.1484)
c             is an upper bound to m(mu) for all mu >= 10 .
c
      l=ifix(mu-1.1484)
c
c     step n. normal sample - snorm(ir) for standard normal deviate
c
   1  g=mu+s*snorm(ir)
      if (g .lt. 0.0) go to 2
      kpoiss=ifix(g)
c
c     step i. immediate acceptance if kpoiss is large enough
c
      if (kpoiss .ge. l) return
c
c     step s. squeeze acceptance - sunif(ir) for (0,1)-sample u
c
      fk=float(kpoiss)
      difmuk=mu-fk
      u=sunif(ir)
      if (d*u .ge. difmuk*difmuk*difmuk) return
c
c     step p. preparations for steps q and h.
c             (recalculations of parameters if necessary)
c             .3989423=(2*pi)**(-.5)  .416667e-1=1./24.  .1428571=1./7.
c             the quantities b1, b2, c3, c2, c1, c0 are for the hermite
c             approximations to the discrete normal probabilities fk.
c             c=.1069/mu guarantees majorization by the 'hat'-function.
c
   2  if (mu .eq. muold) go to 3
      muold=mu
      omega=.3989423/s
      b1=.4166667e-1/mu
      b2=.3*b1*b1
      c3=.1428571*b1*b2
      c2=b2-15.*c3
      c1=b1-6.*b2+45.*c3
      c0=1.-b1+3.*b2-15.*c3
      c=.1069/mu
   3  if (g .lt. 0.0) go to 5
c
c             'subroutine' f is called (kflag=0 for correct return)
c
      kflag=0
      go to 7
c
c     step q. quotient acceptance (rare case)
c
   4  if (fy-u*fy .le. py*exp(px-fx)) return
c
c     step e. exponential sample - sexpo(ir) for standard exponential
c             deviate e and sample t from the laplace 'hat'
c             (if t <= -.6744 then pk < fk for all mu >= 10.)
c
   5  e=sexpo(ir)
      u=sunif(ir)
      u=u+u-1.0
      t=1.8+sign(e,u)
      if (t .le. (-.6744)) go to 5
      kpoiss=ifix(mu+s*t)
      fk=float(kpoiss)
      difmuk=mu-fk
c
c             'subroutine' f is called (kflag=1 for correct return)
c
      kflag=1
      go to 7
c
c     step h. hat acceptance (e is repeated on rejection)
c
   6  if (c*abs(u) .gt. py*exp(px+e)-fy*exp(fx+e)) go to 5
      return
c
c     step f. 'subroutine' f. calculation of px,py,fx,fy.
c             case kpoiss .lt. 10 uses factorials from table fact
c
   7  if (kpoiss .ge. 10) go to 8
      px=-mu
      py=mu**kpoiss/fact(kpoiss+1)
      go to 11
c
c             case kpoiss .ge. 10 uses polynomial approximation
c             a0-a7 for accuracy when advisable
c             .8333333e-1=1./12.  .3989423=(2*pi)**(-.5)
c
   8  del=.8333333e-1/fk
      del=del-4.8*del*del*del
      v=difmuk/fk
      if (abs(v) .le. 0.25) go to 9
      px=fk*alog(1.0+v)-difmuk-del
      go to 10
   9  px=fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)-del
  10  py=.3989423/sqrt(fk)
  11  x=(0.5-difmuk)/s
      xx=x*x
      fx=-0.5*xx
      fy=omega*(((c3*xx+c2)*xx+c1)*xx+c0)
      if (kflag) 4,4,6
c
c     c a s e  b. (start new table and calculate p0 if necessary)
c
  12  muprev=0.0
      if (mu .eq. muold) go to 13
      muold=mu
      m=max0(1,ifix(mu))
      l=0
      p=exp(-mu)
      q=p
      p0=p
c
c     step u. uniform sample for inversion method
c
  13  u=sunif(ir)
      kpoiss=0
      if (u .le. p0) return
c
c     step t. table comparison until the end pp(l) of the
c             pp-table of cumulative poisson probabilities
c             (0.458=pp(9) for mu=10)
c
      if (l .eq. 0) go to 15
      j=1
      if (u .gt. 0.458) j=min0(l,m)
      do 14 k=j,l
      if (u .le. pp(k)) go to 18
  14  continue
      if (l .eq. 35) go to 13
c
c     step c. creation of new poisson probabilities p
c             and their cumulatives q=pp(k)
c
  15  l=l+1
      do 16 k=l,35
      p=p*mu/float(k)
      q=q+p
      pp(k)=q
      if (u .le. q) go to 17
  16  continue
      l=35
      go to 13
  17  l=k
  18  kpoiss=k
      return
      end

