#
# This set of MEM routines solves A*X=D+E where:
#    X is an unknown Image
#    D is some measured Data
#    E is an error term
#    A is a linear operator i.e. a matrix in this discrete representation
#
# The solution is that which:
#    Maximizes the entropy
#       H=-SUM(I) [X(I)*LN(X(I)/M(I))]
#       where M(I) is a default image
#    Fits the data such that
#       SUM(J) [(Dobs(J)-DPredicted(J))/SIGMA(J)]**2 = OMEGA
#    Has total flux
#       SUM(I) X(I) = FLUX
#
# The algorithm used is that due to Wilczek and Drapatz. The dual
# problem of finding the Lagrange multipliers is used instead of the
# primal approach of finding X(I) directly. It is easy to show that the
# MEM image is given by:
#
#       X(I) = FLUX * EXP(-SUM(J) [ A(J,I) * LAMBDA(J)] ) / Z
#
# where the partition function is:
#
#   Z = SUM(I) EXP(-SUM(J) [ A(J,I) * LAMBDA(J)] )
#
# This version calls the Numerical Recipes Conjugate Gradients routine.
#
      subroutine mmsolve (dsz, isz, tol, niter, iflux, omega, defi, 
     $   dobsi, sigmai, memf, lambda)
#
# main routine for mem solution.
#
      real tol, memf(*), defi(*), dobsi(*), sigmai(*), lambda(*),
     $   iflux, omega
      integer dsz, isz, niter
#
      real flux
      integer datasz, imsz, mem, def, dobs, dpred, sigma
      common /memcomi/ flux, datasz, imsz, mem, def, dobs, dpred,
     $   sigma
      real mm(1)
      common /memcom/ mm
#
      real mu
#
#  set up pointers
#
      datasz = dsz
      imsz = isz
      mem = 1
      def = mem + imsz
      dobs = def + imsz + 1
      dpred = dobs + datasz + 1
      sigma = dpred + datasz + 1 
#
      flux = iflux
      do 1 j = 1, datasz
         mm(dobs+j-1) = dobsi(j)
         mm(sigma+j-1) = sigmai(j)
 1    continue
#
      mu = 0.0
      if(omega.gt.0.0) then
         mm(dobs+datasz) = omega
         do 100 j = 1, datasz
            mu = mu + (mm(dobs+j-1)*mm(sigma+j-1))**2
 100     continue
         mu = sqrt(mu/(4.0*omega))
      endif
      lambda(datasz+1) = mu
#
      do 10 i = 1, imsz
         mm(def+i-1) = defi(i)
 10   continue
#
      call frprmn (lambda, datasz, tol, niter, z)
#
      do 20 i = 1, imsz
         memf(i) = mm(mem+i-1) 
 20   continue
#
      end
#
      subroutine f (lambda, z)
      real lambda(*), z
      integer i
      real maxdr, lamtrm, norm
      data maxdr /1e38/
      real flux
      integer datasz, imsz, mem, def, dobs, dpred, sigma
      common /memcomi/ flux, datasz, imsz, mem, def, dobs, dpred,
     $   sigma
      real mm(1)
      common /memcom/ mm
#
# find a^t * lambda
#
      call mmatd (datasz, imsz, lambda, mm(mem))
#
# find new image and partition function
#
      lamtrm = log(maxdr)
      norm = 0.0
      do 10 i = 1, imsz
         mm(mem+i-1) = mm(def+i-1) * exp(min(-mm(mem+i-1),lamtrm))
         norm = norm + mm(mem+i-1)
  10  continue
      do 20 i = 1, imsz
         mm(mem+i-1) = flux * mm(mem+i-1) / norm
  20  continue
#
      z = flux*log(norm)
#
      end
#
      real function func(lambda)
      real lambda(*)
      integer j
      real z, mu, omega
#
      real flux
      integer datasz, imsz, mem, def, dobs, dpred, sigma
      common /memcomi/ flux, datasz, imsz, mem, def, dobs, dpred,
     $   sigma
      real mm(1)
      common /memcom/ mm
#
# find image for this set of lagrange multipliers
#
      call f (lambda, z)
      do 30 j = 1, datasz
         z = z + lambda(j) * mm(dobs+j-1)
   30 continue
      mu = lambda(datasz+1)
      if(mu.gt.0.0) then
         omega = mm(dobs+datasz)
         z = z + mu * omega 
         do 40 j = 1,datasz
            z = z + (lambda(j)*mm(sigma+j-1))**2 / (4.0*mu)
   40    continue
      endif
      func = z
      return
      end
#
      subroutine dfunc (lambda, gradz)
      real lambda(*), gradz(*)
#
      real flux
      integer datasz, imsz, mem, def, dobs, dpred, sigma
      common /memcomi/ flux, datasz, imsz, mem, def, dobs, dpred,
     $   sigma
      real mm(1)
      common /memcom/ mm
#
      real mu, z
      integer j
#
# find current image
#
      call f (lambda, z)
#
# find predicted data for this image using a*x
#
      call mmax (imsz, datasz, mm(mem), mm(dpred))
#
# now calculate gradient
#
      do 10 j = 1, datasz
         gradz(j) = (mm(dobs+j-1) - mm(dpred+j-1))
 10   continue
#
      mu=lambda(datasz+1)
      if(mu.gt.0.0) then
         gradz(datasz+1)=0.0
         do 20 j = 1, datasz
            gradz(datasz+1)=gradz(datasz+1)
     $         + (lambda(j)*mm(sigma+j-1))**2
            gradz(j) = gradz(j)
     $         + lambda(j)*mm(sigma+j-1)**2/(2.0*mu) 
 20      continue
      endif
#
      end
