/*.he mwd -- Maximum Entroy Deconvolution                   @*/
/*.fo                                       - # -*/
/*.ln 10*/
/*.pl 80*/
/*.xr on*/
/*.rw ansi*/
/*.ba Courier*/
/*.aa Courier-Bold*/
/*.us on*/
/*.pa Algorithm Description*/
/*

		A Compact Maximum Entropy Algorithm
		___________________________________

			T.J. Cornwell NRAO/VLA
			J.C. Krueger Alfred U.

			     July 1991


                          Introduction
                         ______________


This maximum entropy program uses an algorithm that is outlined in a
paper published in Astronomy and Astrophysics by Wilczek and Drapatz
(142, 9-12, 1985).  Our program uses a slightly different form of
entropy than outlined, to allow for the use of a priori information.
We also use a different method of optimization than that is outlined
in the paper.

The form of entropy used is: E = -SUM_i f_i ln [f_i/(e.m_i)], where SUM
is a summation from i=1 to the number of pixels in your image.  The a
priori image is m_i, the MEM image is f_i, and e is the base of natural
logarithms i.e. e=2.78..)  The entropy is maximized while being
constrained to constraints on:

	(i) Each data point is fit to within an error comaprable to
	    the noise level.
       (ii) Chi-squared equal to OMEGA.
      (iii) Integral of the image equal to FLUX (this is optional
	    in our algorithm).

Wilczek and Drapatz used a Newton iteration to optimize their function,
but because their images were small, the inversion of the Hessian matrix
was not a problem.  In our application to radio data, our images are much
bigger and the matrix inversion presented a problem, so we choose to use
a conjugate gradient method. Hence this algorithm should be of interest
to those with large datasets. The conjugate gradient algorithm seems to
work reasonably well if no more than about 10 conjugate steps are
performed sequentially before reseting the search direction.

*/
/*.cp*/
/*
                         The core package
                        __________________


This file contains the heart of the algorithm.  This set of procedures
and functions is called for every iteration and does several
calculations.  It is well documented and should be easy to follow,
and is not necessary to change any of the internal routines.  The
interface to the package is via the procedure MwdOneIt:

   MwdOneIt(doinit,datasz,imsz,nconj,flux,mem,def,
            sigma,dobs,dpred,dzdlam,dellam,pdzdla,iomega,tol,length,
            lambda,mu,z,FTrans,LTrans)

The arguments are detailed below:

_______________________________________________________________________________
In or	Name and description				       Type of variable
Out?
_______________________________________________________________________________

I   doinit......True the first call, false otherwise.			logical
I   datasz......The number of data points.				int
I   imsz........The number of pixels in the image.			int
I   nconj.......Number of conjugate steps to take in a row.		int
                We have found that 10 is a good number.
                0 produces steepest ascent
I   flux........The desired integral of final image, if used.		double
                Set flux to 0 if this constraint is not desired.
IO  mem(*)......The current MEM image being reconstructed.	   imsz*double
                The final image will be in this array.
I   def(*)......The a priori image, if used.  Set def to a	   imsz*double
                constant value if an a priori image is not used.
I   sigma(*)....The expected error in the each data point.	 datasz*double
I   dobs(*).....The measured data.				 datasz*double
O   dpred(*)....The data predicted from the current model MEM	 datasz*double
IO  dzdlam(*)...Work array (must be preserved between calls).	 datasz*double
IO  dellam(*)...Work array (must be preserved between calls).	 datasz*double
IO  pdzdla(*)...Work array (must be preserved between calls). 	 datasz*double
I   iomega......The required Chi-squared value.				double
I   tol.........A parameter used to limit the size of the change	double
                in the Lagrange multipliers.  Experimentation
                has shown that 0.01 is a good value.
O   length......GRAD is actually the length of the change in the	double
                Lagrange multiplier.  As the program approaches
                the minimum value of Z, length should get
                smaller.
IO  lambda(*)...Lagrange multipliers for DOBS (the data array).	 datasz*double
O   mu..........Lagrange multipler for flux.				double
O   z...........The potential function (see the Wilczek and		double
                Drapatz paper).
I   FTrans......Procedure to transform "mem".				external
		FTrans (imsz,datasz,mem,dpred)
I   LTrans......Routine for transforming "lambda".			external
		LTrans (datasz,imsz,lambda,mem)
________________________________________________________________________


   Neither the total flux nor chi-squared are calculated explicitly
inside the package. Chi-squared is defined by the loop:

	chisq = 0.0;
	for ( n = 0 ; n < datasz ; n++ )
         {
	  t = (dobs(n) - dpred(n))/sigma(n);
          chisq += t*t;
         }

And so the expected value of Chi-squared (defined by "iomega") is usually
about "datasz".

*/
/*.cp*/
/*

                           Disclaimer
                          ____________


The authors take no responsibility for any errors in this package.  You're
on your own.  However, Tim Cornwell would be interested in any remarks you
may have (tcornwel@nrao.edu).


_________________________________________________________________________

   Maximum Entropy Deconvolution Program

     T.J. Cornwell
     National Radio Astronomy Observatory
     P.O. Box 0,
     Socorro, New Mexico, 87801, USA.

     tcornwel@nrao.edu
     505 835 7333


   C version prepared by

     Robert W. Conley
     PL/LIMI
     3550 Aberdeen Ave SE
     Kirtland AFB, NM 87117-5776
     conley@plk.af.mil
     (505) 846-4844

   Notes to the C version:

     The translation of this code from its original FORTRAN was made
     reflecting the natural way C indexes arrays (which count from 0),
     rather than simply transliterating the FORTRAN (where arrays
     count from 1).

     Comments which begin "/*.xx", where "xx" is a two character code
     such as "cp" or "pa",  are instructions to a C program lister
     and cross reference generator.  Remove these comments if you find
     them objectionable.

_________________________________________________________________________

*/

/*.pa Global Declarations*/
#define	Abs(a)		((a) >= 0   ? (a) : -(a))
#define	Min(a,b)	((a) <= (b) ? (a) :  (b))
#define	Max(a,b)	((a) >= (b) ? (a) :  (b))

/*.pa Supporting procedures and functions*/
/*
   Procedure MwdCDel

     Parameters:  doinit......log   I   True ==> first call to MwdOneIt
                  datasz......int   I   Size of data set
                  nconj.......int   I   Number of successive conjugate steps
                  lambda(*)...real  I   Array of Lagrange multipliers
                  dzdlam(*)...real  I   Gradients of Z w.r.t. lambda
                  dellam(*)...real  IO  Actual step in lambda
                  grad........real  I   Current grad  . Current grad
                  pgrad.......real  I   Previous grad . Previous grad

     Calculate Conjugate gradient step.
*/
void MwdCDel (doinit,datasz,nconj,lambda,dzdlam,dellam,grad,pgrad)
int doinit,datasz,nconj;
double *lambda,*dzdlam,*dellam,grad,pgrad;
{
  static int iconj;
  double gamma;
  int dcntr;


  if (doinit)
   {
    iconj = 0;
    gamma = 0.0;
   }
  else if (iconj >= nconj)
   {
    iconj = 0;
    gamma = 0.0;
   }
  else if (pgrad != 0.0)
   {
    iconj++;
    gamma = grad/pgrad;
    if (Abs(gamma) > 1000.0)
     {
      iconj = 0;
      gamma = 0.0;
     }
   }
  else
   {
    iconj = 0;
    gamma = 0.0;
   }

  for ( dcntr = 0 ; dcntr < datasz ; dcntr++ )
    dellam[dcntr] = dzdlam[dcntr] + gamma*dellam[dcntr];

}

/*.cp*/
/*
   Procedure MwdCF

     Parameters: datasz......int   I   Size of data set
                 imsz........int   I   Size of image
                 flux........int   I   Total flux of image
                 lambda(*)...real  I   Array of Lagrange multipliers
                 sigma(*)....real  I   Array of expected errors
                 dobs(*).....real  I   Observed data
                 mu..........real  I   Mu
                 omega.......real  I   Omega
                 mem(*)......real  O   Current MEM image
                 def(*)......real  I   a prior image
                 z...........real  O   Potential
                 LTrans......ext   I   Name of routine for lambda->F

     Calculate F from transformed Lagrange multipliers.
*/
void MwdCF (datasz,imsz,flux,lambda,sigma,dobs,mu,omega,mem,def,z,LTrans)
int datasz,imsz;
double flux,*lambda,*sigma,*dobs,mu,omega,*mem,*def,*z;
void (*LTrans) ();
{
  static double maxdr = 1e20;
  int dcntr,imcntr;
  double lamtrm,norm,ls;
  double log (),exp ();


  lamtrm = log(maxdr);

  (*LTrans) (datasz,imsz,lambda,mem);

  norm = 0.0;
  for ( imcntr = 0 ; imcntr < imsz ; imcntr++ )
   {
    mem[imcntr] = def[imcntr]*exp(Min(-mem[imcntr],lamtrm));
    norm += mem[imcntr];
   }

  if (flux > 0.0)
   {
    for ( imcntr = 0 ; imcntr < imsz ; imcntr++ )
      mem[imcntr] = flux*mem[imcntr]/norm;

    *z = (-flux)*log(norm);
   }
  else
    *z = flux*norm;

  for ( dcntr = 0 ; dcntr < datasz ; dcntr++ )
    *z -= lambda[dcntr]*dobs[dcntr];

  if (mu > 0.0)
   {
    *z -= mu*omega;

    for ( dcntr = 0 ; dcntr < datasz ; dcntr++ )
     {
      ls = lambda[dcntr] * sigma[dcntr];
      *z -= (ls*ls)/(4.0*mu);
     }
   }

}

/*.cp*/
/*
   Function MwdCGDst

     Parameters:  datasz......int   I   Size of data set
                  lambda(*)...real  I   Array of Lagrange multipliers
                  dellam(*)...real  I   Actual step in lambda
                  mu..........real  I   Mu
                  sigma(*)....real  I   Array of expected errors
                  dobs(*).....real  I   Array of observed data
                  dpred(*)....real  I   Array of predicted data

     Calculate the dot product of the current step with the new gradient.
*/
double MwdCGDst (datasz,lambda,dellam,mu,sigma,dobs,dpred)
int datasz;
double *lambda,*dellam,mu,*sigma,*dobs,*dpred;
{
  int dcntr;
  double t,rc;


  rc = 0.0;
  for ( dcntr = 0 ; dcntr < datasz ; dcntr++ )
    rc += dellam[dcntr]*(dpred[dcntr] - dobs[dcntr]);

  if (mu > 0.0)
   {
    for ( dcntr = 0 ; dcntr < datasz ; dcntr++ )
     {
      t = sigma[dcntr];
      rc -= (dellam[dcntr]*lambda[dcntr]*t*t)/(2.0*mu);
     }
   }

  return (rc);

}

/*.cp*/
/*
   Procedure MwdCGrad

     Parameters:  datasz.......int   I   Size of data set
                  lambda(*)....real  I   Array of Lagrange multipliers
                  dzdlam(*)....real  I   Gradients of Z w.r.t. lambda
                  pdzdlam(*)...real  IO  Previous gradients of Z w.r.t. lambda
                  mu...........real  I   Mu
                  sigma(*).....real  I   Array of expected errors
                  dobs(*)......real  I   Array of observed data
                  dpred(*).....real  I   Array of predicted data
                  grad.........real  I0  Current grad  . Current grad
                  pgrad........real  I0  Previous grad . Previous grad

     Calculate the gradients of Z with respect to the unknowns.
*/
void MwdCGrad (datasz,lambda,dzdlam,pdzdla,mu,sigma,dobs,dpred,grad,pgrad)
int datasz;
double *lambda,*dzdlam,*pdzdla,mu,*sigma,*dobs,*dpred,*grad,*pgrad;
{
  int dcntr;
  double t;


  *pgrad = *grad;
  for ( dcntr = 0 ; dcntr < datasz ; dcntr++ )
   {
    pdzdla[dcntr] = dzdlam[dcntr];
    dzdlam[dcntr] = dpred[dcntr] - dobs[dcntr];
   }

  *grad = 0.0;
  if (mu > 0.0)
   {
    for ( dcntr = 0 ; dcntr < datasz ; dcntr++ )
     {
      t = sigma[dcntr];
      dzdlam[dcntr] -= lambda[dcntr]*(t*t)/(2.0*mu);
     }
   }

  for ( dcntr = 0 ; dcntr < datasz ; dcntr++ )
   {
    t = dzdlam[dcntr];
    *grad += t*t;
   }

}

/*.cp*/
/*
   Function MwdCLen

     Parameters:  datasz......int   I   Size of data set
                  dellam(*)...real  I   Array of changes in Lagrange multipliers

     Calculate length of step, which is returned as the function value.
*/
double MwdCLen (datasz,dellam)
int datasz;
double *dellam;
{
  int dcntr;
  double del,sum;
  double sqrt();


  sum = 0.0;
  for ( dcntr = 0 ; dcntr < datasz ; dcntr++ )
   {
    del = dellam[dcntr];
    sum += del*del;
   }

  return (sqrt(sum/datasz));

}

/*.cp*/
/*
   Function MwdCMu

     Parameters:  datasz......int  I   Size of data set
                  omega.......real I   Desired omega
                  sigma(*)....real I   Array of expected errors
                  lambda(*)...real I   Array of Lagrange multipliers

     Calculate optimum value of mu for a given omega.
*/
double MwdCMu (datasz, omega, sigma, lambda)
int datasz;
double omega,*sigma,*lambda;
{
  int dcntr;
  double dzdmu,t;
  double sqrt();


  dzdmu = 0.0;
  for ( dcntr = 0; dcntr < datasz; dcntr++ )
   {
    t = lambda[dcntr]*sigma[dcntr];
    dzdmu += t*t;
   }

  return (sqrt(dzdmu/(4.0*omega)));

}

/*.cp*/
/*
   Procedure MwdTakSt

     Parameters:  datasz...int  I  Size of data set
                  in1(*)...real I  Array 1
                  in2(*)...real I  Array 2
                  out(*)...real O  Output array
                  tstep....real I  Step

     Take the recommended step.
*/
void MwdTakSt (datasz,in1,in2,out,tstep)
int datasz;
double *in1,*in2,*out,tstep;
{
  int dcntr;


  for ( dcntr = 0 ; dcntr < datasz ; dcntr++ )
    out[dcntr] = in1[dcntr] + tstep*in2[dcntr];

}

/*.pa Procedure MwdOneIt -- Main driver procedure for one iteration*/
/*
   Procedure MwdOneIt

     Parameters:  doinit......log   I   Initialise ? Should be true for first call
                  datasz......int   I   Size of data set
                  imsz........int   I   Size of image
                  nconj.......int   I   Number of successive conjugate steps
                  flux........real  I   Total flux (zero if not constrained)
                  mem(*)......real  IO  Current MEM image
                  def(*)......real  I   Default image
                  sigma(*)....real  I   Array of expected errors
                  dobs(*).....real  I   Array of observed data
                  dpred(*)....real  I   Array of predicted data
                  iomega......real  I   Expected Fit
                  dzdlam(*),
                  dellam(*),
                  pdzdla(*)...real  I   Work space area for gradient arrays
                                        (must be preserved between calls)
                  iomega......real  I   ?
                  tol.........real  I   Tolerance
                  length......real  O   Length of step
                  lambda(*)...real  O   Lagrange multipliers of data
                  mu..........real  O   Lagrange multiplier for fit omega
                  z...........real  O   Potential function
                  FTrans......proc  I   Name of routine for F->dpred
                  LTrans......proc  I   Name of routine for lambda->F
*/
void MwdOneIt(doinit,datasz,imsz,nconj,flux,mem,def,
	      sigma,dobs,dpred,dzdlam,dellam,pdzdla,iomega,tol,length,
	      lambda,mu,z,FTrans,LTrans)
int doinit,datasz,imsz,nconj;
double flux,*mem,*def,*sigma,*dobs,*dpred,*dzdlam,*dellam,
       *pdzdla,iomega,tol,*length,*lambda,*mu,*z;
void (*FTrans) (),(*LTrans) ();
{
  int dcntr;
  static double grad,pgrad,scale;
  double step,stepm,stepp,omega,gds0,gds1,t;


/*
   Do one iteration.
*/
  if (doinit)
    printf ("T");
  else
    printf ("F");
  printf ("%d %d %d %10e %10e %10e %10e %10e %10e %10e %10e\n",
          datasz,imsz,nconj,flux,mem[0],def[0],sigma[0],dobs[0],dpred[0],iomega,tol);

/*
   Initialise.  On the first step we just try for initial values
   of the lambda's so we set omega to zero just for this iteration.
*/
  if (doinit)
   {
    omega = 0.0;
    *mu   = 0.0;

    for ( dcntr = 0 ; dcntr < datasz ; dcntr++ )
     {
      lambda[dcntr] = 0.0;
      dzdlam[dcntr] = 0.0;
      dellam[dcntr] = 0.0;
      pdzdla[dcntr] = 0.0;
     }

    MwdCF (datasz,imsz,flux,lambda,sigma,dobs,*mu,omega,mem,def,z,LTrans);
    (*FTrans) (imsz,datasz,mem,dpred);

    MwdCGrad (datasz,lambda,dzdlam,pdzdla,*mu,sigma,dobs,dpred,&grad,&pgrad);
   }
  else
    omega = iomega;


/*
   Calculate delta lambda and the length of lambda.
*/
  MwdCDel (doinit,datasz,nconj,lambda,dzdlam,dellam,grad,pgrad);
  *length = MwdCLen (datasz,dellam);


/*
   Now set scaling according to initial length.
*/
  if (doinit)
    scale = tol/(*length);
  *length = scale*(*length);

/*
   Find (grad . step) at start.
*/
  gds0 = MwdCGDst (datasz,lambda,dellam,*mu,sigma,dobs,dpred);


/*
   Ensure that we don't step too far.
*/
  if (*length > 0.0)
   {
    t = tol/(*length);			/* step = Min(1,tol/length)		*/
    step = Min(1.0,t);

    t = Max(1.0,t);			/* stepm = Min(10,Max(1,tol/length))	*/
    stepm = Min(10.0,t);
   }
  else
   {
    step = 1.0;
    stepm = 10.0;
   }

  stepp = step;

/*
   Take one step.
*/
  MwdTakSt (datasz,lambda,dellam,lambda,scale*step);


/*
   Calculate (Gradient.Step) at the new position.
*/
  MwdCF (datasz,imsz,flux,lambda,sigma,dobs,*mu,omega,mem,def,z,LTrans);
  (*FTrans) (imsz,datasz,mem,dpred);
  gds1 = MwdCGDst (datasz,lambda,dellam,*mu,sigma,dobs,dpred);


/*
   Calculate the optimum step using the usual formula.
*/
  if (gds0 != gds1)
   {
    t = step*gds0/(gds0 - gds1);	/* step = Max(0,Min(step*gds0/(gds0-gds1),sterm))	*/
    t = Min(t,stepm);
    step = Max(0.0,t);

    if (step <= 0.0)
     {
      printf ("MwdOneIt Error - Incorrect curvature\n");
      step = 1.0;
     }
   }
  else
   {
    printf ("MwdOneIt Error - no progress\n");
    step = 1.0;
   }


/*
   Step to optimum place.
*/
  MwdTakSt (datasz,lambda,dellam,lambda,scale*(step - stepp));


/*
   Now update mu.
*/
  if (omega > 0.0)
    *mu = MwdCMu (datasz,omega,sigma,lambda);


/*
   Calculate the current image resulting from changes in lambda.
*/
  MwdCF (datasz,imsz,flux,lambda,sigma,dobs,*mu,omega,mem,def,z,LTrans);
  (*FTrans) (imsz,datasz,mem,dpred);

  MwdCGrad (datasz,lambda,dzdlam,pdzdla,*mu,sigma,dobs,dpred,&grad,&pgrad);


/*
   Slowly update the scaling factor.
*/
  scale *= (step + 3.0)/4.0;

}

