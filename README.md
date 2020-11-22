# FittingWithOrthoPoly

  Fitting points with Jacobi orthogonal polynomials 
  -------------------------------------------------
(https://en.wikipedia.org/wiki/Jacobi_polynomials)

  Compiler
  --------
  VS 2019 console, can be easily converted into Linux : just std:: and C++ 11/14 are used.

  Parameters
  ----------
  Input : pairs of x-y points, output : Jacobi poly. The code builds ortho poly 
approximation on points f(x). x must be monotonically increasing.

  Parameter U is within [-1.0 .. +1.0] (which maps to [xmin..xmax]) everywhere.
  
  Parameters alpha and beta define a type of Jacobi polynomial; the case 
alpha = beta = 0.0 corresponds to a Legendre poly. This case is common for
approximation of most curves. But if you need to fit a curve y(x) with infinite 
derivatives at the ends, like a half circle, or an aerofoil surface, 
set alpha = beta = 0.5.

  <B>The trick is</B> that actually f(x) / ((1 - x) ^ alpha) * (1 + x) ^ beta), not 
f(x), is approximated.

  Integration
  -----------
  As the poly is orthogonal, it means that approximation coefs are calculated by 
integration. There are two ways to integrate : by numerical Gauss (set fit() 
parameter as GAUSSINT_8 for example) or by trapezoid rule (set OTHER_INTEGRATION). 

Keep in mind that Gauss integration is built upon a single poly itself and it 
cannot well integrate curves of multiple bends. Use trapezoid rule for this, but 
it requires many points.

  Poly degree
  -----------
  Keep in mind that a single poly is built upon the whole region and the poly 
degree must not be very high, say, up to 20; otherwise gamma function / factorial
calculation problems will appear (no checks).

  Tests
  -----
Test 1 : horizontal curve, Legendre poly (0.0,0.0) and Gauss integration

Accuracy 0.000000

Test 2 : horizontal curve, Legendre poly (0.0,0.0) and trapezoid integration, few points

Accuracy 0.000000

Test 3 : horizontal curve, Legendre poly (0.0,0.0) and trapezoid integration, many points

Accuracy 0.000000

Test 4 : half of circle, infinite derivative at ends, Legendre poly (0.0,0.0) and Gauss integration

Accuracy 0.384346

Test 5 : half of circle, infinite derivative at ends, use Jacoby poly (0.5,0.5) and Gauss integration

Accuracy 0.006320

Test 6 : fitting cos curve from 0 to 4 Pi. Legendre poly. Trapezoidal integration on many points

Accuracy 0.052173

Test 7 : fitting cos curve from 0 to 4 Pi. Legendre poly. Gauss integration on many points

Accuracy 0.052173

