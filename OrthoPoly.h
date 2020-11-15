// Copyright 2020 Andrey Kudryavtsev (andrewkoudr@hotmail.com)
//
// Permission to use, copy, modify, and distribute this software and its
// documentation for any purpose and without fee is hereby granted, provided
// that the above copyright notice appears in all copies and that both the
// copyright notice and this permission notice appear in supporting
// documentation, and that the same name not be used in advertising or
// publicity pertaining to distribution of the software without specific,
// written prior permission. 
// We make no representations about the suitability this software for any 
// purpose. It is provided "as is" without express or implied warranty.

#pragma once

#include <vector>
#include "defines.h"
#include "Types.h"
#include "GaussInt.h"


#include <math.h>

using namespace std;

/**
  Fitting points with Jacobi orthogonal polynomials 
  -------------------------------------------------
(https://en.wikipedia.org/wiki/Jacobi_polynomials)

  Compiler
  --------
  VS 2019, can be easily converted into Linux : just std:: and C++ 11/14.

  Parameters
  ----------
  Input : pairs of x-y points, Output : Jacobi poly. The code builds ortho poly 
approximation on points f(x). x must be monotonically increasing.
  Parameter U is within [-1.0 .. +1.0] everywhere.
  Parameters alpha and beta define type of Jacobi polynomial; the case 
alpha = beta = 0.0 corresponds to a Legendre poly. This case is common for
approximation of most curves. But if you need to fit a curve y(x) with infinite 
derivatives at the ends, like a half circle, or an aerofoil surface, 
set alpha = beta = 0.5.
  The trick is that actually f(x) / ((1 - x) ^ alpha) * (1 + x) ^ beta), not 
f(x), is approximated.

  Integration
  -----------
  As the poly is orthogonal, it means that approxomation coefs are calculated by 
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
*/

template <class T> class JacobiPoly {
public:
                            // "alpha" and "beta"
  T alpha = T(0.0);
  T beta = T(0.0);
                            // constructors
  JacobiPoly() = default;
  JacobiPoly(const T palpha, const T pbeta);
                            // fit arbitrary function y(x), use GAUSSINT_... 
                            // for Gauss integration or OTHER_INTEGRATION for 
                            // trapezoid rule (for many points and complicated
                            // curve shape)
  bool fit(const int degree, const std::vector<T> &x, const std::vector<T> &y, 
    int integration);
                            // get value on U [-1..+1]
  T getValue(const T U);
                            // get accuracy of fitting as max y deflection from
                            // original data
  T accuracy(const std::vector<T> &x, const std::vector<T> &y);

private:
                            // coefficients; real degree is one less than the number of coefs
  std::vector<T> coefs;
                            // y min max before scaling
  T ymin = T(0.0);
  T ymax = T(0.0);
                            // get value, degree is one less than the number of coefs,
                            // a - Jacobi poly alpha, b - beta, U [-1..+1]
  static T getPolyValuePrim(const int degree, const T a, const T b, const T U);
                            // get value, degree is one less than the number of coefs,
                            // U [-1..+1]
  T getPolyValue(const int degree, const T U);
                            // get poly values recurrently
  static void getPolyValues(const int n, const T a, const T b, const T U, std::vector<T> &values);
                            // get orthogonality coefs 
  static void getOrthogonalityCoefs(const int n, const T a, const T b, std::vector<T> &values);
};

// Find value inside montonic table
template <class T> static int FindInterval(const std::vector<T> table, const T value)
{
                              // size, signed integer to avoid problems with size() - 2 (last interval)
  int size = static_cast<int>(table.size());
                              // increasing?
  bool increasing = ((table[size - 1] - table[0]) >= 0.0);
                              // table must be monotonic
  #ifdef _DEBUG
    if (increasing)
    {
      for (int i = 0; i < size - 1; i++)
      {
        if (table[i + 1] < table[i])
          assert(false && "Table not monotonic");
      }
    } else
    {
      for (int i = 0; i < size - 1; i++)
      {
        if (table[i + 1] > table[i])
          assert(false && "Table not monotonic");
      }
    }
  #endif
                              // lower and upper, originally out of scope
  int lower = -1;
  int upper = size - 1;
                              // min/max
  T min = increasing ? table[0] : table[size - 1];
  T max = increasing ? table[size - 1] : table[0];
                              // tolerance
  T tolerance = TOLERANCE;
                              // obvious outcome
  if (value < min - tolerance)
    return -1;
  if (value > max + tolerance)
    return -1;

  while ((upper - lower) > 1)
  {
                              // middle point
    int middle = (lower + upper) >> 1;
    if ((value >= table[middle]) == increasing)
    {
      lower = middle;
    } else
    {
      upper = middle;
    }
  }
                              // we are going to return lower, correct it just in case
  if (lower < 0) lower = 0;
  if (lower > size - 2) lower = size - 2;
                              // done!
  return lower;
}

// Scale array of points to the range smin,smax and find its actual min/max 
template <class T> void rescale(std::vector<T> &points, const T smin, const T smax, T &min, T &max)
{
  assert(points.size() > 0);
  assert(smax > smin);
                            // get min/max
  min = max = points[0];
  for (size_t i = 1; i < points.size(); i++)
  {
    if (points[i] < min)
      min = points[i];
    if (points[i] > max)
      max = points[i];
  }
                            // scale to smin,smax
  T d = max - min;
  T sd = smax - smin;

  if (std::abs(d) > TOLERANCE)
  {
    for (size_t i = 0; i < points.size(); i++)
    {
      points[i] = smin + sd * (points[i] - min) / d;
    }
  } else
  {
    std::fill(points.begin(),points.end(),T(0.0));
  }
}

// Get y value for integration, xx must be from xmin to xmax
template <class T> static T getY(const std::vector<T> &x, const std::vector<T> &y, const T xx, T tolerance)
{
                            // find segment by bisection
  int i = FindInterval<T>(x,xx);
  assert(i >= 0);
                            // find vlaue by linear interpolation within segment
  T dx = x[i + 1] - x[i];
  if (std::abs(dx) < tolerance)
  {
    return y[i];
  } else
  {
    return y[i] + (y[i + 1] - y[i]) * (xx - x[i]) / dx;
  }
}

// Factorial
template <class T> static T factorial(int n)
{
  return tgamma(n + 1);
}

// Gamma-function   
template <class T> static T gamma(T x)
{
  return tgamma(x);
}

// Binomial coefficient (z / n)
template <class T> static T binomialCoef(T z, T n)
{
  assert(n >= 0 && z >= n);

  T res = gamma<T>(z + T(1.0)) / (gamma<T>(z - n + T(1.0)) * gamma<T>(n + T(1.0)));

  return res;
}

template <class T> JacobiPoly<T>::JacobiPoly(const T palpha, const T pbeta) : alpha(palpha), beta(pbeta)
{
}

template <class T> T JacobiPoly<T>::getValue(const T U)
{
  T sum = T(0.0);

  for (int i = 0; i < static_cast<int>(coefs.size()); i++)
  {
    sum += coefs[i] * getPolyValue(i,U);
  }
                              // approximate value is without (1 - x)^a * (1 + x)^b
  sum *= pow(T(1.0) - U,alpha) * pow(T(1.0) + U,beta);
                              // y scaled to 0..1
  sum = ymin + (ymax - ymin) * sum;

  return sum;
}

template <class T> T JacobiPoly<T>::getPolyValuePrim(int n, T a, T b, T U)
{
  assert(n >= 0);
  assert(U >= T(-1.0) && U <= T(1.0));

  T sum = T(0.0);

  for (int s = 0; s <= n; s++)
  {
    sum += binomialCoef<T>(n + a, n - s) * binomialCoef<T>(n + b, s) * 
      pow((U - T(1.0)) * T(0.5),s) * pow((U + T(1.0)) * T(0.5),n - s);
  }

  return sum;
}

template <class T> void JacobiPoly<T>::getOrthogonalityCoefs(const int degree, const T a, const T b, 
  std::vector<T> &values)
{
  assert(degree >= 0);
  values.clear();

  T ab1 = a + b + T(1.0);
  for (int n = 0; n <= degree; n++)
  {
    T fn = factorial<T>(n);
    T na1 = n + a + 1;
    T nb1 = n + b + 1;
    T n2ab1 = n * 2 + ab1;
    T nab1 = n + ab1;
    T v = pow(T(2.0),ab1) * gamma<T>(na1) * gamma<T>(nb1) / (n2ab1 * fn * gamma<T>(nab1));
    values.push_back(v);
  }
}

template <class T> void JacobiPoly<T>::getPolyValues(const int n, const T a, const T b, 
  const T U, std::vector<T> &values)
{
  assert(n >= 0);
  assert(U >= -1.0 && U <= 1.0);
  values.clear();

  T ab = a + b;
  int njm = n - 1;
  values.push_back(T(1.0));
  values.push_back((a - b + (ab + 2) * U) * 0.5);

  for (int m = 2; m <= n; m++)
  {
    T nab = m + ab;
    T nab2 = m * 2 + ab;
    T nab21 = nab2 - 1;
    T nab22 = nab2 - 2;
    T na1 = m + a - 1;
    T nb1 = m + b - 1;
    T bracket = nab2 * nab22 * U + a * a - b * b;
    T c0 = 2 * m * nab *  nab22;
    T c1 = nab21 * bracket;
    T c2 = 2 * na1 * nb1 * nab2;
    T v = (c1 * values[m - 1] - c2 * values[m - 2]) / c0;

    values.push_back(v);
  }
}

template <class T> T JacobiPoly<T>::getPolyValue(const int degree, const T U)
{
  assert(degree >= 0);
  assert(U >= -1.0 && U <= 1.0);

  return getPolyValuePrim(degree,alpha,beta,U);
}

template <class T> bool JacobiPoly<T>::fit(const int degree, const std::vector<T> &x, const std::vector<T> &y, 
  int integration)
{
  assert(x.size() > 1);
  assert(x.size() == y.size());

  if (x.size() > 1)
  {
    std::vector<T> yscaled = y;
    rescale(yscaled,T(0.0),T(1.0),ymin,ymax);

    coefs.clear();
    coefs.resize(degree + 1,T(0.0));

    T tolerance = TOLERANCE;

    T xmin = x[0];
    T xmax = x.back();
    assert(xmax > xmin);
    T dx = xmax - xmin;
                              // orthogonality coefs
    std::vector<T> ocoefs;
    getOrthogonalityCoefs(degree,alpha,beta,ocoefs);

                              // integrate by Gauss
    if (integration == GAUSSINT_1 ||
      integration == GAUSSINT_2 ||
      integration == GAUSSINT_4 ||
      integration == GAUSSINT_8)
    {
                          // if y-x function is complex,
                          // do not use Gauss integration here, as
                          // only a few points (max 8 in out case) are 
                          // taken into account
      for (int k = 0; k < GaussInt[integration].numpoints; k++)
      {
        T e = static_cast<T>(GaussInt[integration].knots[k]);
        T w = static_cast<T>(GaussInt[integration].weights[k]);

        T xx = xmin + dx * (e + T(1.0)) * T(0.5);
        T ymiddle = getY(x,yscaled,xx,tolerance);

        std::vector<T> values;
        getPolyValues(degree,alpha,beta,e,values);

        for (int j = 0; j <= degree; j++)
        {
          coefs[j] += ymiddle * values[j] * w / ocoefs[j];
        }
      }
int gsggsgs = 0; //!!!!!!!!!!!!!
    } else
    {
                          // compute coefficients by trapezoidal rule
      for (int j = 1; j < static_cast<int>(x.size()); j++)
      {
        assert(x[j] >= x[j - 1]);

        T xx = (x[j] + x[j - 1]) * T(0.5);
        T u = (xx - xmin) / dx;
        T u1 = (x[j - 1] - xmin) / dx;
        T u2 = (x[j] - xmin) / dx;
        T e = u + u - 1;

        std::vector<T> values;
        getPolyValues(degree,alpha,beta,e,values);

        for (int i = 0; i <= degree; i++)
        {
          coefs[i] += (yscaled[j] + yscaled[j - 1]) * T(0.5) * values[i] * (u2 - u1) * 
            T(2.0) /* -1..+1*/ / ocoefs[i];
        }
      }
    }

    return true;
  } else
  {
    coefs.clear();
    return false;
  }
}

template <class T> T JacobiPoly<T>::accuracy(const std::vector<T> &x, const std::vector<T> &y)
{
  assert(x.size() > 1);
  assert(x.size() == y.size());
  assert(coefs.size() > 0);

  if (coefs.size() > 1)
  {
    T accuracy = T(0.0);

    T xmin = x[0];
    T xmax = x[x.size() - 1];
    assert(xmax > xmin);
    T dx = xmax - xmin;

#ifdef _DEBUG
    std::vector<T> diffs;
    std::vector<T> polys;
#endif

    for (int j = 0; j < static_cast<int>(x.size()); j++)
    {
      T u = (x[j] - xmin) / dx;
                              // -1..+1
      u = u + u - 1;

      T poly = getValue(u);
      T diff = std::abs(poly - y[j]);

#ifdef _DEBUG
      polys.push_back(poly);
      diffs.push_back(diff);
#endif

      if (diff > accuracy)
        accuracy = diff;
    }

    return accuracy;
  } else
  {
    return T(0.0);
  }
}


