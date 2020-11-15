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

#include "OrthoPoly.h"

int main()
{
  printf("Test 1 : horizontal curve, Legendre poly (0.0,0.0) and Gauss integration\n");
  {
    #define T float
    #define DEGREE 4

    JacobiPoly<T> f(0.0,0.0);

    std::vector<T> x,y;
                             
    T xmin = T(-10.0);
    T xmax = T(+25.0);
    int numpoints = 20;
    T dx = (xmax - xmin) / T(numpoints - 1);
    for (int i = 0; i < numpoints; i++)
    {
      T xx = xmin + dx * T(i);
      x.push_back(xx);
      y.push_back(11.5);
    }
                              // fit
#if 1
    bool res = f.fit(DEGREE,x,y,GAUSSINT_8);
#else
    bool res = f.fit(DEGREE,x,y,OTHER_INTEGRATION);
#endif
    T acc = f.accuracy(x,y);
    printf("Accuracy %f\n",acc);

    assert(acc < 0.00001);
  }

  printf("Test 2 : horizontal curve, Legendre poly (0.0,0.0) and trapezoid integration, few points\n");
  {
    #define T float
    #define DEGREE 4

    JacobiPoly<T> f(0.0,0.0);

    std::vector<T> x,y;
                             
    T xmin = T(-10.0);
    T xmax = T(+25.0);
    int numpoints = 20;
    T dx = (xmax - xmin) / T(numpoints - 1);
    for (int i = 0; i < numpoints; i++)
    {
      T xx = xmin + dx * T(i);
      x.push_back(xx);
      y.push_back(11.5);
    }
                              // fit
#if 0
    bool res = f.fit(DEGREE,x,y,GAUSSINT_8);
#else
    bool res = f.fit(DEGREE,x,y,OTHER_INTEGRATION);
#endif
    T acc = f.accuracy(x,y);
    printf("Accuracy %f\n",acc);

                              // integration not accurate, too few points
    assert(acc < 0.6);
  }

  printf("Test 3 : horizontal curve, Legendre poly (0.0,0.0) and trapezoid integration, many points\n");
  {
    #define T float
    #define DEGREE 4

    JacobiPoly<T> f(0.0,0.0);

    std::vector<T> x,y;
                             
    T xmin = T(-10.0);
    T xmax = T(+25.0);
    int numpoints = 2000;
    T dx = (xmax - xmin) / T(numpoints - 1);
    for (int i = 0; i < numpoints; i++)
    {
      T xx = xmin + dx * T(i);
      x.push_back(xx);
      y.push_back(11.5);
    }
                              // fit
#if 0
    bool res = f.fit(DEGREE,x,y,GAUSSINT_8);
#else
    bool res = f.fit(DEGREE,x,y,OTHER_INTEGRATION);
#endif
    T acc = f.accuracy(x,y);
    printf("Accuracy %f\n",acc);

                              // integration not accurate, too few points
    assert(acc < 0.0005);
  }

  printf("Test 4 : half of circle, infinite derivative at ends, Legendre poly (0.0,0.0) and Gauss integration\n");
  {
    #define T float
    #define DEGREE 4

    JacobiPoly<T> f(0.0,0.0);

    std::vector<T> x,y;
                              // cos from 0 to Pi
    T amin = T(0.0);
    T amax = T(PI10);
    int numpoints = 20;
    T da = (amax - amin) / T(numpoints - 1);
    T a = amin;
    T R = 2.1;
    T offset = 0.3;
    for (int i = 0; i < numpoints; i++)
    {
      a = amin + da * T(i);
      x.push_back(-R * cos(a) + offset);
      y.push_back(R * sin(a) + offset);
    }
                              // fit
#if 1
    bool res = f.fit(DEGREE,x,y,GAUSSINT_8);
#else
    bool res = f.fit(DEGREE,x,y,OTHER_INTEGRATION);
#endif
    T acc = f.accuracy(x,y);
    printf("Accuracy %f\n",acc);

    assert(acc < 0.5);
  }

  printf("Test 5 : half of circle, infinite derivative at ends, use Jacoby poly (0.5,0.5) and Gauss integration\n");
  {
    #define T float
    #define DEGREE 4

    JacobiPoly<T> f(0.5,0.5);

    std::vector<T> x,y;
                              // cos from 0 to Pi
    T amin = T(0.0);
    T amax = T(PI10);
    int numpoints = 20;
    T da = (amax - amin) / T(numpoints - 1);
    T a = amin;
    T R = 2.1;
    T offset = 0.3;
    for (int i = 0; i < numpoints; i++)
    {
      a = amin + da * T(i);
      x.push_back(-R * cos(a) + offset);
      y.push_back(R * sin(a) + offset);
    }
                              // fit
#if 1
    bool res = f.fit(DEGREE,x,y,GAUSSINT_8);
#else
    bool res = f.fit(DEGREE,x,y,OTHER_INTEGRATION);
#endif
    T acc = f.accuracy(x,y);
    printf("Accuracy %f\n",acc);

    assert(acc < 0.007);
  }


  printf("Test 6 : fitting cos curve from 0 to 4 Pi. Legendre poly. Trapezoidal integration on many points\n");
  {
    #define T float
    #define DEGREE 16

    JacobiPoly<T> f(0.0,0.0);

    std::vector<T> x,y;
                              // cos from 0 to Pi
    T amin = T(0.0);
    T amax = T(PI10) * T(4.0);
    int numpoints = 200;
    T da = (amax - amin) / T(numpoints - 1);
    T a = amin;
    T R = 1.0;
    T offset = 0.0;
    for (int i = 0; i < numpoints; i++)
    {
      a = amin + da * T(i);
      x.push_back(a + offset);
      y.push_back(R * sin(a) + offset);
    }
                              // fit
#if 0
    bool res = f.fit(DEGREE,x,y,GAUSSINT_8);
#else
    bool res = f.fit(DEGREE,x,y,OTHER_INTEGRATION);
#endif
    T acc = f.accuracy(x,y);
    printf("Accuracy %f\n",acc);

    assert(acc < 0.06);
  }

  printf("Test 7 : fitting cos curve from 0 to 4 Pi. Legendre poly. Gauss integration on many points\n");
  {
    #define T float
    #define DEGREE 16

    JacobiPoly<T> f(0.0,0.0);

    std::vector<T> x,y;
                              // cos from 0 to Pi
    T amin = T(0.0);
    T amax = T(PI10) * T(4.0);
    int numpoints = 200;
    T da = (amax - amin) / T(numpoints - 1);
    T a = amin;
    T R = 1.0;
    T offset = 0.0;
    for (int i = 0; i < numpoints; i++)
    {
      a = amin + da * T(i);
      x.push_back(a + offset);
      y.push_back(R * sin(a) + offset);
    }
                              // fit
#if 1
    bool res = f.fit(DEGREE,x,y,GAUSSINT_20);
#else
    bool res = f.fit(DEGREE,x,y,OTHER_INTEGRATION);
#endif
    T acc = f.accuracy(x,y);
    printf("Accuracy %f\n",acc);

    assert(acc < 0.06);
  }


  printf("\n");
  printf("Press [ENTER]\n");

  getchar();
}
