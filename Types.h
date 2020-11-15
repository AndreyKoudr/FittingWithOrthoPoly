#pragma once

// LINT is always 64-bit
#define LINT long long int

// useful macros
#define LIMIT(x,xmin,xmax) if (x < xmin) x = xmin; if (x > xmax) x = xmax
#define LIMIT_MIN(x,xmin) if (x < xmin) x = xmin
#define LIMIT_MAX(x,xmax) if (x > xmax) x = xmax

// tolerance
#define TOLERANCE std::numeric_limits<T>::epsilon() * static_cast<T>(10.0)

#define ROUND(x)  (int) (floor(x + .5))


//===== PI - associated ========================================================
#ifndef M_PI
  #define M_PI    3.14159265358979323846
#endif

#ifndef PI05
  #define PI05 (M_PI * 0.5)
#endif

#ifndef PI20
  #define PI20 (M_PI * 2.0)
#endif

#ifndef PI10
  #define PI10 M_PI
#endif
