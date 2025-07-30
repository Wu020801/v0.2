/***********************************************************************************
 * This software module was originally developed by Tencent America LLC, Media Lab *
 * in the course of development of dynamic mesh compression. Tencent America LLC,  *
 * Media Lab, retains full right to modify and use the code for its own purpose.   *
 * This copyright notice must be included in all copies or derivative works.       *
 * Copyright (c) Tencent 2023.                                                     *
 ***********************************************************************************/
#pragma once

#include <cstdint>
#include "vector.hpp"

//============================================================================

template<class T>
Vector3D<T>
computeTriangleNormal(const Vector3D<T>& a, const Vector3D<T>& b, const Vector3D<T>& c, const bool normalize = true) {
  auto normal = (b - a) ^ (c - a);
  if (normalize) { normal.normalize(); }
  return normal;
}

//============================================================================

template<class T>
T
computeTriangleArea(const Vector3D<T>& a, const Vector3D<T>& b, const Vector3D<T>& c) {
  return T(0.5) * computeTriangleNormal(a, b, c, false).norm();
}

//============================================================================

inline bool
isDegenerate(const Vector3D<uint32_t>& tri) {
  return tri[0] == tri[1] || tri[0] == tri[2] || tri[1] == tri[2];
}

//============================================================================

template<typename T>
void
getMinEdge02(T const& a11, T const& b1, Vector2D<T>& bc2) {
  T const zero = static_cast<T>(0);
  T const one  = static_cast<T>(1);
  bc2[0]       = zero;
  if (b1 >= zero) {
    bc2[1] = zero;
  } else if (a11 + b1 <= zero) {
    bc2[1] = one;
  } else {
    bc2[1] = -b1 / a11;
  }
}

//============================================================================

template<typename T>
void
getMinEdge12(T const& a01, T const& a11, T const& b1, T const& f10, T const& f01, Vector2D<T>& bc2) {
  T const zero = static_cast<T>(0);
  T const one  = static_cast<T>(1);
  T       h0   = a01 + b1 - f10;
  if (h0 >= zero) {
    bc2[1] = zero;
  } else {
    T h1 = a11 + b1 - f01;
    if (h1 <= zero) {
      bc2[1] = one;
    } else {
      bc2[1] = h0 / (h0 - h1);
    }
  }
  bc2[0] = one - bc2[1];
}

//============================================================================

template<typename T>
void
getMinInterior(Vector2D<T> const& p0, T const& h0, Vector2D<T> const& p1, T const& h1, Vector2D<T>& bc2) {
  T z    = h0 / (h0 - h1);
  T omz  = static_cast<T>(1) - z;
  bc2[0] = omz * p0[0] + z * p1[0];
  bc2[1] = omz * p0[1] + z * p1[1];
}

//============================================================================

template<typename T>
Vector3D<T>
closestPointInTriangle2(const Vector3D<T>& p,
                        const Vector3D<T>& a,
                        const Vector3D<T>& b,
                        const Vector3D<T>& c,
                        Vector3D<T>*       bc = nullptr) {
  T const zero  = static_cast<T>(0);
  T const one   = static_cast<T>(1);
  auto    diff  = p - a;
  auto    edge0 = b - a;
  auto    edge1 = c - a;
  T       a00   = edge0 * edge0;
  T       a01   = edge0 * edge1;
  T       a11   = edge1 * edge1;
  T       b0    = -(diff * edge0);
  T       b1    = -(diff * edge1);

  T f00 = b0;
  T f10 = b0 + a00;
  T f01 = b0 + a01;

  Vector2D<T> p0, p1, bc2;
  T           dt1, h0, h1;

  if (f00 >= zero) {
    if (f01 >= zero) {
      getMinEdge02<T>(a11, b1, bc2);
    } else {
      p0[0] = zero;
      p0[1] = f00 / (f00 - f01);
      p1[0] = f01 / (f01 - f10);
      p1[1] = one - p1[0];
      dt1   = p1[1] - p0[1];
      h0    = dt1 * (a11 * p0[1] + b1);
      if (h0 >= zero) {
        getMinEdge02<T>(a11, b1, bc2);
      } else {
        h1 = dt1 * (a01 * p1[0] + a11 * p1[1] + b1);
        if (h1 <= zero) {
          getMinEdge12<T>(a01, a11, b1, f10, f01, bc2);
        } else {
          getMinInterior<T>(p0, h0, p1, h1, bc2);
        }
      }
    }
  } else if (f01 <= zero) {
    if (f10 <= zero) {
      getMinEdge12<T>(a01, a11, b1, f10, f01, bc2);
    } else {
      p0[0] = f00 / (f00 - f10);
      p0[1] = zero;
      p1[0] = f01 / (f01 - f10);
      p1[1] = one - p1[0];
      h0    = p1[1] * (a01 * p0[0] + b1);
      if (h0 >= zero) {
        bc2 = p0;  // getMinEdge01
      } else {
        h1 = p1[1] * (a01 * p1[0] + a11 * p1[1] + b1);
        if (h1 <= zero) {
          getMinEdge12<T>(a01, a11, b1, f10, f01, bc2);
        } else {
          getMinInterior<T>(p0, h0, p1, h1, bc2);
        }
      }
    }
  } else if (f10 <= zero) {
    p0[0] = zero;
    p0[1] = f00 / (f00 - f01);
    p1[0] = f01 / (f01 - f10);
    p1[1] = one - p1[0];
    dt1   = p1[1] - p0[1];
    h0    = dt1 * (a11 * p0[1] + b1);
    if (h0 >= zero) {
      getMinEdge02<T>(a11, b1, bc2);
    } else {
      h1 = dt1 * (a01 * p1[0] + a11 * p1[1] + b1);
      if (h1 <= zero) {
        getMinEdge12<T>(a01, a11, b1, f10, f01, bc2);
      } else {
        getMinInterior<T>(p0, h0, p1, h1, bc2);
      }
    }
  } else {
    p0[0] = f00 / (f00 - f10);
    p0[1] = zero;
    p1[0] = zero;
    p1[1] = f00 / (f00 - f01);
    h0    = p1[1] * (a01 * p0[0] + b1);
    if (h0 >= zero) {
      bc2 = p0;  // GetMinEdge01
    } else {
      h1 = p1[1] * (a11 * p1[1] + b1);
      if (h1 <= zero) {
        getMinEdge02<T>(a11, b1, bc2);
      } else {
        getMinInterior<T>(p0, h0, p1, h1, bc2);
      }
    }
  }

  if (bc) {
    (*bc)[0] = one - bc2[0] - bc2[1];
    (*bc)[1] = bc2[0];
    (*bc)[2] = bc2[1];
  }

  return a + bc2[0] * edge0 + bc2[1] * edge1;
}

//============================================================================

template<typename T>
Vector3D<T>
closestPointInTriangle(const Vector3D<T>& p,
                       const Vector3D<T>& a,
                       const Vector3D<T>& b,
                       const Vector3D<T>& c,
                       Vector3D<T>*       barycentricCoords = nullptr) {
  const auto ab     = b - a;
  const auto ac     = c - a;
  const auto bc     = c - b;
  const auto ap     = p - a;
  const auto bp     = p - b;
  const auto cp     = p - c;
  const auto snom   = ap * ab;
  const auto sdenom = -(bp * ab);
  const auto tnom   = ap * ac;
  const auto tdenom = -(cp * ac);
  const auto eps    = T(1.0e-10);
  const auto zero   = T(0);
  const auto one    = T(1);

  if (snom <= eps && tnom <= eps) {
    if (barycentricCoords) {
      (*barycentricCoords)[0] = one;
      (*barycentricCoords)[1] = zero;
      (*barycentricCoords)[2] = zero;
    }
    return a;
  }

  const auto unom   = bp * bc;
  const auto udenom = -(cp * bc);

  if (sdenom <= eps && unom <= eps) {
    if (barycentricCoords) {
      (*barycentricCoords)[0] = zero;
      (*barycentricCoords)[1] = one;
      (*barycentricCoords)[2] = zero;
    }
    return b;
  }

  if (tdenom <= eps && udenom <= eps) {
    if (barycentricCoords) {
      (*barycentricCoords)[0] = zero;
      (*barycentricCoords)[1] = zero;
      (*barycentricCoords)[2] = one;
    }
    return c;
  }

  const auto n  = ab ^ ac;
  const auto vc = n * (ap ^ bp);

  if (vc <= zero && snom >= eps && sdenom >= eps) {
    const auto s = snom / (snom + sdenom);
    if (barycentricCoords) {
      (*barycentricCoords)[0] = one - s;
      (*barycentricCoords)[1] = s;
      (*barycentricCoords)[2] = zero;
    }

    return a + s * ab;
  }

  const auto va = n * (bp ^ cp);
  if (va <= zero && unom >= eps && udenom >= eps) {
    const auto s = unom / (unom + udenom);
    if (barycentricCoords) {
      (*barycentricCoords)[0] = zero;
      (*barycentricCoords)[1] = one - s;
      (*barycentricCoords)[2] = s;
    }
    return b + s * bc;
  }

  const auto vb = n * (cp ^ ap);
  if (vb <= zero && tnom >= eps && tdenom >= eps) {
    const auto s = tnom / (tnom + tdenom);
    if (barycentricCoords) {
      (*barycentricCoords)[0] = one - s;
      (*barycentricCoords)[1] = zero;
      (*barycentricCoords)[2] = s;
    }
    return a + s * ac;
  }

  const auto sum = va + vb + vc;
  const auto u   = va / sum;
  const auto v   = vb / sum;
  const auto w   = one - u - v;
  if (barycentricCoords) {
    (*barycentricCoords)[0] = u;
    (*barycentricCoords)[1] = v;
    (*barycentricCoords)[2] = w;
  }
  return u * a + v * b + w * c;
}
