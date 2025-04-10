#pragma once
#include <math.h>

inline bool XMScalarNearEqual
(
    float S1,
    float S2,
    float Epsilon
)
{
    float Delta = S1 - S2;
    return (fabsf(Delta) <= Epsilon);
}

_Use_decl_annotations_
inline void XMScalarSinCos
(
    float* pSin,
    float* pCos,
    float  Value
)
{
    assert(pSin);
    assert(pCos);

    // Map Value to y in [-pi,pi], x = 2*pi*quotient + remainder.
    float quotient = XM_1DIV2PI * Value;
    if (Value >= 0.0f)
    {
        quotient = (float)((int)(quotient + 0.5f));
    }
    else
    {
        quotient = (float)((int)(quotient - 0.5f));
    }
    float y = Value - XM_2PI * quotient;

    // Map y to [-pi/2,pi/2] with sin(y) = sin(Value).
    float sign;
    if (y > XM_PIDIV2)
    {
        y = XM_PI - y;
        sign = -1.0f;
    }
    else if (y < -XM_PIDIV2)
    {
        y = -XM_PI - y;
        sign = -1.0f;
    }
    else
    {
        sign = +1.0f;
    }

    float y2 = y * y;

    // 11-degree minimax approximation
    *pSin = (((((-2.3889859e-08f * y2 + 2.7525562e-06f) * y2 - 0.00019840874f) * y2 + 0.0083333310f) * y2 - 0.16666667f) * y2 + 1.0f) * y;

    // 10-degree minimax approximation
    float p = ((((-2.6051615e-07f * y2 + 2.4760495e-05f) * y2 - 0.0013888378f) * y2 + 0.041666638f) * y2 - 0.5f) * y2 + 1.0f;
    *pCos = sign * p;
}