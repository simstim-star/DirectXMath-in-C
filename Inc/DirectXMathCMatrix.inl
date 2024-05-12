#pragma once

//------------------------------------------------------------------------------
 // Comparison operations
 //------------------------------------------------------------------------------

#if !defined(_XM_NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(push)
#pragma float_control(precise, on)
#endif

// Return true if any entry in the matrix is NaN
inline bool XM_CALLCONV XMMatrixIsNaN(FXMMATRIX M)
{
#if defined(_XM_NO_INTRINSICS_)
    size_t i = 16;
    const uint32_t* pWork = (const uint32_t*)(&M->m[0][0]);
    do {
        // Fetch value into integer unit
        uint32_t uTest = pWork[0];
        // Remove sign
        uTest &= 0x7FFFFFFFU;
        // NaN is 0x7F800001 through 0x7FFFFFFF inclusive
        uTest -= 0x7F800001U;
        if (uTest < 0x007FFFFFU)
        {
            break;      // NaN found
        }
        ++pWork;        // Next entry
    } while (--i);
    return (i != 0);      // i == 0 if nothing matched
#elif defined(_XM_SSE_INTRINSICS_)
    // Load in registers
    XMVECTOR vX = M->r[0];
    XMVECTOR vY = M->r[1];
    XMVECTOR vZ = M->r[2];
    XMVECTOR vW = M->r[3];
    // Test themselves to check for NaN
    vX = _mm_cmpneq_ps(vX, vX);
    vY = _mm_cmpneq_ps(vY, vY);
    vZ = _mm_cmpneq_ps(vZ, vZ);
    vW = _mm_cmpneq_ps(vW, vW);
    // Or all the results
    vX = _mm_or_ps(vX, vZ);
    vY = _mm_or_ps(vY, vW);
    vX = _mm_or_ps(vX, vY);
    // If any tested true, return true
    return (_mm_movemask_ps(vX) != 0);
#else
#endif
}

// Return true if any entry in the matrix is +/-INF
inline bool XM_CALLCONV XMMatrixIsInfinite(FXMMATRIX M)
{
#if defined(_XM_NO_INTRINSICS_)
    size_t i = 16;
    auto pWork = reinterpret_cast<const uint32_t*>(&M->m[0][0]);
    do {
        // Fetch value into integer unit
        uint32_t uTest = pWork[0];
        // Remove sign
        uTest &= 0x7FFFFFFFU;
        // INF is 0x7F800000
        if (uTest == 0x7F800000U)
        {
            break;      // INF found
        }
        ++pWork;        // Next entry
    } while (--i);
    return (i != 0);      // i == 0 if nothing matched
#elif defined(_XM_SSE_INTRINSICS_)
    // Mask off the sign bits
    XMVECTOR vTemp1 = _mm_and_ps(M->r[0], g_XMAbsMask.v);
    XMVECTOR vTemp2 = _mm_and_ps(M->r[1], g_XMAbsMask.v);
    XMVECTOR vTemp3 = _mm_and_ps(M->r[2], g_XMAbsMask.v);
    XMVECTOR vTemp4 = _mm_and_ps(M->r[3], g_XMAbsMask.v);
    // Compare to infinity
    vTemp1 = _mm_cmpeq_ps(vTemp1, g_XMInfinity.v);
    vTemp2 = _mm_cmpeq_ps(vTemp2, g_XMInfinity.v);
    vTemp3 = _mm_cmpeq_ps(vTemp3, g_XMInfinity.v);
    vTemp4 = _mm_cmpeq_ps(vTemp4, g_XMInfinity.v);
    // Or the answers together
    vTemp1 = _mm_or_ps(vTemp1, vTemp2);
    vTemp3 = _mm_or_ps(vTemp3, vTemp4);
    vTemp1 = _mm_or_ps(vTemp1, vTemp3);
    // If any are infinity, the signs are true.
    return (_mm_movemask_ps(vTemp1) != 0);
#endif
}

// Return true if the XMMatrix is equal to identity
inline bool XM_CALLCONV XMMatrixIsIdentity(FXMMATRIX M)
{
#if defined(_XM_NO_INTRINSICS_)
    // Use the integer pipeline to reduce branching to a minimum
    const uint32_t* pWork = (const uint32_t*)(&M->m[0][0]);
    // Convert 1.0f to zero and or them together
    uint32_t uOne = pWork[0] ^ 0x3F800000U;
    // Or all the 0.0f entries together
    uint32_t uZero = pWork[1];
    uZero |= pWork[2];
    uZero |= pWork[3];
    // 2nd row
    uZero |= pWork[4];
    uOne |= pWork[5] ^ 0x3F800000U;
    uZero |= pWork[6];
    uZero |= pWork[7];
    // 3rd row
    uZero |= pWork[8];
    uZero |= pWork[9];
    uOne |= pWork[10] ^ 0x3F800000U;
    uZero |= pWork[11];
    // 4th row
    uZero |= pWork[12];
    uZero |= pWork[13];
    uZero |= pWork[14];
    uOne |= pWork[15] ^ 0x3F800000U;
    // If all zero entries are zero, the uZero==0
    uZero &= 0x7FFFFFFF;    // Allow -0.0f
    // If all 1.0f entries are 1.0f, then uOne==0
    uOne |= uZero;
    return (uOne == 0);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp1 = _mm_cmpeq_ps(M->r[0], g_XMIdentityR0.v);
    XMVECTOR vTemp2 = _mm_cmpeq_ps(M->r[1], g_XMIdentityR1.v);
    XMVECTOR vTemp3 = _mm_cmpeq_ps(M->r[2], g_XMIdentityR2.v);
    XMVECTOR vTemp4 = _mm_cmpeq_ps(M->r[3], g_XMIdentityR3.v);
    vTemp1 = _mm_and_ps(vTemp1, vTemp2);
    vTemp3 = _mm_and_ps(vTemp3, vTemp4);
    vTemp1 = _mm_and_ps(vTemp1, vTemp3);
    return (_mm_movemask_ps(vTemp1) == 0x0f);
#endif
}

//------------------------------------------------------------------------------
// Computation operations
//------------------------------------------------------------------------------

// Perform a 4x4 matrix multiply by a 4x4 matrix
inline XMMATRIX XM_CALLCONV XMMatrixMultiply
(
    FXMMATRIX A,
    CXMMATRIX B
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMMATRIX mResult;
    // Cache the invariants in registers
    float x = A->m[0][0];
    float y = A->m[0][1];
    float z = A->m[0][2];
    float w = A->m[0][3];
    // Perform the operation on the first row
    mResult.m[0][0] = (B->m[0][0] * x) + (B->m[1][0] * y) + (B->m[2][0] * z) + (B->m[3][0] * w);
    mResult.m[0][1] = (B->m[0][1] * x) + (B->m[1][1] * y) + (B->m[2][1] * z) + (B->m[3][1] * w);
    mResult.m[0][2] = (B->m[0][2] * x) + (B->m[1][2] * y) + (B->m[2][2] * z) + (B->m[3][2] * w);
    mResult.m[0][3] = (B->m[0][3] * x) + (B->m[1][3] * y) + (B->m[2][3] * z) + (B->m[3][3] * w);
    // Repeat for all the other rows
    x = A->m[1][0];
    y = A->m[1][1];
    z = A->m[1][2];
    w = A->m[1][3];
    mResult.m[1][0] = (B->m[0][0] * x) + (B->m[1][0] * y) + (B->m[2][0] * z) + (B->m[3][0] * w);
    mResult.m[1][1] = (B->m[0][1] * x) + (B->m[1][1] * y) + (B->m[2][1] * z) + (B->m[3][1] * w);
    mResult.m[1][2] = (B->m[0][2] * x) + (B->m[1][2] * y) + (B->m[2][2] * z) + (B->m[3][2] * w);
    mResult.m[1][3] = (B->m[0][3] * x) + (B->m[1][3] * y) + (B->m[2][3] * z) + (B->m[3][3] * w);
    x = A->m[2][0];
    y = A->m[2][1];
    z = A->m[2][2];
    w = A->m[2][3];
    mResult.m[2][0] = (B->m[0][0] * x) + (B->m[1][0] * y) + (B->m[2][0] * z) + (B->m[3][0] * w);
    mResult.m[2][1] = (B->m[0][1] * x) + (B->m[1][1] * y) + (B->m[2][1] * z) + (B->m[3][1] * w);
    mResult.m[2][2] = (B->m[0][2] * x) + (B->m[1][2] * y) + (B->m[2][2] * z) + (B->m[3][2] * w);
    mResult.m[2][3] = (B->m[0][3] * x) + (B->m[1][3] * y) + (B->m[2][3] * z) + (B->m[3][3] * w);
    x = A->m[3][0];
    y = A->m[3][1];
    z = A->m[3][2];
    w = A->m[3][3];
    mResult.m[3][0] = (B->m[0][0] * x) + (B->m[1][0] * y) + (B->m[2][0] * z) + (B->m[3][0] * w);
    mResult.m[3][1] = (B->m[0][1] * x) + (B->m[1][1] * y) + (B->m[2][1] * z) + (B->m[3][1] * w);
    mResult.m[3][2] = (B->m[0][2] * x) + (B->m[1][2] * y) + (B->m[2][2] * z) + (B->m[3][2] * w);
    mResult.m[3][3] = (B->m[0][3] * x) + (B->m[1][3] * y) + (B->m[2][3] * z) + (B->m[3][3] * w);
    return mResult;
#elif defined(_XM_SSE_INTRINSICS_)
    XMMATRIX mResult;
    // If the first row of A is [X,Y,Z,W]
    // Splat it into registers
    // Use vW to hold the original row
    XMVECTOR vW = A->r[0];
    XMVECTOR vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0)); // X X X X
    XMVECTOR vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1)); // Y Y Y Y
    XMVECTOR vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2)); // Z Z Z Z
    vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3)); // W W W W
    // Perform the operation on the first row
    vX = _mm_mul_ps(vX, B->r[0]); // X*B00 X*B01 X*B02 X*B03
    vY = _mm_mul_ps(vY, B->r[1]); // Y*B10 Y*B11 Y*B12 Y*B13
    vZ = _mm_mul_ps(vZ, B->r[2]); // Z*B20 Z*B21 Z*B22 Z*B23
    vW = _mm_mul_ps(vW, B->r[3]); // W*B30 W*B31 W*B32 W*B33
    // Sum the partials above
    vX = _mm_add_ps(vX, vZ);
    vY = _mm_add_ps(vY, vW);
    vX = _mm_add_ps(vX, vY);
    mResult.r[0] = vX;
    // Repeat for the other 3 rows
    vW = A->r[1];
    vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
    vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
    vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
    vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
    vX = _mm_mul_ps(vX, B->r[0]);
    vY = _mm_mul_ps(vY, B->r[1]);
    vZ = _mm_mul_ps(vZ, B->r[2]);
    vW = _mm_mul_ps(vW, B->r[3]);
    vX = _mm_add_ps(vX, vZ);
    vY = _mm_add_ps(vY, vW);
    vX = _mm_add_ps(vX, vY);
    mResult.r[1] = vX;
    vW = A->r[2];
    vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
    vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
    vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
    vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
    vX = _mm_mul_ps(vX, B->r[0]);
    vY = _mm_mul_ps(vY, B->r[1]);
    vZ = _mm_mul_ps(vZ, B->r[2]);
    vW = _mm_mul_ps(vW, B->r[3]);
    vX = _mm_add_ps(vX, vZ);
    vY = _mm_add_ps(vY, vW);
    vX = _mm_add_ps(vX, vY);
    mResult.r[2] = vX;
    vW = A->r[3];
    vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
    vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
    vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
    vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
    vX = _mm_mul_ps(vX, B->r[0]);
    vY = _mm_mul_ps(vY, B->r[1]);
    vZ = _mm_mul_ps(vZ, B->r[2]);
    vW = _mm_mul_ps(vW, B->r[3]);
    vX = _mm_add_ps(vX, vZ);
    vY = _mm_add_ps(vY, vW);
    vX = _mm_add_ps(vX, vY);
    mResult.r[3] = vX;
    return mResult;
#endif
}