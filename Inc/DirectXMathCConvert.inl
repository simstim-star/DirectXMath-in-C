#pragma once
#include <assert.h>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4701)
// C4701: false positives
#endif

/****************************************************************************
 Load operations
****************************************************************************/

_Use_decl_annotations_
inline XMVECTOR XM_CALLCONV XMLoadFloat3(const XMFLOAT3* pSource)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    V->vector4_f32[0] = pSource->x;
    V->vector4_f32[1] = pSource->y;
    V->vector4_f32[2] = pSource->z;
    V->vector4_f32[3] = 0.f;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    __m128 xy = _mm_castpd_ps(_mm_load_sd((const double*)(pSource)));
    __m128 z = _mm_load_ss(&pSource->z);
    return _mm_movelh_ps(xy, z);
#endif
}

inline XMVECTOR XM_CALLCONV XMLoadInt(_In_ const uint32_t* pSource) {
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    V.vector4_u32[0] = *pSource;
    V.vector4_u32[1] = 0;
    V.vector4_u32[2] = 0;
    V.vector4_u32[3] = 0;
    return V;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_load_ss((const float*)pSource);
#endif
}

inline XMVECTOR XM_CALLCONV XMLoadFloat(_In_ const float* pSource) {
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    V.vector4_f32[0] = *pSource;
    V.vector4_f32[1] = 0.f;
    V.vector4_f32[2] = 0.f;
    V.vector4_f32[3] = 0.f;
    return V;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_load_ss(pSource);
#endif
}


inline XMMATRIX XM_CALLCONV XMLoadFloat3x3(const XMFLOAT3X3* pSource)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_)

    XMMATRIX M;
    M.r[0].vector4_f32[0] = pSource->m[0][0];
    M.r[0].vector4_f32[1] = pSource->m[0][1];
    M.r[0].vector4_f32[2] = pSource->m[0][2];
    M.r[0].vector4_f32[3] = 0.0f;
     
    M.r[1].vector4_f32[0] = pSource->m[1][0];
    M.r[1].vector4_f32[1] = pSource->m[1][1];
    M.r[1].vector4_f32[2] = pSource->m[1][2];
    M.r[1].vector4_f32[3] = 0.0f;
     
    M.r[2].vector4_f32[0] = pSource->m[2][0];
    M.r[2].vector4_f32[1] = pSource->m[2][1];
    M.r[2].vector4_f32[2] = pSource->m[2][2];
    M.r[2].vector4_f32[3] = 0.0f;
    M.r[3].vector4_f32[0] = 0.0f;
    M.r[3].vector4_f32[1] = 0.0f;
    M.r[3].vector4_f32[2] = 0.0f;
    M.r[3].vector4_f32[3] = 1.0f;
    return M;

#elif defined(_XM_SSE_INTRINSICS_)
    __m128 Z = _mm_setzero_ps();

    __m128 V1 = _mm_loadu_ps(&pSource->m[0][0]); // m00, m01, m02, m10
    __m128 V2 = _mm_loadu_ps(&pSource->m[1][1]); // m11, m12, m20, m21
    __m128 V3 = _mm_load_ss(&pSource->m[2][2]);  // m22, 0, 0, 0

    __m128 T1 = _mm_unpackhi_ps(V1, Z); // m02, 0, m10, 0
    __m128 T2 = _mm_unpacklo_ps(V2, Z); // m11, 0, m12, 0
    __m128 T3 = _mm_shuffle_ps(V3, T2, _MM_SHUFFLE(0, 1, 0, 0)); // m22, m22, 0, m11
    __m128 T4 = _mm_movehl_ps(T2, T3); // 0, m11, m12, 0
    __m128 T5 = _mm_movehl_ps(Z, T1);  // m10, 0, 0, 0

    XMMATRIX M;
    M.r[0] = _mm_movelh_ps(V1, T1);                            // m00, m01, m02, 0
    M.r[1] = _mm_add_ps(T4, T5);                               // m10, m11, m12, 0
    M.r[2] = _mm_shuffle_ps(V2, V3, _MM_SHUFFLE(1, 0, 3, 2));  // m20, m21, m22, 0
    M.r[3] = g_XMIdentityR3.v;                                 //   0,   0,   0, 1
    return M;
#endif
}

_Use_decl_annotations_
inline void XM_CALLCONV XMStoreFloat
(
    float* pDestination,
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_)
    *pDestination = XMVectorGetX(V);
#elif defined(_XM_SSE_INTRINSICS_)
    _mm_store_ss(pDestination, XM_PARAM_F(V));
#endif
}


_Use_decl_annotations_
inline void XM_CALLCONV XMStoreFloat3
(
    XMFLOAT3* pDestination,
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_)
    pDestination->x = V->vector4_f32[0];
    pDestination->y = V->vector4_f32[1];
    pDestination->z = V->vector4_f32[2];
#elif defined(_XM_SSE4_INTRINSICS_)
    * (int*)(&pDestination->x) = _mm_extract_ps(*V, 0);
    *(int*)(&pDestination->y) = _mm_extract_ps(*V, 1);
    *(int*)(&pDestination->z) = _mm_extract_ps(*V, 2);
#elif defined(_XM_SSE_INTRINSICS_)
    _mm_store_sd((double*)(pDestination), _mm_castps_pd(XM_PARAM_F(V)));
    __m128 z = XM_PERMUTE_PS(XM_PARAM_F(V), _MM_SHUFFLE(2, 2, 2, 2));
    _mm_store_ss(&pDestination->z, z);
#endif
}


inline void XM_CALLCONV XMStoreFloat4x4
(
    XMFLOAT4X4* pDestination,
    FXMMATRIX M
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_)

    pDestination->m[0][0] = M->r[0].vector4_f32[0];
    pDestination->m[0][1] = M->r[0].vector4_f32[1];
    pDestination->m[0][2] = M->r[0].vector4_f32[2];
    pDestination->m[0][3] = M->r[0].vector4_f32[3];

    pDestination->m[1][0] = M->r[1].vector4_f32[0];
    pDestination->m[1][1] = M->r[1].vector4_f32[1];
    pDestination->m[1][2] = M->r[1].vector4_f32[2];
    pDestination->m[1][3] = M->r[1].vector4_f32[3];

    pDestination->m[2][0] = M->r[2].vector4_f32[0];
    pDestination->m[2][1] = M->r[2].vector4_f32[1];
    pDestination->m[2][2] = M->r[2].vector4_f32[2];
    pDestination->m[2][3] = M->r[2].vector4_f32[3];

    pDestination->m[3][0] = M->r[3].vector4_f32[0];
    pDestination->m[3][1] = M->r[3].vector4_f32[1];
    pDestination->m[3][2] = M->r[3].vector4_f32[2];
    pDestination->m[3][3] = M->r[3].vector4_f32[3];
#elif defined(_XM_SSE_INTRINSICS_)
    _mm_storeu_ps(&pDestination->_11, XM_MATRIX_GET(M,0));
    _mm_storeu_ps(&pDestination->_21, XM_MATRIX_GET(M, 1));
    _mm_storeu_ps(&pDestination->_31, XM_MATRIX_GET(M, 2));
    _mm_storeu_ps(&pDestination->_41, XM_MATRIX_GET(M, 3));
#endif
}

_Use_decl_annotations_
inline void XM_CALLCONV XMStoreFloat4x4A
(
    XMFLOAT4X4A* pDestination,
    FXMMATRIX       M
)
{
    assert(pDestination);
    assert(((uintptr_t)(pDestination) & 0xF) == 0);
#if defined(_XM_NO_INTRINSICS_)

    pDestination->m[0][0] = M->r[0].vector4_f32[0];
    pDestination->m[0][1] = M->r[0].vector4_f32[1];
    pDestination->m[0][2] = M->r[0].vector4_f32[2];
    pDestination->m[0][3] = M->r[0].vector4_f32[3];

    pDestination->m[1][0] = M->r[1].vector4_f32[0];
    pDestination->m[1][1] = M->r[1].vector4_f32[1];
    pDestination->m[1][2] = M->r[1].vector4_f32[2];
    pDestination->m[1][3] = M->r[1].vector4_f32[3];

    pDestination->m[2][0] = M->r[2].vector4_f32[0];
    pDestination->m[2][1] = M->r[2].vector4_f32[1];
    pDestination->m[2][2] = M->r[2].vector4_f32[2];
    pDestination->m[2][3] = M->r[2].vector4_f32[3];

    pDestination->m[3][0] = M->r[3].vector4_f32[0];
    pDestination->m[3][1] = M->r[3].vector4_f32[1];
    pDestination->m[3][2] = M->r[3].vector4_f32[2];
    pDestination->m[3][3] = M->r[3].vector4_f32[3];
#elif defined(_XM_SSE_INTRINSICS_)
    _mm_store_ps(&pDestination->_11, XM_MATRIX_GET(M, 0));
    _mm_store_ps(&pDestination->_21, XM_MATRIX_GET(M, 1));
    _mm_store_ps(&pDestination->_31, XM_MATRIX_GET(M, 2));
    _mm_store_ps(&pDestination->_41, XM_MATRIX_GET(M, 3));
#endif
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMMATRIX XM_CALLCONV XMLoadFloat4x4(const XMFLOAT4X4* pSource)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_)

    XMMATRIX M;
    M.r[0].vector4_f32[0] = pSource->m[0][0];
    M.r[0].vector4_f32[1] = pSource->m[0][1];
    M.r[0].vector4_f32[2] = pSource->m[0][2];
    M.r[0].vector4_f32[3] = pSource->m[0][3];

    M.r[1].vector4_f32[0] = pSource->m[1][0];
    M.r[1].vector4_f32[1] = pSource->m[1][1];
    M.r[1].vector4_f32[2] = pSource->m[1][2];
    M.r[1].vector4_f32[3] = pSource->m[1][3];

    M.r[2].vector4_f32[0] = pSource->m[2][0];
    M.r[2].vector4_f32[1] = pSource->m[2][1];
    M.r[2].vector4_f32[2] = pSource->m[2][2];
    M.r[2].vector4_f32[3] = pSource->m[2][3];

    M.r[3].vector4_f32[0] = pSource->m[3][0];
    M.r[3].vector4_f32[1] = pSource->m[3][1];
    M.r[3].vector4_f32[2] = pSource->m[3][2];
    M.r[3].vector4_f32[3] = pSource->m[3][3];
    return M;
#elif defined(_XM_SSE_INTRINSICS_)
    XMMATRIX M;
    M.r[0] = _mm_loadu_ps(&pSource->_11);
    M.r[1] = _mm_loadu_ps(&pSource->_21);
    M.r[2] = _mm_loadu_ps(&pSource->_31);
    M.r[3] = _mm_loadu_ps(&pSource->_41);
    return M;
#endif
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMMATRIX XM_CALLCONV XMLoadFloat4x4A(const XMFLOAT4X4A* pSource)
{
    assert(pSource);
    assert((((uintptr_t) pSource) & 0xF) == 0);
#if defined(_XM_NO_INTRINSICS_)

    XMMATRIX M;
    M.r[0].vector4_f32[0] = pSource->m[0][0];
    M.r[0].vector4_f32[1] = pSource->m[0][1];
    M.r[0].vector4_f32[2] = pSource->m[0][2];
    M.r[0].vector4_f32[3] = pSource->m[0][3];

    M.r[1].vector4_f32[0] = pSource->m[1][0];
    M.r[1].vector4_f32[1] = pSource->m[1][1];
    M.r[1].vector4_f32[2] = pSource->m[1][2];
    M.r[1].vector4_f32[3] = pSource->m[1][3];

    M.r[2].vector4_f32[0] = pSource->m[2][0];
    M.r[2].vector4_f32[1] = pSource->m[2][1];
    M.r[2].vector4_f32[2] = pSource->m[2][2];
    M.r[2].vector4_f32[3] = pSource->m[2][3];

    M.r[3].vector4_f32[0] = pSource->m[3][0];
    M.r[3].vector4_f32[1] = pSource->m[3][1];
    M.r[3].vector4_f32[2] = pSource->m[3][2];
    M.r[3].vector4_f32[3] = pSource->m[3][3];
    return M;
#elif defined(_XM_SSE_INTRINSICS_)
    XMMATRIX M;
    M.r[0] = _mm_load_ps(&pSource->_11);
    M.r[1] = _mm_load_ps(&pSource->_21);
    M.r[2] = _mm_load_ps(&pSource->_31);
    M.r[3] = _mm_load_ps(&pSource->_41);
    return M;
#endif
}