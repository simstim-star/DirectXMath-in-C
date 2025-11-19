#pragma once

#include "DirectXMathC.h"

typedef enum XMContainmentType
{
    DISJOINT = 0,
    INTERSECTS = 1,
    CONTAINS = 2
} XMContainmentType;

typedef enum XMPlaneIntersectionType
{
    FRONT = 0,
    INTERSECTING = 1,
    BACK = 2
} XMPlaneIntersectionType;

typedef struct XMBoundingSphere
{
    XMFLOAT3 center;
    float r;

} XMBoundingSphere;

XM_INLINE XMContainmentType XM_CALLCONV XMBoundingSphereContainsPoint(_In_ const XMBoundingSphere* const s, _In_ const FXMVECTOR p);

XM_INLINE void XMBoundingSphereMerged(_Out_ XMBoundingSphere* out, _In_ const XMBoundingSphere* s1, _In_ const XMBoundingSphere* s2);

XM_INLINE void XMBoundingSphereFromPoints(_Out_ XMBoundingSphere *Out, _In_ size_t count,
    _In_reads_bytes_(sizeof(XMFLOAT3) + stride * (count - 1)) const XMFLOAT3* points, _In_ size_t stride);


#include "DirectXCollisionC.inl"