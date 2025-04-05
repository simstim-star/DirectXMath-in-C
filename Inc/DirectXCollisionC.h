#pragma once

#include "DirectXMathC.h"

typedef enum DirectX_ContainmentType
{
    DISJOINT = 0,
    INTERSECTS = 1,
    CONTAINS = 2
} DirectX_ContainmentType;

typedef enum DirectX_PlaneIntersectionType
{
    FRONT = 0,
    INTERSECTING = 1,
    BACK = 2
} DirectX_PlaneIntersectionType;

typedef struct DirectX_BoundingSphere
{
    XMFLOAT3 center;
    float r;

} DirectX_BoundingSphere;

inline DirectX_ContainmentType XM_CALLCONV DirectX_BoundingSphere_ContainsPoint(_In_ const DirectX_BoundingSphere* const s, _In_ const FXMVECTOR p);

inline void DirectX_BoundingSphere_Merged(_Out_ DirectX_BoundingSphere* out, _In_ const DirectX_BoundingSphere* s1, _In_ const DirectX_BoundingSphere* s2);

inline void DirectX_BoundingSphere_FromPoints(_Out_ DirectX_BoundingSphere *Out, _In_ size_t count,
    _In_reads_bytes_(sizeof(XMFLOAT3) + stride * (count - 1)) const XMFLOAT3* points, _In_ size_t stride);


#include "DirectXCollisionC.inl"