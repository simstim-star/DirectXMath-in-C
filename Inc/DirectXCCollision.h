#pragma once

#include "DirectXMathC.h"

typedef enum dx_containment_type
{
    DISJOINT = 0,
    INTERSECTS = 1,
    CONTAINS = 2
} dx_containment_type;

typedef enum plane_intersection_type
{
    FRONT = 0,
    INTERSECTING = 1,
    BACK = 2
} plane_intersection_type;

typedef struct dx_bounding_sphere
{
    XMFLOAT3 center;
    float r;

} dx_bounding_sphere;

inline dx_containment_type XM_CALLCONV dx_bounding_sphere_contains_point(_In_ const dx_bounding_sphere* const s, _In_ const FXMVECTOR p);

inline void dx_bounding_sphere_merged(_Out_ dx_bounding_sphere* out, _In_ const dx_bounding_sphere* s1, _In_ const dx_bounding_sphere* s2);

inline void dx_bounding_sphere_from_points(_Out_ dx_bounding_sphere *Out, _In_ size_t count,
    _In_reads_bytes_(sizeof(XMFLOAT3) + stride * (count - 1)) const XMFLOAT3* points, _In_ size_t stride);


#include "DirectXCCollision.inl"