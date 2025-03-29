#pragma once

XMGLOBALCONST XMVECTORF32 g_box_offset[8] =
{
    { { { -1.0f, -1.0f,  1.0f, 0.0f } } },
    { { {  1.0f, -1.0f,  1.0f, 0.0f } } },
    { { {  1.0f,  1.0f,  1.0f, 0.0f } } },
    { { { -1.0f,  1.0f,  1.0f, 0.0f } } },
    { { { -1.0f, -1.0f, -1.0f, 0.0f } } },
    { { {  1.0f, -1.0f, -1.0f, 0.0f } } },
    { { {  1.0f,  1.0f, -1.0f, 0.0f } } },
    { { { -1.0f,  1.0f, -1.0f, 0.0f } } },
};

XMGLOBALCONST XMVECTORF32 g_ray_eps = { { { 1e-20f, 1e-20f, 1e-20f, 1e-20f } } };
XMGLOBALCONST XMVECTORF32 g_ray_neg_eps = { { { -1e-20f, -1e-20f, -1e-20f, -1e-20f } } };
XMGLOBALCONST XMVECTORF32 g_flt_min = { { { -FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX } } };
XMGLOBALCONST XMVECTORF32 g_glt_max = { { { FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX } } };



#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4324 4820)
// C4324: alignment padding warnings
// C4820: Off by default noise
#endif

_Use_decl_annotations_
inline dx_containment_type XM_CALLCONV dx_bounding_sphere_contains_point(const dx_bounding_sphere *const s, const FXMVECTOR p)
{
    XMVECTOR center = XMLoadFloat3(&s->center);
    XMVECTOR r = XMVectorReplicatePtr(&s->r);

    XMVECTOR pt_center = XMVectorSubtract(p, &center);
    XMVECTOR dist_sqr = XMVector3LengthSq(&pt_center);
    XMVECTOR r_sqr = XMVectorMultiply(&r, &r);

    return XMVector3LessOrEqual(&dist_sqr, &r_sqr) ? CONTAINS : DISJOINT;
}

_Use_decl_annotations_
inline void dx_bounding_sphere_merged(_Out_ dx_bounding_sphere *out, _In_ const dx_bounding_sphere *s1, _In_ const dx_bounding_sphere *s2)
{
    XMVECTOR c1 = XMLoadFloat3(&s1->center);
    float r1 = s1->r;

    XMVECTOR c2 = XMLoadFloat3(&s2->center);
    float r2 = s2->r;

    XMVECTOR dist_vec = XMVectorSubtract(&c2, &c1);

    XMVECTOR dist_len = XMVector3Length(&dist_vec);

    float dx = XMVectorGetX(&dist_len);

    // if the distance between centers is less than the sum of the radii, the spheres necessarily overlap or touch
    if (r1 + r2 >= dx)
    {
        if (r1 >= r2 + dx)
        {
            out->r = s1->r;
            out->center = s1->center;
            return;
        }
        else if (r2 >= r1 + dx)
        {
            out->r = s2->r;
            out->center = s2->center;
            return;
        }
    }

    XMVECTOR dir12 = XMVectorDivide(&dist_vec, &dist_len);

    float t1 = fminf(-r1, dx - r2);
    float t2 = fmaxf(r1, dx + r2);
    float new_r = (t2 - t1) * 0.5f;

    XMVECTOR what = XMVectorReplicate(new_r + t1);
    XMVECTOR center_displacement = XMVectorMultiply(&dir12, &what);
    XMVECTOR new_center = XMVectorAdd(&c1, &center_displacement);

    XMStoreFloat3(&out->center, &new_center);
    out->r = new_r;
}


//-----------------------------------------------------------------------------
// Find the approximate smallest enclosing bounding sphere for a set of
// points. Exact computation of the smallest enclosing bounding sphere is
// possible but is slower and requires a more complex algorithm.
// The algorithm is based on  Jack Ritter, "An Efficient Bounding Sphere",
// Graphics Gems.
//-----------------------------------------------------------------------------
_Use_decl_annotations_
inline void dx_bounding_sphere_from_points(dx_bounding_sphere* out, size_t count, const XMFLOAT3* points, size_t stride)
{
    assert(count > 0);
    assert(points);

    XMVECTOR minx, maxx, miny, maxy, minz, maxz;
    minx = maxx = miny = maxy = minz = maxz = XMLoadFloat3(points);

    // Find the points with minimum and maximum x, y, and z
    for (size_t i = 1; i < count; ++i)
    {
        const uint8_t* start = (const uint8_t*) (points) + (i * stride);
        XMVECTOR p_vec = XMLoadFloat3((const XMFLOAT3*) start);

        float px = XMVectorGetX(&p_vec);
        float py = XMVectorGetY(&p_vec);
        float pz = XMVectorGetZ(&p_vec);

        if (px < XMVectorGetX(&minx))
            minx = p_vec;

        if (px > XMVectorGetX(&maxx))
            maxx = p_vec;

        if (py < XMVectorGetY(&miny))
            miny = p_vec;

        if (py > XMVectorGetY(&maxy))
            maxy = p_vec;

        if (pz < XMVectorGetZ(&minz))
            minz = p_vec;

        if (pz > XMVectorGetZ(&maxz))
            maxz = p_vec;
    }

    // Now use the min/max pair that are farthest apart to form the initial sphere
    // with the biggest possible radius. Note that it will not necessarily contain
    // all the points, but it is the optimal starting point to reduce the number of
    // iterations to find the final bounding sphere.

    XMVECTOR deltax = XMVectorSubtract(&maxx, &minx);
    XMVECTOR distx = XMVector3Length(&deltax);

    XMVECTOR deltay = XMVectorSubtract(&maxy, &miny);
    XMVECTOR disty = XMVector3Length(&deltay);

    XMVECTOR deltaz = XMVectorSubtract(&maxz, &minz);
    XMVECTOR distz = XMVector3Length(&deltaz);

    XMVECTOR r;    
    XMVECTOR center;

    if (XMVector3Greater(&distx, &disty))
    {
        if (XMVector3Greater(&distx, &distz))
        {
            // Use min/max x.
            center = XMVectorLerp(&maxx, &minx, 0.5f);
            r = XMVectorScale(&distx, 0.5f);
        }
        else
        {
            // Use min/max z.
            center = XMVectorLerp(&maxz, &minz, 0.5f);
            r = XMVectorScale(&distz, 0.5f);
        }
    }
    else // Y >= X
    {
        if (XMVector3Greater(&disty, &distz))
        {
            // Use min/max y.
            center = XMVectorLerp(&maxy, &miny, 0.5f);
            r = XMVectorScale(&disty, 0.5f);
        }
        else
        {
            // Use min/max z.
            center = XMVectorLerp(&maxz, &minz, 0.5f);
            r = XMVectorScale(&distz, 0.5f);
        }
    }

    // Add any points not inside the sphere and recreate the sphere
    for (size_t i = 0; i < count; ++i)
    {
        const uint8_t* start = (const uint8_t*)(points) + (i * stride);
        XMVECTOR p_vec = XMLoadFloat3((const XMFLOAT3*)start);

        XMVECTOR delta = XMVectorSubtract(&p_vec, &center);

        XMVECTOR dist_center = XMVector3Length(&delta);

        // Adjust sphere to include the new point.
        if (XMVector3Greater(&dist_center, &r))
        {
            // For example, if old radius is 19 and dist to center is 21 we will create a
            // new sphere with radius 20 and with the center moved towards the new point
            const XMVECTOR new_diameter = XMVectorAdd(&r, &dist_center);
            r = XMVectorScale(&new_diameter, 0.5f);
            XMVECTOR ratio = XMVectorDivide(&r, &dist_center);

            XMVECTOR ones = XMVectorReplicate(1.0f);
            XMVECTOR displacement_factor = XMVectorSubtract(&ones, &ratio);

            // In the example above, it will move around 20/21 (around 95%) of the way
            // towards the new point. The closer the new point is, the more we will move
            // towards it. Note that the points can't be much distant, as we already picked
            // the most distant ones from the start.
            XMVECTOR center_displacement = XMVectorMultiply(&delta, &displacement_factor);
            center = XMVectorAdd(&center, &center_displacement);
        }
    }

    XMStoreFloat3(&out->center, &center);
    XMStoreFloat(&out->r, &r);
}