#pragma once

XMGLOBALCONST XMVECTORF32 g_BoxOffset[8] =
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

XMGLOBALCONST XMVECTORF32 g_RayEpsilon = { { { 1e-20f, 1e-20f, 1e-20f, 1e-20f } } };
XMGLOBALCONST XMVECTORF32 g_RayNegEpsilon = { { { -1e-20f, -1e-20f, -1e-20f, -1e-20f } } };
XMGLOBALCONST XMVECTORF32 g_FltMin = { { { -FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX } } };
XMGLOBALCONST XMVECTORF32 g_FltMax = { { { FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX } } };



#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable:4324 4820)
// C4324: alignment padding warnings
// C4820: Off by default noise
#endif

_Use_decl_annotations_
XMContainmentType XM_CALLCONV XMBoundingSphereContainsPoint(const XMBoundingSphere *const s, const FXMVECTOR p)
{
    XMVECTOR center = XMLoadFloat3(&s->center);
    XMVECTOR r = XMVectorReplicatePtr(&s->r);

    XMVECTOR pt_center = XMVectorSubtract(p, XM_REF_1V(center));
    XMVECTOR dist_sqr = XMVector3LengthSq(XM_REF_1V(pt_center));
    XMVECTOR r_sqr = XMVectorMultiply(XM_REF_2V(r, r));

    return XMVector3LessOrEqual(XM_REF_2V(dist_sqr, r_sqr)) ? CONTAINS : DISJOINT;
}

#define XM_BOUNDING_SPHERE_CONTAINS_POINT(SPHERE, P) XMBoundingSphereContainsPoint(SPHERE, XM_REF_1V(P))

_Use_decl_annotations_
void XMBoundingSphereMerged(_Out_ XMBoundingSphere *out, _In_ const XMBoundingSphere *s1, _In_ const XMBoundingSphere *s2)
{
    XMVECTOR c1 = XMLoadFloat3(&s1->center);
    float r1 = s1->r;

    XMVECTOR c2 = XMLoadFloat3(&s2->center);
    float r2 = s2->r;

    XMVECTOR dist_vec = XMVectorSubtract(XM_REF_2V(c2, c1));

    XMVECTOR dist_len = XMVector3Length(XM_REF_1V(dist_vec));

    float dx = XMVectorGetX(XM_REF_1V(dist_len));

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

    XMVECTOR dir12 = XMVectorDivide(XM_REF_2V(dist_vec, dist_len));

    float t1 = fminf(-r1, dx - r2);
    float t2 = fmaxf(r1, dx + r2);
    float new_r = (t2 - t1) * 0.5f;

    XMVECTOR what = XMVectorReplicate(new_r + t1);
    XMVECTOR center_displacement = XMVectorMultiply(XM_REF_2V(dir12, what));
    XMVECTOR new_center = XMVectorAdd(XM_REF_2V(c1, center_displacement));

    XMStoreFloat3(&out->center, XM_REF_1V(new_center));
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
void XMBoundingSphereFromPoints(XMBoundingSphere* out, size_t count, const XMFLOAT3* points, size_t stride)
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
        
        float px = XMVectorGetX(XM_REF_1V(p_vec));
        float py = XMVectorGetY(XM_REF_1V(p_vec));
        float pz = XMVectorGetZ(XM_REF_1V(p_vec));

        if (px < XMVectorGetX(XM_REF_1V(minx)))
            minx = p_vec;

        if (px > XMVectorGetX(XM_REF_1V(maxx)))
            maxx = p_vec;

        if (py < XMVectorGetY(XM_REF_1V(miny)))
            miny = p_vec;

        if (py > XMVectorGetY(XM_REF_1V(maxy)))
            maxy = p_vec;

        if (pz < XMVectorGetZ(XM_REF_1V(minz)))
            minz = p_vec;

        if (pz > XMVectorGetZ(XM_REF_1V(maxz)))
            maxz = p_vec;
    }

    // Now use the min/max pair that are farthest apart to form the initial sphere
    // with the biggest possible radius. Note that it will not necessarily contain
    // all the points, but it is the optimal starting point to reduce the number of
    // iterations to find the final bounding sphere.

    XMVECTOR deltax = XMVectorSubtract(XM_REF_2V(maxx, minx));
    XMVECTOR distx = XMVector3Length(XM_REF_1V(deltax));

    XMVECTOR deltay = XMVectorSubtract(XM_REF_2V(maxy, miny));
    XMVECTOR disty = XMVector3Length(XM_REF_1V(deltay));

    XMVECTOR deltaz = XMVectorSubtract(XM_REF_2V(maxz, minz));
    XMVECTOR distz = XMVector3Length(XM_REF_1V(deltaz));

    XMVECTOR r;
    XMVECTOR center;

    if (XMVector3Greater(XM_REF_2V(distx, disty)))
    {
        if (XMVector3Greater(XM_REF_2V(distx, distz)))
        {
            // Use min/max x.
            center = XMVectorLerp(XM_REF_2V(maxx, minx), 0.5f);
            r = XMVectorScale(XM_REF_1V(distx), 0.5f);
        }
        else
        {
            // Use min/max z.
            center = XMVectorLerp(XM_REF_2V(maxz, minz), 0.5f);
            r = XMVectorScale(XM_REF_1V(distz), 0.5f);
        }
    }
    else // Y >= X
    {
        if (XMVector3Greater(XM_REF_2V(disty, distz)))
        {
            // Use min/max y.
            center = XMVectorLerp(XM_REF_2V(maxy, miny), 0.5f);
            r = XMVectorScale(XM_REF_1V(disty), 0.5f);
        }
        else
        {
            // Use min/max z.
            center = XMVectorLerp(XM_REF_2V(maxz, minz), 0.5f);
            r = XMVectorScale(XM_REF_1V(distz), 0.5f);
        }
    }

    // Add any points not inside the sphere and recreate the sphere
    for (size_t i = 0; i < count; ++i)
    {
        const uint8_t* start = (const uint8_t*)(points) + (i * stride);
        XMVECTOR p_vec = XMLoadFloat3((const XMFLOAT3*)start);

        XMVECTOR delta = XMVectorSubtract(XM_REF_2V(p_vec, center));

        XMVECTOR dist_center = XMVector3Length(XM_REF_1V(delta));

        // Adjust sphere to include the new point.
        if (XMVector3Greater(XM_REF_2V(dist_center, r)))
        {
            // For example, if old radius is 19 and dist to center is 21 we will create a
            // new sphere with radius 20 and with the center moved towards the new point
            const XMVECTOR new_diameter = XMVectorAdd(XM_REF_2V(r, dist_center));
            r = XMVectorScale(XM_REF_1V(new_diameter), 0.5f);
            XMVECTOR ratio = XMVectorDivide(XM_REF_2V(r, dist_center));

            XMVECTOR ones = XMVectorReplicate(1.0f);
            XMVECTOR displacement_factor = XMVectorSubtract(XM_REF_2V(ones, ratio));

            // In the example above, it will move around 20/21 (around 95%) of the way
            // towards the new point. The closer the new point is, the more we will move
            // towards it. Note that the points can't be much distant, as we already picked
            // the most distant ones from the start.
            XMVECTOR center_displacement = XMVectorMultiply(XM_REF_2V(delta, displacement_factor));
            center = XMVectorAdd(XM_REF_2V(center, center_displacement));
        }
    }

    XMStoreFloat3(&out->center, XM_REF_1V(center));
    XMStoreFloat(&out->r, XM_REF_1V(r));
}