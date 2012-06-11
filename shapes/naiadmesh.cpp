
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


// shapes/NaiadMesh.cpp*
#include "stdafx.h"
#include "shapes/naiadmesh.h"
#include "texture.h"
#include "textures/constant.h"
#include "paramset.h"
#include "montecarlo.h"

//Naiad
#include <Nb.h>
#include <NbField.h>
#include <NbBody.h>

//#define VDB
#ifdef VDB
#include <vdb.h>
#endif

// NaiadMesh Method Definitions
NaiadMesh::NaiadMesh(const Transform *o2w, const Transform *w2o,
        bool ro, const ParamSet &params)
    : Shape(o2w, w2o, ro) {
    Nb::begin();

    const char * empStr = params.FindOneString("emp","error").c_str();
    const char * bodyStr = params.FindOneString("body","error").c_str();
    Nb::EmpReader empReader(empStr,"*");
    const Nb::Body * body = empReader.ejectBody(bodyStr);

    //Load shapes from Naiad body
    const Nb::TriangleShape& triangle = body->constTriangleShape();
    const Nb::PointShape& point = body->constPointShape();

    //Get buffers
    const Nb::Buffer3f& posBuf(point.constBuffer3f("position"));
    const Nb::Buffer3i& triIdxBuf(triangle.constBuffer3i("index"));

    ntris = triIdxBuf.size();
    nverts = posBuf.size();

    std::cout << "NaiadMesh:: Found " << nverts << " vertices." << std::endl;
    std::cout << "NaiadMesh:: Found " << ntris << " triangles." << std::endl;

    vertexIndex = new int[3 * ntris];
    for (int i = 0; i < ntris; ++i) {
        vertexIndex[3*i    ] = triIdxBuf[i].v[0];
        vertexIndex[3*i + 1] = triIdxBuf[i].v[1];
        vertexIndex[3*i + 2] = triIdxBuf[i].v[2];
    }
    //memcpy(vertexIndex, &triIdxBuf[0], 3 * ntris * sizeof(int));

    float dispSurface = params.FindOneFloat("disp_surface",0.f);

    p = new Point[nverts];
    for (int i = 0; i < nverts; ++i) {
        Point P(posBuf[i].v[0], posBuf[i].v[1], posBuf[i].v[2]);
        if (dispSurface) {
            const Nb::Buffer3f& normalBuf(point.constBuffer3f("normal"));
            P.x += normalBuf[i].v[0] * dispSurface;
            P.y += normalBuf[i].v[1] * dispSurface;
            P.z += normalBuf[i].v[2] * dispSurface;
        }
        p[i] = (*ObjectToWorld)(P);

#ifdef VDB
        if (ntris > 50000) {
            vdb_sample(0.1);
            vdb_color(0.5f, 0.8f, 1.0f);
        } else {
            vdb_sample(1.0);
            vdb_color(0.5f, 0.5f, 0.5f);
        }

        vdb_point(p[i].x,p[i].y,p[i].z);
#endif
    }

    string uv_empStr = params.FindOneString("uv_emp","error");
    if (uv_empStr != string("error")) {
        Nb::EmpReader empReaderUV(uv_empStr.c_str(),"*");
        string uv_bodyStr = params.FindOneString("uv_body","error");
        const Nb::Body * bodyUV = empReaderUV.ejectBody(uv_bodyStr.c_str());

        const Nb::TriangleShape& tUV = bodyUV->constTriangleShape();
        const Nb::Buffer3i& uvIdx(tUV.constBuffer3i("index"));
        const Nb::Buffer3f& bufU(tUV.constBuffer3f("u"));
        const Nb::Buffer3f& bufV(tUV.constBuffer3f("v"));
        uvs = new float[2*nverts];
        for (int i = 0; i < uvIdx.size(); ++i) {
            int v0 = triIdxBuf[i].v[0];
            int v1 = triIdxBuf[i].v[1];
            int v2 = triIdxBuf[i].v[2];
            uvs[2*v0]     = bufU[i].v[0];
            uvs[2*v0 + 1] = bufV[i].v[0];
            uvs[2*v1]     = bufU[i].v[1];
            uvs[2*v1 + 1] = bufV[i].v[1];
            uvs[2*v2]     = bufU[i].v[2];
            uvs[2*v2 + 1] = bufV[i].v[2];
        }
    } else uvs = NULL;

    //No need for these parameters yet.
    //uvs = NULL;
    n = NULL;
    s = NULL;

    //For foam zbuffer
    foam_block = params.FindOneFloat("foam_block", 0.f);

    Nb::end();
}


NaiadMesh::~NaiadMesh() {
    delete[] vertexIndex;
    delete[] p;
    delete[] s;
    delete[] n;
    delete[] uvs;
}


BBox NaiadMesh::ObjectBound() const {
    BBox objectBounds;
    for (int i = 0; i < nverts; i++)
        objectBounds = Union(objectBounds, (*WorldToObject)(p[i]));
    return objectBounds;
}


BBox NaiadMesh::WorldBound() const {
    BBox worldBounds;
    for (int i = 0; i < nverts; i++)
        worldBounds = Union(worldBounds, p[i]);
    return worldBounds;
}


void NaiadMesh::Refine(vector<Reference<Shape> > &refined) const {
    for (int i = 0; i < ntris; ++i)
        refined.push_back(new NaiadTriangle(ObjectToWorld,
                          WorldToObject, ReverseOrientation,
                          (NaiadMesh *)this, i));
}


BBox NaiadTriangle::ObjectBound() const {
    // Get NaiadTriangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    return Union(BBox((*WorldToObject)(p1), (*WorldToObject)(p2)),
                 (*WorldToObject)(p3));
}


BBox NaiadTriangle::WorldBound() const {
    // Get NaiadTriangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    return Union(BBox(p1, p2), p3);
}


bool NaiadTriangle::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                         DifferentialGeometry *dg) const {
    if (!mesh.GetPtr()->foam_block && ray.zbuffer) return false;

    PBRT_RAY_TRIANGLE_INTERSECTION_TEST(const_cast<Ray *>(&ray), const_cast<NaiadTriangle *>(this));
    // Compute $\VEC{s}_1$

    // Get NaiadTriangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);

    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector d = ray.o - p1;
    float b1 = Dot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    Vector s2 = Cross(d, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;

    // Compute NaiadTriangle partial derivatives
    Vector dpdu, dpdv;
    float uvs[3][2];
    GetUVs(uvs);

    // Compute deltas for NaiadTriangle partial derivatives
    float du1 = uvs[0][0] - uvs[2][0];
    float du2 = uvs[1][0] - uvs[2][0];
    float dv1 = uvs[0][1] - uvs[2][1];
    float dv2 = uvs[1][1] - uvs[2][1];
    Vector dp1 = p1 - p3, dp2 = p2 - p3;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for NaiadTriangle partial derivative matrix
        CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
    }
    else {
        float invdet = 1.f / determinant;
        dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }

    // Interpolate $(u,v)$ NaiadTriangle parametric coordinates
    float b0 = 1 - b1 - b2;
    float tu = b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0];
    float tv = b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1];

    //printf("\nu: %f, v: %f\n",tu,tv);
    // Test intersection against alpha texture, if present
    if (ray.depth != -1) {
    if (mesh->alphaTexture) {
        DifferentialGeometry dgLocal(ray(t), dpdu, dpdv,
                                     Normal(0,0,0), Normal(0,0,0),
                                     tu, tv, this);
        if (mesh->alphaTexture->Evaluate(dgLocal) == 0.f)
            return false;
    }
    }

    // Fill in _DifferentialGeometry_ from NaiadTriangle hit
    *dg = DifferentialGeometry(ray(t), dpdu, dpdv,
                               Normal(0,0,0), Normal(0,0,0),
                               tu, tv, this);
    *tHit = t;

    //If writing to Z buffer
    if (mesh.GetPtr()->foam_block ) {
        ray.zbuffer_depth = t;
    }

    *rayEpsilon = 1e-3f * *tHit;

    PBRT_RAY_TRIANGLE_INTERSECTION_HIT(const_cast<Ray *>(&ray), t);
    return true;
}


bool NaiadTriangle::IntersectP(const Ray &ray) const {
    PBRT_RAY_TRIANGLE_INTERSECTIONP_TEST(const_cast<Ray *>(&ray), const_cast<NaiadTriangle *>(this));
    // Compute $\VEC{s}_1$

    // Get NaiadTriangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);
    
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector d = ray.o - p1;
    float b1 = Dot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    Vector s2 = Cross(d, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;

    // Test shadow ray intersection against alpha texture, if present
    if (ray.depth != -1 && mesh->alphaTexture) {
        // Compute NaiadTriangle partial derivatives
        Vector dpdu, dpdv;
        float uvs[3][2];
        GetUVs(uvs);

        // Compute deltas for NaiadTriangle partial derivatives
        float du1 = uvs[0][0] - uvs[2][0];
        float du2 = uvs[1][0] - uvs[2][0];
        float dv1 = uvs[0][1] - uvs[2][1];
        float dv2 = uvs[1][1] - uvs[2][1];
        Vector dp1 = p1 - p3, dp2 = p2 - p3;
        float determinant = du1 * dv2 - dv1 * du2;
        if (determinant == 0.f) {
            // Handle zero determinant for NaiadTriangle partial derivative matrix
            CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
        }
        else {
            float invdet = 1.f / determinant;
            dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
            dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
        }

        // Interpolate $(u,v)$ NaiadTriangle parametric coordinates
        float b0 = 1 - b1 - b2;
        float tu = b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0];
        float tv = b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1];
        DifferentialGeometry dgLocal(ray(t), dpdu, dpdv,
                                     Normal(0,0,0), Normal(0,0,0),
                                     tu, tv, this);
        if (mesh->alphaTexture->Evaluate(dgLocal) == 0.f)
            return false;
    }
    PBRT_RAY_TRIANGLE_INTERSECTIONP_HIT(const_cast<Ray *>(&ray), t);
    return true;
}


float NaiadTriangle::Area() const {
    // Get NaiadTriangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    return 0.5f * Cross(p2-p1, p3-p1).Length();
}


void NaiadTriangle::GetShadingGeometry(const Transform &obj2world,
        const DifferentialGeometry &dg,
        DifferentialGeometry *dgShading) const {
    if (!mesh->n && !mesh->s) {
        *dgShading = dg;
        return;
    }
    // Initialize _NaiadTriangle_ shading geometry with _n_ and _s_

    // Compute barycentric coordinates for point
    float b[3];

    // Initialize _A_ and _C_ matrices for barycentrics
    float uv[3][2];
    GetUVs(uv);
    float A[2][2] =
        { { uv[1][0] - uv[0][0], uv[2][0] - uv[0][0] },
          { uv[1][1] - uv[0][1], uv[2][1] - uv[0][1] } };
    float C[2] = { dg.u - uv[0][0], dg.v - uv[0][1] };
    if (!SolveLinearSystem2x2(A, C, &b[1], &b[2])) {
        // Handle degenerate parametric mapping
        b[0] = b[1] = b[2] = 1.f/3.f;
    }
    else
        b[0] = 1.f - b[1] - b[2];

    // Use _n_ and _s_ to compute shading tangents for NaiadTriangle, _ss_ and _ts_
    Normal ns;
    Vector ss, ts;
    if (mesh->n) ns = Normalize(obj2world(b[0] * mesh->n[v[0]] +
                                          b[1] * mesh->n[v[1]] +
                                          b[2] * mesh->n[v[2]]));
    else   ns = dg.nn;
    if (mesh->s) ss = Normalize(obj2world(b[0] * mesh->s[v[0]] +
                                          b[1] * mesh->s[v[1]] +
                                          b[2] * mesh->s[v[2]]));
    else   ss = Normalize(dg.dpdu);
    
    ts = Cross(ss, ns);
    if (ts.LengthSquared() > 0.f) {
        ts = Normalize(ts);
        ss = Cross(ts, ns);
    }
    else
        CoordinateSystem((Vector)ns, &ss, &ts);
    Normal dndu, dndv;

    // Compute $\dndu$ and $\dndv$ for NaiadTriangle shading geometry
    if (mesh->n) {
        float uvs[3][2];
        GetUVs(uvs);
        // Compute deltas for NaiadTriangle partial derivatives of normal
        float du1 = uvs[0][0] - uvs[2][0];
        float du2 = uvs[1][0] - uvs[2][0];
        float dv1 = uvs[0][1] - uvs[2][1];
        float dv2 = uvs[1][1] - uvs[2][1];
        Normal dn1 = mesh->n[v[0]] - mesh->n[v[2]];
        Normal dn2 = mesh->n[v[1]] - mesh->n[v[2]];
        float determinant = du1 * dv2 - dv1 * du2;
        if (determinant == 0.f)
            dndu = dndv = Normal(0,0,0);
        else {
            float invdet = 1.f / determinant;
            dndu = ( dv2 * dn1 - dv1 * dn2) * invdet;
            dndv = (-du2 * dn1 + du1 * dn2) * invdet;
        }
    }
    else
        dndu = dndv = Normal(0,0,0);
    *dgShading = DifferentialGeometry(dg.p, ss, ts,
        (*ObjectToWorld)(dndu), (*ObjectToWorld)(dndv),
        dg.u, dg.v, dg.shape);
    dgShading->dudx = dg.dudx;  dgShading->dvdx = dg.dvdx;
    dgShading->dudy = dg.dudy;  dgShading->dvdy = dg.dvdy;
    dgShading->dpdx = dg.dpdx;  dgShading->dpdy = dg.dpdy;
}


NaiadMesh *CreateNaiadMeshShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params,
        map<string, Reference<Texture<float> > > *floatTextures) {
    return new NaiadMesh(o2w, w2o, reverseOrientation, params);
}


Point NaiadTriangle::Sample(float u1, float u2, Normal *Ns) const {
    float b1, b2;
    UniformSampleTriangle(u1, u2, &b1, &b2);
    // Get NaiadTriangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    Point p = b1 * p1 + b2 * p2 + (1.f - b1 - b2) * p3;
    Normal n = Normal(Cross(p2-p1, p3-p1));
    *Ns = Normalize(n);
    if (ReverseOrientation) *Ns *= -1.f;
    return p;
}


