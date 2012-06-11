
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_NaiadMesh_H
#define PBRT_SHAPES_NaiadMesh_H

// shapes/NaiadMesh.h*
#include "shape.h"
#include <map>
using std::map;

// NaiadMesh Declarations
class NaiadMesh : public Shape {
public:
    // NaiadMesh Public Methods
    NaiadMesh(const Transform *o2w, const Transform *w2o, bool ro,
            const ParamSet &params);
    ~NaiadMesh();
    BBox ObjectBound() const;
    BBox WorldBound() const;
    bool CanIntersect() const { return false; }
    void Refine(vector<Reference<Shape> > &refined) const;
    friend class NaiadTriangle;
    template <typename T> friend class VertexTexture;

    float foam_block;
protected:
    // NaiadMesh Protected Data
    int ntris, nverts;
    int *vertexIndex;
    Point *p;
    Normal *n;
    Vector *s;
    float *uvs;
    Reference<Texture<float> > alphaTexture;
};


class NaiadTriangle : public Shape {
public:
    // NaiadTriangle Public Methods
    NaiadTriangle(const Transform *o2w, const Transform *w2o, bool ro,
             NaiadMesh *m, int n)
        : Shape(o2w, w2o, ro) {
        mesh = m;
        v = &mesh->vertexIndex[3*n];
        PBRT_CREATED_TRIANGLE(this);
    }
    BBox ObjectBound() const;
    BBox WorldBound() const;
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
    void GetUVs(float uv[3][2]) const {
        if (mesh->uvs) {
            uv[0][0] = mesh->uvs[2*v[0]];
            uv[0][1] = mesh->uvs[2*v[0]+1];
            uv[1][0] = mesh->uvs[2*v[1]];
            uv[1][1] = mesh->uvs[2*v[1]+1];
            uv[2][0] = mesh->uvs[2*v[2]];
            uv[2][1] = mesh->uvs[2*v[2]+1];
        }
        else {
            uv[0][0] = 0.; uv[0][1] = 0.;
            uv[1][0] = 1.; uv[1][1] = 0.;
            uv[2][0] = 1.; uv[2][1] = 1.;
        }
    }
    float Area() const;
    virtual void GetShadingGeometry(const Transform &obj2world,
            const DifferentialGeometry &dg,
            DifferentialGeometry *dgShading) const;
    Point Sample(float u1, float u2, Normal *Ns) const;
private:
    // NaiadTriangle Private Data
    Reference<NaiadMesh> mesh;
    int *v;
};


NaiadMesh *CreateNaiadMeshShape(const Transform *o2w, const Transform *w2o,
    bool reverseOrientation, const ParamSet &params,
    map<string, Reference<Texture<float> > > *floatTextures = NULL);

#endif // PBRT_SHAPES_NaiadMesh_H
