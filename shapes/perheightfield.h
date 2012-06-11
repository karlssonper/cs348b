
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

#ifndef PBRT_SHAPES_PERHEIGHTFIELD_H
#define PBRT_SHAPES_PERHEIGHTFIELD_H

// shapes/perheightfield.h*
#include "shape.h"
#include "paramset.h"
#include <set>
#include <iostream>

class PerTriangle : public Shape {
public:
    // Triangle Public Methods
    PerTriangle(const Point * op, const Point * wp, const Vector * n,
            const Transform *o2w, const Transform *w2o, bool ro,
            int vertex1, int vertex2, int vertex3)
        : Shape(o2w, w2o, ro), objectPos(op), worldPos(wp), normals(n),
          v1(vertex1), v2(vertex2), v3(vertex3){
        /*bb.pMin = Point(
                min(min(objectPos[v1].x,objectPos[v2].x),objectPos[v3].x),
                min(min(objectPos[v1].y,objectPos[v2].y),objectPos[v3].y),
                min(min(objectPos[v1].z,objectPos[v2].z),objectPos[v3].z)
        );
        bb.pMax = Point(
                max(max(objectPos[v1].x,objectPos[v2].x),objectPos[v3].x),
                max(max(objectPos[v1].y,objectPos[v2].y),objectPos[v3].y),
                max(max(objectPos[v1].z,objectPos[v2].z),objectPos[v3].z)
        );*/
    }

    BBox ObjectBound() const { return bb;};

    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray, float & t) const;
    float hit(const Ray &ray) const;
    void GetUVs(float uv[3][2]) const {
            uv[0][0] = objectPos[v1].x;
            uv[0][1] = objectPos[v1].y;
            uv[1][0] = objectPos[v2].x;
            uv[1][1] = objectPos[v2].y;
            uv[2][0] = objectPos[v3].x;
            uv[2][1] = objectPos[v3].y;
    }

    void GetShadingGeometry(const Transform &obj2world,
            const DifferentialGeometry &dg,
            DifferentialGeometry *dgShading) const;
private:
    // Triangle Private Data
    const Point * objectPos;
    const Point * worldPos;
    const Vector * normals;
    BBox bb;
    int v1,v2,v3;
};

template<int T_CELLS_PER_GRID_AXIS, int T_ORDER>
class GridBox {
public:
    GridBox(){};
    GridBox(const Point * objectPos, const Point * worldPos, const Vector * n,
            const Transform *o2w, const Transform *w2o,
            int x, int y, int size, int nx, int ny) : hasData_(false) {
        for (int i = 0; i < T_CELLS_PER_GRID_AXIS; ++i) {
            for (int j = 0; j < T_CELLS_PER_GRID_AXIS; ++j) {
                int d = size/T_CELLS_PER_GRID_AXIS;
                grids_[j][i] = new GridBox<T_CELLS_PER_GRID_AXIS,T_ORDER-1>(
                        objectPos,worldPos,n,o2w,w2o,x+i*d,y+j*d,d,nx,ny);
            }
        }
        float maxZ = grids_[0][0]->box().pMax.z;
        float minZ = grids_[0][0]->box().pMin.z;
        for (int i = 0; i < T_CELLS_PER_GRID_AXIS; ++i) {
            for (int j = 0; j < T_CELLS_PER_GRID_AXIS; ++j) {
                if (grids_[j][i]->hasData()) {
                    hasData_ = true;
                    maxZ = max(maxZ,grids_[j][i]->box().pMax.z);
                    minZ = min(minZ,grids_[j][i]->box().pMin.z);
                }
            }
        }
        float invNxMinOne = 1.f / (float)(nx-1);
        float invNyMinOne = 1.f / (float)(ny-1);
        int posX = x + size < nx ? x + size : nx - 1;
        int posY = y + size < ny ? y + size : ny - 1;
        bb.pMin = Point((float)x*invNxMinOne, (float)y*invNyMinOne, minZ);
        bb.pMax = Point((float)posX*invNxMinOne,(float)posY*invNyMinOne,maxZ);
    }


    bool hasData() const { return hasData_; } ;

    const BBox& box() const { return bb; } ;

    float hit(const Ray &ray) const {
        float t0,t1;
        if (hasData() && bb.IntersectP(ray,&t0,&t1))
            return t0;
        else
            return 0.f;
    }

    struct IntersectionGrid{
        float t;
        int i,j;
        bool operator <(const IntersectionGrid & ig) const {
            return t < ig.t;
        }
    };

    bool Intersect(const Ray &oRay, const Ray &wRay, float *tHit,
            float *rayEpsilon, DifferentialGeometry *dg) const {
        std::set<IntersectionGrid> igs;
        for (int j = 0; j < T_CELLS_PER_GRID_AXIS; ++j) {
            for (int i = 0; i < T_CELLS_PER_GRID_AXIS; ++i) {
                IntersectionGrid ig = {0.f, i, j };
                if (grids_[j][i]->hasData()) {
                    ig.t = grids_[j][i]->hit(oRay);
                    if (ig.t) {
                        igs.insert(ig);
                    }
                }
            }
        }

        typename std::set<IntersectionGrid>::iterator it;
        for (it = igs.begin(); it != igs.end(); ++it) {
            if (grids_[it->j][it->i]->Intersect(oRay,wRay,
                    tHit, rayEpsilon, dg)) {
                return true;
            }
        }
        return false;

    }


    bool IntersectP(const Ray &oRay, const Ray &wRay) const {
        std::set<IntersectionGrid> igs;
        for (int j = 0; j < T_CELLS_PER_GRID_AXIS; ++j) {
            for (int i = 0; i < T_CELLS_PER_GRID_AXIS; ++i) {
                IntersectionGrid ig = {0.f, i, j };
                if (grids_[j][i]->hasData()) {
                    ig.t = grids_[j][i]->hit(oRay);
                    if (ig.t) {
                        igs.insert(ig);
                    }
                }
            }
        }

        typename std::set<IntersectionGrid>::iterator it;
        for (it = igs.begin(); it != igs.end(); ++it) {
            if (grids_[it->j][it->i]->IntersectP(oRay,wRay)) {
                return true;
            }
        }
        return false;
    }
private:
    GridBox<T_CELLS_PER_GRID_AXIS, T_ORDER-1>* grids_[T_CELLS_PER_GRID_AXIS]
                                                     [T_CELLS_PER_GRID_AXIS];
    bool hasData_;
    BBox bb;
};

template<int T_CELLS_PER_GRID_AXIS>
class GridBox<T_CELLS_PER_GRID_AXIS, 0> {
public:
    GridBox(){};
    GridBox(const Point * objectPos, const Point * worldPos, const Vector * n,
            const Transform *o2w, const Transform *w2o,
            int x, int y, int size, int nx, int ny) {

        hasData_ = (x >= nx - 1 || y >= ny - 1) ? false : true;
        if (!hasData())return;

        float invNxMinOne = 1.f / float(nx-1);
        float invNyMinOne = 1.f / float(ny-1);

        np_x = x + T_CELLS_PER_GRID_AXIS < nx ? T_CELLS_PER_GRID_AXIS : nx-1-x;
        np_y = y + T_CELLS_PER_GRID_AXIS < ny ? T_CELLS_PER_GRID_AXIS : ny-1-y;

        float maxZ = objectPos[x + y * nx].z;
        float minZ = maxZ;
        for (int j = 0; j <= np_y; ++j) {
            for (int i = 0; i <= np_x; ++i) {
                float curZ = objectPos[x+i + (y+j) * nx].z;
                maxZ = max(maxZ, curZ);
                minZ = min(minZ, curZ);
            }
        }

        int posX = x+T_CELLS_PER_GRID_AXIS< nx ? x+T_CELLS_PER_GRID_AXIS : nx-1;
        int posY = y+T_CELLS_PER_GRID_AXIS< ny ? y+T_CELLS_PER_GRID_AXIS : ny-1;
        bb.pMin = Point((float)x*invNxMinOne, (float)y*invNyMinOne, minZ);
        bb.pMax = Point((float)posX*invNxMinOne,(float)posY*invNyMinOne, maxZ);

        for (int j = 0; j < np_y; ++j) {
            for (int i = 0; i < np_x; ++i) {
                triLow_[j][i] = new PerTriangle(objectPos,worldPos,n,o2w,w2o,
                        false, x+i+(y+j)*nx,x+i+1+(y+j)*nx,x+i+1+(y+j+1)*nx);
                triUp_[j][i] = new PerTriangle(objectPos,worldPos,n,o2w,w2o,
                        false, x+i+(y+j)*nx,x+i+1+(y+j+1)*nx,x+i+(y+j+1)*nx);
            }
        }
    };

    const BBox& box() const {
        return bb;
    };

    bool hasData() const {
        return hasData_;
    };

    float hit(const Ray &ray) const {
        float t0,t1;
        if (hasData() && bb.IntersectP(ray,&t0,&t1))
            return t0;
        else
            return 0.f;
    };

    bool Intersect(const Ray &oRay, const Ray &wRay, float *tHit,
            float *rayEpsilon, DifferentialGeometry *dg) const {
        if (!hasData()) return false;
        int ii, jj;
        float tt = 9999999999999999;
        bool lower;
        bool hitTri = false;

        float t;
        for (int j = 0; j < np_y; ++j) {
            for (int i = 0; i < np_x; ++i) {
                //std::cout << "Intersect" <<i<<j<< std::endl;
                //Lower
                if (triLow_[j][i]->IntersectP(wRay,t) && t < tt) {
                    tt = t;
                    ii = i;
                    jj = j;
                    lower = true;
                    hitTri = true;
                }

                //Upper
                if (triUp_[j][i]->IntersectP(wRay,t) && t < tt) {
                    tt = t;
                    ii = i;
                    jj = j;
                    lower = false;
                    hitTri = true;
                }
            }
        }
        if (!hitTri) return false;

        if (lower) {
            triLow_[jj][ii]->Intersect(wRay, tHit, rayEpsilon, dg);
        } else {
            triUp_[jj][ii]->Intersect(wRay, tHit, rayEpsilon, dg);
        }

        return true;
    }

    bool IntersectP(const Ray &oRay, const Ray &wRay) const {
        if (!hasData()) return false;
        float t;
        for (int j = 0; j < np_y; ++j) {
            for (int i = 0; i < np_x; ++i) {
                if (triLow_[j][i]->IntersectP(wRay,t)) {
                    return true;
                }
                if (triUp_[j][i]->IntersectP(wRay,t)) {
                    return true;
                }
            }
        }
        return false;
    }

private:
    bool hasData_;
    PerTriangle* triLow_[T_CELLS_PER_GRID_AXIS][T_CELLS_PER_GRID_AXIS];
    PerTriangle* triUp_[T_CELLS_PER_GRID_AXIS][T_CELLS_PER_GRID_AXIS];
    BBox bb;
    int np_x,np_y;
};

// Heightfield Declarations
template<int T_CELLS_PER_GRID_AXIS, int T_ORDER>
class PerHeightfield : public Shape {
public:
    // Heightfield Public Methods
    PerHeightfield(const Transform *o2, const Transform *w2o, bool ro, int nu,
            int nv, const float *zs);
    ~PerHeightfield();
    bool CanIntersect() const {
        return true;
    };

    BBox WorldBound() const {
        return worldBB;
    }
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;

    BBox ObjectBound() const {
        return BBox(Point(0,0,minZobject), Point(1,1,maxZobject));
    };

    void GetShadingGeometry(const Transform &obj2world,
            const DifferentialGeometry &dg,
            DifferentialGeometry *dgShading) const {
        dg.shape->GetShadingGeometry(obj2world, dg, dgShading);
    }

private:
    // PerHeightfield Private Data
    Point * objectPos;
    Point * worldPos;
    Vector * normals;

    float minZobject, maxZobject;
    BBox worldBB;

    GridBox<T_CELLS_PER_GRID_AXIS,T_ORDER> *uniformGrid;

    int nx, ny;
};



template<int T_CELLS_PER_GRID_AXIS, int T_ORDER>
PerHeightfield<T_CELLS_PER_GRID_AXIS,T_ORDER> *CreatePerHeightfieldShape(
        const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new PerHeightfield<T_CELLS_PER_GRID_AXIS,T_ORDER>(
            o2w, w2o, reverseOrientation, nu, nv, Pz);
};

// PerHeightfield Method Definitions
template<int T_CELLS_PER_GRID_AXIS, int T_ORDER>
PerHeightfield<T_CELLS_PER_GRID_AXIS,T_ORDER>::PerHeightfield(
        const Transform *o2w, const Transform *w2o,
        bool ro, int x, int y, const float *zs)
    : Shape(o2w, w2o, ro) {
    nx = x;
    ny = y;

    int size = pow(2,ceil(log2(nx)));

    std::cout << "nx: " << nx << std::endl;
    std::cout << "ny: " << ny << std::endl;
    //z = new float[nx*ny];
    //memcpy(z, zs, nx*ny*sizeof(float));
    objectPos = new Point[nx*ny];
    worldPos = new Point[nx*ny];




    minZobject = zs[0];
    maxZobject = zs[0];
    worldBB.pMin = Point(0,0,0);
    worldBB.pMax = Point(0,0,0);
    float invNxMinOne = 1.f / float(nx - 1);
    float invNyMinOne = 1.f / float(ny - 1);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (zs[i + j*nx] < minZobject) minZobject = zs[i + j*nx];
            if (zs[i + j*nx] > maxZobject) maxZobject = zs[i + j*nx];
            objectPos[i + j*nx] =
                 Point(float(i)*invNxMinOne,float(j)*invNyMinOne,zs[i + j*nx]);
            worldPos[i + j*nx] = (*o2w)(objectPos[i + j*nx]);
            worldBB = Union(worldBB, worldPos[i + j*nx]);
        }
    }

    Vector  * normalsLow = new Vector[(nx-1)*(ny-1)];
    Vector * normalsUp = new Vector[(nx-1)*(ny-1)];
    for (int j = 0; j < ny-1; ++j) {
        for (int i = 0; i < nx-1; ++i) {
            Vector v1 = objectPos[i + 1 +     j*nx] - objectPos[i + j*nx];
            Vector v2 = objectPos[i + 1 + (j+1)*nx] - objectPos[i + j*nx];
            Vector v3 = objectPos[i +     (j+1)*nx] - objectPos[i + j*nx];
            normalsLow[i + j * (nx-1)] = (Cross (v1,v2));
            normalsUp[ i + j * (nx-1)] = (Cross (v2,v3));
        }
    }

    normals = new Vector[nx *ny];
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            normals[i + j*nx] = Vector(0,0,0);
            if (i > 0 && j < ny-1) {
                normals[i + j*nx] += normalsLow[i-1 + j*(nx-1)];
            }
            if (i > 0 && j > 0) {
                normals[i + j*nx] += normalsLow[i-1 + (j-1)*(nx-1)];
                normals[i + j*nx] += normalsUp[ i-1 + (j-1)*(nx-1)];
            }
            if (i < nx - 1 && j > 0) {
                normals[i + j*nx] += normalsUp[i + (j-1)*(nx-1)];
            }
            if (i < nx - 1 && j < ny-1){
                normals[i + j*nx] += normalsLow[i + j*(nx-1)];
                normals[i + j*nx] += normalsUp[ i + j*(nx-1)];
            }
        }
    }

    uniformGrid = new GridBox<T_CELLS_PER_GRID_AXIS,T_ORDER>(
            objectPos, worldPos,normals, o2w,w2o,0,0,size,nx,ny);
}

template<int T_CELLS_PER_GRID_AXIS, int T_ORDER>
PerHeightfield<T_CELLS_PER_GRID_AXIS,T_ORDER>::~PerHeightfield() {
    //delete[] z;
    delete[] objectPos;
    delete[] worldPos;
    delete uniformGrid;
}

template<int T_CELLS_PER_GRID_AXIS, int T_ORDER>
bool PerHeightfield<T_CELLS_PER_GRID_AXIS,T_ORDER>::Intersect(
        const Ray &ray, float *tHit, float *rayEpsilon,
               DifferentialGeometry *dg) const {
    Ray objectRay = (*WorldToObject)(ray);

    if (uniformGrid->hit(objectRay))
        return uniformGrid->Intersect(objectRay, ray, tHit, rayEpsilon,dg);
    else
        return false;
}

template<int T_CELLS_PER_GRID_AXIS, int T_ORDER>
bool PerHeightfield<T_CELLS_PER_GRID_AXIS,T_ORDER>::IntersectP(
        const Ray &ray) const {
    Ray objectRay = (*WorldToObject)(ray);

    if (uniformGrid->hit(objectRay))
        return uniformGrid->IntersectP(objectRay,ray);
    else
        return false;
}


#endif // PBRT_SHAPES_PERHEIGHTFIELD_H
