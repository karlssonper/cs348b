/*
 * naiadparticle.h
 *
 *  Created on: Jun 5, 2012
 *      Author: per
 */

#ifndef NAIADFOAM_H_
#define NAIADFOAM_H_

#include "shape.h"
#include <Ni.h>
#include <Nb.h>
#include <NbBody.h>
#include <NbBlock.h>

class NaiadFoam;
class NaiadFoamBox;
class NaiadFoamParticle {
public:
    NaiadFoamParticle(Point position, float radius,
            NaiadFoamBox * parent);
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
    NaiadFoamBox * parent;
    Point pos;
    float r;
};

class NaiadFoamBox : public Shape{
public:
    NaiadFoamBox(const Transform *o2w, const Transform *w2o, bool ro,
            Reference<NaiadFoam> p, const Nb::ParticleShape & particle, int blockID);
    BBox ObjectBound() const { return bb; };
    BBox WorldBound() const { return bb; };
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
    float Density() const { return density;};
    void AddToFoamPlane(std::vector<float> &foamPlane, int w, int h,
                        const Transform &w2c,
                        const Transform &c2s,
                        const Transform &s2r,
                        const std::vector<float> & zbuffer);
    Reference<NaiadFoam> parent;
private:
    std::vector<NaiadFoamParticle> particles;
    BBox bb;
    float density;
};

class Scene;
class NaiadFoam : public Shape {
public:
    NaiadFoam(const Transform *o2w, const Transform *w2o,
            bool ro, const ParamSet &params);
    ~NaiadFoam();
    BBox ObjectBound() const { return bb; };
    BBox WorldBound() const { return bb; };
    bool CanIntersect() const { return false; }
    void Refine(vector<Reference<Shape> > &refined) const;
    void CreateFoamPlane(Scene * s);
    static NaiadFoam * cur;
    std::vector<float> & FoamPlane() const {
        return *foamPlane;
    } ;

private:
    std::vector<float> * foamPlane;
    std::vector<NaiadFoamBox *> boxes;
    BBox bb;
    float FoamInLowerBound;
    float FoamInUpperBound;
    float FoamInGamma;
    float FoamOutLowerBound;
    float FoamOutUpperBound;
};

NaiadFoam *CreateNaiadFoam(const Transform *o2w, const Transform *w2o,
    bool reverseOrientation, const ParamSet &params);

#endif /* NAIADFOAM_H_ */
