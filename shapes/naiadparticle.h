/*
 * naiadparticle.h
 *
 *  Created on: Jun 5, 2012
 *      Author: per
 */

#ifndef NAIADPARTICLE_H_
#define NAIADPARTICLE_H_

#include "shape.h"
#include <Ni.h>
#include <Nb.h>
#include <NbBody.h>
#include <NbBlock.h>

class FoamParticle {
public:
    FoamParticle(Point position, float radius,
            Shape * parent);
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
private:
    Shape * parent;
    Point pos;
    float r;
};

class FoamBox : public Shape{
public:
    FoamBox(const Transform *o2w, const Transform *w2o,
            bool ro, const Nb::ParticleShape & particle, int blockID);
    BBox ObjectBound() const { return bb; };
    BBox WorldBound() const { return bb; };
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
    float Density() const { return density;};


    float Area() const {std::cerr << "LOL\n" << std::endl;};

    Point Sample(float u1, float u2, Normal *ns) const{std::cerr << "LOL\n" << std::endl;};
    Point Sample(const Point &p, float u1, float u2, Normal *ns) const{std::cerr << "LOL\n" << std::endl;};
    float Pdf(const Point &p, const Vector &wi) const{std::cerr << "LOL\n" << std::endl;};
private:
    std::vector<FoamParticle> particles;
    BBox bb;
    float density;
};

class NaiadFoamParticles : public Shape {
public:
    NaiadFoamParticles(const Transform *o2w, const Transform *w2o,
            bool ro, const ParamSet &params);
    ~NaiadFoamParticles();
    BBox ObjectBound() const { return bb; };
    BBox WorldBound() const { return bb; };
    bool CanIntersect() const { return false; }
    void Refine(vector<Reference<Shape> > &refined) const;
private:
    std::vector<FoamBox *> boxes;
    BBox bb;
};

NaiadFoamParticles *CreateNaiadFoamParticles(const Transform *o2w, const Transform *w2o,
    bool reverseOrientation, const ParamSet &params);

#endif /* NAIADPARTICLE_H_ */
