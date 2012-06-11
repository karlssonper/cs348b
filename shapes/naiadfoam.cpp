/*
 * naiadparticle.cpp
 *
 *  Created on: Jun 6, 2012
 *      Author: per
 */

#include "naiadfoam.h"
#include "scene.h"
#include "paramset.h"
#include "cameras/perspective.h"
//#include <pngwriter.h>
///#include "/home/per/Downloads/pngwriter-0.5.4/src/pngwriter.h"

//#define VDB
#ifdef VDB
#include <vdb.h>
#endif

NaiadFoam * NaiadFoam::cur = NULL;

NaiadFoamParticle::NaiadFoamParticle(Point position, float radius,
        NaiadFoamBox * p) :  pos(position), r(radius), parent(p){

}

bool NaiadFoamParticle::Intersect(const Ray &rayInc, float *tHit, float *rayEpsilon,
               DifferentialGeometry *dg) const {
    /*const float X0 = ray.o.x;
    const float Y0 = ray.o.y;
    const float Z0 = ray.o.z;

    const float Xd = ray.d.x;
    const float Yd = ray.d.y;
    const float Zd = ray.d.z;

    const float Xc = pos.x;
    const float Yc = pos.y;
    const float Zc = pos.z;

    const float B = 2.0f * (Xd * (X0 - Xc) + Yd * (Y0 - Yc) + Zd * (Z0 - Zc));
    const float C = (X0-Xc)*(X0-Xc) + (Y0-Yc)*(Y0-Yc) + (Z0-Zc)*(Z0-Zc) - r*r;

    const float discriminant = B*B - 4*C;
    if (discriminant < 0.0f)
        return false;
    const float t0 = (- B - sqrt(discriminant)) / 2.0f;
    const float t1 = (- B + sqrt(discriminant)) / 2.0f;
    if (t0 < 0 || t0 != t0) {
        return true;
    }*/

    Matrix4x4 w2o_m;
    w2o_m.m[0][3] = -pos.x;
    w2o_m.m[1][3] = -pos.y;
    w2o_m.m[2][3] = -pos.z;
    Transform w2o = Transform(w2o_m);

    const float thetaMin = -M_PI;
    const float thetaMax = 0;
    const float phiMax = 2.f*M_PI;

    float phi;
    Point phit;
    // Transform _Ray_ to object space
    Ray ray;
    w2o(rayInc, &ray);

    // Compute quadratic sphere coefficients
    float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
    float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z);
    float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y +
              ray.o.z*ray.o.z - r*r;

    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    // Compute intersection distance along ray
    if (t0 > ray.maxt || t1 < ray.mint)
        return false;
    float thit = t0;
    if (t0 < ray.mint) {
        thit = t1;
        if (thit > ray.maxt) return false;
    }

    // Compute sphere hit position and $\phi$
    phit = ray(thit);
    rayInc.hitFoam = true;
    PerspectiveCamera * cam = PerspectiveCamera::cur_cam;
    Transform c2w;
    cam->CameraToWorld.Interpolate(0.f, &c2w);
    Transform w2c = Inverse(c2w);
    Transform c2s = cam->CameraToScreen;
    Transform s2r = cam->ScreenToRaster;

    Point cPos = w2c(phit);
    Point rPos = s2r(c2s(cPos));
    int x = (int)(rPos.x + 0.5);
    int y = (int)(rPos.y + 0.5);
    int w = cam->film->xResolution;
    int h = cam->film->yResolution;

    if (x > 0 && y > 0 && x < w && y < h) {
        rayInc.hitFoam = true;
        /*printf("%i %i %i %i\n", x, y, w ,h);
        std::cout << parent << std::endl;
        std::cout << parent->parent << std::endl;
        std::cout << parent->parent->foamPlane[0] << std::endl;
        std::cout << parent->parent->foamPlane.size() << std::endl;
        std::cout << x + y*w << std::endl;*/
        //std::cout << NaiadFoam::cur->FoamPlane().size() << std::endl;
        rayInc.alphaFoam = NaiadFoam::cur->FoamPlane()[x + y*w];
    }

    /*if (phit.x == 0.f && phit.y == 0.f) phit.x = 1e-5f * r;
    phi = atan2f(phit.y, phit.x);
    if (phi < 0.) phi += 2.f*M_PI;

    // Find parametric representation of sphere hit
    float u = phi / phiMax;
    float theta = acosf(Clamp(phit.z / r, -1.f, 1.f));
    float v = (theta - thetaMin) / (thetaMax - thetaMin);

    // Compute sphere $\dpdu$ and $\dpdv$
    float zradius = sqrtf(phit.x*phit.x + phit.y*phit.y);
    float invzradius = 1.f / zradius;
    float cosphi = phit.x * invzradius;
    float sinphi = phit.y * invzradius;
    Vector dpdu(-phiMax * phit.y, phiMax * phit.x, 0);
    Vector dpdv = (thetaMax-thetaMin) *
        Vector(phit.z * cosphi, phit.z * sinphi,
               -r * sinf(theta));

    // Compute sphere $\dndu$ and $\dndv$
    Vector d2Pduu = -phiMax * phiMax * Vector(phit.x, phit.y, 0);
    Vector d2Pduv = (thetaMax - thetaMin) * phit.z * phiMax *
                    Vector(-sinphi, cosphi, 0.);
    Vector d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) *
                    Vector(phit.x, phit.y, phit.z);

    // Compute coefficients for fundamental forms
    float E = Dot(dpdu, dpdu);
    float F = Dot(dpdu, dpdv);
    float G = Dot(dpdv, dpdv);
    Vector N = Normalize(Cross(dpdu, dpdv));
    float e = Dot(N, d2Pduu);
    float f = Dot(N, d2Pduv);
    float g = Dot(N, d2Pdvv);

    // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
    float invEGF2 = 1.f / (E*G - F*F);
    Normal dndu = Normal((f*F - e*G) * invEGF2 * dpdu +
                         (e*F - f*E) * invEGF2 * dpdv);
    Normal dndv = Normal((g*F - f*G) * invEGF2 * dpdu +
                         (f*F - g*E) * invEGF2 * dpdv);

    //std::cout << "Got here?" << std::endl;

    Matrix4x4 o2w_m;
    o2w_m.m[0][3] = pos.x;
    o2w_m.m[1][3] = pos.y;
    o2w_m.m[2][3] = pos.z;
    Transform o2w = Transform(o2w_m);

    // Initialize _DifferentialGeometry_ from parametric information
    *dg = DifferentialGeometry(o2w(phit), o2w(dpdu), o2w(dpdv),
                               o2w(dndu), o2w(dndv), u, v, parent);

    dg->mult = 1.f;*/

    // Update _tHit_ for quadric intersection
    //*tHit = thit;

    // Compute _rayEpsilon_ for quadric intersection
    //*rayEpsilon = 5e-4f * *tHit;
    return true;
}

bool NaiadFoamParticle::IntersectP(const Ray &rayInc) const {
    Matrix4x4 w2o_m;
    w2o_m.m[0][3] = -pos.x;
    w2o_m.m[1][3] = -pos.y;
    w2o_m.m[2][3] = -pos.z;
    Transform w2o = Transform(w2o_m);

    float phi;
    Point phit;
    // Transform _Ray_ to object space
    Ray ray;
    w2o(rayInc, &ray);

    // Compute quadratic sphere coefficients
    float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
    float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z);
    float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y +
              ray.o.z*ray.o.z - r*r;

    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    // Compute intersection distance along ray
    if (t0 > ray.maxt || t1 < ray.mint)
        return false;
    float thit = t0;
    if (t0 < ray.mint) {
        thit = t1;
        if (thit > ray.maxt) return false;
    }

    return true;
}

NaiadFoamBox::NaiadFoamBox(const Transform *o2w, const Transform *w2o,
        bool ro, Reference<NaiadFoam> par, const Nb::ParticleShape & particle, int blockID)
            : Shape(o2w, w2o, ro), parent(par){
    //Get the blocks of particles
    const Nb::BlockArray3f& blocksPos = particle.constBlocks3f("position");

    const Nb::Block3f& blockPos = blocksPos(blockID);

    density = blockPos.size();

    bb = BBox(Point(blockPos.min().v[0],
                    blockPos.min().v[1],
                    blockPos.min().v[2]),
              Point(blockPos.max().v[0],
                    blockPos.max().v[1],
                    blockPos.max().v[2]));
    //for (int pIdx = 0; pIdx < 1; ++pIdx) {
    for (int pIdx = 0; pIdx < blockPos.size(); pIdx += 20) {
        Point pos(blockPos.data()[pIdx].v[0],
                  blockPos.data()[pIdx].v[1],
                  blockPos.data()[pIdx].v[2]);

        float r = 0.03f;

        particles.push_back(NaiadFoamParticle(pos,r, this));

#ifdef VDB
        vdb_sample(0.01);
        vdb_color(1.0f, 1.0f, 1.0f);
        vdb_point(pos.x, pos.y, pos.z);
#endif
    }
}

bool NaiadFoamBox::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
               DifferentialGeometry *dg) const {
    if (ray.zbuffer) return false;
    if (ray.hitFoam) return false;
    //std::cout << parent.GetPtr()->FoamPlane().size() << std::endl;
    float t;
    if (bb.IntersectP(ray,&t)) {
        Point phit  = ray(t);
        PerspectiveCamera * cam = PerspectiveCamera::cur_cam;
        Transform c2w;
        cam->CameraToWorld.Interpolate(0.f, &c2w);
        Transform w2c = Inverse(c2w);
        Transform c2s = cam->CameraToScreen;
        Transform s2r = cam->ScreenToRaster;

        Point cPos = w2c(phit);
        Point rPos = s2r(c2s(cPos));
        int x = (int)(rPos.x + 0.5);
        int y = (int)(rPos.y + 0.5);
        int w = cam->film->xResolution;
        int h = cam->film->yResolution;

        if (x > 0 && y > 0 && x < w && y < h) {
            ray.hitFoam = true;
            /*printf("%i %i %i %i\n", x, y, w ,h);
            std::cout << parent << std::endl;
            std::cout << parent->parent << std::endl;
            std::cout << parent->parent->foamPlane[0] << std::endl;
            std::cout << parent->parent->foamPlane.size() << std::endl;
            std::cout << x + y*w << std::endl;*/
            //std::cout << NaiadFoam::cur->FoamPlane().size() << std::endl;
            ray.alphaFoam = NaiadFoam::cur->FoamPlane()[x + y*w];
        }
        return false;
/*
        for (int i = 0; i < particles.size(); i+=100) {
            if (particles[i].Intersect(ray, tHit, rayEpsilon, dg)) {
                return false;
            }
        }*/
    }
    return false;
}
bool NaiadFoamBox::IntersectP(const Ray &ray) const {
    return false;
    //if (ray.fromFoam) return false;
    if (bb.IntersectP(ray))
        for (int i = 0; i < particles.size(); ++i)
            if (particles[i].IntersectP(ray))
                return true;
    return false;
}

void NaiadFoamBox::AddToFoamPlane(std::vector<float> &foamPlane, int w, int h,
                                  const Transform &w2c,
                                  const Transform &c2s,
                                  const Transform &s2r,
                                  const std::vector<float> & zbuffer) {
    for (int i = 0; i < particles.size(); ++i) {
        Point cPos = w2c(particles[i].pos);
        Point rPos = s2r(c2s(cPos));
        int x = (int)rPos.x;
        int y = (int)rPos.y;
        bool inFront = zbuffer[x + y * w] ? cPos.z < zbuffer[x + y * w] : true;
        if (inFront && x > 0 && y > 0 && x < w && y < h) {
            ++foamPlane[x + y * w];
        }
    }
}

NaiadFoam::NaiadFoam(const Transform *o2w, const Transform *w2o,
        bool ro, const ParamSet &params): Shape(o2w, w2o, ro){

    Nb::begin();

    int nFoamEmps;
    const string * empStr = params.FindString("emp", &nFoamEmps);
    const string * bodyStr = params.FindString("body", &nFoamEmps);

    for (int i = 0; i < nFoamEmps; ++i) {
        Nb::EmpReader empReader(empStr[i].c_str(),"*");
        const Nb::Body * body = empReader.ejectBody(bodyStr[i].c_str());
        const Nb::ParticleShape & particle = body->constParticleShape();

        bb = BBox(Point(particle.min().v[0],
                        particle.min().v[1],
                        particle.min().v[2]),
                  Point(particle.max().v[0],
                        particle.max().v[1],
                        particle.max().v[2]));

        //Total amount of particles
        const int64_t nParticles = particle.size();
        printf("Found %i particles in %s\n", nParticles, empStr[i].c_str());

        //Number of blocks
        const int nBlocks = particle.constBlocks3f(0).block_count();

        for (int64_t blockIdx = 0; blockIdx < nBlocks; ++blockIdx) {
           boxes.push_back(new NaiadFoamBox(o2w,w2o,ro,this,particle,blockIdx));
        }
    }

    foamPlane = new std::vector<float>(441000,0);

    FoamInLowerBound  = params.FindOneFloat("FoamInLowerBound",0.f);
    FoamInUpperBound  = params.FindOneFloat("FoamInUpperLowerBound",255.f);
    FoamInGamma       = params.FindOneFloat("FoamInGamma",1.f);
    FoamOutLowerBound = params.FindOneFloat("FoamOutUpperLowerBound",0.f);
    FoamOutUpperBound = params.FindOneFloat("FoamOutUpperLowerBound",255.f);

    printf("%f %f %f %f %f\n", FoamInLowerBound, FoamInUpperBound, FoamInGamma, FoamOutLowerBound, FoamOutUpperBound);

    cur = this;
}

NaiadFoam::~NaiadFoam() {
    //for (int i = 0; i < particles.size(); ++i)
    //for (int i = 0; i < 20; ++i)
        //delete particles[i];

    Nb::end();
}

void NaiadFoam::Refine(vector<Reference<Shape> > &refined) const {
    for (int i = 0; i < boxes.size(); ++i)
        refined.push_back(boxes[i]);
}

void NaiadFoam::CreateFoamPlane(Scene * s) {
    PerspectiveCamera * cam = PerspectiveCamera::cur_cam;
    if (cam) {
        std::cerr << "NaiadFoam: Camera found!" << std::endl;

        Transform c2w;
        cam->CameraToWorld.Interpolate(0.f, &c2w);
        Transform w2c = Inverse(c2w);
        Transform c2s = cam->CameraToScreen;
        Transform s2r = cam->ScreenToRaster;

        const int w = cam->film->xResolution;
        const int h = cam->film->yResolution;

        //foamPlane.resize(w * h,0);
        std::cout << "Size of foam plane: " << FoamPlane().size()<<std::endl;

        for (int i = 0; i < boxes.size(); ++i)
            boxes[i]->AddToFoamPlane(FoamPlane(), w, h, w2c, c2s, s2r, s->zBuffer);

        float max_density = 0;
        for (int i = 0; i < w*h; ++i)
            max_density = max(max_density, FoamPlane()[i]);
        for (int i = 0; i < w*h; ++i) {
            FoamPlane()[i] /= max_density;
            FoamPlane()[i] = (powf(((FoamPlane()[i] * 255.0) - FoamInLowerBound) / (FoamInUpperBound - FoamInLowerBound),
                            FoamInGamma) * (FoamOutUpperBound - FoamOutLowerBound) + FoamOutLowerBound) / 255.0;
        }



        //pngwriter png(w,h, 0, "foam.png");

        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                if (FoamPlane()[x + y* w]) {
                    float d = FoamPlane()[x + y*w];
                    //png.plot(x,h-1-y, d, d, d);
                }
            }
        }
        //png.close();
    } else {
        std::cerr << "No camera found, exiting" << std::endl;
        exit(1);
    }
    std::cout << NaiadFoam::cur->FoamPlane()[0] << std::endl;

}

NaiadFoam *CreateNaiadFoam(const Transform *o2w, const Transform *w2o,
    bool reverseOrientation, const ParamSet &params) {
    return new NaiadFoam(o2w, w2o, reverseOrientation, params);
}
