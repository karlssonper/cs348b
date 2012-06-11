#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"
#include <vector>
#include <iostream>
#include "sampler.h"

#define VDB

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#ifdef VDB
#include <vdb.h>
#endif


class RealisticCamera : public Camera {
public:
   RealisticCamera(const AnimatedTransform &cam2world,
      float hither, float yon, float sopen,
      float sclose, float filmdistance, float aperture_diameter,
      const string &specfile,
      float filmdiag,
	  Film *film);
   ~RealisticCamera();
   void CreateLensFlare(Point pointLight, Spectrum pointSpectrum) const;
   float GenerateRay(const CameraSample &sample, Ray *) const;
   struct Lens {
       Vector GetNormal(const Ray * ray, const Point & p) const {
           const float r = fabs(radius);
           const bool posZ = ray->d.z > 0.f;
           const bool convex = radius > 0;
           const float flip = (convex && posZ)||(!convex && !posZ)? -1.f : 1.f;
           return Normalize(flip*Vector(p.x/r,
                                        p.y/r,
                                       (p.z-(zIntersect-radius)) / r));
       }

       bool InsideAperture(const Point &p) const{
           float r = ap*0.4;
           float a = (2.f * r) / (2 + sqrt(2));
           float c = sqrt(2) * a;
           Point vertices[8];
           vertices[0] = Point(-c*0.5, r,  p.z);
           vertices[1] = Point(c*0.5, r,   p.z);
           vertices[2] = Point(r, c*0.5,   p.z);
           vertices[3] = Point(r, -c*0.5,  p.z);
           vertices[4] = Point(c*0.5, -r,  p.z);
           vertices[5] = Point(-c*0.5, -r, p.z);
           vertices[6] = Point(-r, -c*0.5, p.z);
           vertices[7] = Point(-r, c*0.5,  p.z);

           int sign = 0;
           for (int i = 0; i < 8; ++i) {
               Vector seg= vertices[(i+1) % 8] - vertices[i];
               Vector point = p - vertices[i];
               float k = seg.x*point.y - seg.y*point.x;
               if (k == 0.f) return false;
               k = (int)(k/fabs(k));
               if (sign == 0)
                   sign = k;
               else if (k!=sign)
                   return false;
           }
           return true;
       }

       template<int RGB>
       bool GetFresnelTerms(Ray incRay, const Point & p,
               float & R, float & T) const {
           if (radius == 0.f) {
               R = 1.f;
               T = 1.f;
               return true;
           }
           bool from_left = incRay.d.z < 0;

           Vector I = Normalize(incRay.d);
           Vector N = GetNormal(&incRay, p);
           if (!Refract<RGB>(&incRay))
               return false;

           const float cos_theta1 = Dot(-I, N);
           const float cos_theta2 = Dot(-N, incRay.d);

           const float n1 = from_left ? nd_left[RGB] : nd_right[RGB];
           const float n2 = from_left ? nd_right[RGB] : nd_left[RGB];

           const float t1 = (n1 * cos_theta1 - n2 * cos_theta2) /
                            (n1 * cos_theta1 + n2 * cos_theta2);
           const float t2 = (n1 * cos_theta2 - n2 * cos_theta1) /
                            (n1 * cos_theta2 + n2 * cos_theta1);


           R = 0.5 * t1 * t1 + 0.5 * t2 * t2;
           T = 1.f - R;

           if (R != R || T != T)
               return false;

           return true;
       }

       bool Reflect(Ray * ray, bool renderRay = false) const {
           if (ray->d.z == 0.f)
               return false;
           Point p;
           if (!Intersect(*ray, p, renderRay)) return false;

           if (radius == 0.f) {
               *ray = Ray(p, ray->d, 0.f, INFINITY);
               return true;
           }

           Vector I = Normalize(ray->d);
           Vector N =  GetNormal(ray,p);

           *ray = Ray(p, Normalize(I - 2 * Dot(I,N)*N), 0.f, INFINITY);
           return !ray->HasNaNs();
       }

       template <int RGB>
       bool Refract(Ray * ray, bool renderRay = false) const {
           if (ray->d.z == 0.f)
               return false;

           //Get how far to travel until intersection
           Point p;
           if (!Intersect(*ray, p, renderRay)) return false;

           //If the aperture, just continue in the same direction
           if (radius == 0.f) {
               *ray = Ray(p, ray->d, 0.f, INFINITY);
               return true;
           }

           Vector I = Normalize(ray->d);
           Vector N =  GetNormal(ray,p);

           bool from_right = ray->d.z > 0.f;
           float ni = from_right ? nd_right[RGB] : nd_left[RGB];
           float nt = from_right ? nd_left[RGB] : nd_right[RGB];

           float u = ni / nt;
           float sqrt_term = 1.f -u*u*(1.f - (Dot(I,N)*(Dot(I,N))));
           if (sqrt_term < 0.f) return false;
           float y = -u * Dot(I,N) - sqrt(sqrt_term);
           Vector T = I * u + N * y;
           *ray =  Ray(p, Normalize(T), 0.f, INFINITY);

           if (ray->d.z == 0.f)
               return false;

           return true;
       }

       bool Intersect(const Ray & ray, Point & p, bool renderRay = false)const {
           if (radius == 0.f) {
               p = ray((zIntersect - ray.o.z) / ray.d.z);
               if (!InsideAperture(p)) return false;
           } else {
               const float Zc = zIntersect - radius;
               const float Sr = fabs(radius);
               const float X0 = ray.o.x;
               const float Y0 = ray.o.y;
               const float Z0 = ray.o.z;
               const float Xd = ray.d.x;
               const float Yd = ray.d.y;
               const float Zd = ray.d.z;

               const float B = 2 * (Xd * X0 + Yd * Y0 + Zd * (Z0 - Zc));
               const float C = X0*X0 + Y0*Y0 + (Z0 - Zc)*(Z0 - Zc) - Sr*Sr;
               const float discriminant = B*B - 4*C;
               if (discriminant < 0)
                   return false;

               const float t0 = (-B - sqrt(discriminant)) * 0.5f;
               const float t1 = (-B + sqrt(discriminant)) * 0.5f;

               const float & tt = t0 > 0 ? t0 : t1;
               if (tt < 0 || tt != tt)
                   return false;

               p = ray(tt);
               if (p.HasNaNs())
                   return false;

               //Check if aperture blocks ray
               if (sqrt(p.x*p.x+p.y*p.y) > ap*0.5f)
                   return false;
               }



#ifdef VDB
           if (renderRay) {
               vdb_color(r, g, b);
               vdb_line(ray.o.x, ray.o.y, ray.o.z, p.x, p.y, p.z);
           }
#endif

           return true;
       }

       float radius, thick, ap, zIntersect;
       float nd_left[3], nd_right[3];
       int idx;
#ifdef VDB
       float r,g,b;
#endif
   };
private:
   float ShutterOpen;
   float ShutterClose;
   float FilmDistance;
   float FilmDiagonal;
   float Width;
   float Height;
   float WidthHalf;
   float HeightHalf;
   float FilmZ;
   float DiscZ;
   float Hither;
   float SplatFactor;
   int RaysPerPair;
   Film * film;

   std::vector<Lens> LensTable;

   struct LensFlareRay {
       Ray ray;
       unsigned int firstReflectLens, secondReflectLens;
       bool firstReflection;
       Spectrum spec;
   };

   //bool ReflectLensFlareRay(LensFlareRay & lfr, int lens) const;
   //bool RefractLensFlareRay(LensFlareRay & lfr, int lens) const;
   //void TrackLensFlare(LensFlareRay & lfr) const;

   template <int RGB>
   bool ReflectLensFlareRay(LensFlareRay & lfr,
                                             int lens) const {
       float R,T;
       Ray incRay = lfr.ray;
       if (!LensTable[lens].Reflect(&lfr.ray, true)) return false;
       if (!LensTable[lens].GetFresnelTerms<RGB>(incRay,lfr.ray.o, R, T)) return false;
       lfr.spec *= R;
       return true;
   }

   template <int RGB>
   bool RefractLensFlareRay(LensFlareRay & lfr,
                                             int lens) const {
       float R,T;
       Ray incRay = lfr.ray;
       if (!LensTable[lens].Refract<RGB>(&lfr.ray, true)) return false;
       if (!LensTable[lens].GetFresnelTerms<RGB>(incRay,lfr.ray.o, R, T)) return false;
       lfr.spec *= T;
       return true;
   }

   template <int RGB>
   void TrackLensFlare(LensFlareRay & lfr)const{
       //Start with the first lens
       unsigned int lens = 0;
       while (lens < LensTable.size()) {
           bool from_light = lfr.ray.d.z < 0;
           if (from_light) {
               if (lens == lfr.firstReflectLens && !lfr.firstReflection) {
                   //Reflect
                   if (!ReflectLensFlareRay<RGB>(lfr, lens)) return;
                   lfr.firstReflection = true;
                   --lens;
               } else {
                   //Refract
                   if (!RefractLensFlareRay<RGB>(lfr, lens)) return;
                   ++lens;
               }
           } else {
               if (lens == lfr.secondReflectLens) {
                   //Reflect
                   if (!ReflectLensFlareRay<RGB>(lfr, lens)) return;
                   ++lens;
               } else {
                   //Refract
                   if (!RefractLensFlareRay<RGB>(lfr, lens)) return;
                   --lens;
               }
           }
       }

       float t = (FilmZ - lfr.ray.o.z) / (lfr.ray.d.z);
       Point p = lfr.ray.o + lfr.ray.d * t;
       if (p.x > -WidthHalf &&
           p.x < WidthHalf &&
           p.y > -HeightHalf &&
           p.y < HeightHalf) {
           CameraSample cs;
           cs.imageX = (WidthHalf - p.x)  * (float)film->xResolution / Width;
           cs.imageY = (p.y + HeightHalf) * (float)film->yResolution / Height;
           film->Splat(cs, lfr.spec);
       }
   }

#ifdef VDB
template<class VEC3>
void vdb_square(const VEC3 &p0,
                const VEC3 &p1,
                const VEC3 &p2,
                const VEC3 &p3) const {
    vdb_line(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);
    vdb_line(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
    vdb_line(p2.x, p2.y, p2.z, p3.x, p3.y, p3.z);
    vdb_line(p3.x, p3.y, p3.z, p0.x, p0.y, p0.z);
}

void vdb_circle(float r, float z, int numsteps) const {
    float dTheta = 2.0f * M_PI / float(numsteps-1);
    for (int step = 0; step < numsteps; ++step) {
        vdb_line(r * cos(dTheta*step), r * sin(dTheta*step), z,
                 r * cos(dTheta*(step+1)), r * sin(dTheta*(step+1)), z);
    }

}

void vdb_sphere(float r, float z, float ap, int numsteps, bool convex) const {
    float dTheta = M_PI / float(2*numsteps-1);
    float dPhi = 2.0f * M_PI / float(numsteps-1);
    for (int stepTheta = 0; stepTheta < 2*numsteps; ++stepTheta) {
        for (int stepPhi = 0; stepPhi < numsteps; ++stepPhi) {
            Point p(r*cos(dPhi*stepPhi)*sin(dTheta*stepTheta),
                    r*sin(dPhi*stepPhi)*sin(dTheta*stepTheta),
                                  r*cos(dTheta*stepTheta));

            if (sqrt(p.x*p.x + p.y*p.y) > ap*0.5f) continue;
            if (convex && p.z < 0) continue;
            if (!convex && p.z > 0) continue;
            vdb_point(p.x, p.y, z + p.z);
        }
    }
}

#endif
};

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif
