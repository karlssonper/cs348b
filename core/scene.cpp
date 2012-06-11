
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


// core/scene.cpp*
#include "stdafx.h"
#include "scene.h"
#include "camera.h"
#include "film.h"
#include "sampler.h"
#include "volume.h"
#include "parallel.h"
#include "progressreporter.h"
#include "renderer.h"
#include "cameras/perspective.h"
#include "intersection.h"
#include "shapes/naiadfoam.h"

#include <iostream>
//#include <pngwriter.h>

Scene * Scene::cur = NULL;

// Scene Method Definitions
Scene::~Scene() {
    delete aggregate;
    delete volumeRegion;
    for (uint32_t i = 0; i < lights.size(); ++i)
        delete lights[i];
}


Scene::Scene(Primitive *accel, const vector<Light *> &lts,
             VolumeRegion *vr) {
    lights = lts;
    aggregate = accel;
    volumeRegion = vr;
    // Scene Constructor Implementation
    bound = aggregate->WorldBound();
    if (volumeRegion) bound = Union(bound, volumeRegion->WorldBound());
    cur = this;

    PerspectiveCamera * cam = PerspectiveCamera::cur_cam;
    if (cam) {
        std::cerr << "Scene: Camera found for z buffer!" << std::endl;

        const int w = cam->film->xResolution;
        const int h = cam->film->yResolution;

        zBuffer.resize(w * h,0);
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                Point Pras(x, y, 0);
                Point Pcamera;
                cam->RasterToCamera(Pras, &Pcamera);
                Ray ray = Ray(Point(0,0,0), Normalize(Vector(Pcamera)), 0.f, INFINITY);
                Transform c2w;
                cam->CameraToWorld.Interpolate(0.f, &c2w);
                Ray rayWorld = c2w(ray);
                rayWorld.zbuffer = true;
                Intersection isect;
                if (aggregate->Intersect(rayWorld, &isect) && rayWorld.zbuffer_depth) {
                    zBuffer[x+y*w] = rayWorld.zbuffer_depth;
                }
            }
        }

        //pngwriter png(w,h, 0, "zBuffer.png");
        float maxZ = 0.f, minZ = 1 << 20;
        for (int i = 0; i < w*h; ++i) {
            maxZ = max(maxZ,zBuffer[i]);
            if (zBuffer[i])
                minZ = min(minZ,zBuffer[i]);
        }
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                if (zBuffer[x + y*w]) {
                    //std::cout << zBuffer[x + y*w] << std::endl;
                    const float d = (zBuffer[x + y*w] - minZ) / (maxZ - minZ);
                    //png.plot(x,h-1-y, d, d, d);
                }
            }
        }
        //png.close();
    } else {
        std::cerr << "NaiadFoam: Camera not found!" << std::endl;
    }

    NaiadFoam * naiadfoam = NaiadFoam::cur;
    if (naiadfoam) {
        std::cout << "NaiadFoam: Creating foam plane!" << std::endl;
        naiadfoam->CreateFoamPlane(this);
    } else {
        std::cout << "NaiadFoam: Can't write foam plane!" << std::endl;
    }
}


const BBox &Scene::WorldBound() const {
    return bound;
}


