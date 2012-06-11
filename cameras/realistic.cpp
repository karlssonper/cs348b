// cameras/realistic.cpp*
#include "stdafx.h"
#include "cameras/realistic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/image.h"
#include "samplers/stratified.h"
#include "intersection.h"
#include "renderer.h"

#include <fstream>
#include <algorithm>
#include <vector>
#include <string>


using namespace std;

std::ostream& operator<<(std::ostream& os,const RealisticCamera::Lens & tr) {
    os << "Lens properties:" <<
            "\n\t Radius:    " << tr.radius <<
            "\n\t Thickness: " << tr.thick <<
            "\n\t nd:        " << tr.nd_right <<
            "\n\t Aperture : " << tr.ap << std::endl;
    return os;
}

std::istream& operator>>(std::istream& is, RealisticCamera::Lens & tr) {
#ifdef VDB
    tr.r = (double) rand() / RAND_MAX;
    tr.g = (double) rand() / RAND_MAX;
    tr.b = (double) rand() / RAND_MAX;
#endif
    return is >> tr.radius >> tr.thick >> tr.nd_right[1] >> tr.ap;
}

void RealisticCamera::CreateLensFlare(Point pointLight,
        Spectrum pointSpectrum) const {
    printf("Lens flare from pos (%f, %f, %f)\n",
            pointLight.x, pointLight.y, pointLight.z);

    //The disc infront of the first lens
    //Get where on the Z axis the sample disc should be placed
    bool convex = LensTable.front().radius > 0;
    float r = fabs(LensTable.front().radius);
    float dz = r * cos(asin(-LensTable.front().ap/(2.0f*r)));
    float discZ = convex ? dz-r : r-dz;

    printf("Lens flare sampling disc at z: %f\n", dz);

#ifdef VDB
    vdb_color(1.0f, 0.f, 1.f);
    vdb_circle(LensTable.front().ap * 0.5, discZ, 40);
    vdb_point(pointLight.x, pointLight.y, pointLight.z);
    Point prev = pointLight;
#endif

    float lensU, lensV;
    const float discRadius = LensTable.front().ap * 0.5;
    Vector Z(0.f, 0.f, -1.f);

    //Every pixel pair
    for (unsigned int i = 1; i < LensTable.size()+1; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            for (int nRays = 0; nRays < RaysPerPair; ++ nRays) {
                //Can't reflect with itself, and it also need a surface
                if (LensTable[i].radius == 0.f || LensTable[j].radius == 0.f)
                    continue;

                //Generate ray
                float u = (float) rand() / (float) RAND_MAX;
                float v = (float) rand() / (float) RAND_MAX;
                ConcentricSampleDisk(u, v, &lensU, &lensV);
                Point discPos(lensU*discRadius,lensV*discRadius, discZ);
                Ray ray = Ray(pointLight,
                              Normalize(discPos-pointLight),
                              0.f,
                              INFINITY);

                LensFlareRay lfr[3];
                for (int RGB = 0; RGB < 3; ++ RGB) {
                    lfr[RGB].firstReflectLens = i;
                    lfr[RGB].secondReflectLens = j;
                    lfr[RGB].firstReflection = false;
                    lfr[RGB].ray = ray;
                }

                float rgb[3];
                pointSpectrum.ToRGB(rgb);

                float r[3] = { rgb[0], 0.f ,0.f };
                float g[3] = { 0.f, rgb[1] ,0.f };
                float b[3] = { 0.f, 0.f , rgb[2] };

                lfr[0].spec = pointSpectrum.FromRGB(r);
                lfr[1].spec = pointSpectrum.FromRGB(g);
                lfr[2].spec = pointSpectrum.FromRGB(b);

                //Traverse
                TrackLensFlare<0>(lfr[0]);
                TrackLensFlare<1>(lfr[1]);
                TrackLensFlare<2>(lfr[2]);
            }
        }
    }
    film->WriteImage(SplatFactor);
}

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	   // Extract common camera parameters from \use{ParamSet}
	   float hither = params.FindOneFloat("hither", -1);
	   float yon = params.FindOneFloat("yon", -1);
	   float shutteropen = params.FindOneFloat("shutteropen", -1);
	   float shutterclose = params.FindOneFloat("shutterclose", -1);

	   // Realistic camera-specific parameters
	   string specfile = params.FindOneString("specfile", "");
	   // about 70 mm default to film
	   float filmdistance = params.FindOneFloat("filmdistance", 70.0);
	   float fstop = params.FindOneFloat("aperture_diameter", 1.0);
	   float filmdiag = params.FindOneFloat("filmdiag", 35.0);
	   assert(hither != -1 && yon != -1 && shutteropen != -1 &&
	      shutterclose != -1 && filmdistance!= -1);
	   if (specfile == "") {
	       Severe( "No lens spec file supplied!\n" );
	   }
	   return new RealisticCamera(cam2world, hither, yon,
	      shutteropen, shutterclose, filmdistance, fstop,
	      specfile, filmdiag, film);
}

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
                                 float hither, float yon,
                                 float sopen, float sclose,
                                 float filmdistance, float aperture_diameter_,
                                 const string &specfile,
                                 float filmdiag,
								 Film *f)
                                 : Camera(cam2world, sopen, sclose, f),
								   ShutterOpen(sopen),
								   ShutterClose(sclose),
                                   FilmDistance(filmdistance),
                                   FilmDiagonal(filmdiag),
                                   Hither(hither),
								   film(f)

{
    Width = filmdiag * cos(M_PI/4);
    Height = filmdiag * sin(M_PI/4);
    WidthHalf  = Width * 0.5;
    HeightHalf = Height * 0.5;
    FilmZ = 0;

    ifstream datFile(specfile.c_str());
    if (datFile.is_open()) {
        string comments;
        while (datFile.peek() == '#') {
            getline(datFile,comments);
        };
        Lens read;
        cout << "Reading Lens data from " << specfile << endl;
        float z = 0.f;
        float nd_left[3] = { 1.f, 1.f, 1.f};

        int idx = 0;
        while (datFile >> read) {
            if (read.nd_right[1] == 0.f){
                read.nd_right[0] = 1.f;
                read.nd_right[1] = 1.f;
                read.nd_right[2] = 1.f;
            }
            read.zIntersect = z;
            read.idx = idx++;

            read.nd_right[0] = read.nd_right[1] - 0.02;
            read.nd_right[2] = read.nd_right[1] + 0.035;

            read.nd_left[0] = nd_left[0];
            read.nd_left[1] = nd_left[1];
            read.nd_left[2] = nd_left[2];

            nd_left[0] = read.nd_right[0];
            nd_left[1] = read.nd_right[1];
            nd_left[2] = read.nd_right[2];

            z-= read.thick;
            LensTable.push_back(read);
            FilmZ -= LensTable.back().thick;
            cout << read << endl;
        }
    } else {
        cerr << "Could not open " << specfile << endl;
        exit(1);
    }
    FilmZ -= FilmDistance;

    int nPairs = 0;
    for (unsigned int i = 1; i < LensTable.size(); ++i)
        for (unsigned int j = 0; j < i; ++j)
            ++nPairs;
    RaysPerPair = 1000000;
    SplatFactor = 1.f / (RaysPerPair * nPairs);

    //Get where on the Z axis the sample disc should be placed
    bool convex = LensTable.back().radius > 0;
    float r = fabs(LensTable.back().radius);
    float dz = r * cos(asin(LensTable.back().ap/(2.0f*r)));
    if (convex)
        DiscZ = FilmZ + FilmDistance - r + dz;
    else
        DiscZ = FilmZ + FilmDistance + r - dz;

#ifdef VDB
    vdb_frame();
    Point p0(-WidthHalf, -HeightHalf, FilmZ);
    Point p1( WidthHalf, -HeightHalf, FilmZ);
    Point p3(-WidthHalf,  HeightHalf, FilmZ);
    Point p2( WidthHalf,  HeightHalf, FilmZ);
    vdb_square(p0, p1, p2, p3);
    vdb_color(1.0f, 0.f, 1.f);
    vdb_circle(LensTable.back().ap * 0.5, DiscZ, 40);
    vdb_color(1.0f, 1.f, 0.f);
    vdb_line(0.f, 0.f, FilmZ, 0.f, 0.f, 2.f);
    vdb_color(1.0f, 0.f, 0.f);
    float t = 0.f;
    for (unsigned int i = 0; i < LensTable.size(); ++i) {
        if (LensTable[i].radius != 0.f) {
            vdb_sphere(fabs(LensTable[i].radius), t - LensTable[i].radius,
                    LensTable[i].ap, 90, LensTable[i].radius > 0);
        } else {
            /*const int rings = 10;
            vdb_color(0.3f, 0.3f, 0.3f);
            for (int k = 0; k < rings; ++k)
                vdb_circle(LensTable[i].ap*(0.5+(float)k/float(3*rings)),t, 40);
            vdb_color(1.f, 0.f, 0.f);*/

            float r = LensTable[i].ap*0.5;
            float a = (2.f * r) / (2 + sqrt(2));
            float c = sqrt(2) * a;


            Point vertices[8];
            vertices[0] = Point(-c*0.5, r, t);
            vertices[1] = Point(c*0.5, r, t);
            vertices[2] = Point(r, c*0.5, t);
            vertices[3] = Point(r, -c*0.5, t);
            vertices[4] = Point(c*0.5, -r, t);
            vertices[5] = Point(-c*0.5, -r, t);
            vertices[6] = Point(-r, -c*0.5, t);
            vertices[7] = Point(-r, c*0.5, t);

            vdb_color(0.3f, 0.3f, 0.3f);
            for (int i = 0; i < 8; ++i) {
                vdb_line(vertices[i].x, vertices[i].y, vertices[i].z,
                         vertices[(i+1) % 8].x, vertices[(i+1) % 8].y, vertices[(i+1) % 8].z);
            }
            vdb_color(1.f, 0.f, 0.f);
        }
        t -= LensTable[i].thick;
    }

#endif
}


RealisticCamera::~RealisticCamera()
{

}


float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
    int x_start, x_end, y_start, y_end;
    film->GetPixelExtent(&x_start, &x_end, &y_start, &y_end);
    const Point filmplanePos(WidthHalf  -sample.imageX/float(x_end)*Width,
                             -HeightHalf +sample.imageY/float(y_end)*Height,
                             FilmZ);

    //Get points on the disc infront of the first lens
    float lensU, lensV;
    ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
    const float discRadius = LensTable.back().ap * 0.5;
    Point discPos(lensU*discRadius,lensV*discRadius, DiscZ);

    *ray = Ray(filmplanePos, Normalize(discPos-filmplanePos),0.f, INFINITY);
	bool render = false;
#ifdef VDB
    char * env = getenv ("VDB_RAYS_LIMIT");
    render = rand() % atoi(env) >= atoi(env)-1;
#endif

    //Traverse ray
    for (int i = LensTable.size()-1; i >= 0; --i){
        if (!LensTable[i].Refract<1>(ray,render)) return 0.f;
    }

#ifdef VDB
        //Ray to scene
        if (render) {
            vdb_color(1.f, 1.f, 1.f);
            Point scene = (*ray)(10);
            vdb_line(ray->o.x, ray->o.y, ray->o.z, scene.x, scene.y, scene.z);
        }
#endif

    // GenerateRay() should return the weight of the generated ray
    CameraToWorld(*ray, ray);
    ray->d = Normalize(ray->d);

    Vector Z(0.f, 0.f, 1.f);
    Vector ToDisc = Normalize(discPos - filmplanePos);
    float costheta = Dot(Z,ToDisc);
    float costheta4 = costheta*costheta*costheta*costheta;
    float A_over_Z2 = (discRadius*discRadius*M_PI) /(FilmDistance*FilmDistance);
    return costheta4 * A_over_Z2;
}
