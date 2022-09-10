#ifndef TRACE_CPP_PHONG_H
#define TRACE_CPP_PHONG_H

#include "slVector.h"
#include "trace.h"
#include <cmath>

static const double ka = 0.1;

class Shading {
public:
    // phong shading
    static SlVector3 phong(const Light &pl, const HitRecord &hr) {
        SlVector3 color{};
        SlVector3 lDir = pl.p - hr.p;
        normalize(lDir);
        double cosTheta = dot(hr.n, lDir);
        // diffuse
        SlVector3 Kd = hr.f.kd * hr.f.color;
        color += Kd * pl.c * std::max(cosTheta, 0.0);

        // specular
        SlVector3 revLDir = -lDir;
        SlVector3 Ks = hr.f.ks * hr.f.color;

        SlVector3 refLDir = revLDir + 2.0 * cosTheta * hr.n;
        color += Ks * pl.c * pow(std::max(dot(refLDir, hr.v), 0.0), hr.f.shine);

        // ambient
        SlVector3 Ka = ka * hr.f.color;
        color += Ka * Kd * pl.c;

        return color;
    }
};


#endif //TRACE_CPP_PHONG_H
