#include "trace.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <cfloat>
#include <vector>
#include <chrono>
#include <cassert>
#include <stack>

#ifdef __APPLE__
#define MAX std::numeric_limits<double>::max()
#else

#include "values.h"
#include "Shading.h"

#define MAX DBL_MAX
#endif

// return the determinant of the matrix with columns a, b, c.
double det(const SlVector3 &a, const SlVector3 &b, const SlVector3 &c) {
    return a[0] * (b[1] * c[2] - c[1] * b[2]) +
           b[0] * (c[1] * a[2] - a[1] * c[2]) +
           c[0] * (a[1] * b[2] - b[1] * a[2]);
}

inline double sqr(double x) { return x * x; }

Node* BuildNode(const std::vector<std::pair<Surface *, Fill> >& surfaces) {
    // Create a node with a bounding box and all surfaces within it, then recursively divide the bounding box
    // in two of the half size as the child nodes.
    if (surfaces.empty()) return nullptr;
    Node *node = new Node();

    // surfaces
    node->sfs = surfaces;

    SlVector3 minCenter = surfaces[0].first->center();
    SlVector3 maxCenter = minCenter;
    // bounding box
    AABB *box = new AABB();
    box->minCorner = surfaces[0].first->minPos();
    box->maxCorner = surfaces[0].first->maxPos();
    for (int i = 1; i < surfaces.size(); i++) {
        Surface *sf = surfaces[i].first;
        minCenter = min(minCenter, sf->center());
        maxCenter = max(maxCenter, sf->center());
        box->minCorner = min(box->minCorner, sf->minPos());
        box->maxCorner = max(box->maxCorner, sf->maxPos());
    }
    node->box = box;

    // stop recursion if is a leaf node
    if (surfaces.size() == 1) {
        return node;
    }

    // left node and right node
    std::vector<std::pair<Surface *, Fill> > lSfs;
    std::vector<std::pair<Surface *, Fill> >rSfs;
    SlVector3 dif = maxCenter - minCenter;
    int longestAxis = 0;
    for (int i = 1 ; i < 3; i++) {
        if (dif[i] > dif[longestAxis])
            longestAxis = i;
    }
    double mid = (minCenter[longestAxis] + maxCenter[longestAxis]) / 2;
    for (const auto & p : surfaces) {
        if (p.first->center()[longestAxis] < mid) {
            lSfs.push_back(p);
        } else {
            rSfs.push_back(p);
        }
    }
    assert(!lSfs.empty() && !rSfs.empty());
    node->lNode = BuildNode(lSfs);
    node->rNode = BuildNode(rSfs);

    return node;
}

bool AABB::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    this->minCorner;
    this->maxCorner;
    // x
    double tmin = (minCorner.x() - r.e.x()) / r.d.x();
    double tmax = (maxCorner.x() - r.e.x()) / r.d.x();
    if (tmin > tmax) std::swap(tmin, tmax);

    // y
    double tymin = (minCorner.y() - r.e.y()) / r.d.y();
    double tymax = (maxCorner.y() - r.e.y()) / r.d.y();
    if (tymin > tymax) std::swap(tymin, tymax);
    if (tymax < tmin || tymin > tmax) return false;
    tmin = std::min(tmin, tymin);
    tmax = std::max(tmax, tymax);

    // z
    double tzmin = (minCorner.z() - r.e.z()) / r.d.z();
    double tzmax = (maxCorner.z() - r.e.z()) / r.d.z();
    if (tzmin > tzmax) std::swap(tzmin, tzmax);
    if (tzmax < tmin || tzmin > tmax) return false;
    tmin = std::min(tmin, tzmin);
    tmax = std::max(tmax, tzmax);

    if (tmax < t0 || tmin > t1) return false;
    hr.t = tmin;
    return true;
}


bool Triangle::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    // Step 1 Ray-triangle test
    // (V2-V1)*beta + (V3-V1)*gamma - D*t = A - V1
    SlVector3 V2minusV1 = this->b - this->a;
    SlVector3 V3minusV1 = this->c - this->a;
    SlVector3 AminusV1 = r.e - this->a;

    // solve system of equations using matrix
    double adj[3][3];
    double det;
    double beta{}, gamma{}, t{};

    adj[0][0] = V3minusV1[1] * -r.d[2] - -r.d[1] * V3minusV1[2];
    adj[1][0] = -r.d[1] * V2minusV1[2] - V2minusV1[1] * -r.d[2];
    adj[2][0] = V2minusV1[1] * V3minusV1[2] - V3minusV1[1] * V2minusV1[2];

    det = V2minusV1[0] * adj[0][0] + V3minusV1[0] * adj[1][0] + -r.d[0] * adj[2][0];

    if (det == 0) {
        return false;
    }

    adj[0][1] = -r.d[0] * V3minusV1[2] - V3minusV1[0] * -r.d[2];
    adj[1][1] = V2minusV1[0] * -r.d[2] - V2minusV1[2] * -r.d[0];
    adj[2][1] = V3minusV1[0] * V2minusV1[2] - V2minusV1[0] * V3minusV1[2];
    adj[0][2] = V3minusV1[0] * -r.d[1] - -r.d[0] * V3minusV1[1];
    adj[1][2] = -r.d[0] * V2minusV1[1] - V2minusV1[0] * -r.d[1];
    adj[2][2] = V2minusV1[0] * V3minusV1[1] - V3minusV1[0] * V2minusV1[1];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            adj[i][j] /= det;
        }
    }

    for (int i = 0; i < 3; i++) {
        beta += adj[0][i] * AminusV1[i];
        gamma += adj[1][i] * AminusV1[i];
        t += adj[2][i] * AminusV1[i];
    }

    if (gamma < 0 || beta < 0 || beta + gamma > 1) return false;
    if (t0 < t && t < t1) {
        // if hit, update the hit record
        hr.beta = beta;
        hr.gamma = gamma;
        hr.alpha = 1 - beta - gamma;
        hr.t = t;
        hr.p = r.e + r.d * hr.t;
        hr.n = cross(V2minusV1, V3minusV1);
        normalize(hr.n);
        return true;
    }
    return false;
}

SlVector3 Triangle::minPos() {
    return {
            std::min(std::min(a.x(), b.x()), c.x()),
            std::min(std::min(a.y(), b.y()), c.y()),
            std::min(std::min(a.z(), b.z()), c.z())
    };
}

SlVector3 Triangle::maxPos() {
    return {
            std::max(std::max(a.x(), b.x()), c.x()),
            std::max(std::max(a.y(), b.y()), c.y()),
            std::max(std::max(a.z(), b.z()), c.z())
    };
}

SlVector3 Triangle::center() {
    return (a + b + c) / 3;
}


bool TrianglePatch::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    bool temp = Triangle::intersect(r, t0, t1, hr);
    if (temp) {
        hr.n = hr.alpha * n1 + hr.beta * n2 + hr.gamma * n3;
        normalize(hr.n);
        return true;
    }
    return false;
}


bool Sphere::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    // Step 1 Ray-sphere test
    // D^2t^2 + 2ADt - 2DCt + (A-C)^2 - r^2
    SlVector3 EminusC = r.e - this->c;
    double a = dot(r.d, r.d);
    double b = 2 * dot(r.d, EminusC);
    double c = dot(EminusC, EminusC) - sqr(this->rad);
    double delta = sqr(b) - 4 * a * c;
    if (delta < 0) return false;
    double sqrtDelta = sqrt(delta);

    double tA = (-b - sqrtDelta) / (2 * a);
    double tB = (-b + sqrtDelta) / (2 * a);
    double tLess, tGreater;
    if (tA < tB) {
        tLess = tA;
        tGreater = tB;
    } else {
        tLess = tB;
        tGreater = tA;
    }

    // check hit record and save the nearest valid hit(update hr)
    if (t0 < tLess && tLess < t1) {
        hr.t = tLess;
        hr.p = r.e + r.d * hr.t;
        hr.n = hr.p - this->c;
        normalize(hr.n);
        hr.v = r.e - hr.p;
        normalize(hr.v);
        return true;
    }
    if (t0 < tGreater && tGreater < t1) {
        hr.t = tGreater;
        hr.p = r.e + r.d * hr.t;
        hr.n = hr.p - this->c;
        normalize(hr.n);
        hr.v = r.e - hr.p;
        normalize(hr.v);
        return true;
    }
    return false;
}

SlVector3 Sphere::minPos() {
    return {
            c.x() - rad,
            c.y() - rad,
            c.z() - rad
    };
}

SlVector3 Sphere::maxPos() {
    return {
            c.x() + rad,
            c.y() + rad,
            c.z() + rad
    };
}

SlVector3 Sphere::center() {
    return this->c;
}



Tracer::Tracer(const std::string &fname) {
    std::ifstream in(fname.c_str(), std::ios_base::in);
    std::string line;
    char ch;
    Fill fill;
    bool coloredlights = false;
    while (in) {
        if (in.eof()) break;
        getline(in, line);
        switch (line[0]) {
            case 'b': {
                std::stringstream ss(line);
                ss >> ch >> bcolor[0] >> bcolor[1] >> bcolor[2];
                break;
            }

            case 'v': {
                getline(in, line);
                std::string junk;
                std::stringstream fromss(line);
                fromss >> junk >> eye[0] >> eye[1] >> eye[2];

                getline(in, line);
                std::stringstream atss(line);
                atss >> junk >> at[0] >> at[1] >> at[2];

                getline(in, line);
                std::stringstream upss(line);
                upss >> junk >> up[0] >> up[1] >> up[2];

                getline(in, line);
                std::stringstream angless(line);
                angless >> junk >> angle;

                getline(in, line);
                std::stringstream hitherss(line);
                hitherss >> junk >> hither;

                getline(in, line);
                std::stringstream resolutionss(line);
                resolutionss >> junk >> res[0] >> res[1];
                break;
            }

            case 'p': {
                bool patch = false;
                std::stringstream ssn(line);
                unsigned int nverts;
                if (line[1] == 'p') {
                    patch = true;
                    ssn >> ch;
                }
                ssn >> ch >> nverts;
                std::vector<SlVector3> vertices;
                std::vector<SlVector3> normals;
                for (unsigned int i = 0; i < nverts; i++) {
                    getline(in, line);
                    std::stringstream ss(line);
                    SlVector3 v, n;
                    if (patch) ss >> v[0] >> v[1] >> v[2] >> n[0] >> n[1] >> n[2];
                    else ss >> v[0] >> v[1] >> v[2];
                    vertices.push_back(v);
                    normals.push_back(n);
                }
                bool makeTriangles = false;
                if (vertices.size() == 3) {
                    if (patch) {
                        surfaces.push_back(
                                std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2],
                                                                             normals[0], normals[1], normals[2]),
                                                           fill));
                    } else {
                        surfaces.push_back(
                                std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                    }
                } else if (vertices.size() == 4) {
                    SlVector3 n0 = cross(vertices[1] - vertices[0], vertices[2] - vertices[0]);
                    SlVector3 n1 = cross(vertices[2] - vertices[1], vertices[3] - vertices[1]);
                    SlVector3 n2 = cross(vertices[3] - vertices[2], vertices[0] - vertices[2]);
                    SlVector3 n3 = cross(vertices[0] - vertices[3], vertices[1] - vertices[3]);
                    if (dot(n0, n1) > 0 && dot(n0, n2) > 0 && dot(n0, n3) > 0) {
                        makeTriangles = true;
                        if (patch) {
                            surfaces.push_back(
                                    std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2],
                                                                                 normals[0], normals[1], normals[2]),
                                                               fill));
                            surfaces.push_back(
                                    std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[2], vertices[3],
                                                                                 normals[0], normals[2], normals[3]),
                                                               fill));
                        } else {
                            surfaces.push_back(
                                    std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]),
                                                               fill));
                            surfaces.push_back(
                                    std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[2], vertices[3]),
                                                               fill));
                        }
                    }
                    if (!makeTriangles) {
                        std::cerr << "I didn't make triangles.  Poly not flat or more than quad.\n";
                    }
                }
                break;
            }

            case 's' : {
                std::stringstream ss(line);
                SlVector3 c;
                double r;
                ss >> ch >> c[0] >> c[1] >> c[2] >> r;
                surfaces.push_back(std::pair<Surface *, Fill>(new Sphere(c, r), fill));
                break;
            }
            case 'f' : {
                std::stringstream ss(line);
                ss >> ch >> fill.color[0] >> fill.color[1] >> fill.color[2] >> fill.kd >> fill.ks >> fill.shine
                   >> fill.t >> fill.ior;
                break;
            }

            case 'l' : {
                std::stringstream ss(line);
                Light l;
                ss >> ch >> l.p[0] >> l.p[1] >> l.p[2];
                if (!ss.eof()) {
                    ss >> l.c[0] >> l.c[1] >> l.c[2];
                    coloredlights = true;
                }
                lights.push_back(l);
                break;
            }

            default:
                break;
        }
    }
    if (!coloredlights) for (unsigned int i = 0; i < lights.size(); i++) lights[i].c = 1.0 / sqrt(lights.size());
    im = new SlVector3[res[0] * res[1]];
    shadowbias = 1e-6;
    samples = 1;
    aperture = 0.0;
}


Tracer::~Tracer() {
    if (im) delete[] im;
    for (unsigned int i = 0; i < surfaces.size(); i++) delete surfaces[i].first;
    FreeNode(this->BVHTreeRoot);
}


SlVector3 Tracer::shade(HitRecord &hr) const {
    if (color) return hr.f.color;
    if (hr.raydepth > maxraydepth) return {};

    SlVector3 color(0.0);
    HitRecord dummy;

    for (unsigned int i = 0; i < lights.size(); i++) {
        const Light &light = lights[i];
        bool shadow = false;

        // Step 3 Check for shadows here
        SlVector3 SrftoLT = light.p - hr.p;
        Ray shadowRay = Ray(hr.p, SrftoLT);
        normalize(shadowRay.d);
        shadow = BVHIntersection(BVHTreeRoot, shadowRay, shadowbias, mag(SrftoLT), dummy,
                                 true);

        if (!shadow) {
            // Step 2 do shading here
            color += Shading::phong(light, hr);
        }
    }

    // Step 4 Add code for computing reflection color here
    hr.raydepth += 1;
    Light refLight;
    SlVector3 refDir = -hr.v + 2.0 * hr.n * dot(hr.n, hr.v);
    Ray reflectRay = Ray(hr.p, refDir); // normalize?
    refLight.c = trace(reflectRay, shadowbias, MAX);
    refLight.p = hr.p + refDir;
    color += hr.f.ks * Shading::phong(refLight, hr);

    // Step 5 Add code for computing refraction color here

    return color;
}


SlVector3 Tracer::trace(const Ray &r, double t0, double t1) const {
    HitRecord hr{};
    hr.t = MAX;
    HitRecord debugHr; debugHr.t = MAX;
    SlVector3 color(bcolor);

    bool hit = false;

    // Step 1 See what a ray hits
    hit = BVHIntersection(BVHTreeRoot, r, t0, t1, hr);

    if (hit) {
        color = shade(hr);
    }
    return color;
}


void Tracer::traceImage() {
    // set up coordinate system
    SlVector3 w = eye - at;
    w /= mag(w);
    SlVector3 u = cross(up, w);
    normalize(u);
    SlVector3 v = cross(w, u);
    normalize(v);

    double d = mag(eye - at);
    double h = tan((M_PI / 180.0) * (angle / 2.0)) * d;
    double l = -h;
    double r = h;
    double b = h;
    double t = -h;

    SlVector3 *pixel = im;

    for (unsigned int j = 0; j < res[1]; j++) {
        for (unsigned int i = 0; i < res[0]; i++, pixel++) {

            SlVector3 result(0.0, 0.0, 0.0);

            for (int k = 0; k < samples; k++) {

                double rx = 1.1 * rand() / RAND_MAX;
                double ry = 1.1 * rand() / RAND_MAX;

                double x = l + (r - l) * (i + rx) / res[0];
                double y = b + (t - b) * (j + ry) / res[1];
                SlVector3 dir = -d * w + x * u + y * v;

                Ray r(eye, dir);
                normalize(r.d);

                result += trace(r, hither, MAX);

            }
            (*pixel) = result / samples;
        }
    }
}


void Tracer::writeImage(const std::string &fname) {
#ifdef __APPLE__
    std::ofstream out(fname, std::ios::out | std::ios::binary);
#else
    std::ofstream out(fname.c_str(), std::ios_base::binary);
#endif
    out << "P6" << "\n" << res[0] << " " << res[1] << "\n" << 255 << "\n";
    SlVector3 *pixel = im;
    char val;
    for (unsigned int i = 0; i < res[0] * res[1]; i++, pixel++) {
        val = (unsigned char) (std::min(1.0, std::max(0.0, (*pixel)[0])) * 255.0);
        out.write(&val, sizeof(unsigned char));
        val = (unsigned char) (std::min(1.0, std::max(0.0, (*pixel)[1])) * 255.0);
        out.write(&val, sizeof(unsigned char));
        val = (unsigned char) (std::min(1.0, std::max(0.0, (*pixel)[2])) * 255.0);
        out.write(&val, sizeof(unsigned char));
    }
    out.close();
}

void Tracer::buildBVHTree() {
    BVHTreeRoot = BuildNode(surfaces);
}


bool Tracer::BVHIntersection(Node *nd, const Ray &r, double t0, double t1, HitRecord &hr, bool shadowTest) const {
    bool hit = false;
    HitRecord dummy;
    std::stack<Node*> callStack;
    callStack.push(nd);

    while (!callStack.empty()) {
        if (shadowTest && hit) return true;
        Node *currNode = callStack.top();
        callStack.pop();
        if (currNode->box->intersect(r, t0, t1, dummy)) {
            // in leaf node, check intersection with surface
            if (currNode->sfs.size() == 1) {
                auto p = currNode->sfs[0];
                if (p.first->intersect(r, t0, t1, hr)) {
                    hit = true;
                    t1 = hr.t;
                    hr.f = p.second;
                    hr.v = r.e - hr.p;
                    normalize(hr.v);
                }
            }
            if (currNode->rNode) callStack.push(currNode->rNode);
            if (currNode->lNode) callStack.push(currNode->lNode);
        }
    }
    return hit;
}




int main(int argc, char *argv[]) {
    int c;
    double aperture = 0.0;
    int samples = 1;
    int maxraydepth = 5;
    bool color = false;
    while ((c = getopt(argc, argv, "a:s:d:c")) != -1) {
        switch (c) {
            case 'a':
                aperture = atof(optarg);
                break;
            case 's':
                samples = atoi(optarg);
                break;
            case 'c':
                color = true;
                break;
            case 'd':
                maxraydepth = atoi(optarg);
                break;
            default:
                abort();
        }
    }
    if (argc - optind != 2) {
        std::cout << "usage: trace [opts] input.nff output.ppm" << std::endl;
        for (unsigned int i = 0; i < argc; i++) std::cout << argv[i] << std::endl;
        exit(0);
    }
    Tracer tracer(argv[optind++]);
    tracer.aperture = aperture;
    tracer.samples = samples;
    tracer.color = color;
    tracer.maxraydepth = maxraydepth;

    auto beginOfBuildingTree = std::chrono::steady_clock::now();
    tracer.buildBVHTree();
    auto endOfBuildingTree = std::chrono::steady_clock::now();
    std::cout << "BVH Tree Creating Time:"<< std::chrono::duration_cast<std::chrono::seconds>(endOfBuildingTree - beginOfBuildingTree).count() << "s\n";

    auto beginOfRendering = std::chrono::steady_clock::now();
    tracer.traceImage();
    auto endOfRendering = std::chrono::steady_clock::now();
    std::cout << "Rendering Time:"<< std::chrono::duration_cast<std::chrono::seconds>(endOfRendering - beginOfRendering).count() << "s\n";

    tracer.writeImage(argv[optind++]);
};