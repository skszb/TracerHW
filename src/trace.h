#ifndef TRACE_H
#define TRACE_H

#include "slVector.h"

#include <vector>

class Ray {
public:
    SlVector3 e;
    SlVector3 d;
    int depth;
    Ray(const SlVector3 &_e, const SlVector3 &_d, int _depth = 0) : e(_e), d(_d), depth(_depth) {};
};

// Material
class Fill {
public:
    SlVector3 color;
    // shine is for specular
    double kd, ks, shine, t, ior;
};

class HitRecord {
public:
    double t, alpha, beta, gamma;
    // inters pos, normal, view
    SlVector3 p, n, v;
    Fill f;
    int raydepth;
};

class Light {
public:
    // position and color
    SlVector3 p, c;
};

class Surface {
public:
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const = 0;
    virtual SlVector3 minPos() = 0;
    virtual SlVector3 maxPos() = 0;
    virtual SlVector3 center() = 0;
    virtual ~Surface() {};
};

class Triangle : public Surface {
    SlVector3 a,b,c;
public:
    Triangle(const SlVector3 &_a, const SlVector3 &_b, const SlVector3 &_c) : a(_a), b(_b), c(_c) {};
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
    SlVector3 minPos();
    SlVector3 maxPos();
    SlVector3 center();
};

class TrianglePatch : public Triangle {
    SlVector3 n1, n2, n3;
public:
    TrianglePatch(const SlVector3 &_a, const SlVector3 &_b, const SlVector3 &_c,
                  const SlVector3 &_n1, const SlVector3 &_n2, const SlVector3 &_n3)
            : Triangle(_a,_b,_c), n1(_n1), n2(_n2), n3(_n3) {};
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
};

class Sphere : public Surface {
    SlVector3 c;
    double rad;
public:
    Sphere(const SlVector3 &_c, double _r) : c(_c), rad(_r) {};
    bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
    SlVector3 minPos();
    SlVector3 maxPos();
    SlVector3 center();
};

class AABB {
public:
    AABB() = default;
    bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
    SlVector3 minCorner;
    SlVector3 maxCorner;

};

class Node {
public:
    Node() = default;
    Node(Node &nd) {
        box = nd.box;
        lNode = nd.lNode;
        rNode = nd.rNode;
        sfs = nd.sfs;
    }
    AABB *box;
    Node* lNode = nullptr;
    Node* rNode = nullptr;
    std::vector<std::pair<Surface *, Fill> > sfs;
};

Node* BuildNode(const std::vector<std::pair<Surface *, Fill> >& surfaces);

void FreeNode(Node* nd) {
    if (!nd) return;
    free(nd->box);
    FreeNode(nd->lNode);
    FreeNode(nd->rNode);
}


class Tracer {
    SlVector3 bcolor, eye, at, up;
    double angle, hither;
    unsigned int res[2];
    std::vector<std::pair<Surface *, Fill> > surfaces;
    std::vector<Light> lights;
    double shadowbias;

    SlVector3 *im;
public:
    Tracer(const std::string &fname);
    ~Tracer();
    void traceImage();
    SlVector3 trace(const Ray &ray, double t0, double t1) const;
    SlVector3 shade(HitRecord &hr) const;
    void writeImage(const std::string &fname);
    void buildBVHTree();
    bool BVHIntersection(Node* nd, const Ray &r, double t0, double t1, HitRecord &hr, bool shadowTest= false) const;

    bool color;
    int samples;
    double aperture;
    int maxraydepth;
    Node* BVHTreeRoot = nullptr;

private:
};

#endif
