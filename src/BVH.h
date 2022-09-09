#include "trace.h"
#include <vector>
#include <cassert>

class BVH {
public:
    virtual bool intersect() const = 0;
    virtual ~BVH() = default;
};

class AABB : BVH {
public:
    SlVector3 min[3];
    SlVector3 max[3];
};

class Node {
public:
    BVH* box;
    Node* l;
    Node* r;
    std::vector<std::pair<Surface *, Fill> > surfaces;
    bool isLeaf;
};


Node* BuildTree(const std::vector<std::pair<Surface *, Fill> >& surfaces);

Node* buildNode(std::vector<std::pair<Surface *, Fill> >* surfaces, int leftIdx, int rightIdx);



