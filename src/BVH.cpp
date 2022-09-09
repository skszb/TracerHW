#include "BVH.h"

Node* BuildTree(const std::vector<std::pair<Surface *, Fill> >& surfaces) {
    assert(!surfaces.empty());
    SlVector3 min, max;
    for (int i = 1; i < surfaces.size(); i++) {

    }
}

Node* buildNode(std::vector<std::pair<Surface *, Fill> >* surfaces, int leftIdx, int rightIdx) {

}