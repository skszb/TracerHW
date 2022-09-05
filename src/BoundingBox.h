//
// Created by zb on 2022/8/31.
//

#ifndef TRACE_CPP_BOUNDINGBOX_H
#define TRACE_CPP_BOUNDINGBOX_H

#endif //TRACE_CPP_BOUNDINGBOX_H

#include "trace.h"
class BoundingBox {
    virtual bool intersect() const = 0;
    virtual ~BoundingBox() {};
};