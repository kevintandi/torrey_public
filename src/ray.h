#ifndef RAY_H
#define RAY_H

#include "vector.h"

class ray {
    public:
        ray() {}
        ray(const Vector3& origin, const Vector3& direction)
            : orig(origin), dir(direction)
        {}

        Vector3 origin() const  { return orig; }
        Vector3 direction() const { return dir; }

        Vector3 at(double t) const {
            return orig + t*dir;
        }

    public:
        Vector3 orig;
        Vector3 dir;
};

#endif