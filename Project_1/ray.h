//
//  ray.h
//  Project_1
//
//  Created by John Levy on 07/05/2023.
//

#ifndef ray_h
#define ray_h

#include "vector.h"

class Ray{
public:
    Vector O;
    Vector dir;
    Ray(const Vector &origin, const Vector& direction): O(origin), dir(direction){};
};

#endif /* ray_h */
