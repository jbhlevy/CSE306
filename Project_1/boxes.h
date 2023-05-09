//
//  boxes.h
//  Project_1
//
//  Created by John Levy on 09/05/2023.
//

#ifndef boxes_h
#define boxes_h

#include "vector.h"
#include "ray.h"

class BoundingBox{
public:
    Vector m, M;
    
    BoundingBox(){
        m = Vector(1E10, 1E10, 1E10);
        M = Vector(-1E10, -1E10, -1E10);
    }
    
    bool intersect(const Ray& ray, double& tmin){
        double tx1 = (m[0] - ray.O[0]) / ray.dir[0];
        double tx2 = (M[0] - ray.O[0]) / ray.dir[0];
        
        double txmin = std::min(tx1, tx2);
        double txmax = std::max(tx1, tx2);
        
        double ty1 = (m[1] - ray.O[1]) / ray.dir[1];
        double ty2 = (M[1] - ray.O[1]) / ray.dir[1];
        
        double tymin = std::min(ty1, ty2);
        double tymax = std::max(ty1, ty2);
        
        double tz1 = (m[2] - ray.O[2]) / ray.dir[2];
        double tz2 = (M[2] - ray.O[2]) / ray.dir[2];
        
        double tzmin = std::min(tz1, tz2);
        double tzmax = std::max(tz1, tz2);
        
        tmin = std::max(txmin, std::max(tymin, tzmin));
        if (std::min(txmax, std::min(tymax, tzmax)) > std::max(txmin, std::max(tymin, tzmin))) {
            tmin = std::max(txmin, std::max(tymin, tzmin)); // updating inter_distance
            return true;
        }
        else{
            return false;
        }
    }
};

class BVH{
public:
    BoundingBox bbox;
    size_t starting_triangle, ending_triangle;
    BVH *left, *right;
    
    BVH(int start, size_t end){
        starting_triangle = start;
        ending_triangle = end;
        left = NULL;
        right = NULL;
    }
    
};

#endif /* boxes_h */
