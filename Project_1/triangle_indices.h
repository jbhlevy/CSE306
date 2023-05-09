//
//  triangle_indices.h
//  Project_1
//
//  Created by John Levy on 09/05/2023.
//

#ifndef triangle_indices_h
#define triangle_indices_h

#include "vector.h"

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

    bool intersect_triangle(const Ray& ray, const Vector &A, const Vector &B, const Vector &C,  const Vector &NA, const Vector &NB, const Vector &NC, Vector &P, Vector &N, double &t){
        
        double alpha; 

        Vector AB = B - A;
        Vector AC = C - A;
        N = cross(AB, AC);
        double _UN = 1 / dot(ray.dir, N);
        
        Vector _P = cross(A-ray.O, ray.dir);
        double gamma = -dot(AB, _P) * _UN;
        double beta = dot(AC, _P) * _UN;
        t = dot(A - ray.O, N) * _UN;
        
        if(t<0){
            return false;
        }
        if(beta < 0 || beta > 1){
            return false;
        }
        if(gamma < 0 || gamma > 1){
            return false;
        }
        alpha = 1 - beta - gamma;
        if (alpha < 0){
            return false;
        }
        P = ray.O + (t *  ray.dir); // origini + t*direction
        N = alpha * NA + beta * NB + gamma * NC;
        N.normalize();
        return true;
    }

#endif /* triangle_indices_h */
