//
//  random.h
//  Project_1
//
//  Created by John Levy on 07/05/2023.
//

#ifndef random_h
#define random_h

#include <random>
#include "vector.h"

static std::default_random_engine engine[8] ;
static std::uniform_real_distribution<double> uniform(0, 1);

namespace Random{

std::vector<double>box_muller(int thread_id){
    double u1 = uniform(engine[thread_id]);
    double u2 = uniform(engine[thread_id]);
    
    double randomX = sqrt(-2 * log(u1)) * cos(2* M_PI * u2) * 0.5;
    double randomY = sqrt(-2 * log(u1)) * sin(2* M_PI * u2) * 0.5;
    
    std::vector<double> res;
    res.push_back(randomX);
    res.push_back(randomY);
    
    return res;
}

Vector random(const Vector& N, int thread_id){
    double n1 = uniform(engine[thread_id]);
    double n2 = uniform(engine[thread_id]);
    
    Vector T1;
    Vector T2;
    
    double x = sqrt(1 - n2) * cos(2 * M_PI * n1);
    double y = sqrt(1 - n2) * sin(2 * M_PI * n1);
    double z = sqrt(n2);
    
    if ((std::abs(N[0]) <= std::abs(N[1])) && std::abs(N[0]) <= std::abs(N[2])){
        T1 = Vector(0, N[2], -N[1]);
    }
    else{
        if ((std::abs(N[1]) <= std::abs(N[2])) && (std::abs(N[1]) <= std::abs(N[0]))){
            T1 = Vector(N[2], 0, -N[0]);
        }
        else{
            T1 = Vector(-N[1], N[0], 0);
        }
    }
    T1.normalize();
    T2 = cross(N, T1);
    
    return x * T1 + y * T2 + z * N;
}
}
#endif /* random_h */
