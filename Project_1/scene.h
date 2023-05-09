//
//  scene.h
//  Project_1
//
//  Created by John Levy on 07/05/2023.
//

#ifndef scene_h
#define scene_h

#include "objects.h"
#include "random.h"

class Scene{
public:
    std::vector<Objects::Object*> objects;
    double Intensity;
    Vector L; // Light direction
    
    Scene(double i, const Vector& l): Intensity(i), L(l) {};
    
    ~Scene(){
        objects.clear();
    }
    
    void addSphere(Objects::Sphere *s){
        objects.push_back(s);
    }
    
    void addMesh(Objects::TriangleMesh *m){
        objects.push_back(m);
    }
    
    //Scene intersect: call the intersect function for each object in the scene
    bool const intersect(const Ray& ray, Vector& P, Vector& N, double& t, int& _id){
        Vector localP, localN;
        double localt = 0;
        bool has_inter = false;
        
        t = std::numeric_limits<double>::max();
        for (int i = 0; i < objects.size(); i++) {
            if(objects[i]->intersect(ray, localP, localN, localt)){
                if(localt < t){
                    t = localt;
                    P = localP;
                    N = localN;
                    _id = i;
                    has_inter = true;
                }
            }
        }
        return has_inter;
    }
    
    Vector get_color(const Ray &ray, int max_bounces){
        if(max_bounces <= 0)
            return Vector(0, 0, 0);
        
        Vector P, N;
        double t;
        int _id;
        Vector color(0, 0, 0);
        
        if(intersect(ray, P, N, t, _id)){ //2.Looking for intersections with objects
            //Hollowness
            if(objects[_id]->is_hollow){
                N = -N;
            }
            //Mirrors
            if(objects[_id]->is_mirror){
                Vector mirror_dir = ray.dir - 2 * dot(ray.dir, N) * N;
                Ray mirror_ray(P + 0.001*N, mirror_dir);
                return get_color(mirror_ray, max_bounces-1);
            }
            //Transparence
            if(objects[_id]->is_transparent){
                double n1 = 1;
                double n2 = 1.4;
                double rad;
                Vector N_transparency = N;
                Vector t_tangent, t_normal;
                
                if(dot(ray.dir, N) > 0){
                    std::swap(n1, n2);
                    N_transparency = -N_transparency;
                }
                t_tangent = (n1/n2) * (ray.dir - (dot(ray.dir, N_transparency)*N_transparency));
                rad = 1 - (sqr(n1/n2) * (1 - sqr(dot(ray.dir, N_transparency))));
                
                if(rad < 0){
                    Vector mirror_dir = ray.dir - 2*dot(ray.dir, N)*N;
                    Ray mirror_ray(P - 0.001*N, mirror_dir);
                    return get_color(mirror_ray, max_bounces - 1);
                }
                else{
                    t_normal = -sqrt(rad) * N_transparency;
                    Ray refracted_ray(P - N_transparency*0.001, t_tangent + t_normal);
                    return get_color(refracted_ray, max_bounces - 1);
                }
                
            }
            //Regular objects
            Vector shadowP, shadowN;
            double shadow_t;
            int shadow_id;
            bool in_shadow = false;
            
            double d2 = (L - P).norm2();
            
            Vector light_dir = (L - P);
            light_dir.normalize();
            Ray shadow_ray(P + 0.001*N, light_dir);
            
            if(intersect(shadow_ray, shadowP, shadowN, shadow_t, shadow_id)){
                if(shadow_t * shadow_t < d2){
                    in_shadow = true;
                }
            }
            if(!in_shadow){
                color = (Intensity/(4.* M_PI * d2)) * (objects[_id]->Albedo/M_PI) * std::max(0., dot(N, light_dir));
                //3. Computing the color of the object
            }
            //Adding indirect lighting
            int thread_id = omp_get_thread_num();
            Ray indirect_ray(P + 0.0001*N, Random::random(N, thread_id));
            color+= objects[_id]->Albedo * get_color(indirect_ray, max_bounces - 1);
            return color;
        }
        return color;
    }
};

#endif /* scene_h */
