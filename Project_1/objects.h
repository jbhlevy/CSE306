//
//  objects.h
//  Project_1
//
//  Created by John Levy on 07/05/2023.
//

#ifndef objects_h
#define objects_h

#include "triangle_indices.h"

#include <list>
#include "vector.h"
#include "ray.h"
#include "boxes.h"

namespace Objects {

class Object{
public:
    Vector Albedo;
    bool is_mirror;
    bool is_transparent;
    bool is_hollow;
    
    virtual bool const intersect(const Ray& ray, Vector& P, Vector& N, double& t) = 0;
    
    Object(){};
    Object(Vector const &albedo, bool mirror=false, bool transparent=false, bool hollow=false): Albedo(albedo), is_mirror(mirror), is_transparent(transparent), is_hollow(hollow){};
};

class Sphere : public Objects::Object{
public:
    Vector C;
    double r;
    
    Sphere(const Vector& center, double radius, const Vector &albedo,  bool mirror = false, bool transparency = false, bool hollow = false): ::Objects::Object(albedo, mirror, transparency, hollow), C(center), r(radius) {};
    
    const bool intersect(const Ray& ray, Vector& P, Vector &N, double& t){
        double delta;
        double positive_root;
        double negative_root;
        
        delta = sqr(dot(ray.dir, ray.O - C)) - ((ray.O - C).norm2() - r*r);
                
        if(delta >= 0){
            positive_root =  dot(ray.dir, C - ray.O) - sqrt(delta);
            negative_root = dot(ray.dir, C - ray.O) + sqrt(delta);
            
            if(negative_root < 0)
                return false;
            
            t = (positive_root > 0) ? positive_root : negative_root;
            P = ray.O + (t * ray.dir);
            N = P - C;
            N.normalize();
            return true;
        }
        return (delta >= 0);
    }
};

class TriangleMesh : public Objects::Object{
public:
  ~TriangleMesh() {}
    TriangleMesh() {};
    TriangleMesh(const Vector & rho, bool is_mirror=false, bool is_transparent=false, bool is_hollow=false):
    ::Objects::Object(rho, is_mirror, is_transparent, is_hollow){};
    
    bool const intersect(const Ray& r, Vector& P, Vector &N, double& t){
        double inter_distance;
        double tmp;
        bool has_inter=false;
        
        if(!root->bbox.intersect(r, tmp))
            return false;

        std::list<BVH*> nodes_to_visit;
        nodes_to_visit.push_front(root);

        t = std::numeric_limits<double>::max();
        while (!nodes_to_visit.empty()){
            BVH* curNode = nodes_to_visit.back();
            nodes_to_visit.pop_back();

            if(curNode->left){
                if(curNode->left->bbox.intersect(r, inter_distance)){
                    if(inter_distance < t){
                        nodes_to_visit.push_back(curNode->left);
                    }
                }
                if(curNode->right->bbox.intersect(r, inter_distance)){
                    if(inter_distance < t){
                        nodes_to_visit.push_back(curNode->right);
                    }
                }
            }
            else{
                for(size_t i = curNode->starting_triangle; i < curNode->ending_triangle; i++){
                    
                    const Vector& A = vertices[indices[i].vtxi];
                    const Vector& B = vertices[indices[i].vtxj];
                    const Vector& C = vertices[indices[i].vtxk];
                    
                    const Vector& NA = normals[indices[i].ni];
                    const Vector& NB = normals[indices[i].nj];
                    const Vector& NC = normals[indices[i].nk];

                    Vector local_P, local_N;
                    double local_t;
                    
                    bool local_has_inter = intersect_triangle(r, A, B, C, NA, NB, NC, local_P, local_N, local_t);
                    if(local_has_inter && local_t < t){
                        t = local_t;
                        has_inter = true;
                        N = local_N;
                        P = local_P;
                    }
                }
            }
        }
        return has_inter;
    }
    
    void scale_and_translate(double scalling, const Vector& translation){
        for(int i = 0; i< vertices.size(); i++){
            vertices[i] = vertices[i] * scalling + translation;
        }
    }
    
    BoundingBox compute_bounding_box(int starting_triangle, size_t ending_triangle){
        
        BoundingBox box = BoundingBox();
        for (int i = starting_triangle ; i < ending_triangle; i++){
            for (int j = 0; j<3; j++){
                box.m[j] = std::min(box.m[j], vertices[i][j]);
                box.M[j] = std::max(box.M[j], vertices[i][j]);
            }
        }
        return box;
    }
    
    void compute_bvh(BVH* node, int starting_triangle, size_t ending_triangle){
        //Computation of current node bounding box
        node->bbox = compute_bounding_box(starting_triangle, ending_triangle);
        
        //Diagonal and longest axis
        Vector diag = node->bbox.M - node->bbox.m;
        Vector middle_diag = node->bbox.m + diag * 0.5;
        int longest_axis = -1;
        if(diag[0] > diag[1] && diag[0] > diag[2]){
            longest_axis = 0;
        }
        else if (diag[1] > diag[0] && diag[1] > diag[2]){
            longest_axis = 1;
        }
        else if(diag[2] > diag[0] && diag[2] > diag[1]){
            longest_axis = 2;
        }
        int pivot = starting_triangle;
        //Dividing triangles in two sets
        for (size_t i = starting_triangle; i<ending_triangle; i++){
            //Compute barycenter
            const Vector& A = vertices[indices[i].vtxi];
            const Vector& B = vertices[indices[i].vtxj];
            const Vector& C = vertices[indices[i].vtxk];

            Vector Barycenter = findBarycenter(A, B, C);
            
            //Longest axis criterion
            if(Barycenter[longest_axis] < middle_diag[longest_axis]){
                std::swap(indices[i], indices[pivot]);
                pivot++;
            }
        }
        //Exiting condition
        if(pivot <= starting_triangle || ending_triangle-pivot <=1 || ending_triangle - starting_triangle <5){
            return;
        }
        //Recursive call
        node->left = new BVH(starting_triangle, pivot);
        node->right = new BVH(pivot, ending_triangle);
        compute_bvh(node->left, starting_triangle, pivot);
        compute_bvh(node->right, pivot, ending_triangle);
        return;
    }
    
    
    void readOBJ(const char* obj) {
 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
 
            }
 
        }
        fclose(f);
    }
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    //BoundingBox box;
    BVH *root;
};


}

#endif /* objects_h */
