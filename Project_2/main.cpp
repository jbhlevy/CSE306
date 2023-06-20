//
//  main.cpp
//  TD2
//
//  Created by John Levy on 10/05/2023.
//

#include <iostream>
#include "vector.h"
#include "omp.h"
#include "lbfgs.h"

#include <sstream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define VOLUME_AIR 0.7
#define VOLUME_FLUID 0.3

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {
public:
    std::vector<Vector> vertices;
    
    double area() const{
        if(vertices.size() < 3)
            return 0;
        
        double A = 0;
        for (int i = 0; i < vertices.size(); i++) {
            A += vertices[i][0] * vertices[(i+1) %  vertices.size()][1] - vertices[i][1] * vertices[(i+1) % vertices.size()][0];
        }
        return std::fabs(A/2);
    }
    
    double integrate_square_distance(const Vector& P){
        if(vertices.size() < 3)
            return 0;
        
        
        double value = 0;
        for (int i = 1; i < vertices.size() - 1; i++) {
            double local_value = 0;
            Vector triangle[3] = {vertices[0], vertices[i], vertices[i+1]};
            for (int k = 0 ; k < 3; k++) {
                for (int l = k; l < 3; l++) {
                    local_value += dot(triangle[k] - P, triangle[l] - P);
                }
            }
            
            Vector e1 = triangle[1] - triangle[0];
            Vector e2 = triangle[2] - triangle[0];
            double area_triangle =  0.5 * std::fabs(e1[1] * e2[0] - e1[0] * e2[1]);
            
            value += area_triangle / 6. * local_value;
        }
        return value;
    }
    
    Vector centroid() const {
        Vector c = Vector(0., 0., 0.);
        int N = (int)vertices.size();
        double a =  area();
        for (int i = 0; i < N; i++) {
            int ip1 = (i+1) % N;
            c[0] += (vertices[i][0] + vertices[ip1][0]) * (vertices[i][0] * vertices[ip1][1] - vertices[ip1][0] * vertices[i][1]);
            c[1] += (vertices[i][1] + vertices[ip1][1]) * (vertices[i][0] * vertices[ip1][1] - vertices[ip1][0] * vertices[i][1]);
        }
        return c / (6. * a);
    }
}; // End of Polyogn Class
 
// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
        FILE* f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
            fprintf(f, "<g>\n");
            fprintf(f, "<polygon points = \"");
            for (int j = 0; j < polygons[i].vertices.size(); j++) {
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
            fprintf(f, "</g>\n");
        }
        fprintf(f, "</svg>\n");
        fclose(f);
    }

int sgn(double x){
    if(x < 0) return -1;
    if(x > 0) return 1;
    return 0;
}


void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W*H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W-1., std::max(0., W * bminx));
        bminy = std::min(H-1., std::max(0., H * bminy));
        bmaxx = std::max(W-1., std::max(0., W * bmaxx));
        bmaxy = std::max(H-1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
                    int sign = sgn(det);
                    if (prevSign == 0) prevSign = sign; else
                        if (sign == 0) sign = prevSign; else
                        if (sign != prevSign) {
                            isInside = false;
                            break;
                        }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                    double distEdge = std::abs(det)/ edgeLen;
                    double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                    if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {
                    {   // the N first particles may represent fluid, displayed in blue
                      image[((H - y - 1)*W + x) * 3] = 0;
                      image[((H - y - 1)*W + x) * 3 + 1] = 0;
                      image[((H - y - 1)*W + x) * 3 + 2] = 255;
                    }
                    if (mindistEdge <= 2) {
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 0;
                    }

                }
                
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}


class PowerDiagram{
public:
    std::vector<Vector> points;
    std::vector<double> weights;
    std::vector<Polygon> voronoi_diagram;
    Polygon disk;
    
    PowerDiagram(){
        const int N_disk = 30;
        disk.vertices.resize(N_disk);
        
        for (int i = 0; i < N_disk; i++) {
            double t = i /(double)N_disk * M_PI * 2;
            disk.vertices[i][0] = cos(t);
            disk.vertices[i][1] = sin(t);
            disk.vertices[i][2] = 0;
        }
        
    };
    
    PowerDiagram(const std::vector<Vector>& pts, const std::vector<double> &ws){
        points = pts;
        weights = ws;
        
        const int N_disk = 30;
        disk.vertices.resize(N_disk);
        
        for (int i = 0; i < N_disk; i++) {
            double t = i /(double)N_disk * M_PI * 2;
            disk.vertices[i][0] = cos(t);
            disk.vertices[i][1] = sin(t);
            disk.vertices[i][2] = 0;
        }
    }
    
    
    Polygon clip_by_bissector(const Polygon &poly, int index_0, int index_i, const Vector& P_0, const Vector& P_I) const{
        //Sutherland Hodgeman alg
        
        Polygon result;
        
        Vector M = (P_0 + P_I) /2;
        
        Vector Mprime = M + (weights[index_0]-weights[index_i])/(2*((P_0 - P_I).norm2())) * (P_I - P_0);
        
        for (int i = 0; i<poly.vertices.size(); i++) {
            //i corresponds to j in the code/slide conversion
            
            const Vector& A = (i == 0) ? poly.vertices.back() : poly.vertices[i - 1];
            
            const Vector& B = poly.vertices[i];
            
            double t = dot(Mprime-A, P_I-P_0) / dot(B - A, P_I - P_0);
            Vector P = A + t * (B -A);
            
            if ((B - P_0).norm2() - weights[index_0] < (B - P_I).norm2() - weights[index_i]) { // B is inside
                if((A - P_0).norm2() - weights[index_0]>= (A - P_I).norm2() - weights[index_i]){
                    //Add point of intersection to polygon
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else {
                if((A - P_0).norm2() - weights[index_0] < (A - P_I).norm2()- weights[index_i]){ // B outside and A inside
                    result.vertices.push_back(P); //Prof had else {if {}}
                }
            }
        }
        return result;
    }
    
    Polygon intersect_with_disk(const Polygon& polygon, const Vector& center, double radius) const {
        
        Polygon result(polygon);
        for (int i = 0; i < disk.vertices.size(); i++) {
            const Vector& u = disk.vertices[i] * radius + center;
            const Vector& v = disk.vertices[(i+1) % disk.vertices.size()] * radius + center;
            result = clip_by_edge(result, u, v);
        }
        return result;
    }
    
    Polygon clip_by_edge(const Polygon &poly, const Vector& u, const Vector& v) const{
        //Sutherland Hodgeman alg
        
        Polygon result;
        result.vertices.reserve(poly.vertices.size() + 1);
        Vector N(v[1] - u[1], -v[0] + u[0], 0);
        
        for (int i = 0; i<poly.vertices.size(); i++) {
            //i corresponds to j in the code/slide conversion
            
            const Vector& A = (i == 0) ? poly.vertices.back() /*[poly.vertices.size() - 1]*/ : poly.vertices[i - 1];
            
            const Vector& B = poly.vertices[i];
            
            double t = dot(u-A, N) / dot(B - A, N);
            Vector P = A + t * (B - A);
            
            if (dot(B - u, N) < 0) { // B is inside
                if(dot(A - u, N) > 0){ // A outside
                    //Add point of intersection to polygon
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else if(dot(A - u, N) < 0){ // B outside and A inside
                result.vertices.push_back(P); //Prof had else {if {}}
            }
        }
        return result;
    }
    
    
    Polygon compute_voronoi_cell(int index){
        
        Polygon polygon;
        
        polygon.vertices.resize(4);
        polygon.vertices[0] = Vector(0, 0, 0); // we compute the diagram restricted t unit square
        polygon.vertices[1] = Vector(1, 0, 0);
        polygon.vertices[2] = Vector(1, 1, 0);
        polygon.vertices[3] = Vector(0, 1, 0);
        
        for (int i=0; i<points.size(); i++) {
            if(index == i)
                continue;
            polygon = clip_by_bissector(polygon, index ,i, points[index], points[i]);
        }
        polygon = intersect_with_disk(polygon, points[index], sqrt(weights[index] - weights[weights.size() - 1]));
        return polygon;
    }
    
    void compute(){
        
        voronoi_diagram.resize(points.size());
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i< points.size(); i++) {
            voronoi_diagram[i] = compute_voronoi_cell(i);
        }
        
    }
    void save(std::string filename){
        save_svg(voronoi_diagram, filename, "blue");
    }
    
}; // End of Power diagram class

class OT{
public:
    std::vector<Vector> points;
    std::vector<double> lambdas;
    PowerDiagram solution;
    
    OT() {};
    
    OT(std::vector<Vector>& pts, const std::vector<double> &lambdas): points(pts), lambdas(lambdas) {};
    
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<OT*>(instance)->evaluate(x, g, n, step);
    }
    
    lbfgsfloatval_t evaluate(
              const lbfgsfloatval_t *x,
              lbfgsfloatval_t *g,
              const int n,
              const lbfgsfloatval_t step
              )
    {
        lbfgsfloatval_t fx = 0.0;
        for (int i=0; i<n; i++) {
            solution.weights[i] = x[i];
        }
        solution.compute();
        
        
        
        
        double s1 = 0;
        double s2 = 0;
        double s3 = 0;
        double estimated_volume_fluid = 0;
        for (int i = 0; i < n-1;i ++) {
            double cell_area = solution.voronoi_diagram[i].area();
            g[i] = - (lambdas[i] - cell_area);
            s1 += solution.voronoi_diagram[i].integrate_square_distance(solution.points[i]);
            s2 -= solution.weights[i] * cell_area;
            s3 += solution.weights[i] * lambdas[i];
            estimated_volume_fluid += cell_area;
        }
        fx = s1 + s2 + s3;
        
        double estimated_volume_air = 1.0 - estimated_volume_fluid;
        //energy
        fx += x[n-1]*(VOLUME_AIR - estimated_volume_air);
        g[n-1] = -(VOLUME_AIR - estimated_volume_air);
        
        return -fx;
    }

          static int _progress(
              void *instance,
              const lbfgsfloatval_t *x,
              const lbfgsfloatval_t *g,
              const lbfgsfloatval_t fx,
              const lbfgsfloatval_t xnorm,
              const lbfgsfloatval_t gnorm,
              const lbfgsfloatval_t step,
              int n,
              int k,
              int ls
              )
          {
              return reinterpret_cast<OT*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
          }

          int progress(
              const lbfgsfloatval_t *x,
              const lbfgsfloatval_t *g,
              const lbfgsfloatval_t fx,
              const lbfgsfloatval_t xnorm,
              const lbfgsfloatval_t gnorm,
              const lbfgsfloatval_t step,
              int n,
              int k,
              int ls
              )
          {
//              printf("Iteration %d:\n", k);
//              printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
//              printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
//              printf("\n");
              return 0;
          }
        
        
    
    void solve(){
        solution.points = points;
        solution.weights.resize(points.size()+1);
        std::fill(solution.weights.begin(), solution.weights.end(), 1.0);
        solution.weights[solution.weights.size() - 1] = 1 - 0.0001;
        
        double fx = 0;
        //LBFGS...
        int ret = lbfgs(points.size() + 1, &solution.weights[0], &fx, _evaluate, _progress, this, NULL);
        solution.compute();
    }
        

}; // End of OT Class

class Fluid{
public:
    OT otsolver;
    std::vector<Vector> particles;
    std::vector<Vector> velocities;
    
    Fluid(int N) {
        particles.resize(N);
        for (int i = 0; i < N; i++) {
            particles[i] = Vector(rand()/(double)RAND_MAX, rand()/(double)RAND_MAX);
        }
        velocities.resize(N, Vector(0., 0., 0.));
    };
    
    void stepFluid() {
        
        otsolver.points = particles;
        otsolver.lambdas = std::vector<double>(particles.size(), 1./(double)particles.size() * VOLUME_FLUID);
        otsolver.solve();
        
        const double mass_particle = 200;
        const double epsilon2 = sqr(0.004);
        const double dt = 0.002;
    
        for (int i = 0; i < particles.size(); i++) {
            
            Vector gravity = Vector(0, -9.81, 0) * mass_particle;
            Vector centroid = otsolver.solution.voronoi_diagram[i].centroid();
            Vector otForce = 1./epsilon2 * (centroid - particles[i]);
            Vector forces = gravity + otForce;
            velocities[i] += dt/mass_particle * forces;
            particles[i] += dt* velocities[i];
        }
    }
    
    void runFluid(){
        
        for (int i = 0; i < 1000 ; i++) {
            stepFluid();
            save_frame(otsolver.solution.voronoi_diagram, "animation", i);
        }
    }
};



int main(int argc, const char * argv[]) {

    std::cout << "Hello, World!\n";
    
    
    Fluid fluid(100);
    std::cout << "Running fluid simulation..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    fluid.runFluid();
    
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = finish - start;
    double time_duration = duration.count();
    std::cout << "Finished in " << time_duration << " seconds!" << std::endl;
    
    return 0;
}
