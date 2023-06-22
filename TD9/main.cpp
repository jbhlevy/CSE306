#include <iostream>
#include <map>
#include <vector>

#include "vector.h"
#include "mesh_reader.hpp"

//Disclaimer: worked on this code with Philippe Guyard as it was not a graded assignement

class Edge {
public:
    int u, v;
    Edge(int u = -1, int v = -1): u(u), v(v) {}
};

bool operator<(const Edge& e1, const Edge& e2) {
    if (std::min(e1.u, e1.v) < std::min(e2.u, e2.v))
        return true; 
    
    return std::max(e1.u, e1.v) < std::max(e1.u, e2.v);
}


int main() {
    // Read goethe.obj with the mesh reader
    TriangleMeshDescriptor mesh;
    mesh.readOBJ("goethe.obj");

    std::map<Edge, std::vector<size_t>> edge_to_tri;
    std::map<int, std::vector<int>> vertex_to_tri;
    for(int m = 0; m < mesh.indices.size(); ++m) {
        int i = mesh.indices[m].vtxi;
        int j = mesh.indices[m].vtxj;
        int k = mesh.indices[m].vtxk;

        edge_to_tri[Edge(i, j)].push_back(m);
        edge_to_tri[Edge(j, k)].push_back(m);
        edge_to_tri[Edge(k, i)].push_back(m);

        vertex_to_tri[i].push_back(m);
        vertex_to_tri[j].push_back(m);
        vertex_to_tri[k].push_back(m);
    }


    std::vector<Edge> boundary_edges;
    std::vector<bool> is_boundary(mesh.vertices.size(), false);
    for(auto it = edge_to_tri.begin(); it != edge_to_tri.end(); it++) {
        if (it->second.size() == 1) {
            boundary_edges.push_back(it->first);
            is_boundary[it->first.u] = true;
            is_boundary[it->first.v] = true;
        }
    }

    std::vector<Edge> ordered_boundary_edges(boundary_edges.size());
    ordered_boundary_edges[0] = boundary_edges[0];
    for(int i = 1; i < boundary_edges.size(); ++i) {
        for(int j = 0; j < boundary_edges.size(); ++j) {
            if (boundary_edges[j].a == ordered_boundary_edges[i - 1].b) {
                ordered_boundary_edges[i] = boundary_edges[j];
                break;
            }
        }
    }


    for(int i = 0; i < ordered_boundary_edges.size(); ++i) {
        double theta = i / (double)ordered_boundary_edges.size() * 2 * M_PI;
        Vector3 circle_vtx(0.5 * std::cos(theta), 0.5 + 0.5 * std::sin(theta), 0.);        
        mesh.vertices[ordered_boundary_edges[i].a] = circle_vtx;
    }

    for(int iter = 0; iter < 10; ++iter) {
        std::vector<Vector3> updated_vertices = mesh.vertices;
        for(int i = 0; i < mesh.vertices.size(); ++i) {
            if (is_boundary[i])
                continue;

            Vector3 sum_neighbors(0., 0., 0);
            size_t total_neighbors = 0;
            for(int j = 0; j < vertex_to_tri[i].size(); ++j) {
                size_t tri_idx = vertex_to_tri[i][j];
                if (mesh.indices[tri_idx].vtxi != i)
                    sum_neighbors += mesh.vertices[mesh.indices[tri_idx].vtxi];
                if (mesh.indices[tri_idx].vtxj != i)
                    sum_neighbors += mesh.vertices[mesh.indices[tri_idx].vtxj];
                if (mesh.indices[tri_idx].vtxk != i)
                    sum_neighbors += mesh.vertices[mesh.indices[tri_idx].vtxk];

                total_neighbors += 2;
            }
            updated_vertices[i] = sum_neighbors / total_neighbors;
        }

        mesh.vertices = updated_vertices;
    }

    mesh.writeOBJ("goethe_new.obj");

    return 0;
}
