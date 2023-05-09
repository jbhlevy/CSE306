//
//  main.cpp
//  Project_1
//
//  Created by John Levy on 05/05/2023.
//

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
#include <fstream>

#include "omp.h"

#include "boxes.h"
#include "vector.h"
#include "scene.h"
void print_caption(double time, int ray, int size, int bounce, std::vector<std::string> names, std::vector<Vector> albedos){
    
    std::string caption = "\'insert_title_here\' " + std::to_string(time) + " seconds of computation. " + std::to_string(ray) + " rays per pixel, " + std::to_string(size) + " x " + std::to_string(size) + " image, " + std::to_string(bounce) + " bounces.";
    for (int i = 0; i<names.size(); i++) {
        caption += " Albedo of " + names[i] + ": Vector(" + std::to_string(albedos[i][0]) + " " +  std::to_string(albedos[i][1]) + " " +  std::to_string(albedos[i][2]) + "). ";
    }
    caption += "light intensity: 3E10, f.o.v = 60 degrees.";
    
    std::ofstream stream;
    stream.open("caption.txt");
    stream << caption << std::endl;
    stream.close();
    return;
};

int main(int argc, const char * argv[]) {
    std::cout << "Hello, World!\n";
    
    int W = 512;
    int H = 512;
    
    Scene scene(3E10, Vector(-10, 20, 40));
    Vector camera_center(0, 0, 55);
    double alpha = M_PI / 3;
    int max_bounces = 5;
    
    //Scene set-up
    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    Objects::Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7));
    Objects::Sphere ceiling(Vector(0, 1000,0), 940, Vector(0.3, 0.5, 0.3));
    Objects::Sphere front_wall(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
    Objects::Sphere left_wall(Vector(-1000, 0, 0), 940 , Vector(0.5, 0.8, 0.1));
    Objects::Sphere right_wall(Vector(1000, 0, 0), 940,  Vector(0.9, 0.2 , 0.3));
    Objects::Sphere behind_wall(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));

    scene.addSphere(&floor);
    scene.addSphere(&ceiling);
    scene.addSphere(&left_wall);
    scene.addSphere(&right_wall);
    scene.addSphere(&front_wall);
    scene.addSphere(&behind_wall);
    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    
    std::vector<std::string> names;
    std::vector<Vector> albs;
    
    //Spheres
    
    //Uncomment these lines for a simple white sphere
//    Objects::Sphere S(Vector(0, 0, 0), 10, Vector(0.5, 0.5, 0.5));
//    scene.addSphere(&S);
//    names.push_back("sphere");
//    albs.push_back(Vector(0.5, 0.5, 0.5));
        
    //Uncomment these lines for a mirror sphere
//    Objects::Sphere S(Vector(0, 0, 0), 10, Vector(0.5, 0.5, 0.5), true);
//    scene.addSphere(&S);
//    names.push_back("mirror sphere");
//    albs.push_back(Vector(0.5, 0.5, 0.5));
    
    //Uncomment these linse for a transparent sphere
//    Objects::Sphere S(Vector(0, 0, 0), 12, Vector(0.5, 0.5, 0.5), false, true);
//    scene.addSphere(&S);
//    names.push_back("transparent sphere");
//    albs.push_back(Vector(0.5, 0.5, 0.5));
        
    //Uncomment these lines for a hollow sphere
//    Objects::Sphere S1(Vector(0, 0, 0), 12, Vector(0.5, 0.5, 0.5), false, true);
//    Objects::Sphere S2(Vector(0, 0, 0), 11, Vector(0.5, 0.5, 0.5), false, true, true);
//    scene.addSphere(&S1);
//    scene.addSphere(&S2);
//    names.push_back("hollow sphere");
//    albs.push_back(Vector(0.5, 0.5, 0.5));
    
    
    //Mesh handling
    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    //Cat mesh
    Objects::TriangleMesh mesh(Vector(0.3, 0.2, 0.25));
    mesh.readOBJ("cat.obj");
    scene.addMesh(&mesh);
    mesh.scale_and_translate(0.36, Vector(0, -10, 0));
    mesh.root = new BVH(0, mesh.indices.size());
    mesh.compute_bvh(mesh.root, 0, mesh.indices.size());
    names.push_back("cat");
    albs.push_back(Vector(0.3, 0.2, 0.25));
    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    
    
    int nb_of_paths = 60;
    std::cout << "Finished setting initial conditions" << std::endl;

    //Image creation loop
    std::vector<unsigned char> image(W * H * 3, 0);
    
    std::cout << "Creating image..." << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++){
        //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
        //Uncomment these lines to visualize the work of each thread during the execution
//        std::cout << "Hello from thread " << omp_get_thread_num() << " of "
//                          << omp_get_num_threads() << " on iteration " << i << "  of "
//                          << H << std::endl;

        //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
        int thread_id = omp_get_thread_num();
        for (int j = 0; j < W; j++){
            Vector color(0, 0, 0);
            for (int k=0; k<nb_of_paths; k++){
                
                std::vector<double>box_muller = Random::box_muller(thread_id);
                double randomX = box_muller[0];
                double randomY = box_muller[1];
                
                Vector ray_dir;
                ray_dir[0] = j - W / 2. +0.5 + randomX;
                ray_dir[1] = -i + H / 2. +0.5 + randomY;
                ray_dir[2] = - W / (2. * tan(alpha/2));
                ray_dir.normalize();
                Ray ray(camera_center, ray_dir);
                
                color += scene.get_color(ray, max_bounces); //1. Entry point of the program
            }
            color = color / nb_of_paths;
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    //Uncomment these lines to activate the timer
    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = finish - start;
    double time_duration = duration.count();
    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    print_caption(time_duration, nb_of_paths, W, max_bounces, names, albs); //// COMENT THIS LINE AND THE ABOVE SHIT
    std::cout << "Process ended sucessfully :) in " << std::to_string(time_duration) << std::endl << "Image was saved in \'image.png\' " << std::endl;
    return 0;
}
