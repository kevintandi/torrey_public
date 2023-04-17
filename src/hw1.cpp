#include "hw1.h"
#include "hw1_scenes.h"
#include "ray.h" 

using namespace hw1;

Image3 hw_1_1(const std::vector<std::string> &/*params*/) {
    // Homework 1.1: generate camera rays and output the ray directions
    // The camera is positioned at (0, 0, 0), facing towards (0, 0, -1),
    // with an up vector (0, 1, 0) and a vertical field of view of 90 degree.

    Image3 img(640 /* width */, 480 /* height */);

    // Camera
    auto viewport_height = 2.0;
    auto viewport_width = 640/480.0 * viewport_height; // Make sure its a float
    auto focal_length = 1.0;

    auto origin = Vector3(0., 0., 0.);
    auto horizontal = Vector3(viewport_width, 0., 0.);
    auto vertical = Vector3(0., viewport_height, 0.);
    auto upper_left_corner = origin - horizontal/2. + vertical/2. - Vector3(0., 0., focal_length);

    // Render
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            // Initialize ray 
            auto u = double(x+0.5) / (img.width);
            auto v = double(y+0.5) / (img.height);
            ray r(origin, upper_left_corner + u*horizontal - v*vertical - origin);

            img(x,y) = normalize(r.dir);  

            //img(x,y) = normalize(Vector3(0., 2*v-1.0, -1.0));
        }
    }
    return img;
}

double hit_sphere(const Vector3& center, double radius, const ray& r){
    Vector3 oc = r.origin() - center; 
    auto a = dot(r.direction(), r.direction()); 
    auto b = 2.0 * dot(oc, r.direction()); 
    auto c = dot(oc, oc) - radius*radius; 
    auto discriminant = b*b -4*a*c; 

    auto epsilon = 0.0001; 

    // Move this out of the function! 
    auto t_min = epsilon;
    auto t_max = 9999999999999999999.; 

    if(discriminant < 0){
        return -1.0; 
    } else{
        auto root = (-b - sqrt(discriminant))/(2.0*a); 
        if (root < t_min || t_max < root) {
            root = (-b + sqrt(discriminant))/(2.0*a);
            if (root < t_min || t_max < root){return -1.0;}
        }
        return root; 
    }
}

Image3 hw_1_2(const std::vector<std::string> &/*params*/) {
    // Homework 1.2: intersect the rays generated from hw_1_1
    // with a unit sphere located at (0, 0, -2)

    Image3 img(640 /* width */, 480 /* height */);

    // Camera
    auto viewport_height = 2.0;
    auto viewport_width = 640/480.0 * viewport_height; // Make sure its a float
    auto focal_length = 1.0;

    auto origin = Vector3(0., 0., 0.);
    auto horizontal = Vector3(viewport_width, 0., 0.);
    auto vertical = Vector3(0., viewport_height, 0.);
    auto upper_left_corner = origin - horizontal/2. + vertical/2. - Vector3(0., 0., focal_length);

    // Render
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            // Initialize ray 
            auto u = double(x+0.5) / (img.width);
            auto v = double(y+0.5) / (img.height);
            ray r(origin, upper_left_corner + u*horizontal - v*vertical - origin);

            // Intersect ray with surface normal
            auto sphere_c = Vector3(0.0, 0.0, -2.0);
            auto t = hit_sphere(sphere_c, 1.0, r);
            if(t > 0.0){ 
                Vector3 N = (r.at(t) - sphere_c); // Normal vector
                img(x,y) = normalize((N+1.0)/2.0); //need normalize(N)?  
            } else{ 
                img(x, y) = Vector3(0.5, 0.5, 0.5);
            }
            //img(x,y) = normalize(Vector3(0., 2*v-1.0, -1.0));
        }
    }
    return img;
}

class camera {
    public:
        camera(
            Vector3d lookfrom,
            Vector3d lookat,
            Vector3d   vup,
            double vfov, // vertical field-of-view in degrees
            double aspect_ratio
        ) {
            auto theta = radians(vfov);
            auto h = tan(theta/2);
            auto viewport_height = 2.0 * h;
            auto viewport_width = aspect_ratio * viewport_height;

            auto w = normalize(lookfrom - lookat);
            auto u = normalize(cross(vup, w));
            auto v = cross(w, u);

            origin = lookfrom;
            horizontal = viewport_width * u;
            vertical = viewport_height * v;
            upper_left_corner = origin - horizontal/2. + vertical/2. - w;
        }

        ray get_ray(double s, double t) const {
            return ray(origin, upper_left_corner + s*horizontal - t*vertical - origin);
        }

    private:
        Vector3 origin;
        Vector3 upper_left_corner;
        Vector3 horizontal;
        Vector3 vertical;
};



Image3 hw_1_3(const std::vector<std::string> &params) {
    // Homework 1.3: add camera control to hw_1_2. 
    // We will use a look at transform:
    // The inputs are "lookfrom" (camera position),
    //                "lookat" (target),
    //                and the up vector
    // and the vertical field of view (in degrees).
    // If the user did not specify, fall back to the default
    // values below.
    // If you use the default values, it should render
    // the same image as hw_1_2.

    Vector3 lookfrom = Vector3{0, 0,  0};
    Vector3 lookat   = Vector3{0, 0, -2};
    Vector3 up       = Vector3{0, 1,  0};
    Real    vfov     = 90;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-lookfrom") {
            Real x = std::stof(params[++i]);
            Real y = std::stof(params[++i]);
            Real z = std::stof(params[++i]);
            lookfrom = Vector3{x, y, z};
        } else if (params[i] == "-lookat") {
            Real x = std::stof(params[++i]);
            Real y = std::stof(params[++i]);
            Real z = std::stof(params[++i]);
            lookat = Vector3{x, y, z};
        } else if (params[i] == "-up") {
            Real x = std::stof(params[++i]);
            Real y = std::stof(params[++i]);
            Real z = std::stof(params[++i]);
            up = Vector3{x, y, z};
        } else if (params[i] == "-vfov") {
            vfov = std::stof(params[++i]);
        }
    }

    // avoid unused warnings
    UNUSED(lookfrom);
    UNUSED(lookat);
    UNUSED(up);
    UNUSED(vfov);

    Image3 img(640 /* width */, 480 /* height */);

    // Initialize camera 
    double aspect_ratio = 640./480.; 
    camera cam(lookfrom, lookat, up, vfov, aspect_ratio);

    // Render
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            // Initialize ray 
            auto u = double(x+0.5) / (img.width);
            auto v = double(y+0.5) / (img.height);

            ray r = cam.get_ray(u,v); 

            // Intersect ray with surface normal
            auto sphere_c = Vector3(0.0, 0.0, -2.0);
            auto t = hit_sphere(sphere_c, 1.0, r);
            if(t > 0.0){ 
                Vector3 N = normalize(r.at(t) - sphere_c); // Normal vector
                img(x,y) = ((N+1.0)/2.0); //need normalize(N)?  
            } else{ 
                img(x, y) = Vector3(0.5, 0.5, 0.5);
            }
            //img(x,y) = normalize(Vector3(0., 2*v-1.0, -1.0));
        }
    }

    return img;
}

Image3 hw_1_4(const std::vector<std::string> &params) {
    // Homework 1.4: render the scenes defined in hw1_scenes.h
    // output their diffuse color directly.
    if (params.size() == 0) {
        return Image3(0, 0);
    }

    int scene_id = std::stoi(params[0]);
    UNUSED(scene_id); // avoid unused warning
    // Your scene is hw1_scenes[scene_id]
    Scene scene = hw1_scenes[scene_id]; 

    Image3 img(640 /* width */, 480 /* height */);

    // Initialize camera 
    double aspect_ratio = 640./480.; 
    camera cam(scene.camera.lookfrom,scene.camera.lookat, scene.camera.up, scene.camera.vfov, aspect_ratio);

    // Render
    // Loop through scene
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {

        // Initialize ray 
        auto u = double(x+0.5) / (img.width);
        auto v = double(y+0.5) / (img.height);
        ray r = cam.get_ray(u,v); 

        // loop through objects 
        auto i = 0; 
        bool hit = false; 
        double min_d = 0; 
        for(auto sphere: scene.shapes){
            // Intersect ray with surface normal
            auto t = hit_sphere(sphere.center, sphere.radius, r);
            if(t > 0.0){ 
                hit = true; 
                auto distance_i = distance(r.at(t), r.at(0.));
                if(i==0){
                    min_d = distance_i;
                }
                // If object closer, overwrite 
                if(distance_i <= min_d){
                    min_d = distance_i;
                    img(x,y) = scene.materials[sphere.material_id].color; //need normalize(N)?                            
                }
                i++; 
            } 
        }
        if(!hit){
            img(x, y) = Vector3(0.5, 0.5, 0.5);
        }
                  
        }
    }
    return img;
}



// Just for debugging: refactor later if it works 
double hit_sphere_shadow(const Vector3& center, double radius, const ray& r, double t_min, double t_max){
    Vector3 oc = r.origin() - center; 
    auto a = dot(r.direction(), r.direction()); 
    auto b = 2.0 * dot(oc, r.direction()); 
    auto c = dot(oc, oc) - radius*radius; 
    auto discriminant = b*b -4*a*c; 

    if(discriminant < 0){
        return -1.0; 
    } else{
        auto root = (-b - sqrt(discriminant))/(2.0*a); 
        if (root < t_min || t_max < root) {
            root = (-b + sqrt(discriminant))/(2.0*a);
            if (root < t_min || t_max < root){return -1.0;}
        }
        return root; 
    }
}

Image3 hw_1_5(const std::vector<std::string> &params) {
    // Homework 1.5: render the scenes defined in hw1_scenes.h,
    // light them using the point lights in the scene.
    if (params.size() == 0) {
        return Image3(0, 0);
    }

    int scene_id = std::stoi(params[0]);
    UNUSED(scene_id); // avoid unused warning
    // Your scene is hw1_scenes[scene_id]
    Scene scene = hw1_scenes[scene_id]; 

    Image3 img(640 /* width */, 480 /* height */);
    // Initialize camera 
    double aspect_ratio = 640./480.; 
    camera cam(scene.camera.lookfrom,scene.camera.lookat, scene.camera.up, scene.camera.vfov, aspect_ratio);

    // Render
    // Loop through scene
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {

        // Initialize ray 
        auto u = double(x+0.5) / (img.width);
        auto v = double(y+0.5) / (img.height);
        ray r = cam.get_ray(u,v); 

        // loop through objects 
        auto i = 0; 
        bool hit = false; 
        double min_d = 0; 
        for(auto sphere: scene.shapes){
            // Intersect ray with surface normal
            auto t = hit_sphere(sphere.center, sphere.radius, r);

            Vector3 shade_p = r.at(t);
            auto distance_i = distance(shade_p, r.orig);

            if(t > 0.0){ 
                hit = true; 
                if(i==0){
                    min_d = distance_i;
                } 
                // If object closer, overwrite 
                if(distance_i <= min_d){
                    min_d = distance_i;
                    Vector3 N = normalize(shade_p - sphere.center); // Normal vector
                    
                    
                    // Calculate color using response...
                    img(x,y) = Vector3(0,0,0); // Init at 0
                    
                    for(auto light: scene.lights){                        
                        // l: normalized vector between shading point & point light


                        // Fixed the bug! Visibility needs to be reset to 1 each time!
                        double visibility = 1.;
                        Vector3 l = normalize(light.position - shade_p);
                        ray shadow(light.position, -l);   
                        //auto v = hit_sphere(shade_p, 0.0, shadow); // point is sphere radius 0 
                        double shade_d = distance(shade_p, light.position); 

                        // Check if shadow ray gets blocked by any spheres...
                        for(auto sphere_l: scene.shapes){
                            double eps = 0.00001; 
                            double k = hit_sphere_shadow(sphere_l.center, sphere_l.radius, shadow, eps, (1.-eps)*shade_d); 
                            if(k == -1.0){continue;}
                            if(k < (1.-eps)*shade_d){visibility = 0.; break;}
                            // auto distance_sphere_l = distance(shadow.at(k), light.position); 

                            // //shadow ray test here. Remove the one in the hit sphere function!
                            // if(distance_sphere_l > eps && distance_sphere_l < (1.-eps)*shade_d){visibility = 0.; break;}
                        }
                        
                        //why do we need this?
                        double dot_nl = dot(N,l); 
                        if(dot_nl < 0.){dot_nl = -dot_nl;} 

                        img(x,y) += scene.materials[sphere.material_id].color*std::max(dot_nl,0.0)*light.intensity
                                        * visibility
                                        /(c_PI*distance_squared(shade_p, light.position));
                    }                    
                }
                i++; 
            } 
        }
        if(!hit){
            img(x, y) = Vector3(0.5, 0.5, 0.5);
        }           
        }
    }
    return img;
}

#include <random> 

inline double random_double(std::mt19937 &rng) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}


inline double clamp(double x, double min, double max){ 
    if(x < min) return min;  
    if(x > max) return max;
    return x; 
}

Image3 hw_1_6(const std::vector<std::string> &params) {
    // Homework 1.6: add antialiasing to homework 1.5
    if (params.size() == 0) {
        return Image3(0, 0);
    }

    int scene_id = 0;
    int spp = 64;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        } else {
            scene_id = std::stoi(params[i]);
        }
    }

    UNUSED(scene_id); // avoid unused warning
    UNUSED(spp); // avoid unused warning
    // Your scene is hw1_scenes[scene_id]

    Image3 img(160 /* width */, 120 /* height */);
    Scene scene = hw1_scenes[scene_id]; 

    // Initialize camera 
    double aspect_ratio = 160./120.; 
    camera cam(scene.camera.lookfrom,scene.camera.lookat, scene.camera.up, scene.camera.vfov, aspect_ratio);
    std::mt19937 rng(time(nullptr)); 

    // Render
    // Loop through scene
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {

        // Shoot multiple rays and average 
        for(int s=0; s < spp; ++s){
            auto u = double(x+random_double(rng)) / (img.width);
            auto v = double(y+random_double(rng)) / (img.height);
            ray r = cam.get_ray(u,v); 
            Vector3 tempcolor(0.,0.,0.);

            // loop through objects 
            auto i = 0; 
            bool hit = false; 
            double min_d = 0; 
            for(auto sphere: scene.shapes){
                // Intersect ray with surface normal
                auto t = hit_sphere(sphere.center, sphere.radius, r);

                Vector3 shade_p = r.at(t);
                auto distance_i = distance(shade_p, r.orig);

                if(t > 0.0){ 
                    std::cout << 'i'; 
                    hit = true; 
                    if(i==0){
                        min_d = distance_i;
                    }
                    // If object closer, overwrite 
                    if(distance_i <= min_d){
                        min_d = distance_i;
                        Vector3 N = normalize (shade_p - sphere.center); // Normal vector

                        // Calculate color using response...
                        tempcolor = Vector3(0,0,0); // Init at 0
                        
                        for(auto light: scene.lights){               

                            // Visibility needs to be reset inside this loop!
                            double visibility = 1.;         
                            // l: normalized vector between shading point & point light
                            Vector3 l = normalize(light.position - shade_p);
                            ray shadow(light.position, -l);   
                            //auto v = hit_sphere(shade_p, 0.0, shadow); // point is sphere radius 0 
                            double shade_d = distance(shade_p, light.position); 

                            // Check if shadow ray gets blocked by any spheres...
                            for(auto sphere_l: scene.shapes){
                                auto eps = 0.00001; 
                                auto k = hit_sphere_shadow(sphere_l.center, sphere_l.radius, shadow, eps, (1.-eps)*shade_d); 
                                if(k == -1.0){continue;}
                                if(k < shade_d){visibility = 0.; break;}
                            }
                            
                            //why do we need this?
                            double dot_nl = dot(N,l); 
                            if(dot_nl < 0.){dot_nl = -dot_nl;} 

                            tempcolor += scene.materials[sphere.material_id].color*std::max(dot_nl,0.)*light.intensity
                                            * visibility
                                            /(c_PI*distance_squared(shade_p, light.position));
                        }                    
                    }
                    i++; 
                } 
            }
            if(!hit){
                tempcolor += Vector3(0.5, 0.5, 0.5);
            }  

            img(x,y) += tempcolor;         
        }        
        
        // Normalize for anti aliasing
        img(x,y) = img(x,y)/double(spp); 
        img(x,y).x = clamp(img(x,y).x, 0, 0.999);
        img(x,y).y = clamp(img(x,y).y, 0, 0.999);
        img(x,y).z = clamp(img(x,y).z, 0, 0.999);
        }
    }
    return img;
}

Vector3 reflect(const Vector3& v, const Vector3& n){ 
    return (v - 2*dot(v,n)*n); 
}

struct hit_record{
    Vector3 point; 
    Vector3 normal; 
    double t; 
};


// Only call if mirror surface
Vector3 raycolor(ray r, const Scene& world, Vector3 color, int depth){

    // loop through objects 
    auto i = 0; 
    bool hit = false; 
    double min_d = 0; 

    // Init variables
    Vector3 surface_color = Vector3(0.,0.,0.);
    Vector3 shade_p = Vector3(0.,0.,0.);
    Vector3 N = Vector3(0.,0.,0.); 

    for(auto sphere: world.shapes){ 
        // Intersect ray with surface normal
        auto t = hit_sphere(sphere.center, sphere.radius, r);
        auto distance_i = distance(r.at(t), r.orig);

        if(t > 0.0){ 
            hit = true; 
            if(i==0){
                min_d = distance_i;
            }
            // If object closer, overwrite 
            if(distance_i <= min_d){
                min_d = distance_i;
                shade_p = r.at(t); 
                
                N = normalize(shade_p - sphere.center); // Normal vector
                surface_color = world.materials[sphere.material_id].color;  

                // std::cout << "current color is " << surface_color << '\n'; 
                // std::cout << "ray color is " << color << '\n'; 

                //If mirror, reflect light ray and accumulate color...
                if(world.materials[sphere.material_id].type == MaterialType::Mirror){


                    // // Make sure N and r is in same direction!
                    // double dot_nr = dot(N,normalize(r.dir)); 
                    // if(dot_nr > 0.){N = -N;} 

                    Vector3 dir_reflect = reflect(normalize(r.dir), N); 
                    ray r_reflect(shade_p, dir_reflect); 
                    surface_color *= color;  
                    depth--; 
                    //std::cout << "reflect!, depth = "<< depth << "\n"; 
                    return surface_color*raycolor(r_reflect, world, surface_color, depth); 
                }else{ // Lambertian
                    // Don't accumulate! 
                    //std::cout << "Lambertian! \n "; 
                } 

                if(depth == 0){
                    std::cout << "end recurse" << "\n" ; 
                    Vector3 lambertian_color(0.,0.,0.); 

                    ///// Light the shading points ////
                    for(auto light: world.lights){               

                        // Visibility needs to be reset inside this loop!
                        double visibility = 1.;         
                        // l: normalized vector between shading point & point light
                        Vector3 l = normalize(light.position - shade_p);
                        ray shadow(light.position, -l);   
                        double shade_d = distance(shade_p, light.position); 

                        // Check if shadow ray gets blocked by any spheres...
                        for(auto sphere_l: world.shapes){
                            auto eps = 0.00001; 
                            auto k = hit_sphere_shadow(sphere_l.center, sphere_l.radius, shadow, eps, (1.-eps)*shade_d); 
                            if(k == -1.0){continue;}
                            if(k < shade_d){visibility = 0.; break;}
                        }
                        
                        //why do we need this?
                        double dot_nl = dot(N,l); 
                        if(dot_nl < 0.){dot_nl = -dot_nl;} 

                        lambertian_color += surface_color*std::max(dot_nl,0.)*light.intensity
                                        * visibility
                                        /(c_PI*distance_squared(shade_p, light.position));
                        
                    }
                    //// Light the shading points ////
                    
                    return lambertian_color;
                }
            }
            i++; 
        }
    }
    if(hit){
        //std::cout << "Color is: " << color*surface_color <<"\n"; 

        Vector3 lambertian_color(0.,0.,0.); 

        ///// Light the shading points ////
        for(auto light: world.lights){               

            // Visibility needs to be reset inside this loop!
            double visibility = 1.;         
            // l: normalized vector between shading point & point light
            Vector3 l = normalize(light.position - shade_p);
            ray shadow(light.position, -l);   
            double shade_d = distance(shade_p, light.position); 

            // Check if shadow ray gets blocked by any spheres...
            for(auto sphere_l: world.shapes){
                auto eps = 0.00001; 
                auto k = hit_sphere_shadow(sphere_l.center, sphere_l.radius, shadow, eps, (1.-eps)*shade_d); 
                if(k == -1.0){continue;}
                if(k < shade_d){visibility = 0.; break;}
            }
            
            //why do we need this?
            double dot_nl = dot(N,l); 
            if(dot_nl < 0.){dot_nl = -dot_nl;} 

            lambertian_color += surface_color*std::max(dot_nl,0.)*light.intensity
                            * visibility
                            /(c_PI*distance_squared(shade_p, light.position));
            
        }
        return lambertian_color;
        //// Light the shading points ////

    }else{ // Light lambertian point
        //std::cout << "hit background! color is: " << color*Vector3(0.5, 0.5, 0.5) <<"\n"; 
        return Vector3(0.5, 0.5, 0.5);        
    }

}

Image3 hw_1_7(const std::vector<std::string> &params) {
    // Homework 1.7: add mirror materials to homework 1.6
    if (params.size() == 0) {
        return Image3(0, 0);
    }

    int scene_id = 0;
    int spp = 64;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        } else {
            scene_id = std::stoi(params[i]);
        }
    }

    UNUSED(scene_id); // avoid unused warning
    UNUSED(spp); // avoid unused warning
    // Your scene is hw1_scenes[scene_id]

    Image3 img(640 /* width */, 480 /* height */);
    Scene scene = hw1_scenes[scene_id]; 

    // Initialize camera 
    double aspect_ratio = 640./480.; 
    camera cam(scene.camera.lookfrom,scene.camera.lookat, scene.camera.up, scene.camera.vfov, aspect_ratio);
    std::mt19937 rng(time(nullptr)); 

    // Render
    // Loop through scene
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            // Shoot multiple rays and average 
            for(int s=0; s < spp; ++s){
                auto u = double(x+random_double(rng)) / (img.width);
                auto v = double(y+random_double(rng)) / (img.height);
                ray r = cam.get_ray(u,v); 
                Vector3 temp_color(0.,0.,0.);
                
                // loop through objects to find shading point
                auto i = 0; 
                bool hit = false; 
                // Vector3 shade_p(0.,0.,0.);
                double min_d = 0; 
                Vector3 shade_p(0.,0.,0.); 
                Vector3 N(0.,0.,0.); 
                MaterialType hit_mat = MaterialType::Mirror; 


                // Find hit point and properties
                for(auto sphere: scene.shapes){
                    // Intersect ray with surface normal
                    auto t = hit_sphere(sphere.center, sphere.radius, r);
                    Vector3 int_p = r.at(t);
                    auto distance_i = distance(int_p, r.orig);

                    if(t > 0.0){ 
                        hit = true; 
                        if(i==0){
                            min_d = distance_i;
                        }
                        // If object closer, overwrite 
                        if(distance_i <= min_d){
                            min_d = distance_i;
                            shade_p = int_p; 
                            N = normalize (shade_p - sphere.center); // Normal vector
                            hit_mat = scene.materials[sphere.material_id].type; //Update Material 
                        }
                        i++; 
                    }
                }


                if(hit){
                    // Reflect rays and accumulate color! 
                    int maxdepth = 200; 
                    Vector3 color = raycolor(r, scene, Vector3(1.,1.,1.), maxdepth); // Any color works, here init to zero, raycolor will handle it!      
                    // ///// Light the shading points ////
                    // for(auto light: scene.lights){               

                    //     // Visibility needs to be reset inside this loop!
                    //     double visibility = 1.;         
                    //     // l: normalized vector between shading point & point light
                    //     Vector3 l = normalize(light.position - shade_p);
                    //     ray shadow(light.position, -l);   
                    //     double shade_d = distance(shade_p, light.position); 

                    //     // Check if shadow ray gets blocked by any spheres...
                    //     for(auto sphere_l: scene.shapes){
                    //         auto eps = 0.00001; 
                    //         auto k = hit_sphere_shadow(sphere_l.center, sphere_l.radius, shadow, eps, (1.-eps)*shade_d); 
                    //         if(k == -1.0){continue;}
                    //         if(k < shade_d){visibility = 0.; break;}
                    //     }
                        
                    //     //why do we need this?
                    //     double dot_nl = dot(N,l); 
                    //     if(dot_nl < 0.){dot_nl = -dot_nl;} 

                    //     if(hit_mat == MaterialType::Diffuse){
                    //         temp_color += color*std::max(dot_nl,0.)*light.intensity
                    //                         * visibility
                    //                         /(c_PI*distance_squared(shade_p, light.position));
                    //     }else{ // Mirror
                    //         temp_color += color*std::max(dot_nl,0.)*light.intensity
                    //                         * visibility
                    //                         /(c_PI*distance_squared(shade_p, light.position)); 
                    //     }

                        
                    // }
                    // //// Light the shading points ////

                    temp_color = color; 

                }else{ 
                    temp_color = Vector3(0.5,0.5, 0.5);  
                }
                img(x,y) += temp_color;         
            }
            // Normalize for anti aliasing
            img(x,y) = img(x,y)/double(spp); 
            img(x,y).x = clamp(img(x,y).x, 0, 0.999);
            img(x,y).y = clamp(img(x,y).y, 0, 0.999);
            img(x,y).z = clamp(img(x,y).z, 0, 0.999);             
    } 
}
return img;

}

#include "parallel.h"

Image3 hw_1_8(const std::vector<std::string> &params) {
    // Homework 1.8: parallelize HW 1.7
    if (params.size() == 0) {
        return Image3(0, 0);
    }

    int scene_id = 0;
    int spp = 64;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-spp") {
            spp = std::stoi(params[++i]);
        } else {
            scene_id = std::stoi(params[i]);
        }
    }

    UNUSED(scene_id); // avoid unused warning
    UNUSED(spp); // avoid unused warning
    // Your scene is hw1_scenes[scene_id]

    Image3 img(1280 /* width */, 960 /* height */);
    Scene scene = hw1_scenes[scene_id]; 

    // Initialize camera 
    double aspect_ratio = 640./480.; 
    camera cam(scene.camera.lookfrom,scene.camera.lookat, scene.camera.up, scene.camera.vfov, aspect_ratio);

    int w = 1280, h = 960;

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;    

    parallel_for([&](const Vector2i &tile) {
        std::mt19937 RNG(tile[1] * num_tiles_x + tile[0] /* seed */); 
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
            // render pixel (x, y)
                for(int s=0; s < spp; ++s){
                    auto u = double(x+random_double(RNG)) / (img.width);
                    auto v = double(y+random_double(RNG)) / (img.height);
                    ray r = cam.get_ray(u,v); 
                    Vector3 temp_color(0.,0.,0.);
                    
                    // loop through objects to find shading point
                    auto i = 0; 
                    bool hit = false; 
                    // Vector3 shade_p(0.,0.,0.);
                    double min_d = 0; 
                    Vector3 shade_p(0.,0.,0.); 
                    Vector3 N(0.,0.,0.); 
                    MaterialType hit_mat = MaterialType::Mirror; 


                    // Find hit point and properties
                    for(auto sphere: scene.shapes){
                        // Intersect ray with surface normal
                        auto t = hit_sphere(sphere.center, sphere.radius, r);
                        Vector3 int_p = r.at(t);
                        auto distance_i = distance(int_p, r.orig);

                        if(t > 0.0){ 
                            hit = true; 
                            if(i==0){
                                min_d = distance_i;
                            }
                            // If object closer, overwrite 
                            if(distance_i <= min_d){
                                min_d = distance_i;
                                shade_p = int_p; 
                                N = normalize (shade_p - sphere.center); // Normal vector
                                hit_mat = scene.materials[sphere.material_id].type; //Update Material 
                            }
                            i++; 
                        }
                    }

                    if(hit){

                        // Reflect rays and accumulate color! 
                        int maxdepth = 200; 
                        Vector3 color = raycolor(r, scene, Vector3(1.,1.,1.), maxdepth); 
                        // Any color works, here init to zero, raycolor will handle it!      
                        temp_color = color; 

                    }else{ 
                        temp_color = Vector3(0.5,0.5, 0.5);  
                    }
                    img(x,y) += temp_color;         
                }
                // Normalize for anti aliasing
                img(x,y) = img(x,y)/double(spp); 
                img(x,y).x = clamp(img(x,y).x, 0, 0.999);
                img(x,y).y = clamp(img(x,y).y, 0, 0.999);
                img(x,y).z = clamp(img(x,y).z, 0, 0.999);

            }
        }
    }, Vector2i(num_tiles_x, num_tiles_y));    

    return img;
}
