#include "opt_algorithms.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <random>

// int main(int argc, char *argv[]) {
//     std::string input_filename = "tests/data.json";     // default input
//     std::string output_filename = "output.json";        // default name for output

//     bool flip = false;

//     for (int i = 1; i < argc; i++) {
//         std::string arg = argv[i];
//         if (arg == "-i" && i + 1 < argc) {
//             input_filename = argv[++i];
//         } else if (arg == "-o" && i + 1 < argc) {
//             output_filename = argv[++i];
//         } else if (arg == "-f") {                       // in case of flipping first
//             flip = true;
//         } else {
//             std::cout << "Invalid argument! " << arg << std::endl;
//             return -1;
//         }
//     }

//     FILE* file = fopen(input_filename.c_str(), "r");
//     if (file == NULL) {
//         std::cout << "Input file not found, ending program" << std::endl;
//         return -1;
//     }
//     fclose(file);

//     PSLG *graph = new PSLG(input_filename);
//     CDT_C cdt;
//     CDT_C cdt_copy;

//     bool delaunay_flag = graph->get_delaunay_flag();
//     std::string method = graph->get_method();
//     std::map<std::string, double> params = graph->get_parameters();

//     // Only do delaunay if specified in parameters
//     if (delaunay_flag) {
//         graph->delaunay_passer(&cdt);
//         graph->delaunay_passer(&cdt_copy);
//     } else {
//         std::cout << "Starting without delaunay." << std::endl;
//     }

//     std::cout << "Before processing:" << std::endl;
//     std::cout << "Obtuse angles are " << graph->is_obtuse_gen(&cdt) << std::endl;
//     std::cout << std::endl;

//     if (flip) {
//         graph->flip_edges(&cdt);
//         std::cout << "Flipped!" << std::endl;
//     }

//     double p = 0.0; // initial convergence rate

//     // Find which method to execute and get parameters from .json file, then execute that method
//     if (method == "local") {
//         int max_iterations = params["L"];
//         local_search(graph, &cdt, max_iterations);
//     } else if (method == "sa") {
//         int alpha = params["alpha"];
//         int beta = params["beta"];
//         int max_iterations = params["L"];
//         simulated_annealing(graph, &cdt, alpha, beta, max_iterations);
//     } else if (method == "ant") {
//         int w1 = params["alpha"];
//         int w2 = params["beta"];
//         int x = params["xi"];
//         int y = params["psi"];
//         int lambda = params["lambda"];
//         int K = params["kappa"];
//         int max_iterations = params["L"];
//         ant_colony(graph, &cdt, w1, w2, x, y, K, max_iterations, lambda);
//     } else {
//         std::cerr << "Invalid method! " << std::endl;
//         return -1;
//     }


//     // New part added

//     p = graph->compute_convergence_rate(&cdt_copy);   // convergence rate
//     std::cout << "Average convergence rate p̅ is : " << p << std::endl;

//     // Checking convergence
//     if (p < 0.1) {  // small rate -> better convergence
//         std::cout << "Converged, using p " << std::endl;


//     } else {    // big rate -> worse convergence
//         std::cout << "Did not converge, using randomization " << std::endl;


//     }


//     // End of new part 


//     std::cout << "After processing:" << std::endl;
//     std::cout << "Obtuse angles are " << graph->is_obtuse_gen(&cdt) << std::endl;

//     CGAL::draw(cdt);
//     graph->produce_output(cdt);
//     std::cout << "Output written! File is : " << output_filename << std::endl;

//     return 0;
// }

Point random_point(Point centroid) {
    // Using gaussian normal distribution for two coordinates
    double deviation = 0.5; // random choice, can change

    std::random_device rd;
    std::mt19937 gen(rd());

    std::normal_distribution<> x_distribution(CGAL::to_double(centroid.x()), deviation);
    std::normal_distribution<> y_distribution(CGAL::to_double(centroid.y()), deviation);

    double random_x = x_distribution(gen);
    double random_y = y_distribution(gen);

    return Point(random_x, random_y);   // the new random point
}





int main(int argc, char *argv[]) {
    std::string input_filename = "tests/data.json";     // default input
    std::string output_filename = "output.json";        // default name for output

    bool flip = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            input_filename = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            output_filename = argv[++i];
        } else if (arg == "-f") {                       // in case of flipping first
            flip = true;
        } else {
            std::cout << "Invalid argument! " << arg << std::endl;
            return -1;
        }
    }

    FILE* file = fopen(input_filename.c_str(), "r");
    if (file == NULL) {
        std::cout << "Input file not found, ending program" << std::endl;
        return -1;
    }
    fclose(file);

    PSLG *graph = new PSLG(input_filename);
    CDT_C cdt;

    bool delaunay_flag = graph->get_delaunay_flag();
    std::string method = graph->get_method();
    std::map<std::string, double> params = graph->get_parameters();

    // Only do delaunay if specified in parameters
    if (delaunay_flag) {
        graph->delaunay_passer(&cdt);
    } else {
        std::cout << "Starting without delaunay." << std::endl;
    }

    std::cout << "Before processing:" << std::endl;
    std::cout << "Obtuse angles are " << graph->is_obtuse_gen(&cdt) << std::endl;
    std::cout << std::endl;

    if (flip) {
        graph->flip_edges(&cdt);
        std::cout << "Flipped!" << std::endl;
    }

    double p = 0.0; // initial convergence rate

    // Find which method to execute and get parameters from .json file, then execute that method
    if (method == "local") {
        int max_iterations = params["L"];
        p = local_search(graph, &cdt, max_iterations);
    } else if (method == "sa") {
        int alpha = params["alpha"];
        int beta = params["beta"];
        int max_iterations = params["L"];
        p = simulated_annealing(graph, &cdt, alpha, beta, max_iterations);
    } else if (method == "ant") {
        int w1 = params["alpha"];
        int w2 = params["beta"];
        int x = params["xi"];
        int y = params["psi"];
        int lambda = params["lambda"];
        int K = params["kappa"];
        int max_iterations = params["L"];
        p = ant_colony(graph, &cdt, w1, w2, x, y, K, max_iterations, lambda);
    } else {
        std::cerr << "Invalid method! " << std::endl;
        return -1;
    }


    // New part added

    std::cout << "Average convergence rate p̅ is : " << p << std::endl;

    // Checking convergence
    if (p < 0.1) {  // small rate -> better convergence
        std::cout << "Converged, using p " << std::endl;


    } else {    // big rate -> worse convergence
        std::cout << "Did not converge, using randomization " << std::endl;
        // CDT::Finite_faces_iterator it;
        // std::vector<CDT::Face_handle> obtuse_faces;

        // for (it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); it++) {
        //     if (graph->is_obtuse_face(it)) {
        //         obtuse_faces.push_back(it);
        //     }
        // }

        // for (CDT::Face_handle face : obtuse_faces) {
        //     Point a = it->vertex(0)->point();
        //     Point b = it->vertex(1)->point();
        //     Point c = it->vertex(2)->point();

        //     Point centroid = CGAL::centroid(a, b, c);

        //     Point point_to_insert = random_point(centroid);
        //     cdt.insert(point_to_insert);
        // }
    }


    // End of new part 


    std::cout << "After processing:" << std::endl;
    std::cout << "Obtuse angles are " << graph->is_obtuse_gen(&cdt) << std::endl;

    CGAL::draw(cdt);
    graph->produce_output(cdt);
    std::cout << "Output written! File is : " << output_filename << std::endl;

    return 0;
}