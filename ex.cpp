#include "opt_algorithms.hpp"
#include <fstream>
#include <iostream>
#include <map>

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
    CDT cdt;

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

    // Find which method to execute and get parameters from .json file, then execute that method
    if (method == "local") {
        int max_iterations = params["L"];
        local_search(graph, &cdt, max_iterations);
    } else if (method == "sa") {
        int alpha = params["alpha"];
        int beta = params["beta"];
        int max_iterations = params["L"];
        simulated_annealing(graph, &cdt, alpha, beta, max_iterations);
    } else if (method == "ant") {
        int w1 = params["alpha"];
        int w2 = params["beta"];
        int x = params["xi"];
        int y = params["psi"];
        int lambda = params["lambda"];
        int K = params["kappa"];
        int max_iterations = params["L"];
        ant_colony(graph, &cdt, w1, w2, x, y, K, max_iterations, lambda);
    } else {
        std::cerr << "Invalid method! " << std::endl;
        return -1;
    }

    std::cout << "After processing:" << std::endl;
    std::cout << "Obtuse angles are " << graph->is_obtuse_gen(&cdt) << std::endl;

    CGAL::draw(cdt);
    graph->produce_output(cdt);
    std::cout << "Output written! File is : " << output_filename << std::endl;

    return 0;
}