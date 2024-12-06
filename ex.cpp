#include "opt_algorithms.hpp"
#include <fstream>
#include <iostream>

// MAIN TESTING LOCAL SEARCH, SIMULATED ANNEALING
int main(int argc, char *argv[]) {
    std::string filename = "tests/data.json";
    bool flip = false;
    int max_iterations = 10; // this is l, default is 10. Can be modified by command line argument (currently, before parameters)
    int alpha, beta;
    int w1, w2, x, y, lambda, K;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-f") {
            flip = true;
        // local search argument
        // using stoi instead of atoi because of weird cgal errors
        } else if (arg == "-l" && i + 1 < argc) {
            std::cout << "chose local" << std::endl;
            max_iterations = std::stoi(argv[++i]); // this is L's value
        // simulated annealing argument
        } else if (arg == "-s" && i + 3 < argc) {
            std::cout << "chose sa" << std::endl;
            // supposing the user runs this as : -s alpha_value beta_value L_value
            alpha = std::stoi(argv[++i]);   // this is alpha's value
            beta = std::stoi(argv[++i]);     // this is beta's value
            max_iterations = std::stoi(argv[++i]); // this is L's value
        } else if (arg == "-a" && i + 7 < argc) {
            std::cout << "chose ant colony" << std::endl;
            // supposing the user runs this as : -a weight 1 weight 2 x y lambda K L
            w1 = std::stoi(argv[++i]);
            w2 = std::stoi(argv[++i]);
            x = std::stoi(argv[++i]);
            y = std::stoi(argv[++i]);
            lambda = std::stoi(argv[++i]);
            K = std::stoi(argv[++i]);
            max_iterations = std::stoi(argv[++i]); // this is L's value
        } else if (arg[0] != '-') { // no flag for method given (current implementation)
            filename = arg;
        } else {
            std::cout << "invalid argument :" << arg << std::endl;
            return -1;
        }
    }

    // debugging
    std::cout << "argument count: " << argc << std::endl;
    for (int i = 0; i < argc; i++) {
        std::cout << "argv[" << i << "]: " << argv[i] << std::endl;
    }
    // debugging end


    FILE* file = fopen(filename.c_str(), "r");
    if (file == NULL) {
        std::cout << "File not found, ending program" << std::endl;
        return 0;
    }
    fclose(file);

    PSLG *graph = new PSLG(filename);
    CDT cdt;
    graph->delaunay_passer(&cdt);

    std::cout << "Before processing:" << std::endl;
    std::cout << "Obtuse angles are " << graph->is_obtuse_gen(&cdt) << std::endl;
    std::cout << std::endl;

    if (flip) {
        graph->flip_edges(&cdt);
    }

    // testing for local search
    // local_search(graph, &cdt, max_iterations);

    // testing for simulated annealing
    // simulated_annealing(graph, &cdt, alpha, beta, max_iterations);

    // testing for ant colony 
    ant_colony(graph, &cdt, w1, w2, x, y, K, max_iterations, lambda);

    std::cout << "After processing:" << std::endl;
    std::cout << "Obtuse angles are " << graph->is_obtuse_gen(&cdt) << std::endl;
    std::cout << std::endl;

    CGAL::draw(cdt);
    graph->produce_output(cdt);

    return 0;
}


// // MAIN TESTING NEW INPUT AND OUTPUT
// int main(int argc, char* argv[]) {
//     std::string filename = "tests/data.json";   // the default, if no input specified
//     std::string output_file = "output.json";
//     bool flip = false;

//     // not sure if this is right
//     for (int i = 1; i < argc; i++) {
//         std::string arg = argv[i];
//         if (arg == "-f") {
//             flip = true;
//         } else if (arg == "-i" && i + 1 < argc) {
//             filename = argv[++i];
//         } else if (arg == "-o" && i + 1 < argc) {
//             output_file = argv[++i];
//         }
//     }

//     FILE* file = fopen(filename.c_str(), "r");
//     if (file == NULL) {
//         std::cout << "File not found: " << filename << ", ending program" << std::endl;
//         return 0;
//     }
//     fclose(file);

//     PSLG* graph = new PSLG(filename);
//     CDT cdt;

//     // maybe make a func to use delaunay?
//     // if (graph->use_delaunay) {
//     //     std::cout << "Using Delaunay to start." << std::endl;
//     //     graph->delaunay_passer(&cdt);
//     // }

//     std::cout << std::endl;
//     std::cout << "Before processing:" << std::endl;
//     std::cout << "Obtuse angles are " << graph->is_obtuse_gen(&cdt) << std::endl;
//     std::cout << std::endl;

//     // this is the different part
//     // std::string method = graph->get_method();
//     // std::cout << "Processing begins. Method chosen is: " << method << std::endl;
    
//     // if (method == "local") {
//     //    // int L = somehow get L from parameters list
//     //     // std::cout << "Parameters: L chosen is " << L << std::endl;
//     //     // local_search(graph, &cdt, L);
//     // // } else if (method == "sa"), else if (method == "ant") ... {
//     // }
      
//     //process_instance(graph, &cdt); -> not needed
//     std::cout << "Processing ended." << std::endl;

//     std::cout << std::endl;
//     std::cout << "After processing:" << std::endl;
//     std::cout << "Obtuse angles are " << graph->is_obtuse_gen(&cdt) << std::endl;
//     std::cout << std::endl;

//     CGAL::draw(cdt);
//     graph->produce_output(cdt); // here, also pass output_file 

//     return 0;
// }


// // OLD MAIN
// int main(int  argc, char *argv[]) {
//     std::string filename = "tests/data.json";
//     bool flip = false;

//     for (int i = 1; i < argc; i++) {
//         std::string arg = argv[i];
//         if (arg == "-f") {
//             flip = true;
//         } else {
//             filename = arg;
//         }
//     }

//     FILE* file = fopen(filename.c_str(), "r");
//     if (file == NULL) {
//         std::cout << "File not found, ending program" << std::endl;
//         return 0;
//     }
//     fclose(file);
    

//     PSLG * graph = new PSLG(filename);
//     CDT cdt;
//     graph->delaunay_passer(&cdt);

//     std::cout << std::endl;
//     std::cout << "Before processing:" << std::endl;
//     std:: cout << "Obtuse angles are "<< graph->is_obtuse_gen(&cdt) << std::endl;
//     std::cout << std::endl;

//     // Processing
//     std::cout << "Processing begins : " << std::endl; 

//     if (flip) {
//         graph->flip_edges(&cdt);
//     }

//     process_instance(graph, &cdt);
//     std::cout << "Processing ended. " << std::endl; 
//     // End

//     std::cout << std::endl;
//     std::cout << "After processing:" << std::endl;
//     std:: cout << "Obtuse angles are "<< graph->is_obtuse_gen(&cdt) << std::endl;
//     std::cout << std::endl;

//     CGAL::draw(cdt);
//     // std::cout << "before writing" << std::endl;
//     graph->produce_output(cdt);
//     // std::cout << "wrote to json" << std::endl;


//     return 0;
// }