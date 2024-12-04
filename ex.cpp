#include "pslg.hpp"
#include <fstream>
#include <iostream>


// Optimization algorithm 1 : local search
void local_search(PSLG *graph, CDT *cdt, int L) {
    int iterations = 0;

    while (iterations < L) {
        int inserted = 0;
        CDT::Finite_faces_iterator it;
        int num_obtuse = graph->is_obtuse_gen(cdt);
        int max_improvement = -1;
        std::pair<Point, int> max_pair;
        std::pair<Point, int> result;

        for (it = cdt->finite_faces_begin(); it != cdt->finite_faces_end(); it++) {
            max_improvement = -1;

            // Check steiner methods and compare improvement value
            result = graph->insert_steiner_center(*cdt, it, num_obtuse);
            if (max_improvement < result.second) {
                max_improvement = result.second;
                max_pair = result;
            }

            result = graph->insert_steiner_mid(*cdt, it, num_obtuse);
            if (max_improvement < result.second) {
                max_improvement = result.second;
                max_pair = result;
            }

            result = graph->insert_steiner_bisection(*cdt, it, num_obtuse);
            if (max_improvement < result.second) {
                max_improvement = result.second;
                max_pair = result;
            }
        
            result = graph->insert_steiner_projection(*cdt, it, num_obtuse);
            if (max_improvement < result.second) {
                max_improvement = result.second;
                max_pair = result;
            }

            result = graph->insert_steiner_circumcenter(*cdt, it, num_obtuse);
            if (max_improvement < result.second) {
                max_improvement = result.second;
                max_pair = result;
            }

            if (max_improvement < 0) {
                continue;
            }

            graph->insert_steiner_point(max_pair.first);    // insert steiner point picked
            inserted++;
        } 

        if (inserted <= 0) {
            std::cout << "No more improvements can be made" << std::endl;
            break;
        }

        graph->insert_all_steiner(cdt);
        iterations++;
    }

    std::cout << "Iterations were : " << iterations << std::endl;
    std::cout << "Final number of obtuse angles is currently: " << graph->is_obtuse_gen(cdt) << std::endl;
}


// Optimization algorithm 2 : simulated annealing
// Energy function to use in simulated annealing
double energy(PSLG *graph, CDT *cdt, double alpha, double beta) {
    int num_obtuse = graph->is_obtuse_gen(cdt);
    int num_steiner = graph->get_num_steiner_points();
    
    // energy is the weighted sum of obtuse triangles and steiner points
    return alpha * num_obtuse + beta * num_steiner;
}


void simulated_annealing(PSLG *graph, CDT *cdt, double alpha, double beta, int max_iterations) {
    double E = energy(graph, cdt, alpha, beta); // starter energy
    double T = 1.0;                             // starter temperature
    std::pair<Point, int> max_pair;
    std::pair<Point, int> result;

    std::vector<Point> new_steiner_points; // vector to keep all steiner points inserted, and add them all in the end 

    while(T >= 0) {
        CDT::Finite_faces_iterator it;
        int num_obtuse = graph->is_obtuse_gen(cdt);

        // For each obtuse triangle
        for (it = cdt->finite_faces_begin(); it != cdt->finite_faces_end(); it++) {

            if (graph->is_obtuse_face(it)) {    // if it's obtuse

                // randomly select a steiner method, and insert steiner point.
                int choice = rand() % 5;
                switch (choice) {
                    case 0:
                        result = graph->insert_steiner_center(*cdt, it, num_obtuse);
                        break;
                    case 1:
                        result = graph->insert_steiner_mid(*cdt, it, num_obtuse);
                        break;
                    case 2:
                        result = graph->insert_steiner_bisection(*cdt, it, num_obtuse);
                        break;
                    case 3:
                        result = graph->insert_steiner_projection(*cdt, it, num_obtuse);
                        break;
                    case 4: 
                        result = graph->insert_steiner_circumcenter(*cdt, it, num_obtuse);
                        break;
                    default:
                        result = std::make_pair(Point(NAN, NAN), -1);
                        break;
                }


                // Save the inserted steiner point to check for improvement
                Point steiner = result.first;
                std::cout << "Steiner point: (" << steiner.x() << ", " << steiner.y() << ")" << std::endl;

                // Skip nan (invalid) points
                if (std::isnan(steiner.x()) || std::isnan(steiner.y())) {
                    continue;
                }

                // Calculate ΔΕ after insertion
                CDT new_cdt = graph->return_cdt(*cdt, steiner);

                // // debbuging 2
                // if (cdt->is_infinite(cdt->locate(steiner))) {   // checking if it's outside the triangulation
                //     std::cerr << "Steiner point is outside the triangulation" << std::endl;
                // }
                // // end

                // // added to debug
                // if (!new_cdt.is_valid(false)) {
                //     std::cerr << "New cdt is invalid after inserting steiner" << std::endl;
                //     continue;
                // }
                // // done debugging

                double new_energy = energy(graph, &new_cdt, alpha, beta);
                double delta = new_energy - E;

                std::cout << "-----------------------" << std::endl;
                std::cout << "delta is " << delta << std::endl;


                if (delta < 0) {                // if difference is negative, accept new configuration
                    std::cout << "difference in energy is negative " << std::endl;
                    // cdt = &new_cdt;
                    new_steiner_points.push_back(steiner);
                    E = new_energy;
                    std::cout << "-----------------------" << std::endl;
                    std::cout << "difference is negative and found new energy " << E << std::endl;
                } else {                    // accept with probability 
                    double probability = exp(-delta / T);
                    // calculate random number and compare, accepting or not. maybe this is not the right way to do it, though ...
                    double rand_value = (double)rand() / RAND_MAX;
                    std::cout << "value (random) is " << rand_value << std::endl;
                    std::cout << "probability is " << probability << std::endl;
        
                    if (rand_value < probability) {
                        std::cout << "accepted the new cdt with probability." << std::endl;
                        // cdt = &new_cdt;
                        new_steiner_points.push_back(steiner);
                        E = new_energy;
                    } else { 
                        std::cout << "did not accept the new cdt with probability" << std::endl; 
                    }
                }

            }

        }

        if (!new_steiner_points.empty()) {
            for (const Point& p: new_steiner_points) {
                cdt->insert(p);
            }
        } 
        new_steiner_points.clear();             // reset for next iteration


        // recalculate steiner number
        num_obtuse = graph->is_obtuse_gen(cdt);

        // decrease temperature with formula given
        T = T - (1.0/max_iterations);           // assume they're positive, here
        std::cout << "new temperature is : " << T << std::endl;

        if (T < 0 ) break;                      // possibly redundant
    }
}


// Optimization algorithm 3 : ant colony
// Helper function for ant colony : heuristic, which uses radius-to-height ratio calculation
// followed the class slides
double ant_evaluate(double r) {
    double h; // heuristic
    // projection
    if (r > 2.0) {
        h = std::max(0.0, (r - 1) / r);
    // circumcenter
    } else if (r >= 1.0) {
        h = (r/ 2) + r;
    // midpoint
    } else if (r < 1.0) {
       h = std::max(0.0, (3 - 2 * r) / 3);
    } else h = 1.0;    // keep or not 
    // + mean of adjacent obtuse ?
    return h;
}



// Helper function for ant colony : pheromone updates
// differences vector is Δτsp​=1/(1+α⋅obtuse count+β⋅Steiner count) ​if the option improves the triangulation. 0, otherwise.​
// they will be calculated beforehand
// evaporation rate is L, given by user 

// use reference here?
void update_pheromones(std::vector<double> pheromones, const std::vector<double> differences, double evaporation_rate) {
    for (int i = 0; i < pheromones.size(); i++) {
        pheromones[i] = (1.0 - evaporation_rate) * pheromones[i] + differences[i];
    }
}


// Helper function for ant colony : improve triangulation
// (evaluate through heuristic)


// lambda : evaporation rate. K = number of ants, L = number of loops / cycles
void ant_colony(PSLG *graph, CDT *cdt, double w1, double w2, int K, int L, double lambda) {
    std::vector<double> pheromones;
    std::vector<double> differences; // ?

    // start from 1 ?
    for (int c = 0; c <= L; c++) {
        std::vector<std::vector<Point>> solutions(K); // for every ant in a certain cycle, a vector of solutions

        for (int ant = 0; ant <= K; ant++) {

        }

    }


}


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
    local_search(graph, &cdt, max_iterations);

    // testing for simulated annealing
    // simulated_annealing(graph, &cdt, alpha, beta, max_iterations);

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







// void process_instance(PSLG *graph, CDT *cdt) {
//     int iterations = 0, max_iterations = 20;    // "default" max iterations

//     while (true) {
//         int inserted = 0;
//         CDT::Finite_faces_iterator it;
//         int num_obtuse = graph->is_obtuse_gen(cdt);
//         int max_improvement = -1;
//         std::pair<Point, int> max_pair;
//         std::pair<Point, int> result;

//         for (it = cdt->finite_faces_begin(); it != cdt->finite_faces_end(); it++) {
//             max_improvement = -1;

//             // Check steiner methods and compare improvement value
//             result = graph->insert_steiner_center(*cdt, it, num_obtuse);
//             if (max_improvement < result.second) {
//                 max_improvement = result.second;
//                 max_pair = result;
//             }

//             result = graph->insert_steiner_mid(*cdt, it, num_obtuse);
//             if (max_improvement < result.second) {
//                 max_improvement = result.second;
//                 max_pair = result;
//             }

//             result = graph->insert_steiner_bisection(*cdt, it, num_obtuse);
//             if (max_improvement < result.second) {
//                 max_improvement = result.second;
//                 max_pair = result;
//             }
        
//             result = graph->insert_steiner_projection(*cdt, it, num_obtuse);
//             if (max_improvement < result.second) {
//                 max_improvement = result.second;
//                 max_pair = result;
//             }

//             result = graph->insert_steiner_circumcenter(*cdt, it, num_obtuse);
//             if (max_improvement < result.second) {
//                 max_improvement = result.second;
//                 max_pair = result;
//             }

//             if (max_improvement < 0) {
//                 continue;
//             }

//             graph->insert_steiner_point(max_pair.first);    // insert steiner point picked
//             inserted++;
//         } 

//         if (inserted <= 0) {
//             break;
//         }
        
//         graph->insert_all_steiner(cdt);

//         std::cout << "Final number of obtuse angles is currently: " << graph->is_obtuse_gen(cdt) << std::endl;

//         iterations++;
//         // for instances that do not converge
//         if (iterations >= max_iterations) {
//             break;
//         }
//     }
// }