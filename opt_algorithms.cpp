#include "opt_algorithms.hpp"

// Optimization algorithm 1 : local search
void local_search(PSLG *graph, CDT_C *cdt, int L) {
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
            result = graph->insert_steiner_center_single(*cdt, it, num_obtuse);
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
}


// Optimization algorithm 2 : simulated annealing
// Energy function to use in simulated annealing
double energy(PSLG *graph, CDT *cdt, double alpha, double beta) {
    int num_obtuse = graph->is_obtuse_gen(cdt);
    int num_steiner = graph->get_num_steiner_points();
    
    // energy is the weighted sum of obtuse triangles and steiner points
    return alpha * num_obtuse + beta * num_steiner;
}


void simulated_annealing(PSLG *graph, CDT_C *cdt, double alpha, double beta, int max_iterations) {
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

            if (graph->is_obtuse_face(it)) {

                // randomly select a steiner method, and insert steiner point.
                int choice = rand() % 5;
                switch (choice) {
                    case 0:
                        result = graph->insert_steiner_center_single(*cdt, it, num_obtuse);
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

                // Skip nan (invalid) points
                double new_x = CGAL::to_double(steiner.x());
                double new_y = CGAL::to_double(steiner.y());
                if (steiner == Point(NULL, NULL)) {
                    continue;
                }

                // Calculate ΔΕ after insertion
                CDT new_cdt = graph->return_cdt(*cdt, steiner);
                double new_energy = energy(graph, &new_cdt, alpha, beta);
                double delta = new_energy - E;



                if (delta < 0) {                // if difference is negative, accept new configuration
                    new_steiner_points.push_back(steiner);
                    graph->insert_steiner_point(steiner);
                    E = new_energy;
                } else {                    // accept with probability 
                    double probability = exp(-delta / T);
                    // calculate random number and compare, accepting or not
                    double rand_value = (double)rand() / RAND_MAX;
        
                    if (rand_value < probability) {
                        new_steiner_points.push_back(steiner);
                        graph->insert_steiner_point(steiner);
                        E = new_energy;
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

        if (T < 0 ) break;  // just in case
    }
}


// Optimization algorithm 3 : ant colony
// Helper function for ant colony : heuristic, which uses radius-to-height ratio calculation
// followed the class slides
// using this, we get a probability for the next step
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
    } else h = 1.0;
    return h;
}

// Helper function : Print pheromones
void print_pheromones(std::vector<double> pheromones) {
    std::cout << "PHEROMONES : " << std::endl;
    for (int i = 0; i < pheromones.size(); i++) {
        std::cout << pheromones[i] << std::endl;
    }
}

// Helper function for ant colony : compute radius-to-height ratio
double compute_r(CDT::Finite_faces_iterator it) {
    Point p1 = it->vertex(0)->point();
    Point p2 = it->vertex(1)->point();
    Point p3 = it->vertex(2)->point();

    // find every length, to find the longest side, like in the pslg functions
    double length_a = std::sqrt(CGAL::to_double(CGAL::squared_distance(p2, p3))); // p2-p3 edge
    double length_b = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p3))); // p1 - p3 edge
    double length_c = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p2))); // p1 - p2 edge

    // calculate the circumradius R
    // formula is R = (abc) / sqrt((a + b + c)(b + c - a)(c + a - b)(a + b - c))
    double calc = (length_a + length_b + length_c)*(length_b + length_c - length_a)*(length_c + length_a - length_b)*(length_a + length_b - length_c);
    double R = (length_a * length_b * length_c) / std::sqrt(calc);
    
    // find which side is the longest (max of all the lengths)
    double longest = std::max({length_a, length_b, length_c});  

    // find the height from the longest side, according to the formula height = (2*area) / longest height
    double s = (length_a + length_b + length_c) / 2.0;                              // the semi perimeter
    double area = std::sqrt(s * (s - length_a) * (s - length_b) * (s - length_c));  // calculated with heron's formula A=sqrt(s⋅(s−a)⋅(s−b)⋅(s−c))
    if (area == 0) {
        return -1.0;
    }

    double height = (2.0 * area) / longest;

    // this is the final r
    return R / height;
}



// Helper function for ant colony : pheromone updates
// differences vector is Δτsp​=1/(1+α⋅obtuse count+β⋅Steiner count) ​if the option improves the triangulation. 0, otherwise.​
// evaporation rate is L, given by user 
void update_pheromones(PSLG *graph, CDT_C *cdt, std::vector<double>& pheromones, double evaporation_rate, int num_obtuse, int num_steiner, double alpha, double beta, std::vector<std::vector<Point>>solutions) {
    std::vector<double> differences(5, 0.0);        // this is the Δτsp vector, which will be used to update pheromones

    // Before updating the pheromones, check if improvement was made
    for (int i = 0; i < 5; i++) {
        bool was_improved = false;

        // Iterate over all solutions (ant_sol) of each ant
        for (std::vector<Point> ant_sol : solutions) {
            // Get status after ant's solution
            int obtuse_count_after = graph->is_obtuse_gen(cdt); 
            int steiner_count_after = graph->get_num_steiner_points();
            
            // Check if there's an improvement
            if (obtuse_count_after < num_obtuse || steiner_count_after < num_steiner) {
                was_improved = true;
                break;
            }
        }

        // Calculate the Δsp accordingly
        if (was_improved) {
            differences[i] = 1.0 / (1 + alpha * num_obtuse + beta * num_steiner);
        } else { 
            differences[i] = 0.0;
        }
    }

    // Update pheromones using Δsp
    for (int i = 0; i < pheromones.size(); i++) {
        pheromones[i] = (1.0 - evaporation_rate) * pheromones[i] + differences[i];
    }
}

// lambda : evaporation rate. K = number of ants, L = number of loops / cycles
void ant_colony(PSLG *graph, CDT_C *cdt, double alpha, double beta, double x, double y, int K, int L, double lambda) {
    std::vector<double> pheromones(5, 1.0);         // 1.0 is the starter value t_0, which must be higher than 0, initializing for every steiner option

    // Keeping the initial steiner point number to check for improvements later
    int num_obtuse = graph->is_obtuse_gen(cdt);
    int num_steiner = graph->get_num_steiner_points();

    // So we can find the best solution after each ant's solution
    int max_num_obtuse = std::numeric_limits<int>::max(); 
    int max_num_steiner = std::numeric_limits<int>::max();
    std::vector<Point> new_steiner_vec;     // pretty much like in simulated annealing, best steiner will be stored here
    
    for (int c = 1; c <= L; c++) {
        std::vector<std::vector<Point>> solutions(K); // for every ant in a certain cycle, a vector of solutions

        std::vector<std::pair<int, int>> ant_values(K); // obtuse count and steiner count for every ant

        for (int ant = 1; ant <= K; ant++) {
            std::vector<Point> ant_sol;               // The current ant's solution, which is the steiner point insertion it'll choose
            // Improve triangulation and Evaluate triangulation using heuristic

            // Each ant identifies obtuse triangles in the current triangulation and selects
            // randomly an obtuse triangle and a Steiner point location.
            // Each Steiner point option is chosen with probability as in slides.
            CDT::Finite_faces_iterator it;
            for (it = cdt->finite_faces_begin(); it != cdt->finite_faces_end(); it++) {

                if (graph->is_obtuse_face(it)) {

                    // We need the information to insert in the heuristic, to evaluate and choose.
                    double r = compute_r(it);

                    // Evaluating the triangulation using heuristic, creating a vector for the heuristics for all the options
                    std::vector<double> heuristics;
                    for (int i = 0; i < 5; i++) {
                        double heuristic = ant_evaluate(r); 
                        heuristics.push_back(heuristic);
                    }
                    
                    // Randomly select the steiner option and the point, according to probability.

                    // Calculate the probability for the current choice, which needs the heuristic scores too.
                    std::vector<double> probabilities(5, 0.0);
                    for (int i = 0; i < 5; i++) {
                        probabilities[i] = pow(pheromones[i], x) * pow(heuristics[i], y);
                    }

                    // Now, normalize the probabilities
                    // accumulation = total weight of the probabilities, to normalize
                    long double accumulation = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
                    if (accumulation > 0) {
                        for (double& prob : probabilities) {
                            prob /= accumulation;
                        }
                    }

                    // Using discrete distribution to randomly generate a number to choose from
                    std::random_device rd;          // seed for random number   
                    std::mt19937 generator(rd());
                    std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
                    int method_choice = distribution(generator);    // the steiner method for the ant

                    Point point_choice;                             // the steiner point for the ant
                    // Taking the point directly instead of using "result" variable because we don't need the improvement
                    switch (method_choice) {
                        case 0: 
                            point_choice = graph->insert_steiner_center_single(*cdt, it, num_obtuse).first; 
                            break;
                        case 1: 
                            point_choice = graph->insert_steiner_mid(*cdt, it, num_obtuse).first; 
                            break;
                        case 2: 
                            point_choice = graph->insert_steiner_bisection(*cdt, it, num_obtuse).first; 
                            break;
                        case 3: 
                            point_choice = graph->insert_steiner_projection(*cdt, it, num_obtuse).first; 
                            break;
                        case 4: 
                            point_choice = graph->insert_steiner_circumcenter(*cdt, it, num_obtuse).first; 
                            break;
                        default: 
                            point_choice = std::make_pair(Point(NULL, NULL), -1).first;
                            break;
                    }

                    // If the steiner point chosen is valid
                    if (!(point_choice == Point(NULL, NULL))) {
                        ant_sol.push_back(point_choice);            
                    }           
                }
            }

            solutions[ant - 1] = ant_sol;   // doesn't work with push_back since we have a constant number of elements, K 

        }

        // Find the best ant
        // Iterate over the solutions of the ants
       for (std::vector<Point> ant_sol : solutions) {
            int num_obtuse_ant = graph->is_obtuse_gen(cdt);
            int num_steiner_ant= ant_sol.size();

            // Compare new values with the values before, to see if new solution is better
            if (num_obtuse_ant < max_num_obtuse || (num_obtuse_ant == max_num_obtuse && num_steiner_ant < max_num_steiner)) {
                max_num_obtuse = num_obtuse_ant;
                max_num_steiner = num_steiner_ant;
                new_steiner_vec = ant_sol;  // this is the steiner point of the best current solution
                // So if it is, we keep it (save best triangulation)
            }
        }

        // Finally, insert the vector into the cdt, accepting the current configuration as the best one
        if (!new_steiner_vec.empty()) {
            for (Point& point : new_steiner_vec) {
                graph->insert_steiner_point(point);
                cdt->insert(point);
            }
        }

        // Update pheromones for next cycle
        update_pheromones(graph, cdt, pheromones, lambda, num_obtuse, num_steiner, alpha, beta, solutions);

    }   
}