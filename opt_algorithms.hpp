#include "pslg.hpp"
#include <iostream>
#include <cmath>
#include <random>
#include <numeric>
#include <vector>

// The three optimization algorithms required, to improve the triangulation
// 1. local search
void local_search(PSLG *graph, CDT *cdt, int L);

// 2. simulated annealing
double energy(PSLG *graph, CDT *cdt, double alpha, double beta);
void simulated_annealing(PSLG *graph, CDT *cdt, double alpha, double beta, int max_iterations);

// 3. ant colony
double ant_evaluate(double r);
void print_pheromones(std::vector<double> pheromones);
double compute_r(CDT::Finite_faces_iterator it);
void update_pheromones(PSLG *graph, CDT *cdt, std::vector<double>& pheromones, double evaporation_rate, int num_obtuse, int num_steiner, double alpha, double beta, std::vector<std::vector<Point>>solutions);
void ant_colony(PSLG *graph, CDT *cdt, double alpha, double beta, double x, double y, int K, int L, double lambda);