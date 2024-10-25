#include "pslg.hpp"
#include <fstream>
#include <iostream>

void process_instance(PSLG *graph, CDT *cdt) {
    int iterations = 0, max_iterations = 20;

    while (true) {
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

            if (max_improvement < 0) {
                continue;
            }

            graph->insert_steiner_point(max_pair.first);    // insert steiner point picked
            inserted++;
        } 

        if (inserted <= 0) {
            break;
        }
        
        graph->insert_all_steiner(cdt);

        std::cout << "Final number of obtuse angles is currently: " << graph->is_obtuse_gen(cdt) << std::endl;

        iterations++;
        // for instances that do not converge
        if (iterations >= max_iterations) {
            break;
        }
    }
}


int main(int  argc, char *argv[]) {
    std::string filename = "tests/data.json";
    bool flip = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-f") {
            flip = true;
        } else {
            filename = arg;
        }
    }

    FILE* file = fopen(filename.c_str(), "r");
    if (file == NULL) {
        std::cout << "File not found, ending program" << std::endl;
        return 0;
    }
    fclose(file);
    

    PSLG * graph = new PSLG(filename);
    CDT cdt;
    graph->delaunay_passer(&cdt);

    std::cout << std::endl;
    std::cout << "Before processing:" << std::endl;
    std:: cout << "Obtuse angles are "<< graph->is_obtuse_gen(&cdt) << std::endl;
    std::cout << std::endl;

    // Processing
    std::cout << "Processing begins : " << std::endl; 

    if (flip) {
        graph->flip_edges(&cdt);
    }

    process_instance(graph, &cdt);
    std::cout << "Processing ended. " << std::endl; 
    // End

    std::cout << std::endl;
    std::cout << "After processing:" << std::endl;
    std:: cout << "Obtuse angles are "<< graph->is_obtuse_gen(&cdt) << std::endl;
    std::cout << std::endl;

    CGAL::draw(cdt);
    graph->produce_output();

    return 0;
}