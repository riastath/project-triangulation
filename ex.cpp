#include "pslg.hpp"
#include <fstream>
#include <iostream>

enum action {
    none,
    flip,
    center_point,
    mid_point,
    bisector_point,
    projection_point
};

void test(PSLG *graph, CDT *cdt) {
    CDT::Finite_faces_iterator it;
    int num_obtuse = graph->is_obtuse_gen(cdt);
    int max_improvement = -1;
    std::pair<Point, int> max_pair;
    std::pair<Point, int> result;

    for (it = cdt->finite_faces_begin(); it != cdt->finite_faces_end(); it++) {
        max_improvement = -1;

        result = graph->insert_steiner_center(*cdt, it, num_obtuse);
        std::cout << "centroid is" << result.first << std::endl;
        std::cout << "improvement is" << result.second << std::endl;
        if (max_improvement < result.second) {
            max_improvement = result.second;
            max_pair = result;
        }

        result = graph->insert_steiner_mid(*cdt, it, num_obtuse);
        std::cout << "midpoint is" << result.first << std::endl;
        std::cout << "improvement is" << result.second << std::endl;
        if (max_improvement < result.second) {
            max_improvement = result.second;
            max_pair = result;
        }

        result = graph->insert_steiner_bisection(*cdt, it, num_obtuse);
        std::cout << "bisection point is" << result.first << std::endl;
        std::cout << "improvement is" << result.second << std::endl;
        if (max_improvement < result.second) {
            max_improvement = result.second;
            max_pair = result;
        }

    
        result = graph->insert_steiner_projection(*cdt, it, num_obtuse);
        std::cout << "projection point is" << result.first << std::endl;
        std::cout << "improvement is" << result.second << std::endl;
        if (max_improvement < result.second) {
            max_improvement = result.second;
            max_pair = result;
        }


        std::cout << "improvement is currently !!!!!!!!!!!!!! " << max_improvement << std::endl;
        if (max_improvement < 0) {
            continue;
        }
        std::cout << "point to insert is !!!!!!!!!!!! " << result.first << std::endl;
        graph->insert_steiner_point(max_pair.first);
    }
    graph->insert_all_steiner(cdt);
    std::cout << " ----------  inserted ------------- " << std::endl;
}


void argument_handler(int argc, char* argv[], enum action* actions, std::string*filename) {
    int a = 0;
    for (int i = 1;  i < argc; i++) {
        if (argv[i][0] != '-') {
            continue;
        }
        if (argv[i][1] == 'i') { // get filename
            if (i+1 == argc) {
                std::cout << "filename not given, reverting to default" << std::endl;
                return;
            }
            if (argv[i+1][0] == '-') {
                std::cout << "filename not given, reverting to default" << std::endl;
                return;
            }
            *filename = argv[i+1];
            continue;
        }
        switch (argv[i][1]) {
            case 'f':
                actions[a] = flip;
                a++;
                break;
            case 'c':
                actions[a] = center_point;
                a++;
                break;
            case 'm':
                actions[a] = mid_point;
                a++;
                break;
            case 'b':
                actions[a] = bisector_point;
                a++;
                break;
            case 'p':
                actions[a] = projection_point;
                a++;
                break;
            default: 
                std::cout << "unrecognized command: " << argv[i] << std::endl;
                break;
        }
    }
}

int main(int  argc, char *argv[]) {
    std::string filename = "tests/data.json";
    
    enum action actions[10];
    for (int i = 0; i < 10; i++) {
        actions[i] = none;
    }

    argument_handler(argc, argv, actions, &filename);

    // std::cout << "filename is: " << filename << std::endl;
    FILE* file = fopen(filename.c_str(), "r");
    if (file == NULL) {
        std::cout << "file not found, ending program" << std::endl;
        return 0;
    }
    fclose(file);
    

    PSLG * graph = new PSLG(filename);
    graph->printer();

    CDT cdt;
    graph->delaunay_passer(&cdt);

    std::cout << std::endl;
    std::cout << "Before processing:" << std::endl;
    graph->is_obtuse_gen(&cdt);
    std::cout << std::endl;

    int a = 0;
    while (actions[a] != none) {
        switch(actions[a]){
            case flip:
                graph->flip_edges(&cdt);
                break;
            case center_point:
                test(graph, &cdt);
                break;
            case mid_point:
                test(graph, &cdt);

                break;
            case bisector_point:
                test(graph, &cdt);
                break;
            case projection_point:
                test(graph, &cdt);
                break;
        }
        a++;
    }

    std::cout << std::endl;
    std::cout << "After processing:" << std::endl;
    graph->is_obtuse_gen(&cdt);
    std::cout << std::endl;
    CGAL::draw(cdt);

    graph->produce_output();

    return 0;
}