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
    // CGAL::draw(cdt);

    // for (CDT::Face_handle fh: cdt.finite_face_handles()) {
    //     std::cout << "flipable: " << cdt.is_flipable(fh,0) << std::endl;
    //     cdt.flip(fh,0);  // from face fh, flip edge "facing" the corner 0
    //     break;
    // }
    // CGAL::draw(cdt);

    // bool res = graph->is_obtuse(&cdt);
    // CGAL::draw(cdt);
    // std::cout << res << std::endl;
    // if (res) std::cout << "Obtuse triangles exist in the instance" << std::endl;
    // else std::cout << "No obtuse triangles in the instance!" << std::endl;

    // graph->insert_steiner_center(&cdt);
    // graph->insert_steiner_mid(&cdt);
    // graph->is_obtuse_gen(&cdt);

    // graph->insert_steiner_bisection(&cdt);
    // graph->is_obtuse_gen(&cdt);

    std::cout << std::endl;
    std::cout << "Before processing:" << std::endl;
    graph->is_obtuse_gen(&cdt);
    std::cout << std::endl;

    int a = 0;
    while (actions[a] != none) {
        switch(actions[a]){
            case flip:
                // std::cout << "flip" << std::endl;
                graph->flip_edges(&cdt);
                break;
            case center_point:
                // std::cout << "center" << std::endl;
                graph->insert_steiner_center(&cdt);
                break;
            case mid_point:
                // std::cout << "midpoint" << std::endl;
                graph->insert_steiner_mid(&cdt);
                break;
            case bisector_point:
                // std::cout << "bisector" << std::endl;
                graph->insert_steiner_bisection(&cdt);
                break;
            case projection_point:
                graph->insert_steiner_projection(&cdt);
                break;
        }
        a++;
    }


    // graph->flip_edges(&cdt);
    // graph->insert_steiner_center(&cdt);
    // graph->insert_steiner_mid(&cdt);
    // graph->insert_steiner_bisection(&cdt);

    std::cout << std::endl;
    std::cout << "After processing:" << std::endl;
    graph->is_obtuse_gen(&cdt);
    std::cout << std::endl;

    // int faxes = 0;
    // for (CDT::Face_handle fh: cdt.finite_face_handles()) {
    //     std::cout << "face:" << std::endl;
    //     std::cout << "validity: " << fh->is_valid() << std::endl;  
    //     std::cout << fh->vertex(0)->point() << " " << fh->vertex(1)->point() << " " << fh->vertex(2)->point() << std::endl;
    //     faxes++;
    // }
    // std::cout << "num of faces: " << faxes << std::endl;

    CGAL::draw(cdt);

    // -- Starting here, testing for the output json using property tree -- //
    graph->produce_output();

    return 0;
}


// for looping : have an iteration limit and flag for obtuse's existence, so as long as there's an obtuse left, keep looping