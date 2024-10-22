#include "pslg.h"

int main(void) {
    std::string filename = "tests/data.json";
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

    std::cout << "Before fliping" << std::endl;
    graph->is_obtuse_gen(&cdt);
    graph->flipper_not_0(&cdt);
    std::cout << "After fliping" << std::endl;
    // graph->insert_steiner_center(&cdt);
    // graph->insert_steiner_mid(&cdt);
    graph->insert_steiner_bisection(&cdt);
    graph->is_obtuse_gen(&cdt);
    CGAL::draw(cdt);

    // -- Starting here, testing for the output json using property tree -- //
    graph->produce_output();

    return 0;
}