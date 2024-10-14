// #define CGAL_USE_BASIC_VIEWER
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <iostream>

//  #define CGAL_USE_BASIC_VIEWER

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

int main() {
    // Initialize the Constrained Delaunay Triangulation (CDT)
    CDT cdt;

    // Define the points from the PSLG (x, y coordinates)
    std::vector<Point> points = {
        Point(632, 1588), Point(1330, 1097), Point(3051, 470), Point(5040, 1077),
        Point(5883, 2766), Point(8130,3629), Point(9280, 2836), Point(9613, 4963),
        Point(9422, 6363), Point(8996, 7327), Point(8020, 7611), Point(8467, 9720),
        Point(6735, 9183), Point(4674, 7865), Point(2519, 7692), Point(973, 9797),
        Point(1205, 6005), Point(1929, 5812), Point(3203, 6301), Point(5345, 2923)
    };

    // Insert points into the triangulation
    for (const Point& p: points) {
        cdt.insert(p);
    }

    // Define and add the constrained edges (from additional_constraints)
    std::vector<std::pair<int, int>> constraints = {
        {3, 4}, {5, 6}, {9, 10}, {10, 11}, {11, 12}, {12, 13}, {13, 14},
        {14, 15}, {15, 16}, {18,19}, {19,0}
    };

    // Insert constrainde edges based on the provided indices
    for (const auto& constraint: constraints) {
        cdt.insert_constraint(points[constraint.first], points[constraint.second]);
    }

    CGAL::draw(cdt);

    return 0;
}