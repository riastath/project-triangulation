#include <iostream>
#include <cmath>
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"
#include <string.h>
#include <vector>
#include <sstream>
#include <typeinfo>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Vector_2.h>                                      // used for angle calculation

#include <CGAL/Triangle_2.h>

#include <numeric> // for std::gcd -> output requirements

// for json parsing, reading and writing
namespace pt = boost::property_tree;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Edge Edge;
typedef K::Vector_2 Vector;

typedef CGAL::Triangle_2<K> Triangle;

// Class for the triangulation process, using delaunay and steiner methods as requested
class PSLG {
private:
    std::string uid;
    std::vector<int> bounds;
    std::vector<int> points_x;
    std::vector<int> points_y;
    std::vector<Point> point_vec;
    std::vector<std::pair<Point, Point>> additional_constraints;
    std::vector<Point> steiner_points;  // to insert the steiner points
    int num_points;
    int num_constraints;

public:
    PSLG(std::string filename);         // constructor

    // Helper functions
    void delaunay_passer(CDT* delaunay_instance);
    void produce_output();

    // Obtuse checking functions
    double angle(Point a, Point b, Point c);
    bool is_obtuse(Point a, Point b, Point c);
    int is_obtuse_gen(CDT* instance);

    // Steiner point functions
    void insert_steiner_point(Point point);
    void insert_all_steiner(CDT *cdt);
    std::pair <Point, int> insert_steiner_center(CDT instance, CDT::Face_handle face, int num_obtuse);
    std::pair <Point, int> insert_steiner_mid(CDT instance, CDT::Face_handle face, int num_obtuse);
    std::pair <Point, int> insert_steiner_bisection(CDT instance, CDT::Face_handle face, int num_obtuse);
    std::pair <Point, int> insert_steiner_projection(CDT instance, CDT::Face_handle face, int num_obtuse);
    std::pair <Point, int> insert_steiner_circumcenter(CDT instance, CDT::Face_handle face, int num_obtuse);

    // Edge flipping functions
    bool face_is_infinite(CDT::Face_handle face, CDT *instance);
    void flip_edges(CDT *cdt);

};