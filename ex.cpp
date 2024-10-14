#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"
#include <string.h>
#include <sstream>
#include <typeinfo>
#include <vector>
#include <iostream>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cmath>
#include <CGAL/Vector_2.h> // used for angle calculation

namespace pt = boost::property_tree;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Edge Edge;
typedef K::Vector_2 Vector;

class PSLG {
private:
    std::string uid;
    std::vector<int> bounds;
    // points are ints 
    std::vector<int> points_x;
    std::vector<int> points_y;
    std::vector<Point> point_vec;
    std::vector<std::pair<Point, Point>> additional_constraints;
    // std::vector<std::pair<int, int>> del_constr;
    int num_points;
    int num_constraints;

public:
    // constructor
    PSLG(std::string filename) {
        pt::ptree root;
        pt:read_json(filename, root);

        // Insert .json data into the corresponding class members
        uid = root.get<std::string>("instance_uid");
        num_points = root.get<int>("num_points");
        num_constraints = root.get<int>("num_constraints");

        for (pt::ptree::value_type &bound : root.get_child("region_boundary")) {
            bounds.push_back(bound.second.get_value<int>());
        }

        // Get coordinates
        for (pt::ptree::value_type &x : root.get_child("points_x")) {
            points_x.push_back(x.second.get_value<int>());
        }
        
        for (pt::ptree::value_type &y : root.get_child("points_y")) {
            points_y.push_back(y.second.get_value<int>());
        }

        for (int i = 0; i < points_x.size(); i++) {
            point_vec.push_back(Point(points_x.at(i), points_y.at(i)));
        }

        // Find out why this doesn't work 
        for (pt::ptree::value_type &add : root.get_child("additional_constraints")) {


            // what type is this ?? 
            // std::cout << typeid(add.first).name() << std::endl;  // this is a string
            // if only it were that easy XD. this is a ptree (lol)
            // yes, every single pair of points is a ptree (why? idk, i guess it would be too easy otherwise?)
            int graph_edge[2];
            int x = 0;
            for (pt::ptree::value_type &pedge : add.second) {
                graph_edge[x] = pedge.second.get_value<int>();
                x++;
            }

            // del_constr.push_back(std::make_pair(graph_edge[0], graph_edge[1]));

            std::pair<Point, Point> constr = std::make_pair(point_vec.at(graph_edge[0]), point_vec.at(graph_edge[1]));
            additional_constraints.push_back(constr);
        }

    }

    void printer() {
        std::cout << "uid is: " << uid << std::endl;
        std::cout << "bounds: " << std::endl;
        for (int i = 0; i < bounds.size(); i++) {
            std::cout << bounds.at(i) << " ";
        }
        std::cout << std::endl;
    }
    
    void delaunay_passer(CDT* delaunay_instance) {
        
        // pass the points delaunay instance
        for (const Point& p: point_vec){
            delaunay_instance->insert(p);
        }

        // pass modified constraints
        // for (const auto& constraint: del_constr) {
        //     delaunay_instance->insert_constraint(point_vec[constraint.first], point_vec[constraint.second]);
        // }
        for (const auto& constraint: additional_constraints) {
            delaunay_instance->insert_constraint(constraint.first, constraint.second);
        }
    }

    // Function to calculate angles using the dot product method -> used to check for non-obtuse triangles later
    double angle(Point a, Point b, Point c) {
        // calculating the vectors needed for the dot product formula
        Vector u = b - a;
        Vector v = c - a;

        double product = u * v;
        // |u| and |v|. no length function exists, so using square length and taking its root
        double ulength = std::sqrt(u.squared_length());
        double vlength = std::sqrt(v.squared_length());

        double cosine = product/(ulength * vlength);
    
        // formula is dot(v,u) = |v| * |u| * cos(angle) 
        double angle = std::acos(cosine);   // arc cosine

        // check angle's degrees. formula given an angle is degrees = angle * 180/pi
        double degrees = angle * (180/M_PI);

        return degrees;
    }

    // Function to check if triangle is non-obtuse by calculating angle degrees, using the function above
    bool is_obtuse(Point a, Point b, Point c) {
        double angle_A = angle(a, b, c);
        double angle_B = angle(c, a, b);
        double angle_C = angle(b, a, c);

        if (angle_A > 90.0 || angle_B > 90.0 || angle_C > 90.0) {
            return true;    // triangle is obtuse
        }
        return false;
    }

    // to check for the actual triangulation, modify given a cdt instance 
    // that iterates with finite_faces methods and checks the vertices, after studying how it works
    // bool is_obtuse(CDT* instance ) {

    //     return false;
    // }
};


int main(void) {
    std::string filename = "tests/data.json";
    PSLG * graph = new PSLG(filename);
    graph->printer();

    CDT cdt;
    graph->delaunay_passer(&cdt);
    // CGAL::draw(cdt);

    // testing with a triangle that should be obtuse
    Point a(0, 0);
    Point b(2, 0);
    Point c(1, 3);
    bool res = graph->is_obtuse(a, b, c);
    std::cout << res << std::endl;

    // this is useful for iterating over faces, so i'm keeping it here
    // Delaunay::Finite_faces_iterator it;
    // for (it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++) {
    //     std::cout << dt.triangle(it) << std::endl;
    // }

    return 0;
}

