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

namespace pt = boost::property_tree;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

class PSLG {
private:
    std::string uid;
    std::vector<int> bounds;
    // points are ints 
    std::vector<int> points_x;
    std::vector<int> points_y;
    std::vector<Point> point_vec;
    std::vector<std::pair<Point, Point>> additional_constraints;
    int num_points;
    int num_constraints;

public:
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


        // Find out why this doesn't work 
        for (pt::ptree::value_type &add : root.get_child("additional_constraints")) {


            // what type is this ?? 
            std::cout << typeid(add.first).name() << std::endl;  // this is a string
            std::cout << typeid(add.second).name() << std::endl;
            // additional_constraints.push_back(std::make_pair(px, py));
            // additional_constraints.push_back(add.second.get_value<std::pair<Point, Point>>());
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
    
};


int main(void) {
    std::string filename = "tests/data.json";
    PSLG * graph = new PSLG(filename);
    graph->printer();

    return 0;
}

