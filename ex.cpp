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


        for (pt::ptree::value_type &add : root.get_child("additional_constraints")) {

            // second -> the node's value, which is an [i, j] ptree represenation, so we can use get_child to get the child element
            // and with front and back, we get the first and second elements in the ptree, so for example 5 and 6 in [5,6]
            int i = add.second.get_child("").front().second.get_value<int>(); //  add.second.get_child("").front() is the first child, i 
            int j = add.second.get_child("").back().second.get_value<int>();  //  add.second.get_child("").back() is the second child, j 

            // Construct the points using the ints we got, i and j
            Point px(points_x[i], points_y[i]);
            Point py(points_x[j], points_y[j]);

            // then make them a pair, and push to constraints
            additional_constraints.push_back(std::make_pair(px, py));
        }
    }

    void printer() {
        std::cout << "uid is: " << uid << std::endl;
        for (int i = 0; i < bounds.size(); i++) {
            std::cout << bounds.at(i) << " ";
        }
        std::cout << std::endl;

        // test to print point pairs in constraints
        // with iterator?
        std::vector<std::pair<Point, Point>>::iterator it;
        for (it = additional_constraints.begin(); it != additional_constraints.end(); it++) {
            Point px = it->first;
            Point py = it->second;

            std::cout << "first point (px) in pair is : " << px.x() << ", " << px.y() << std::endl;
            std::cout << "second point (py) in pair is : " << py.x() << ", " << py.y() << std::endl;
        }
    }
};


int main(void) {
    std::string filename = "tests/data.json";
    PSLG * graph = new PSLG(filename);
    graph->printer();

    return 0;
}
