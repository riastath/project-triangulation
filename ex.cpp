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
    // std::vector<std::pair<int, int>> del_constr;
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

            // std::cout << typeid(add.second).name() << std::endl;
            std::pair<Point, Point> constr = std::make_pair(point_vec.at(graph_edge[0]), point_vec.at(graph_edge[1]));
            additional_constraints.push_back(constr);
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
};


int main(void) {
    std::string filename = "tests/data.json";
    PSLG * graph = new PSLG(filename);
    graph->printer();

    CDT cdt;
    graph->delaunay_passer(&cdt);
    CGAL::draw(cdt);

    return 0;
}

