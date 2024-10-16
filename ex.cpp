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
    // std::vector<Point> steiner_points;  // to insert the steiner points
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

        for (pt::ptree::value_type &add : root.get_child("additional_constraints")) {
            int graph_edge[2];
            int x = 0;
            for (pt::ptree::value_type &pedge : add.second) {
                graph_edge[x] = pedge.second.get_value<int>();
                x++;
            }

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
    bool is_obtuse_gen(CDT* instance) {
        CDT::Finite_faces_iterator it;  // initialize iterator
        for (it = instance->finite_faces_begin(); it != instance->finite_faces_end(); it++) {
            // need to examine every vertex using the iterator, and find the 3 points of each triangle
            // these examine the 1st, 2nd and 3rd vertex respectively.
            // vertex() apparently returns a handle, which is akin to a pointer to an object.
            // so, using this we get the point with its coordinates (x,y) from the corresponding vertex.
            Point a = it->vertex(0)->point();
            Point b = it->vertex(1)->point();
            Point c = it->vertex(2)->point();

            // std::cout << "Points of the current triangle :" << std::endl;
            // std::cout << a << std::endl;
            // std::cout << b << std::endl;
            // std::cout << c << std::endl;

            double angle_A = angle(a, b, c);
            double angle_B = angle(c, a, b);
            double angle_C = angle(b, a, c);

            // std::cout << "Angles of the current triangle :" << std::endl;
            // std:: cout << angle_A << std::endl;
            // std::cout << angle_B << std::endl;
            // std::cout << angle_C << std:: endl;

            if (angle_A > 90.0 || angle_B > 90.0 || angle_C > 90.0) {
                // return true;    // an obtuse triangle is found in the instance
                std::cout << "found" << std::endl;  // to test
            }
        }

        return false;   // no obtuse triangle found in the instance
    }



    void insert_steiner(CDT *instance) {
        CDT::Finite_faces_iterator it;
        std::vector<Point> steiner_points;  // to insert the steiner points

        for (it = instance->finite_faces_begin(); it != instance->finite_faces_end(); it++) {
            Point a = it->vertex(0)->point();
            Point b = it->vertex(1)->point();
            Point c = it->vertex(2)->point();

            if (is_obtuse(a, b, c)) {
                Point center = CGAL::centroid(a, b, c);
                steiner_points.push_back(center);
            }
        }
        for (const Point& p: steiner_points) {
            instance->insert(p);
        }
        std::cout << "steiner points inserted are" << steiner_points.size() << std::endl;

    }
};


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

    graph->insert_steiner(&cdt);
    graph->is_obtuse_gen(&cdt);

    return 0;
}

