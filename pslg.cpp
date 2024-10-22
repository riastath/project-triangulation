#include "pslg.hpp"

// Constructor
PSLG::PSLG(std::string filename) {
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

// get uid -> probably redundant
std::string PSLG::get_instance_uid() {
    return uid;
}
    
// get steiner points
std::vector<Point> PSLG::get_steiner() {
    return steiner_points;
}

// To print class members -> for testing
void PSLG::printer() {
    std::cout << "uid is: " << uid << std::endl;
    std::cout << "bounds: " << std::endl;
    for (int i = 0; i < bounds.size(); i++) {
        std::cout << bounds.at(i) << " ";
    }
    std::cout << std::endl;
}

void PSLG::delaunay_passer(CDT* delaunay_instance) {
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


// Function to calculate angles using the dot product method. Used to check for non-obtuse triangles later.
double PSLG::angle(Point a, Point b, Point c) {
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
bool PSLG::is_obtuse(Point a, Point b, Point c) {
    double angle_A = angle(a, b, c);
    double angle_B = angle(b, a, c);
    double angle_C = angle(c, a, b);
    if (angle_A > 90.0 || angle_B > 90.0 || angle_C > 90.0) {
        return true;    // triangle is obtuse
    }
    return false;
}


// Version 2 of the obtuse checking function, to check for the actual triangulation (not just a triangle given)
bool PSLG::is_obtuse_gen(CDT* instance) {
    int number_of_obtuce = 0;
    CDT::Finite_faces_iterator it;  // initialize iterator
    for (it = instance->finite_faces_begin(); it != instance->finite_faces_end(); it++) {
        // we need to examine every vertex using the iterator, and find the 3 points of each triangle
        // these examine the 1st, 2nd and 3rd vertex respectively.
        // vertex() returns a handle, which is akin to a pointer to an object,
        // so, using this we get the point with its coordinates (x,y) from the corresponding vertex.
        Point a = it->vertex(0)->point();
        Point b = it->vertex(1)->point();
        Point c = it->vertex(2)->point();

        // std::cout << "Points of the current triangle :" << std::endl;
        // std::cout << a << std::endl;
        // std::cout << b << std::endl;
        // std::cout << c << std::endl;

        double angle_A = angle(a, b, c);
        double angle_B = angle(b, a, c);
        double angle_C = angle(c, a, b);

        // std::cout << "Angles of the current triangle :" << std::endl;
        // std:: cout << angle_A << std::endl;
        // std::cout << angle_B << std::endl;
        // std::cout << angle_C << std:: endl;

        if (angle_A > 90.0 || angle_B > 90.0 || angle_C > 90.0) {
            // return true;    // an obtuse triangle is found in the instance
            std::cout << "found" << std::endl;  // to test
            number_of_obtuce++;
        }
    }
    std::cout << "///// obtuce triangles foud: " << number_of_obtuce << " \\\\\\\\\\" << std::endl; 
    return false;   // no obtuse triangle found in the instance
}

// Simply insert the point in the steiner vector 
void PSLG::insert_steiner_point(Point point) {
    steiner_points.push_back(point);
}

// Steiner point method 1 : insert in center 
void PSLG::insert_steiner_center(CDT *instance) {
    CDT::Finite_faces_iterator it;
    // std::vector<Point> steiner_points;  // to insert the steiner points

    // if (angle(fh->vertex(j)->point(), fh->vertex((j+1)%3)->point(), fh->vertex((j+2)%3)->point()) <= 90) {
                // continue;
            // }

    for (it = instance->finite_faces_begin(); it != instance->finite_faces_end(); it++) {
        // std::cout << "entered loop" << std::endl;
        // Point a = it->vertex(0)->point();
        // Point b = it->vertex(1)->point();
        // Point c = it->vertex(2)->point();
        Point extra;
        int i;
        bool found = false;
        for (i = 0; i < 3; i++) {
            // std::cout << "help" << std::endl;
            if (angle(it->vertex(i)->point(), it->vertex((i+1)%3)->point(), it->vertex((i+2)%3)->point()) <= 90) {
                // std::cout << "it not obtuce" << std::endl;
                continue;
            }
            found = true;
            // std::cout << "found obtuce" << std::endl;
            CDT::Face_handle neigh = it->neighbor(i);
            for (int j = 0; j < 3; j++) {
                if (neigh->vertex(j)->point() == it->vertex(i)->point() || neigh->vertex(j)->point() == it->vertex((i+1)%3)->point() || neigh->vertex(j)->point() == it->vertex((i+2)%3)->point()) {
                    continue;
                }
                extra = neigh->vertex(j)->point();
                break;
            }
            break;
        }
        // std::cout << "exit loop" << std::endl;
        if (found == false) {
            continue;
        }
        Point a = it->vertex(i)->point();
        Point b = it->vertex((i+1)%3)->point();
        Point c = it->vertex((i+2)%3)->point();
        // std::cout << "Points: " << std::endl;
        // std::cout << a << std::endl << b << std::endl << c << std::endl << extra << std::endl;
        Point center = CGAL::centroid(a, b, c, extra);
        // std::cout << "center is: " << center << std::endl;
        insert_steiner_point(center);

        // if (is_obtuse(a, b, c)) {
        //     Point center = CGAL::centroid(a, b, c);
        //     std::cout << "center is: " << center << std::endl;
        //     insert_steiner_point(center);
        // }
    }

    for (const Point& p: steiner_points) {
        instance->insert(p);
    }
    std::cout << "steiner points inserted are" << steiner_points.size() << std::endl;
}

// Steiner point method 2 : insert midpoint
void PSLG::insert_steiner_mid(CDT *instance) {
    CDT::Finite_faces_iterator it;

    for (it = instance->finite_faces_begin(); it != instance->finite_faces_end(); it++) {
        Point a = it->vertex(0)->point();
        Point b = it->vertex(1)->point();
        Point c = it->vertex(2)->point();

        double angle_A = angle(a, b, c);
        double angle_B = angle(b, a, c);
        double angle_C = angle(c, a, b);

        if (angle_A > 90.0) {
            Point mid = CGAL::midpoint(b, c);
            insert_steiner_point(mid);
        } else if (angle_B > 90.0) {
            Point mid = CGAL::midpoint(a, c);
            insert_steiner_point(mid);

        } else if (angle_C > 90.0) {
            Point mid = CGAL::midpoint(a, b);
            insert_steiner_point(mid);
        }
    }

    for (const Point& p: steiner_points) {
        instance->insert(p);
    }
    std::cout << "steiner points inserted are" << steiner_points.size() << std::endl;
}

// Steiner point method 3 : insert a steiner point so that the obtuse angle is bisected
void PSLG::insert_steiner_bisection(CDT *instance) {
    CDT::Finite_faces_iterator it;

    for (it = instance->finite_faces_begin(); it != instance->finite_faces_end(); it++) {
        Point a = it->vertex(0)->point();
        Point b = it->vertex(1)->point();
        Point c = it->vertex(2)->point();

        double angle_A = angle(a, b, c);
        double angle_B = angle(b, a, c);
        double angle_C = angle(c, a, b);

        Point obtuse_vertex;
        Point p_a, p_b;
        double length_a = 0.0, length_b = 0.0;

        if (angle_A > 90.0) {
            // I need ab and ac
            obtuse_vertex = a;
            p_a = b;
            p_b = c;             
        } else if (angle_B > 90.0) {
            // I need ab and bc
            obtuse_vertex = b;
            p_a = a;
            p_b = c; 
        } else if (angle_C > 90.0) {
            // I need ac and bc
            obtuse_vertex = c;
            p_a = a;
            p_b = b;
        } else continue;    // no obtuse angle -> check next face

        // Use CGAL'S squared distance, instead of squared_length for vectors, because now we're using the vertex and points 
        length_a = std::sqrt(CGAL::squared_distance(obtuse_vertex, p_a));
        length_b = std::sqrt(CGAL::squared_distance(obtuse_vertex, p_b));

        // Calculate ratio to find the point, using the math formula / equation with the coordinates of the points
        // (Couldn't find the intersection using CGAL vectors)
        double ratio_a = length_a / (length_a + length_b);
        double ratio_b = length_b / (length_a + length_b);

        Point bisection_point {
            ratio_a * p_a.x() + ratio_b * p_b.x(),
            ratio_a * p_a.y() + ratio_b * p_b.y(),
        };

        insert_steiner_point(bisection_point);
    }

    // Insert all the steiner points at the end 
    for (const Point& p: steiner_points) {
        instance->insert(p);
    }
    std::cout << "steiner points inserted are" << steiner_points.size() << std::endl;
}


bool PSLG::face_is_infinite(CDT::Face_handle face, CDT *instance) {
    for (int i = 0; i < 3; i++) {
        // std::cout << face->vertex(i)->point() << std::endl;
        if (face->vertex(i) == instance->infinite_vertex()) {
            // std::cout << "hello?" << std::endl;
            return true;
        }
    }
    return false;
}

void PSLG::flipper_not_0(CDT *cdt) {
    std::vector<std::pair<CDT::Face_handle, int>> flip_vec; 
    std::vector<CDT::Face_handle> fliped_neighbors;
    // std::pair<Point, Point> constr = std::make_pair(point_vec.at(graph_edge[0]), point_vec.at(graph_edge[1]));

    int i = 1;
    for (CDT::Face_handle fh: cdt->finite_face_handles()) {

        // std::cout << "face " << i << " has neighbors: " << std::endl;
        // for (int j = 0; j < 3; j++) {
        //     std::string res = face_is_infinite(fh->neighbor(j), &cdt) ? "invalid" : "valid";
        //     std::cout << j << ": " << res << std::endl;
        // }
        for (int j = 0; j < 3; j++) {
            // check if neighbor exists to flip
                
            if (!is_obtuse(fh->vertex(j)->point(), fh->vertex((j+1)%3)->point(), fh->vertex((j+2)%3)->point())) {
                continue;
            }

            if (face_is_infinite(fh->neighbor(j), cdt)){
                continue;
            }
                
            // has the face been in a flip
            auto it = find(fliped_neighbors.begin(), fliped_neighbors.end(), fh);
            if (it != fliped_neighbors.end()) {
                std::cout << "face handle allready in vector" << std::endl;
                continue;
            }

            // has the neighbor initiated a flip
            CDT::Face_handle key = fh->neighbor(j);
            auto it2 = find_if(flip_vec.begin(), flip_vec.end(),[key](const auto& p) { return p.first == key; });
            if (it2 != flip_vec.end()) {
                std::cout << "neighbor face handle allready in vector" << std::endl;
                continue;
            }

            // has the neighbor been in a flip
            auto it3 = find(fliped_neighbors.begin(), fliped_neighbors.end(), fh->neighbor(j));
            if (it3 != fliped_neighbors.end()) {
                std::cout << "neighbor face handle has been fliped" << std::endl;
                continue;
            }

            // CDT::Face_handle key = fh->neighbor(j);
            // auto it = find_if(flip_vec.begin(), flip_vec.end(),[key](const auto& p) { return p.first == key; });
            // if (it != flip_vec.end()) {
            //     std::cout << "neighbor face handle allready in vector" << std::endl;
            //     continue;
            // }

            std::pair<CDT::Face_handle, int> flip_item = std::make_pair(fh, j);
            flip_vec.push_back(flip_item);
            fliped_neighbors.push_back(fh->neighbor(j));
        }
        // std::cout << "### end ###" << std::endl;
    }

    std::cout << "to flip:" << flip_vec.size() << std::endl;
    for (int i = 0; i < flip_vec.size(); i++) {
        std::cout << "flip" << std::endl;
        std::pair<CDT::Face_handle, int> flip_item = flip_vec.at(i);
        std::cout << flip_item.first->vertex(0)->point() << "  " 
        << flip_item.first->vertex(1)->point() << "  " 
        << flip_item.first->vertex(2)->point() << std::endl;
        cdt->flip(flip_item.first, flip_item.second);
    }
    std::cout << "##### end of flips #####" << std::endl;
}

// Function that produces the final .json output with the data as requested
void PSLG::produce_output() {
    pt::ptree root;
    root.put("content_type", "CG_SHOP_2025_Solution");
    root.put("instance_uid", uid);
        
    // trees to insert the rest of the data
    pt::ptree steiner_points_x, steiner_points_y, edges;

    for (Point point : steiner_points) {

        pt::ptree x_node, y_node;       // different nodes for each coordinate, to make the pairs with the empty string to insert
        x_node.put("",point.x());
        y_node.put("", point.y());

        // this way, the points are shown as in the input .json files
        steiner_points_x.push_back(std::make_pair("", x_node));
        steiner_points_y.push_back(std::make_pair("", y_node));
    }

    // add these points to the root of the tree 
    root.add_child("steiner_points_x", steiner_points_x);
    root.add_child("steiner_points_y", steiner_points_y);

    // now, create the output edges
    for (std::pair<Point, Point> edge : additional_constraints) {
        pt::ptree edge_node;
        pt::ptree p1_tree, p2_tree; // p1 and p2 are basically the positions of the two points

        // we need to find where the two points are in the point vector 
        auto it1 = std::find(point_vec.begin(), point_vec.end(), edge.first);
        auto it2 = std::find(point_vec.begin(), point_vec.end(), edge.second);
        int p1_pos = std::distance(point_vec.begin(), it1);
        int p2_pos = std::distance(point_vec.begin(), it2);

        p1_tree.put("", p1_pos);
        p2_tree.put("", p2_pos);
        edge_node.push_back(std::make_pair("", p1_tree));
        edge_node.push_back(std::make_pair("", p2_tree));

        edges.push_back(std::make_pair("", edge_node));
    }

    root.add_child("edges", edges);

    write_json(std::cout, root);    // with this, we print the output (can also redirect to an output file)
}
