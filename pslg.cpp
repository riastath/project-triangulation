#include "pslg.hpp"
#include <gmp.h>

// Constructor
PSLG::PSLG(std::string filename) {
    pt::ptree root;
    pt::read_json(filename, root);

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

    // Create vector with graph points
    for (int i = 0; i < points_x.size(); i++) {
        point_vec.push_back(Point(points_x.at(i), points_y.at(i)));
    }


    // ====== convex hull start ======
    std::vector<Point_2> hprev;

    for (int &point: bounds) {
        // std::cout << "point is:" << point_vec[point] << std::endl;
        hprev.push_back(Point_2(point_vec[point].x(), point_vec[point].y()));
    }

    std::vector<Point_2> hull;
    CGAL::convex_hull_2(hprev.begin(), hprev.end(), std::back_inserter(hull));
    // std::cout << "convex hull points: " << hull.size() << std::endl;

    if (hull.size() != bounds.size()) {
        std::cout << "PSLG is not convex" << std::endl;
    }
    else {
        bool conv = true;
        for (int i = 0; i < hull.size(); i++) {
            if (hull[i] != point_vec[bounds[i]]) {
                std::cout << hull[i] << " != " << point_vec[bounds[i]] << std::endl;
                conv = false;
                break;
            }
        }
        if (conv) {
            std::cout << "is convex" << std::endl;
        }
    }
    // ====== convex hull end ======



    // Add constraints
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

    // New parameters added for part 2 : method, delaunay flag, parameters of algorithm
    method = root.get<std::string>("method");
    delaunay = root.get<bool>("delaunay");

    // Parameters will be different for each algorithm
    for(pt::ptree::value_type &param : root.get_child("parameters")) {
        parameters[param.first] = param.second.get_value<double>();
    }

}


bool PSLG::has_circles() {
    std::vector<Point> visited;
    return false;
}

// Returns the cdt after inserting (steiner) point
CDT PSLG::return_cdt(CDT cdt, Point point) {
    cdt.insert(point);
    return cdt;
}

// Getter functions
int PSLG::get_num_steiner_points() {
    return steiner_points.size();
}

std::string PSLG::get_method() {
    return method;
}

std::map<std::string, double> PSLG::get_parameters() {
    return parameters;
}

bool PSLG::get_delaunay_flag() {
    return delaunay;
}

// To pass all points and edges to the delaunay triangulation instance
void PSLG::delaunay_passer(CDT* delaunay_instance) {
    // pass the points delaunay instance
    for (const Point& p: point_vec){
        delaunay_instance->insert(p);
    }

    // pass the constraints to delaunay istance
    for (const auto& constraint: additional_constraints) {
        delaunay_instance->insert_constraint(constraint.first, constraint.second);
    }
}


// Function to calculate angles using the dot product method. Used to check for non-obtuse triangles later.
double PSLG::angle(Point a, Point b, Point c) {
    // calculating the vectors needed for the dot product formula
    Vector u = b - a;
    Vector v = c - a;
    
    K::FT prod = u * v;
    double product = CGAL::to_double(prod);
    // |u| and |v|. no length function exists, so using square length and taking its root
    double ulength = std::sqrt(CGAL::to_double(u.squared_length()));
    double vlength = std::sqrt(CGAL::to_double(v.squared_length()));

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
int PSLG::is_obtuse_gen(CDT* instance) {
    int number_of_obtuse = 0;
    CDT::Finite_faces_iterator it;
    for (it = instance->finite_faces_begin(); it != instance->finite_faces_end(); it++) {
        // we need to examine every vertex using the iterator, and find the 3 points of each triangle. these examine the 1st, 2nd and 3rd vertex respectively.
        // vertex() returns a handle, which is akin to a pointer to an object,
        // so, using this we get the point with its coordinates (x,y) from the corresponding vertex.
        Point a = it->vertex(0)->point();
        Point b = it->vertex(1)->point();
        Point c = it->vertex(2)->point();

        double angle_A = angle(a, b, c);
        double angle_B = angle(b, a, c);
        double angle_C = angle(c, a, b);

        if (angle_A >= 180.0 || angle_B >= 180.0 || angle_C >= 180.0) {
            continue;
        }

        if (angle_A > 90.0 || angle_B > 90.0 || angle_C > 90.0) {
            number_of_obtuse++;
        }
    }
    return number_of_obtuse;
}

// Function to check for obtuse triangles using face handle
bool PSLG::is_obtuse_face(CDT::Face_handle f) {

    Point point_a = (Point)f->vertex(0)->point();
    Point point_b = (Point)f->vertex(1)->point();
    Point point_c = (Point)f->vertex(2)->point();

    if (is_obtuse(point_a, point_b, point_c)) return true;
    else return false;
}


// Simply insert the point in the steiner vector 
void PSLG::insert_steiner_point(Point point) {
    steiner_points.push_back(point);
}

// Steiner point method 1 : insert in center 
std::pair<Point, int> PSLG::insert_steiner_center(CDT_C instance, CDT::Face_handle face, int num_obtuse) {
    Point extra;
    int i;
    bool found = false;
    for (i = 0; i < 3; i++) {
        // reject non obtuce angles of face
        if (angle(face->vertex(i)->point(), face->vertex((i+1)%3)->point(), face->vertex((i+2)%3)->point()) <= 90) {
            continue;
        }

        found = true;
        CDT::Face_handle neigh = face->neighbor(i);
        for (int j = 0; j < 3; j++) { // find the neighboring face's vertex that is not common with the current face
            if (neigh->vertex(j)->point() == face->vertex(i)->point() || neigh->vertex(j)->point() == face->vertex((i+1)%3)->point() || neigh->vertex(j)->point() == face->vertex((i+2)%3)->point()) {
                continue;
            }
            extra = neigh->vertex(j)->point();
            break;
        }
        break;
    }
    if (found == false) {
        return std::make_pair(Point(NULL, NULL), -1); // return this value if the function was unable to find the steiner point
    }
    Point a = face->vertex(i)->point();
    Point b = face->vertex((i+1)%3)->point();
    Point c = face->vertex((i+2)%3)->point();

    Point center = CGAL::centroid(a, b, c, extra);
    insert_all_steiner(&instance);
    instance.insert(center);
    int num_after = is_obtuse_gen(&instance);

    int improvement = num_obtuse - num_after;
    return std::make_pair(center, improvement); // return steiner point and the improvement it will make if inserted

}

// Steiner point method 1.1 : insert in center, single triangle
std::pair<Point, int> PSLG::insert_steiner_center_single(CDT_C instance, CDT::Face_handle face, int num_obtuce) {
    int i;
    bool found = false;
    for (i = 0; i < 3; i++) {
        // reject non obtuce angles of face
        if (angle(face->vertex(i)->point(), face->vertex((i+1)%3)->point(), face->vertex((i+2)%3)->point()) <= 90) {
            continue;
        }

        found = true;
        break;
    }
    if (found == false) {
        return std::make_pair(Point(NULL, NULL), -1); // return this value if the function was unable to find the steiner point
    }
    Point a = face->vertex(i)->point();
    Point b = face->vertex((i+1)%3)->point();
    Point c = face->vertex((i+2)%3)->point();

    Point center = CGAL::centroid(a, b, c);
    insert_all_steiner(&instance);
    instance.insert(center);
    int num_after = is_obtuse_gen(&instance);

    int improvement = num_obtuce - num_after;
    return std::make_pair(center, improvement);
}

Point centroid_custom(std::vector<Point> points) {
    int n = points.size();
    long double A = 0;
    Point Pi,Pi1;
    for (int i = 0; i < n; i++) {
        Pi = points.at(i);
        Pi1 = points.at((i+1)%n);
        A += (CGAL::to_double(Pi.x())*CGAL::to_double(Pi1.y()) - CGAL::to_double(Pi1.x())*CGAL::to_double(Pi.y()))/2; 
    }
    long double x = 0, y = 0;
    long double mult;
    for (int i = 0; i < n; i++) {
        Pi = points.at(i);
        Pi1 = points.at((i+1)%n);
        mult = (CGAL::to_double(Pi.x())*CGAL::to_double(Pi1.y()) - CGAL::to_double(Pi1.x())*CGAL::to_double(Pi.y()));
        x += ((CGAL::to_double(Pi.x()) + CGAL::to_double(Pi1.x())) * mult)/(6*A);
        y += ((CGAL::to_double(Pi.y()) + CGAL::to_double(Pi1.y())) * mult)/(6*A);
    }

    double x_final, y_final;
    x_final = x;
    y_final = y;
    return Point(x_final, y_final);
}

// Steiner point method 1.2 : insert in center, check and include neighbors 
std::pair<CDT_C, Point> PSLG::insert_steiner_center_neighbors(CDT_C instance, CDT::Face_handle face, int num_obtuse, int* res) {

    Point polygon_points[6]; // 3 points for the triangle, and 3 more if we include all neighbors.
    for (int i = 0; i < 6; i++) {
        polygon_points[i] = Point(NULL, NULL);
    }

    int i;
    bool found = false;
    for (i = 0; i < 3; i++) {
        // reject non obtuce angles of face
        if (angle(face->vertex(i)->point(), face->vertex((i+1)%3)->point(), face->vertex((i+2)%3)->point()) <= 90) {
            continue;
        }

        found = true;
        break;
    }
    if (found == false) {
        *res = -1;
        return std::make_pair(instance, Point(NULL, NULL)); // return if triangle is not obtuce
    }
    // insert the obtuce triangle's points to polygon points
    polygon_points[0] = face->vertex(i)->point(); // 0 is the obtuce angle
    polygon_points[1] = face->vertex((i+1)%3)->point();
    polygon_points[2] = face->vertex((i+2)%3)->point();


    // check the 3 neighbors
    CDT::Face_handle neigh;
    for (int j = 0; j < 3; j++) {
        neigh = face->neighbor(j);

        // check the neighbor's points
        for (int k = 0; k < 3; k++) {
            // if common point: ignore
            if (neigh->vertex(k)->point() == face->vertex(i)->point() || neigh->vertex(k)->point() == face->vertex((i+1)%3)->point() || neigh->vertex(k)->point() == face->vertex((i+2)%3)->point()) {
                continue;
            }
            // if non-obtuce: ignore
            if (angle(neigh->vertex(k)->point(), neigh->vertex((k+1)%3)->point(), neigh->vertex((k+2)%3)->point()) <= 90) {
                continue;
            }
            polygon_points[3+j] = neigh->vertex(k)->point();
            break;
        }
    }

    
    // constraints logic
    /* if 3 not exists: 1-2, else: 1-3, 3-2
       if 4 not exists: 2-0, else: 2-4, 4-0
       if 5 not exists: 0-1, else: 0-5, 5-1*/
    // insert constraints
    Point a, b, current;
    std::vector<std::pair<Point, Point>> constr_vec;
    for (int j = 3; j < 6; j++) {
        current = polygon_points[j];
        a = polygon_points[(j-2)%3];
        b = polygon_points[((j-2)%3 +1)%3];
        if (current == Point(NULL, NULL)) {
            // constraint a - b
            constr_vec.push_back(std::make_pair(a, b));
        }
        else {
            // constraint a - current
            // constraint current - b
            constr_vec.push_back(std::make_pair(a, current));
            constr_vec.push_back(std::make_pair(current, b));
        }
    }


    CDT::Face_handle it = face;

    // identify points to remove
    bool to_remove[6] = {false, false, false, false, false, false};
    for (int j = 3; j < 6; j++) {
        current = polygon_points[j];
        int a = (j-2)%3;
        int b = (a+1)%3;
        if (current != Point(NULL, NULL)) {
            to_remove[a] = true;
            to_remove[b] = true;
            to_remove[j] = true;
        }
    }
    // remove points
    int number_to_remove = 0;
    for (int j = 0; j < 6; j++) {
        if (to_remove[j] == true) {
            number_to_remove++;
        }
    }
    std::vector<CDT::Vertex_handle> remove_vect; // we need vertex handles to remove points, and faces change after every operation
    CDT::Finite_faces_iterator temp_face;
    for (temp_face = instance.finite_faces_begin(); temp_face != instance.finite_faces_end(); temp_face++) {
        bool cont = false;
        for (int k = 0; k < 6; k++) {
            if (!to_remove[k]) {
                continue;
            }
            Point temp_point = polygon_points[k];
            if (temp_face->vertex(0)->point() == temp_point) {
                if (instance.are_there_incident_constraints(temp_face->vertex(0)) == false) {
                    remove_vect.push_back(temp_face->vertex(0));
                }
                cont = true;
                to_remove[k] = false;
            }
            else if (temp_face->vertex(1)->point() == temp_point) {
                if (instance.are_there_incident_constraints(temp_face->vertex(1)) == false) {
                    remove_vect.push_back(temp_face->vertex(1));
                }
                cont = true;
                to_remove[k] = false;
            }
            else if (temp_face->vertex(2)->point() == temp_point) {
                if (instance.are_there_incident_constraints(temp_face->vertex(2)) == false) {
                    remove_vect.push_back(temp_face->vertex(2));
                }
                cont = true;
                to_remove[k] = false;
            }
        }
        if (cont == false) {
            continue;
        }
        if (remove_vect.size() == number_to_remove) {
            break;
        }
    }
    for (CDT::Vertex_handle vh : remove_vect) {
        instance.remove(vh);
    }
    

    // calculate centroid
    std::vector<Point> points;
    for (int j = 0; j < 6; j++) {
        if (polygon_points[j] == Point(NULL, NULL)) {
            continue;
        }
        points.push_back(polygon_points[j]);
    }
    Point center = centroid_custom(points);

    // insert centroid    
    instance.insert_no_flip(center);


    // re-insert removed points in constraints
    for (std::pair<Point, Point> pconstr : constr_vec) {
        instance.insert_constraint(pconstr.first, pconstr.second);
    }


    // remove inserted constraints without removing their points
    for (temp_face = instance.finite_faces_begin(); temp_face != instance.finite_faces_end(); temp_face++) {
        for (int k = 0; k < 3; k++) {
            for (int ve = 0; ve < constr_vec.size(); ve++) {
                std::pair<Point, Point> pconstr = constr_vec.at(ve);
                if ((pconstr.first == temp_face->vertex((k+1)%3)->point() && pconstr.second == temp_face->vertex((k+2)%3)->point())||(pconstr.first == temp_face->vertex((k+2)%3)->point() && pconstr.second == temp_face->vertex((k+1)%3)->point())) {
                    instance.remove_constrained_edge(temp_face, k);
                    constr_vec.erase(constr_vec.begin() + ve);
                }
            }
            if (constr_vec.size() == 0) {
                break;
            }
        }
        if (constr_vec.size() == 0) {
            break;
        }
    }

    int improvement = num_obtuse - is_obtuse_gen(&instance);
    *res = improvement;
    return std::make_pair(instance, center);
}

// Steiner point method 2 : insert midpoint
std::pair<Point, int> PSLG::insert_steiner_mid(CDT_C instance, CDT::Face_handle face, int num_obtuse) {
    Point a = face->vertex(0)->point();
    Point b = face->vertex(1)->point();
    Point c = face->vertex(2)->point();
    Point mid = Point(NULL, NULL);

    double angle_A = angle(a, b, c);
    double angle_B = angle(b, a, c);
    double angle_C = angle(c, a, b);

    insert_all_steiner(&instance);
    if (angle_A > 90.0) {
        mid = CGAL::midpoint(b, c);
        instance.insert(mid);
    } else if (angle_B > 90.0) {
        mid = CGAL::midpoint(a, c);
        instance.insert(mid);
    } else if (angle_C > 90.0) {
        mid = CGAL::midpoint(a, b);
        instance.insert(mid);
    } else {
        return std::make_pair(Point(NULL, NULL), -1); // no obtuse angle
    }

    int num_after = is_obtuse_gen(&instance);
    int improvement = num_obtuse - num_after;

    return std::make_pair(mid, improvement);
}

// Steiner point method 3 : insert a steiner point so that the obtuse angle is bisected
std::pair<Point, int> PSLG::insert_steiner_bisection(CDT_C instance, CDT::Face_handle face, int num_obtuse) {
    Point a = face->vertex(0)->point();
    Point b = face->vertex(1)->point();
    Point c = face->vertex(2)->point();

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
    } else {
        return std::make_pair(Point(NULL, NULL), -1);
    }

    // Use CGAL'S squared distance, instead of squared_length for vectors, because now we're using the vertex and points 
    length_a = std::sqrt(CGAL::to_double(CGAL::squared_distance(obtuse_vertex, p_a)));
    length_b = std::sqrt(CGAL::to_double(CGAL::squared_distance(obtuse_vertex, p_b)));

    // Calculate ratio to find the point, using the math formula / equation with the coordinates of the points
    // (Couldn't find the intersection using CGAL vectors)
    // Added this check later : just in case, to avoid division by zero
    if (length_a + length_b == 0) {
        return std::make_pair(Point(NAN, NAN), -1);
    }
    double ratio_a = length_a / (length_a + length_b);
    double ratio_b = length_b / (length_a + length_b);

    Point bisection_point {
        ratio_a * p_a.x() + ratio_b * p_b.x(),
        ratio_a * p_a.y() + ratio_b * p_b.y(),
    };
    insert_all_steiner(&instance);
    instance.insert(bisection_point);
    int num_after = is_obtuse_gen(&instance);

    int improvement = num_obtuse - num_after;
    return std::make_pair(bisection_point, improvement);
}

std::pair<Point, int> PSLG::insert_steiner_projection(CDT_C instance, CDT::Face_handle face, int num_obtuse) {
    Point projection_point;
    int improvement;
    bool found = false;
    for (int j = 0; j < 3; j++) {
        Point p0 = face->vertex(j)->point();
        Point p1 = face->vertex((j+1)%3)->point();
        Point p2 = face->vertex((j+2)%3)->point();

        if (angle(p0,p1,p2) <= 90) {
            continue;
        }

        double mult = 0.000000001; // offset because cgal is cgal
        int dirx = 1;
        int diry = 1;
        double x,y;
        if (p2.x() == p1.x()) { // if line equation is: x=b
            x = CGAL::to_double(p1.x());
            y = CGAL::to_double(p0.y());
        }
        else if (p2.y() == p1.y()) { // if line equation in y=b
            x = CGAL::to_double(p0.x());
            y = CGAL::to_double(p1.y());
        }
        else {
            // calculate line equations
            double slope = (CGAL::to_double(p2.y()) - CGAL::to_double(p1.y()))/CGAL::to_double((p2.x()) - CGAL::to_double(p1.x()));
            double slope_per = -1/slope;

            
            double b = CGAL::to_double(p1.y()) - slope * CGAL::to_double(p1.x());
            double b_per = CGAL::to_double(p0.y()) - slope_per * CGAL::to_double(p0.x());

            // find intersection point
            x = (b_per - b) / (slope - slope_per);
            y = slope_per * x + b_per;

            // modifying points (extend slightly away from the obtuce point, while on the perpendicular line)
            // find direction for each coordinate:
            dirx = (x - CGAL::to_double(p0.x())) / abs(x - CGAL::to_double(p0.x()));
            diry = (y - CGAL::to_double(p0.y())) / abs(y - CGAL::to_double(p0.y()));
            // add apropriate number to extend
            x += dirx * 1 * mult;
            y += diry * abs(slope_per) * mult;
        }

        projection_point = {x,y};
        found = true;
    }

    if (found == false) {
        return std::make_pair(Point(NULL,NULL), -1);
    }

    insert_all_steiner(&instance);
    instance.insert(projection_point);
    int num_after = is_obtuse_gen(&instance);

    improvement = num_obtuse - num_after;

    return std::make_pair(projection_point, improvement);
}

std::pair<Point, int> PSLG::insert_steiner_circumcenter(CDT_C instance, CDT::Face_handle face, int num_obtuse) {
    Point a = face->vertex(0)->point();
    Point b = face->vertex(1)->point();
    Point c = face->vertex(2)->point();

    Point circumcenter = CGAL::circumcenter(a, b, c);

    Triangle triangle(a, b, c);
    if (!triangle.has_on_bounded_side(circumcenter)) {  // checking if it's within boundaries
        return std::make_pair(Point(NULL, NULL), -1);
    }

    // circumcenter is valid
    instance.insert_no_flip(circumcenter);

    int num_after = is_obtuse_gen(&instance);
    int improvement = num_obtuse - num_after;

    return std::make_pair(circumcenter, improvement);
}

void PSLG::insert_all_steiner(CDT_C *cdt) {
    for (const Point& p: steiner_points) {
        cdt->insert(p);
    }
}

// Check if a triangle is infinite (through face handle)
bool PSLG::face_is_infinite(CDT::Face_handle face, CDT *instance) {
    for (int i = 0; i < 3; i++) {
        // if a triangle includes atleast one infinite vertex, then it is infinite
        if (face->vertex(i) == instance->infinite_vertex()) {
            return true;
        }
    }
    return false;
}

void PSLG::flip_edges(CDT *cdt) {
    // faces that initiated a flip
    std::vector<std::pair<CDT::Face_handle, int>> flip_vec; 
    // faces that were fliped by another face
    std::vector<CDT::Face_handle> fliped_neighbors;

    int i = 1;
    for (CDT::Face_handle fh: cdt->finite_face_handles()) {

        for (int j = 0; j < 3; j++) {
            
            // flip only if obtuse angle exists    
            if (!is_obtuse(fh->vertex(j)->point(), fh->vertex((j+1)%3)->point(), fh->vertex((j+2)%3)->point())) {
                continue;
            }

            // check if neighbor exists to flip (neghbor != infinite face)
            if (face_is_infinite(fh->neighbor(j), cdt)){
                continue;
            }
                
            // check if the face has been involved in a flip
            auto it = find(fliped_neighbors.begin(), fliped_neighbors.end(), fh);
            if (it != fliped_neighbors.end()) {
                continue;
            }

            // check if the neighbor has initiated a flip
            CDT::Face_handle key = fh->neighbor(j);
            auto it2 = find_if(flip_vec.begin(), flip_vec.end(),[key](const auto& p) { return p.first == key; });
            if (it2 != flip_vec.end()) {
                continue;
            }

            // check if the neighbor has been involved in a flip
            auto it3 = find(fliped_neighbors.begin(), fliped_neighbors.end(), fh->neighbor(j));
            if (it3 != fliped_neighbors.end()) {
                continue;
            }

            // insert faces to apropriate vectors to flip
            std::pair<CDT::Face_handle, int> flip_item = std::make_pair(fh, j);
            flip_vec.push_back(flip_item);
            fliped_neighbors.push_back(fh->neighbor(j));
        }

    }

    // flip faces 
    for (int i = 0; i < flip_vec.size(); i++) {
        std::pair<CDT::Face_handle, int> flip_item = flip_vec.at(i);
        cdt->flip(flip_item.first, flip_item.second);
    }
    
    std::cout << "Edges fliped are: " << flip_vec.size() << std::endl;
}

void print_rational(const K::FT& coord, unsigned long* a, unsigned long* b) {
    const auto exact_coord = CGAL::exact(coord);

    // Convert the exact coordinate to a GMP rational (mpq_t)
    const mpq_t* gmpq_ptr = reinterpret_cast<const mpq_t*>(&exact_coord);

    // Declare GMP integers to hold the numerator and denominator
    mpz_t num, den;
    mpz_init(num);
    mpz_init(den);

    // Extract the numerator and denominator using GMP functions
    mpq_get_num(num, *gmpq_ptr);  // Get the numerator
    mpq_get_den(den, *gmpq_ptr);  // Get the denominator

    // Print the numerator and denominator
    *a = *num->_mp_d;
    *b = *den->_mp_d;

    // Clear GMP integers
    mpz_clear(num);
    mpz_clear(den);
} 

// Function that produces the final .json output with the data as requested
void PSLG::produce_output(CDT instance) {
    pt::ptree root;
    root.put("content_type", "CG_SHOP_2025_Solution");
    root.put("instance_uid", uid);
        
    // trees to insert the rest of the data
    pt::ptree steiner_points_x, steiner_points_y, edges;
    unsigned long a,b;
    char ch[50];

    for (Point point : steiner_points) {
        
        pt::ptree x_node, y_node;       // different nodes for each coordinate, to make the pairs with the empty string to insert
        
        print_rational(point.x(), &a, &b);
        if (b == 1) {
            x_node.put("", point.x());
        }
        else {
            sprintf(ch, "%lu/%lu", a, b);
            x_node.put("", ch);
        }
        
        print_rational(point.y(), &a, &b);
        if (b == 1) {
            y_node.put("", point.y());
        }
        else {
            sprintf(ch, "%lu/%lu", a, b);
            y_node.put("", ch);
        }


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

    // New class members, adding to output
    int obtuse_count = is_obtuse_gen(&instance);
    root.put("obtuse_count", obtuse_count);
    root.put("method", method);;

    pt::ptree param_tree;                  // for the parameters
    for (std::pair<const std::string, double> &param : parameters) {
        param_tree.put(param.first, param.second);
    }
    root.add_child("parameters", param_tree);

    write_json("output.json", root);    // redirect output to a file 
}