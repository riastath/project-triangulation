#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"
#include <string.h>
#include <vector>
#include <iostream>
using namespace std;
namespace pt = boost::property_tree;

class PSLG {
    string uid;
    std::vector<int> bounds;

public:
    PSLG(string name) {
        pt::ptree root;
        pt:read_json(name, root);
        uid = root.get<string>("instance_uid");

        for (pt::ptree::value_type &bound : root.get_child("region_boundary")) {
            bounds.push_back(bound.second.get_value<int>());
        }
    }

    void printer() {
        cout << "uid is: " << uid << endl;
        cout << "bounds: " << endl;
        for (int i = 0; i < bounds.size(); i++) {
            std::cout << bounds.at(i) << " ";
        }
        std::cout << std::endl;
    }
};



int main(void) {
    PSLG * graph = new PSLG("tests/data.json");
    graph->printer();

    return 0;
}