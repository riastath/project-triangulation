#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"
#include <iostream>
#include <vector>
// using namespace boost::property_tree;

// Short alias for this namespace
namespace pt = boost::property_tree;


int main(void) {
    // Create a root
    pt::ptree root;

    // Load the json file in this ptree
    pt::read_json("data.json", root);


    std::string id = root.get<std::string>("instance_uid");
    std::cout << id << std::endl;

    std::vector<int> bounds;
    for (pt::ptree::value_type &bound : root.get_child("region_boundary")) {
        bounds.push_back(bound.second.get_value<int>());
    }

    for (int i = 0; i < bounds.size(); i++) {
        std::cout << bounds.at(i) << " ";
    }
    std::cout << std::endl;
    // for (pt::ptree::value_type &fruit : root.get_child("fruits")){
    //     // fruit.first contain the string ""
    //     bounds.push_back(fruit.second.data());
    // }

}
