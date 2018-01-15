#include <Rcpp.h>
#include <string>
#include <vector>
//#include <map>

// typedef Rcpp::LogicalVector lv_type;
//typedef Rcpp::LogicalVector::const_iterator::difference_type lv_ptrdiff_t;

//typedef std::vector<double> c_type;
////typedef std::string o_key_type;
//typedef std::vector<double> o_value_type;
////typedef std::pair<o_key_type, o_value_type> o_allocator_type;
////typedef std::unordered_map<o_key_type, o_value_type> o_type;
//typedef std::vector<o_value_type> o_type;

// GCC-4.7 doesn't have std::unordered_map from C++11 so we have to use tr1's unordered map instead
// GCC-4.7:
/* 
#include <tr1/unordered_map>
typedef std::tr1::unordered_map<std::string, std::size_t> my_universal_unordered_map;
*/

// Others:
// #include <unordered_map>
// typedef std::unordered_map my_universal_unordered_map;

// Or maybe with boost?
#include <boost/unordered_map.hpp> 
typedef boost::unordered_map<std::string, std::size_t> my_universal_unordered_map;
