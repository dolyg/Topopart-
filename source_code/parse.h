#ifndef PARSE_H_
#define PARSE_H_


#include <fstream>
#include <filesystem>

#include "db.h"

using std::fstream;

namespace topart {

void parse_input(DB &db, fstream &input);

namespace fs = std::filesystem;
int parse_input_eda3(DB &db,fs::path &directory,std::fstream &input_are,std::fstream &input_info,
                    std::fstream &input_net,std::fstream &input_topo);

}  // namespace topart

#endif  // PARSE_H_
