#include <fstream>
#include <iostream>
#include <filesystem>

#include "db.h"
#include "parse.h"

using namespace std;

int main(int argc, char *argv[]) {
    ios_base::sync_with_stdio(0);
    cin.tie(0);

    if (argc != 5) {
        cout << "********************** [Error] *******************" << endl;
        cout << "*  [Usage] ./partitioner -t <input file> -s <output file> *" << endl;
        cout << "********************** [Error] *******************" << endl;
        return 1;
    }
    fstream input_are,input_info,input_net,input_topo;
    // fstream input(argv[1], std::fstream::in);
    filesystem::path input=argv[2];
    fstream output(argv[4], std::fstream::out);

    topart::DB db;
    //节点的代号减1
    intg hop_limit=topart::parse_input_eda3(db, input, input_are,input_info,input_net,input_topo);
    if(hop_limit==-1){
        cout<<"failed to parser the input file\n";
        return 0;
    }
    db.timer.output_time("FINISH Parsing");

    //生成固定节点
    db.find_root();
    db.timer.output_time("FINISH find_root");    

    db.calculate_candidate_fpga();
    db.timer.output_time("FINISH calculate candidate fpga");
    db.solve_src();
    db.partition();
    db.timer.output_time("FINISH Initial Partition");

    // db.refine_CV();
    // db.mixed_refine(50,0,1);

    // db.refine_CV();
    // db.mixed_refine(50,0,1);

    db.refine_CV();
    db.mixed_refine(50,1);
    // db.multi_refine(4);
    db.timer.output_time("ALL FINISH");
    db.output(output);

    // input.close();
    input_are.close();
    input_info.close();
    input_net.close();
    input_topo.close();
    output.close();
    return 0;
}
