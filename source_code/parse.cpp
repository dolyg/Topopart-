#include "parse.h"

#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <filesystem>
#include <regex>

namespace topart {

using std::getline;
using std::stringstream;
using std::istringstream;
using std::pair;

//#define LOG 0
//解析EDA3的输入片段
namespace fs = std::filesystem;
void parser_are(std::fstream &input_are,std::map<intg,std::array<intg,8>> &node_capacity,intg &total_node){
    std::string line;
    while (getline(input_are, line)) { // 逐行读取
        istringstream iss(line);
        intg index;
        char ttmp;
        iss>>ttmp>>index;
        index-=1;           //编号缩减
        total_node++;
        for(int i=0;i<8;++i){
            intg source_cnt=0;
            iss>>source_cnt;
            node_capacity[index][i]=source_cnt;
        }
    }   
}
void parser_net(std::fstream &input_net,vector<vector<intg>> &g_circuit){
    vector<intg> tmp;
    std::string line;
    while (getline(input_net, line)) { // 逐行读取
        tmp.clear();
        istringstream iss(line);
        char ttmp;              //用来暂存g的
        iss>>ttmp;
        intg node;              //源点编号
        iss>>node;
        node-=1;                //编号缩减
        tmp.emplace_back(node);     //塞入源点

        intg weight; //获取超边权重
        iss>>weight;

        while(iss>>ttmp&&iss>>node)
            tmp.emplace_back(node-1);

        tmp.emplace_back(weight);
        g_circuit.emplace_back(tmp);  
    }
}
void parser_topo(std::fstream &input_topo,vector<pair<intg, intg>> &g_fpga,intg &hop_limit){
    //获取hop_limit
    input_topo>>hop_limit;
    // NOTE: weird
    string s;
    getline(input_topo, s);

    std::string line;
    while (getline(input_topo, line)) { // 逐行读取
        istringstream iss(line);
        intg u,v;
        char ttmp;
        iss>>ttmp>>ttmp>>ttmp>>ttmp>>u>>ttmp>>ttmp>>ttmp>>ttmp>>v;
        u-=1;v-=1;            //编号缩减
        g_fpga.emplace_back(std::make_pair(u, v)); 
    }
}
void parser_info(std::fstream &input_info,std::map<intg,std::array<intg,8>> &fpga_capacity,intg &total_fpga,
                std::map<intg,intg> &wire_capacity){
    std::string line;
    while (getline(input_info, line)) { // 逐行读取
        istringstream iss(line);
        char ttmp;
        intg index;
        iss>>ttmp>>ttmp>>ttmp>>ttmp>>index;
        index-=1;                //编号缩减
        total_fpga++;
        intg source_cnt=0;
        iss>>source_cnt;
        wire_capacity[index]=source_cnt;
        for(int i=0;i<8;i++){
            iss>>source_cnt;
            fpga_capacity[index][i]=source_cnt;
        }
    }   
}
void debug_are(const std::map<intg,std::array<intg,8>> &node_capacity,const intg & total_node);
void debug_info(const std::map<intg,std::array<intg,8>> &fpga_capacity,const intg &total_fpga,
                const std::map<intg,intg> &wire_capacity);
void debug_topo(const vector<pair<intg, intg>> &g_fpga,const intg &hop_limit);
void debug_net(const vector<vector<intg>> &g_circuit);


int parse_input_eda3(DB &db,fs::path &directory,std::fstream &input_are,std::fstream &input_info,
                    std::fstream &input_net,std::fstream &input_topo){ 
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            if(entry.path().extension()==".are")
                input_are.open(entry.path(), std::fstream::in);
            if(entry.path().extension()==".info")
                input_info.open(entry.path(), std::fstream::in);
            if(entry.path().extension()==".net")
                input_net.open(entry.path(), std::fstream::in);
            if(entry.path().extension()==".topo")
                input_topo.open(entry.path(), std::fstream::in);
        }
    }
    if (!(input_are.is_open()&&input_info.is_open()&&input_net.is_open()&&input_topo.is_open())){
        std::cout<<"one file is not opened";
        return -1;
    }

    //are  节点gx 八个数字代表八种资源
    std::map<intg,std::array<intg,8>> node_capacity;
    intg total_node=0;
    parser_are(input_are,node_capacity,total_node);
    //debug_are(node_capacity,total_node);

    //net   源节点gx 超边权重数字 被驱动点gx 暂且把超图存成简单的图
    vector<vector<intg>> g_circuit;     //最后一个位置存储的是权重
    parser_net(input_net,g_circuit);
    //debug_net(g_circuit);    

    //info  节点FPGAx 八个数字代表八种资源
    std::map<intg,intg> wire_capacity;      //最大联通量
    std::map<intg,std::array<intg,8>> fpga_capacity;
    intg total_fpga=0;
    parser_info(input_info,fpga_capacity,total_fpga,wire_capacity);
    //debug_info(fpga_capacity,total_fpga,wire_capacity);
    //topo 开头一个数字代表最大hop,然后FPGAx FPGAx代表谁和谁相连 
    vector<pair<intg, intg>> g_fpga;
    intg hop_limit;
    parser_topo(input_topo,g_fpga,hop_limit);
    //debug_topo(g_fpga,hop_limit);


    db.build_fpga_graph(total_fpga, fpga_capacity, g_fpga,wire_capacity,hop_limit);
    db.build_circuit_graph(total_node, g_circuit, node_capacity,total_fpga);
    //结束
    return hop_limit;
}

//用来调试的函数 
void debug_net(const vector<vector<intg>> &g_circuit){
    #ifdef LOG_PARSE
    std::cout << "=================Circuit===========================\n" ;
    for(auto circuit:g_circuit){
        std::cout<<"g"<<circuit[0]+1<<" ";
        std::cout<<circuit[circuit.size()-1]<<" ";
        for(int i=1;i<circuit.size()-1;i++){
            //std::cout<<"g"<<circuit[i]+1<<(i==circuit.size()-2?"":" ");
            std::cout<<"g"<<circuit[i]+1<<" ";
        }
        std::cout<<"\n";
    }
    std::cout << "=============================================\n" ;
    #endif   
}
void debug_info(const std::map<intg,std::array<intg,8>> &fpga_capacity,const intg &total_fpga,
                const std::map<intg,intg> &wire_capacity){
    #ifdef LOG_PARSE
    std::cout << "=============================================\n" ;
    std::cout<<"fpga's number is "<<total_fpga<<'\n';
    for(int i=0;i<total_fpga;i++){
        std::cout<<"fpga"<<i+1<<" ";
        std::cout<<wire_capacity.at(i)<<" ";
        for(int j=0;j<8;j++){
            std::cout<<fpga_capacity.at(i).at(j)<<" ";
        }
        std::cout<<"\n";
    }
    std::cout << "=============================================\n" ;
    #endif
}
void debug_topo(const vector<pair<intg, intg>> &g_fpga,const intg &hop_limit){
    #ifdef LOG_PARSE
    std::cout << "=============================================\n" ;
        std::cout<<"hop_limit:"<<hop_limit<<'\n';
        for(auto fpga:g_fpga)
            cout<<fpga.first+1<<" "<<fpga.second+1<<"\n";
    std::cout << "=============================================\n" ;
    #endif
}
void debug_are(const std::map<intg,std::array<intg,8>> &node_capacity,const intg & total_node){
#ifdef LOG_PARSE
    std::cout << "=============================================\n" ;
    std::cout<<"node's number is "<<total_node<<'\n';
    for(auto node:node_capacity){
        std::cout<<"node"<<node.first+1<<" ";
        for(int i=0;i<8;i++){
            std::cout<<node.second[i]<<" ";
        }
        std::cout<<"\n";
    }
    std::cout << "=============================================\n" ;
#endif
}
}  // namespace topart
