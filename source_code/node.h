#ifndef NODE_H_
#define NODE_H_

#include <limits>
#include <iostream>
#include <algorithm>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <vector>
#include <utility>
#include <functional>
#include <cassert>
#include "tensor.h"
#include "type.h"


namespace topart {

using std::endl;
using std::set;
using std::string;
using std::stringstream;
using std::unordered_map;
using std::unordered_set;


class Net;
class CircuitNode;
class FPGANode;

//用于FPGA之间判断各自的关�?,�?��矩阵，在bd.cpp的build_hop_fpga�?��始化
extern unordered_map<FPGANode *, intg> fpga_to_index;
extern vector<FPGANode*> index_to_fpga;       
extern vector<vector<intg>> fpga_weight;      //合法的有距�?，不合法�?��-1,对自己的距�?�?0
extern vector<vector<intg>> fpga_parents;     //fpga[i]�?��i为起始点,所有点的父�?

class Node {

public:
    intg name;

public:
    Node(intg s) : name(s){};
    virtual string status() = 0;
};


class FPGANode : public Node {
public:
    std::array<intg,8> fpga_capacity;
    std::array<intg,8> usage;
    unordered_map<intg, set<FPGANode *>> S_hat;
    intg usage_wire;
    intg wire_capacity;
    unordered_set<CircuitNode*> placed_circuit;
    unordered_set<CircuitNode *> copy_circuit;
public:
    FPGANode(intg s, std::array<intg,8> &c,intg wire_number) : Node(s), fpga_capacity(c), wire_capacity(wire_number) {
        std::fill(usage.begin(), usage.end(), 0);
        usage_wire=0;
    }
    virtual string status() override {      //输出S hat
        stringstream ss;
        for (auto &p : S_hat) {
            auto &dist = p.first;
            auto &ff = p.second;
            ss << "hat{S}(F" << name << ", " << dist << ") = ";
            for (auto &s : ff) {
                ss << s->name << ", ";
            }
            ss << endl;
        }
        return ss.str();
    }
    void add_circuit(CircuitNode* Node);
    void remove_circuit(CircuitNode* Node);
    void add_copy_circuit(CircuitNode* c);
    void remove_copy_circuit(CircuitNode* c);
    void add_circuit_undercopy(CircuitNode *c);
    void remove_circuit_undercopy(CircuitNode *c);
    //intg valid(CircuitNode *c);
    //intg exchange_valid(CircuitNode *origin,CircuitNode *new_come,FPGANode *new_fpga);
    intg exchange_valid_test(CircuitNode *origin,CircuitNode *new_come);
    std::pair<double,double> free_space(CircuitNode *base);
    intg free_wire();
    void get_demand(CircuitNode* c,std::array<intg,8> &resource_demand,unordered_map<FPGANode*,intg>& wire_demand);
};

class CircuitNode : public Node {
public:
    FPGANode *fpga_node;

    set<FPGANode *> cddt;  // candidate fpga for this node to be assigned to.
    Tensor<intg> tsr_cddt;

    bool fixed;
    set<CircuitNode *> S;  // only for fixed node.

    unordered_set<Net *> nets;              //存储所参与的超�?
    unordered_set<Net *> src_nets;              //�?��作为src参与的超边是
    unordered_set<Net *> drain_nets;            //�?��作为drain参与的超�?
    unordered_map<intg, intg> cut_increment_map;
    unordered_map<intg, intg> cut_increment_undercopy_map;
    unordered_map<intg, intg> connect_increment_map;//这里�??外互联数�?/�?��的capacity
    unordered_map<intg, intg> copy_increment_map;//这里�??外互联数�?/�?��的capacity
    bool should_defer;
    intg neighbor_assign_cnt;
    intg neighbor_sum_cnt;
    std::array<intg, 8> node_capacity;
    unordered_set<FPGANode *> copy_position;

public:
    CircuitNode(intg s, intg num_f,std::array<intg,8> &capacity) : Node(s),node_capacity(capacity) {
        fixed = false;
        fpga_node = nullptr;
        tsr_cddt = Tensor<intg>(num_f);
        should_defer = false;
        cut_increment_map.clear();
        connect_increment_map.clear();
        copy_position.clear();
        copy_increment_map.clear();
        cut_increment_undercopy_map.clear();
        neighbor_assign_cnt=0;
        neighbor_sum_cnt=0;
    }
    void set_fixed(FPGANode *fn);
    void add_net(Net *n,bool is_src) {
        nets.emplace(n);   
        if(is_src){
            src_nets.emplace(n);
        }else{
            drain_nets.emplace(n);
        }
    }
    void add_fpga(FPGANode *fn);
    void remove_fpga(FPGANode *fn);
    void add_copy_fpga(FPGANode *fn);
    void remove_copy_fpga(FPGANode *fn);
    void reset_cut_increment() { cut_increment_map.clear(); }
    void reset_connect_increment() { connect_increment_map.clear(); }
    void reset_copy_increment() { copy_increment_map.clear(); }
    void calculate_cut_increment(FPGANode *if_f);
    void calculate_connect_increment(FPGANode *if_f);
    void calculate_copy_increment(FPGANode *if_f);
    void calculate_cut_increment_undercopy(FPGANode *if_f) ;
    bool is_fixed() { return fixed; }
    virtual string status() override {
        stringstream ss;
        ss << "S(" << name << ") = ";
        for (auto &p : S) {
            ss << p->name << ", ";
        }
        ss << endl;
        return ss.str();
    }

    void flush_tsr_to_cddt(vector<FPGANode *> &mapping);
    void flush_cddt_to_tsr();

    bool assigned() { return is_fixed() || fpga_node != nullptr; }
    void defer() { should_defer = true; }

    // intg try_move(Tensor<intg> &move_cddt);
    intg try_move(FPGANode *f);

};

class MergeNode : public Node {
public:
    FPGANode *fpga_node;
    unordered_set<CircuitNode *> belong_nodes;
    CircuitNode * behalf_c;
    unordered_set<Net *> nets;              //存储所参与的超�?
    unordered_set<Net *> src_nets;              //�?��作为src参与的超边是
    unordered_set<Net *> drain_nets;            //�?��作为drain参与的超�?
    unordered_set<Net *> common_nets;
    std::array<intg, 8> node_capacity;
    unordered_map<FPGANode *,unordered_set<CircuitNode*>> copy_position;
public:
    MergeNode(intg s,CircuitNode *behalf_c,unordered_set<CircuitNode *> &belong_nodes,FPGANode* f) :
     Node(s),behalf_c(behalf_c),belong_nodes(belong_nodes),fpga_node(f) {
        for(auto node:belong_nodes){
            for(int i=0;i<8;i++){
                node_capacity[i]+=node->node_capacity[i];
            }
            nets.insert(node->nets.begin(),node->nets.end());
            src_nets.insert(node->src_nets.begin(),node->src_nets.end());
            drain_nets.insert(node->drain_nets.begin(),node->drain_nets.end());
            for(auto copy_f:node->copy_position){
                copy_position[copy_f].insert(node);
            }
        }
        for(auto net_:src_nets){
            if(drain_nets.count(net_)){
                //belong_node��һ����src_node����һЩ��drain_node
                drain_nets.erase(net_);
                common_nets.insert(net_);
            }
        }
    
    }
    virtual string status() override {
        stringstream ss;
        ss << "MergeNode(" << name << ") = {";
        for (auto &node : belong_nodes) {
            ss << node->name << ", ";
        }
        ss << "} on FPGA " << (fpga_node ? fpga_node->name : -1) << endl;
        return ss.str();
    }
    void add_fpga_circuit(FPGANode *fn);
    void remove_fpga_circuit(FPGANode *fn);
    void add_copy_fpga(FPGANode *fn);
    void remove_copy_fpga(FPGANode *fn);
    void calculate_copy_increment(FPGANode *if_f);
    void calculate_cut_increment_undercopy(FPGANode *if_f) ;
};

class Net {
public:
    unordered_set<CircuitNode *> net_cell;      //drain_node and src
    unordered_set<FPGANode *> used_fpga_node;
    unordered_map<intg, intg> used_fpga_node_count;
    CircuitNode* source;           // source node
    intg weight;           // 超边权重 
    unordered_set<CircuitNode *> copy_circuit;
    //如果复制src到fpga_a，�?于当前fpga要考虑消除src对于drain_node(used_fpga_node和used_copy_fpga)产生的割
    //如果复制drain到fpga-a，�?考虑当前fpga在fpga_a上产生新的割
    unordered_set<FPGANode *> used_copy_fpga;
    unordered_map<intg, intg> used_copy_fpga_count;

    unordered_set<FPGANode *> used_all_fpga;
public:
    Net(unordered_set<CircuitNode *> nc,CircuitNode *src,intg wweight) : net_cell(nc),source(src),weight(wweight) {
    };
    void add_fpga(FPGANode *f) {
        used_fpga_node.emplace(f);
        used_all_fpga.emplace(f);
        if (used_fpga_node_count.count(f->name) <= 0) {
            used_fpga_node_count[f->name] = 1;
        } else {
            used_fpga_node_count[f->name]++;
        }
    }
    void remove_fpga(FPGANode *f) {
        if (used_fpga_node_count[f->name] == 1) {
            used_fpga_node_count.erase(f->name);
            used_fpga_node.erase(f);
            if(!used_copy_fpga.count(f))
                used_all_fpga.erase(f);
        } else {
            used_fpga_node_count[f->name] -= 1;
        }
    }
    void add_copy_fpga(FPGANode *f);
    void remove_copy_fpga(FPGANode *f);
    void add_copy_circuit(CircuitNode *c);
    void remove_copy_circuit(CircuitNode *c);  
    intg estimate_increase_cut_size(FPGANode *if_f);
    intg estimate_increase_cut_size_refine(FPGANode *if_f,FPGANode* origin,bool src);
    intg cost() {
        if (used_fpga_node.size() <= 1)
            return 0;
        return used_fpga_node.size();
    }
    bool try_move(FPGANode *f);

};


}  // namespace topart


#endif  // NODE_H_
