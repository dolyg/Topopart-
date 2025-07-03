#ifndef DB_H
#define DB_H

#include <fstream>
#include <queue>
#include <map>

#include "graph.h"
#include "node.h"
#include "timer.h"

namespace topart {

using std::fstream;
using std::queue;

class DB {
public:
    //find_root的时候是否除去fixed_node的CV
    bool button_CV_fixed_check=0;
    //connection_function的a
    //hp1,hp2,hp3,hp4

    Graph<FPGANode> fpga;
    Graph<FPGANode> hop_fpga;   
    Graph<CircuitNode> circuit;
    Graph<CircuitNode> dir_circuit;
    vector<Net *> nets;
    vector<pair<double,double>> fpga_neighbor_free_space;        
    vector<intg> fpga_neighbor_free_wire;
    unordered_set<CircuitNode *> no_indegree_nodes;
    intg queue_capacity=10;
    set<CircuitNode *> circuit_defer;
    set<CircuitNode *> fixed_defer;
    vector<pair<intg,intg>> fixed_nodes;
    queue<pair<intg, intg>> Q_fixed;//calculate_cddt_fpga使用
    std::array<double,8> ave_source;
    std::array<double,8> sum_source;
    std::array<double,8> rate_source;
    intg hop_after_parition=0;
    intg all_copy_cnt=0;
    vector<intg> bfs_order;
    Timer timer;
    class MonotonicQueue;

    // topart procedure - algorithm 1

    //优化后的cal_circuit_candidate
    void cal_circuit_candidate();
    intg estimate_cut_increment(intg node_id, intg to_fpga_id);
    bool enough_space_for_neighbor(CircuitNode *c, FPGANode *f);

    void calculate_fpga_neighbor_free_space(intg fpga_id,CircuitNode *base);
    intg output_loss();
    intg topo_vio(FPGANode *f, intg c);

    // void try_legalize(intg threshold);
    void debugnet();
    void get_bfs_order();
public:
    DB() { timer.set_start(); }

    // Build DB
    void build_fpga_graph(intg num_vertex,
                            std::map<intg,std::array<intg,8>> &fpga_capacity,
                            vector<pair<intg, intg>> &e,
                            std::map<intg,intg> &wire_capacity,
                            intg hop);
    void build_circuit_graph(intg num_vertex, vector<vector<intg>> &e,std::map<intg,std::array<intg,8>> &node_capacity,intg total_fpga);
    void set_fixed_circuit(vector<intg> &fixed_node,vector<intg> &fixed_fpga,const vector<CircuitNode*> &vv,intg fix_number);

    // topart procedure
    void calculate_candidate_fpga();
    void partition();
    void partition_after(set<intg> &Q_before);
    void refine();
    void output(fstream &out);
    void build_hop_fpga(intg hop,intg num_vertex);
    //寻根算法
    bool topo_valid(FPGANode *f,CircuitNode *c);
    bool topo_valid_undercopy(FPGANode *f,CircuitNode *c);
    bool copy_valid(FPGANode *f,CircuitNode *c);
    bool try_copy(CircuitNode *c,Tensor<intg> &cddt_copy);
    void find_root();
    void solve_CV();
    void add_source(CircuitNode *c);
    void remove_source(CircuitNode *c);
    void solve_src();
    intg try_move(CircuitNode* c,Tensor<intg> &cddt_move);
    intg try_move_undercopy(CircuitNode* c,Tensor<intg> &cddt_move);
    void debug_resouce();
    void get_fixed_CV(vector<intg> &fixed_node,vector<intg> &fixed_fpga,intg fix_number);
    void free_fixed_CV();
    intg try_exchange(CircuitNode *origin,CircuitNode *new_come);
    intg exchange_valid(CircuitNode *origin,FPGANode *origin_fpga,CircuitNode *new_come,FPGANode *new_fpga);
    intg valid(FPGANode *f,CircuitNode *c);
    intg valid_undercopy(FPGANode *f,CircuitNode *c);
    intg free_space_function(FPGANode *f,CircuitNode *c);
    double quantizi_src(CircuitNode *c);
    void refine_CV();
    void refine_space_wire();
    void force_CV(double argu_conn,vector<CircuitNode*> &re_q);
    void mixed_refine(int loop_time,int copy_button);
    void refine_copy();
    void refine_under_copy();
    intg try_move_CV_space_increment(CircuitNode *c,FPGANode *f);
    intg move_originf_connect_increment(CircuitNode *c);
    intg reject_place(CircuitNode *c);
    void force_place(intg c,intg f);
    void force_remove(intg c,intg f);
    void force_exchange(intg a,intg b);
    intg try_move_increment(CircuitNode *c,FPGANode *f);
    intg try_move_increment_undercopy(CircuitNode *c,FPGANode *f);
    intg try_move_copy_node_increment(CircuitNode *c,FPGANode *f);
    double get_score(int total_hop);
    double get_node_powerful_capacity(CircuitNode *c);
    void output();
    bool source_include(CircuitNode *removed_node,CircuitNode *placing_node);
    bool is_src_or_drain(CircuitNode *c,CircuitNode *verified);
    int remove_copy_valid(CircuitNode *c,FPGANode *f);
    void remove_some_node(vector<CircuitNode *> &re_q);
    //多级方法
    vector<MergeNode *> merge_nodes;
    Graph<MergeNode > merge_graph;
    void build_merge_graph();
    void merge_vertex(MergeNode *c,MergeNode *new_node);
    void multi_refine(int loop_time);
    intg try_multimove_undercopy(unordered_set<CircuitNode *> &union_c,Tensor<intg> &cddt_move);
    intg try_multimove_increment_undercopy(MergeNode &c, FPGANode *if_f);
    int multi_valid_undercopy(FPGANode *if_f,MergeNode &c);
    int topo_multi_valid_undercopy(FPGANode *if_f,MergeNode &c);
    int free_space_function(MergeNode &c,FPGANode *if_f);
};


}  // namespace topart

#endif /* DB_H */
