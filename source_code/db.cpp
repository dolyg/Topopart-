#include "db.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <numeric>
#include <set>

#include "log.h"
#include "node.h"
#include "tensor.h"
// #include "km.h"

namespace topart {

using std::cout;
using std::endl;
using std::fstream;
using std::make_pair;
using std::min;
using std::pair;
using std::queue;
using std::map;

template<class T,class I,class D>
void SPFA(
        const unordered_map<T *, I> &vv_map,
        const vector<T*> &vv,
        std::vector<D> &dist,
        std::vector<intg> &parents,
        const std::vector<std::vector<T*>> &gg,
        const std::vector<std::vector<D>> &wweight,
        intg root_index);

template<class T,class I>
void print_G(
            unordered_map<T *, intg> &v_map,
            vector<T *> &vv,
            vector<vector<T *>> &gg,
            vector<vector<I>> &wweight,
            std::string s);

void DB::build_fpga_graph(intg num_vertex,
                          std::map<intg,std::array<intg,8>> &fpga_capacity,
                          vector<pair<intg, intg>> &e,
                          std::map<intg,intg> &wire_capacity,
                          intg hop) {
    fpga_neighbor_free_space.resize(num_vertex,make_pair(0,0));
    fpga_neighbor_free_wire.resize(num_vertex,0);
    vector<vector<FPGANode *>> gg(num_vertex);
    vector<vector<intg>> e_weight(num_vertex);
    vector<FPGANode *> fpga_node;
    for (int i = 0; i < num_vertex; ++i) {
        fpga_node.emplace_back(new FPGANode(i, fpga_capacity[i],wire_capacity[i]));
    }

    for (auto &p : e) {
        auto &v = p.first;
        auto &w = p.second;
        gg[v].emplace_back(fpga_node[w]);
       e_weight[v].emplace_back(1);
        gg[w].emplace_back(fpga_node[v]);
        e_weight[w].emplace_back(1);
    }
    fpga.init(fpga_node, gg,e_weight);

    //print_G(fpga.v_map,fpga.v,fpga.g,fpga.weight,"fpga");
    build_hop_fpga(hop,num_vertex);

    log("HOP FPGA is built");
}
void DB::build_hop_fpga(intg hop,intg num_vertex){
    //初始化node_FPGA里的矩阵。
    index_to_fpga=fpga.v;
    fpga_to_index=fpga.v_map;
    fpga_weight.resize(num_vertex,vector<intg>(num_vertex,-1));
    //hop_graph是由graph变化的，对每个点跑一次最短路，把最短路小于等于hop的点连接起来，边权代表的是两点之间的最短距离

    vector<vector<FPGANode *>> gg(num_vertex);
    vector<vector<intg>> e_weight(num_vertex);
    // vector<vector<intg>> parents;
    for(int i=0;i<num_vertex;++i){
        vector<intg> dist(num_vertex, std::numeric_limits<intg>::max());       //鍒濓拷?鍖栨墍鏈夎窛绂讳负INT64_MAX
        vector<intg> parent(num_vertex, -1);
        SPFA(fpga.v_map, fpga.v, dist,parent, fpga.g, fpga.weight, i);     //鍦╢pga涓婅窇鏈€锟�?锟斤拷
        for(int j=0;j<num_vertex;++j){
            // if(j==i)    continue;
            if(j==i){
                fpga_weight[i][j]=dist[j];
                continue;   
            }
            if(dist[j]<=hop){
                gg[i].emplace_back(fpga.v[j]);
                e_weight[i].emplace_back(dist[j]);
                fpga_weight[i][j]=dist[j];
            }
        }
        fpga_parents.emplace_back(parent);
    }
    hop_fpga.init(fpga.v,gg,e_weight);

    //更新node.h里的
    
    //print_G(hop_fpga.v_map,hop_fpga.v,hop_fpga.g,hop_fpga.weight,"hop_fpga");
}

void DB::build_circuit_graph(intg num_vertex, vector<vector<intg>> &e,std::map<intg,std::array<intg,8>> &node_capacity,intg total_fpga) {
    vector<CircuitNode *> circuit_node;                 //椤剁偣鍚戦噺

    vector<vector<CircuitNode *>> gg(num_vertex);       //鏃犲悜鍥鹃偦鎺ヨ〃
    vector<vector<intg>> e_weight(num_vertex);
    vector<intg> indegree(num_vertex,0);            //璁板綍姣忎釜鑺傜偣鐨勫叆锟�?
    vector<vector<CircuitNode *>> dir_gg(num_vertex);       //鏈夊悜鍥鹃偦鎺ヨ〃
    vector<vector<intg>> dir_e_weight(num_vertex);

    for (int i = 0; i < num_vertex; ++i) {              //鍒涘缓Circuit鑺傜偣鐭╅樀 
        circuit_node.emplace_back(new CircuitNode(i, total_fpga ,node_capacity[i]));
    }
    // Build Net and graph
    nets.clear();                                       //瀛樺偍鎵€鏈夌殑瓒呰竟
    unordered_set<CircuitNode *> us;                    //瓒呰竟,琛ㄧず涓轰竴锟�?锟斤拷锟�?

    for (auto &n : e) {
        intg start_v = n[0];
        intg weight = n[n.size()-1];
        unordered_set<CircuitNode *> us;
        us.clear();
        us.emplace(circuit_node[start_v]);
        for (int i = 1; i < n.size()-1; ++i) {
            auto &neighbor = n[i];
            us.emplace(circuit_node[neighbor]);
            gg[start_v].emplace_back(circuit_node[neighbor]);
            e_weight[start_v].emplace_back(weight);
            gg[neighbor].emplace_back(circuit_node[start_v]);
            e_weight[neighbor].emplace_back(weight);
            //有向图
            dir_gg[start_v].emplace_back(circuit_node[neighbor]);
            dir_e_weight[start_v].emplace_back(weight);
        }

        //每个顶点存储对应的超边 
        nets.push_back(new Net(us,circuit_node[start_v],weight));
        circuit_node[start_v]->add_net(nets.back(),1);
        for (int i = 1; i < n.size()-1; ++i) {
            auto &neighbor = n[i];
            circuit_node[neighbor]->add_net(nets.back(),0);
            indegree[neighbor]++;
        }
    }
    //处理度数
    for(int i=0;i<num_vertex;++i)
        if(indegree[i]==0)
            no_indegree_nodes.emplace(circuit_node[i]);

    circuit.init(circuit_node, gg, e_weight);
    dir_circuit.init(circuit_node, dir_gg, dir_e_weight);

    //debug
    //print_G(circuit.v_map,circuit_node,gg,e_weight,"circuit");
    //print_G(dir_circuit.v_map,circuit_node,dir_gg,dir_e_weight,"dir_circuit");
}


void DB::find_root(){
    //根据circuit创建图G'
    intg fix_number=hop_fpga.num_vertex;
    vector<vector<CircuitNode *>> &gg(circuit.g);

    intg root_index=0;          //Initial root node with minimal degree node
    vector<intg> fixed_node;        //store all the fixed node
    vector<intg> fixed_fpga(fpga.num_vertex,0);         //store all the fixed fpga node
    std::iota(fixed_fpga.begin(), fixed_fpga.end(), 0);
    for(intg i=0;i<gg.size();i++){           //i is the node`s index ,g[i][j]是nodeb的指针，需要找到i和j同时存在的超边
        if(gg[i].size()<gg[root_index].size())
            root_index=i;
    }

    // //从root_index开始对vv,gg,wweight这个图进行SPFA
    fixed_node.emplace_back(root_index);
    int all_root_number=min(fix_number,circuit.num_vertex);       //还需要寻找root数量
    // std::vector<double> dist(gg.size(), std::numeric_limits<intg>::max());       //all dist be -1
    std::vector<intg> dist(gg.size(), std::numeric_limits<intg>::max());      
    std::vector<intg> parents(gg.size(), -1);
    while(all_root_number){
        // SPFA(circuit.v_map, vv, dist, parents,gg, wweight, root_index);
        SPFA(circuit.v_map, circuit.v, dist, parents,circuit.g, circuit.weight, root_index);
        root_index=max_element(dist.begin(),dist.end())-dist.begin();   //find the max distance node
        fixed_node.emplace_back(root_index);
        dist[root_index]=0;
        all_root_number--;
    }

    //solve the priority
    auto sum_degree=[this](intg a){
            double sum_weight_a=0;
            for(int i=0;i<circuit.weight[a].size();i++)
                sum_weight_a+=circuit.weight[a][i];
            sum_weight_a=circuit.g[a].size()/sum_weight_a;
            return sum_weight_a;
        };

    //第一优先级是SPFA上的权重，第二优先级才是sum_degree,因此只排列前半段，不排列后半段
    auto Priority_end=fixed_node.begin()+fix_number;
    //    //对fixed_node里的节点按照权重之和乘以度数从大到小排序
    sort(fixed_node.begin(),Priority_end,[&](intg a,intg b){
        //对于两个固定节点，比较权重和/度数。
        return sum_degree(a)<sum_degree(b);
    });
    //对于fpga按照度数排序，fpga度数大的放前面,跟重要的fixed_node对应
    sort(fixed_fpga.begin(),fixed_fpga.end(),[this](intg a,intg b){
        return hop_fpga.g[a].size()>hop_fpga.g[b].size();
    });

    //处理违例fixed_fpga;
    int fixed_node_size=fixed_node.size();
    int CV_num=0;bool force_place=0;
    if(button_CV_fixed_check){
        for(int i=0;i<fix_number;i++){  
            for(int j=0;j<i;j++){
                intg node_i=fixed_node[i],node_j=fixed_node[j];
                if(circuit.g_set[node_i].count(circuit.v[node_j])){       //濡傛灉i鍜宩锟�?锟斤拷锟�?
                    intg fpga_i=fixed_fpga[i],fpga_j=fixed_fpga[j];
                    if(fpga_weight[fpga_i][fpga_j]==-1){                 //瀵瑰簲鐨刦pga杩樿繚渚嬩簡
                        int temp=fixed_node[i];
                        for(int k=i;k<fixed_node.size()-1;++k)          //閭ｅ氨鎶奿闄ゅ幓
                            fixed_node[k]=fixed_node[k+1];
                        fixed_node[fixed_node_size-1]=temp;
                        //对新的候选节点重新排序
                        for(int i=fix_number-1;i>0;i--)
                            if(sum_degree(fixed_node[i])>sum_degree(fixed_node[i-1]))
                                std::swap(fixed_node[i],fixed_node[i-1]);
                            else
                                break;
                        //重新匹配fixed_node
                        i--;
                        fixed_node_size--;
                        CV_num++;
                        if(CV_num<fpga.num_vertex)
                            break;
                        else{
                            #ifdef LOG
                            cout<<"Can't guarantee no CV fixed_Node\n";
                            #endif
                            force_place=1;
                        }
                    }
                }
            }// for(int j=0;j<i;j++)
            if(force_place){
                get_fixed_CV(fixed_node,fixed_fpga,fix_number);
                break;
            }
        }//for(int i=0;i<fixed_fpga.size();i++)
    }//if(button_CV_fixed_check)
    else{
        get_fixed_CV(fixed_node,fixed_fpga,fix_number);
    }

    if(fix_number>fixed_node_size){
        #ifdef LOG
        cout<<"The fixed_node is not equal to fixed_fpga\n";
        #endif
        return ;
    }
    else{    
        //设置固定节点。
        set_fixed_circuit(fixed_node,fixed_fpga,circuit.v,fix_number);
    }
}
void DB::get_fixed_CV(vector<intg> &fixed_node,vector<intg> &fixed_fpga,intg fix_number){
    //鎶奀V鐨刦ixed_node娣诲姞鍒癴ixed_defer
    int fixed_node_size=fixed_node.size();
    int CV_num=0;bool force_place=0;
    for(int i=0;i<fix_number;i++){   
        for(int j=0;j<i;j++){
            intg node_i=fixed_node[i],node_j=fixed_node[j];
            if(circuit.g_set[node_i].count(circuit.v[node_j])){     
                intg fpga_i=fixed_fpga[i],fpga_j=fixed_fpga[j];
                if(fpga_weight[fpga_i][fpga_j]==-1)                 
                    fixed_defer.emplace(circuit.v[node_i]);
            }
        }
    }
}
void DB::set_fixed_circuit(vector<intg> &fixed_node,vector<intg> &fixed_fpga,const vector<CircuitNode*> &vv,intg fix_number) {
    const auto &circuit_vertex = circuit.get_all_vertex();
    const auto &fpga_vertex = fpga.get_all_vertex();
    #ifdef LOG
    cout<<"fixed node:"<<endl;
    #endif
    for(int i=0;i<fix_number;i++){
        fpga_vertex[fixed_fpga[i]]->add_circuit(vv[fixed_node[i]]);
        circuit_vertex[fixed_node[i]]->set_fixed(fpga_vertex[fixed_fpga[i]]);
        fixed_nodes.push_back(make_pair(fixed_node[i],fixed_fpga[i]));
        Q_fixed.push(make_pair(fixed_node[i],fixed_fpga[i]));
        #ifdef LOG
        cout<<fixed_node[i]+1<<" ";
        if((i-1)%10==0&&i!=0)
            cout<<endl;
        #endif
    }
    #ifdef LOG
    cout<<endl;
    #endif
}

void DB::calculate_candidate_fpga() {
    // calculate fpga maxDist and S-hat
    /*
        求FPGA的S_hat
        处理cddt
    */

    const auto &fpga_list = hop_fpga.get_all_vertex();
    auto insert_S_hat = [&](intg src, intg dst, intg dist) {
        for (int i = 1; i < dist; ++i) {          
            for (auto &s : fpga_list[src]->S_hat[i]) {
                fpga_list[src]->S_hat[dist].emplace(s);
            }
        }
        fpga_list[src]->S_hat[dist].emplace(fpga_list[dst]);    
        fpga_list[src]->S_hat[dist].emplace(fpga_list[src]);
    };
    hop_fpga.calculate_max_dist(insert_S_hat);
    //hop_fpga.get_status();
    log("FINISH fpga pre-calculation");

    for (auto &c_node : circuit.get_all_vertex()) {
        if (c_node->is_fixed()) {
            c_node->cddt.emplace(c_node->fpga_node);
        } else {
            for (auto &f : fpga_list) {
                c_node->cddt.emplace(f);
            }
        }
    }
    log("FINISH circuit node cddt pre-calculation");
    //circuit.get_status();

    cal_circuit_candidate();
}

void DB::cal_circuit_candidate(){
    vector<FPGANode *> intersection_result;     //store cddt and S_hat
    unordered_set<intg> propagation_node;       //The propagation node
    unordered_map<intg, unordered_map<intg, intg>> circuit_node_s_dist;     
    const auto &circuit_list = circuit.get_all_vertex();
    while (!Q_fixed.empty()) {
        //我居然对一个变化的容器使用了&去引用，太傻逼了
        auto top = Q_fixed.front();
        Q_fixed.pop();

        auto node_id = top.first;
        auto fpga_id = top.second;
        auto c_node = circuit.get_vertex(node_id);
        auto fpga_node = fpga.get_vertex(fpga_id);
        //优化求S
        auto d_max = hop_fpga.max_dist[fpga_id];
        intg c_node_maxd=circuit.opt_get_max_dist(c_node->name,
                                [&](intg src, intg dst, intg dist) {     
                                    circuit_node_s_dist[src][dst] = dist;
                                    if (dist < d_max) {
                                        c_node->S.emplace(circuit_list[dst]);
                                    }
                                    c_node->S.emplace(circuit_list[src]);
                                },propagation_node);
        //已经传播过的节点记录下来
        propagation_node.emplace(node_id);
        //枚举moveable节点
        for (auto &neighbor : c_node->S) {     
            if (neighbor->is_fixed())
                continue;

            // propagation_node.emplace(neighbor->name);
            //求取fixed_node到moveable_node的距离,小于等于这个距离的FPGA才能被更新候选使用thory 1
            intg k = -1;
            if (circuit_node_s_dist.count(node_id) > 0 &&
                circuit_node_s_dist[node_id].count(neighbor->name) > 0) {
                k = circuit_node_s_dist[node_id][neighbor->name];
            } 
            else {
                throw std::logic_error("Didn't calculate the dist");
            }

            // 求交集尺寸 Cddt_fpga(moveable_node) and S(fixed_node, k)
            intg desired_size =
                std::min(neighbor->cddt.size(), fpga_node->S_hat[k].size());
            if (intersection_result.size() < desired_size)
                intersection_result.resize(desired_size);

            auto it = std::set_intersection(      
                neighbor->cddt.begin(), neighbor->cddt.end(),
                fpga_node->S_hat[k].begin(), fpga_node->S_hat[k].end(),
                intersection_result.begin());

            //求出交集更新Cddt_fpga(moveable_node)
            neighbor->cddt = set<FPGANode *>(intersection_result.begin(), it);

            if (neighbor->cddt.size() == 1) {
                auto vj_hat = *(neighbor->cddt.begin());
                Q_fixed.push(make_pair(neighbor->name, vj_hat->name));
            }

            if (neighbor->cddt.size() == 0) {
                throw std::logic_error("Unreachable");
            }
        }
    }

    //create cddt
    for (auto &c : circuit.get_all_vertex()) {
        c->flush_cddt_to_tsr();
    }

#ifdef LOG_DB
    log("Cddt result: ");
    for (auto &c : circuit.get_all_vertex()) {
        cout << "Cddt(" << c->name << ") = ";
        for (auto &s : c->cddt) {
            cout << s->name << ", ";
        }
        cout << endl;
    }
#endif
}

void DB::calculate_fpga_neighbor_free_space(intg fpga_id,CircuitNode *base) {     
                                                                        //仅用于FPGA排序
    intg fpga_num=1;
    auto &f_usage = fpga_neighbor_free_space[fpga_id];
    auto &f_usage_wire=fpga_neighbor_free_wire[fpga_id];
    const auto &f = fpga.get_vertex(fpga_id);
    auto space_f=f->free_space(base);                   //计算f自身的剩余资源与base所需资源的关系
    auto wire_f=f->free_wire();               //计算base后，有关所有fpga的剩余的线的情况{最小线剩余，平均线剩余 }
    f_usage.first=space_f.first;  f_usage.second=space_f.second;
    f_usage_wire=wire_f;
    for (auto &f_neighbor : hop_fpga.g[f->name]) {
        auto tmp=f_neighbor->free_space(base);
        f_usage.first =f_usage.first + tmp.first;
        f_usage.second =f_usage.second + tmp.second;
        fpga_num++;
    }
}
intg DB::try_move(CircuitNode* c,Tensor<intg> &cddt_move){
    //评价是否存在一个fpga，可以跟c的邻居的fpga不产生违例
    for(auto neibor:circuit.g[c->name]){
        if(neibor->fpga_node==nullptr) continue;
        Tensor<intg> tmp_(hop_fpga.num_vertex,1);
        for(auto fpga_neibor:hop_fpga.g[neibor->fpga_node->name])
            tmp_.at_is(fpga_neibor->name,0);
        tmp_.at_is(neibor->fpga_node->name,0);
        cddt_move-=tmp_;

    }
    return cddt_move.all_zero();
}
void DB::free_fixed_CV(){
    for(auto CV_node:fixed_defer){
        auto remove_fpga=CV_node->fpga_node;
        remove_fpga->remove_circuit(CV_node);
        CV_node->remove_fpga(remove_fpga);
        circuit_defer.emplace(CV_node);
        #ifdef LOG
        cout<<"CV fixed node is "<<CV_node->name+1<<endl;
        #endif
    }
}
bool DB::topo_valid(FPGANode *f,CircuitNode *c){
    for(auto neibor:circuit.g[c->name]){
        if(neibor->fpga_node==nullptr)continue;
        if(fpga_weight[f->name][neibor->fpga_node->name]==-1)
            return 0;
    }
    return 1;
}
intg DB::valid(FPGANode *f,CircuitNode *c){
    //c没有fpga
    //topo_valid 是-2是最严重的问题
    if(topo_valid(f,c)==0){
        return -2;
    }
    //return (usage < fpga_capacity)
    #ifdef LOG
    assert(c->fpga_node==nullptr);
    assert(f!=nullptr);
    #endif
    for(int i=0;i<8;i++)
        if(f->usage[i]+c->node_capacity[i]>f->fpga_capacity[i]){
            return 0;
        }
    vector<intg> used_wire(hop_fpga.num_vertex,0);
    for(auto &net:c->src_nets){   
        //对作为src的点来说，一旦放置后会产生割，那所有fpga都要检验   

        //net里面还没有fpga，或者只有自身那肯定不产生割
        if(net->used_fpga_node.size()==0||net->used_fpga_node.size()==1&&net->used_fpga_node.count(f))
            continue;

        //剩下的情况下都会产生割
        if(f->usage_wire+net->weight+used_wire[f->name]>f->wire_capacity)
            return -1;
        else
            used_wire[f->name]+=net->weight;
        for(auto &fpga_:net->used_fpga_node){
            if(fpga_==f) continue;
            if(fpga_->usage_wire+net->weight+used_wire[fpga_->name]>fpga_->wire_capacity)
                return -1;
            else
                used_wire[fpga_->name]+=net->weight;
        }
    }
    for(auto &net:c->drain_nets){ 
        auto source_c=net->source;
        auto source_f=source_c->fpga_node;
                //对于drain点来说
        //如果source未放置那无所谓
        //如果src已经放置，当前点放置后产生了第二个fpga，那要给自身和src都添加割
        //如果src已经放置，当前点放置后没有产生新的fpga，那就不影响
        //如果src已经放置，当前点放置后产生了新的fpga那就只处理自己
        if(net->used_fpga_node.count(f))
            continue;
        //如果src_f不存在，那就要保证drain_f的放置可以保证有一个fpga与所有drain_f相邻
        if(source_f==nullptr){
            continue;
        }
        //如果src_f不存在，那就要保证drain_f的放置可以保证有一个fpga与所有drain_f相邻
        if(!net->used_fpga_node.count(f)){
            if(f->usage_wire+net->weight+used_wire[f->name]>f->wire_capacity)
                return -1;
            else
                used_wire[f->name]+=net->weight;
            //如果只有src_fpga一个fpga，那src_fpga也处理
            if(net->used_fpga_node.size()==1)
                if(source_f->usage_wire+net->weight+used_wire[source_f->name]>source_f->wire_capacity)
                    return -1;
                else    
                    used_wire[source_f->name]+=net->weight;
        }
    }
    return 1;
}

void DB::solve_src(){
    //计算每个fpga分布的平均资源
    for(int i=0;i<8;i++)
        ave_source[i]=0;
    for(auto node_:circuit.get_all_vertex())
        for(int i=0;i<8;i++)
            ave_source[i]+=node_->node_capacity[i];
    for(int i=0;i<8;i++){
        sum_source[i]=ave_source[i];
        // ave_source[i]/=circuit.num_vertex;
        ave_source[i]=(double)ave_source[i]/hop_fpga.num_vertex;
        rate_source[i]=(double)ave_source[i]/fpga.get_vertex(0)->fpga_capacity[i];
    }
}

void DB::partition() {  //浣跨敤hop_FPGA
    // NOTE: for algorithm 2, Q with size of cddt of node id, node id
    // It need to be sorted.
    //Priority queue Q in order of the number of candidate FPGAs
    
    intg fpga_capacity_sum=0;
    double powerful_capacity=0;
    for(int i=0;i<8;++i)
        fpga_capacity_sum+=fpga.get_vertex(0)->fpga_capacity[i];
    for(int i=0;i<8;++i)
        powerful_capacity+=(double)fpga.get_vertex(0)->fpga_capacity[i]*(double)fpga.get_vertex(0)->fpga_capacity[i]/fpga_capacity_sum;

    int capacity_digit=log10((intg)powerful_capacity) + 1;
    int wire_digit=log10(fpga.get_vertex(0)->wire_capacity)+1;
    int circuit_digit=log10(circuit.num_vertex)+1;
    int net_digit=log10(nets.size())+1;
    int node_threshold=6;
    // NOTE: for algorithm 2, R with element pair:
    double hp_1,hp_3,hp_4,hp_5,hp_2;
    //0.1 2 10 0.1
    //0.00001 0.0002 0.00004
    hp_1=pow(10,capacity_digit-wire_digit)*0.1;       //对外互联数目增长数目
    hp_2=pow(10,capacity_digit-wire_digit)*2;        //hop增长
    hp_3=10*1;                      //circuit 的资源对于FPGA剩余资源的比例
    hp_4=pow(10,capacity_digit)*0.00001;        //剩余的线空间决定了分散程度
    //瓒筹拷?澶ч偅灏卞彲浠ラ€傚綋鏇磋拷姹俬op
    if(circuit_digit>=node_threshold){
        hp_2*=100000;
    }

    if(0){
        hp_1=pow(10,capacity_digit-wire_digit)*0;       
        hp_2=pow(10,capacity_digit-wire_digit)*100;        
        hp_3=10*0;     
        hp_4=pow(10,capacity_digit)*0;      \
    }
    #ifdef LOG
    cout<<hp_1<<" "<<hp_2<<" "<<hp_4<<endl;
    #endif
    auto Rcmp = [&](const pair<intg, intg> &lhs, const pair<intg, intg> &rhs) {
        //优先hop较小的
        const auto &cl = circuit.get_vertex(lhs.first);
        const auto &cr = circuit.get_vertex(rhs.first);
        const auto &fl = hop_fpga.get_vertex(lhs.second);
        const auto &fr = hop_fpga.get_vertex(rhs.second);
        assert(fl!=fr);
        //对外互联数目增长(可扩展性)越小越好
        double fl_connect_incre= cl->connect_increment_map[fl->name];
        double fr_connect_incre= cr->connect_increment_map[fr->name];
        //hop增长(可扩展性)越小越好
        double l_cost = cl->cut_increment_map[fl->name];
        double r_cost = cr->cut_increment_map[fr->name];
        //空间最小值和空间总值的倍数(合法性，可扩展性)越大越好
        pair<double,double> fl_usage_space = fpga_neighbor_free_space[fl->name];
        pair<double,double> fr_usage_space = fpga_neighbor_free_space[fr->name];
        //剩余的对外互联数目(合法性)越大越好
        double fl_usage_wire = fpga_neighbor_free_wire[fl->name];
        double fr_usage_wire= fpga_neighbor_free_wire[fr->name];

        //按照重量级排序
        //整体越大越好
        double l_ans=hp_1*(-fl_connect_incre)-hp_2*l_cost-hp_3*fl_usage_space.first+
                    hp_4*fl_usage_wire;
        double r_ans=hp_1*(-fr_connect_incre)-hp_2*r_cost-hp_3*fr_usage_space.first+
                    hp_4*fr_usage_wire;
        if(l_ans==r_ans)
            return lhs.second < rhs.second;
        return l_ans>r_ans;
    };
    double qp_1=(pow(10,capacity_digit)*0.001);
    double qp_2=(pow(10,net_digit-capacity_digit-circuit_digit)/fpga.num_vertex*100);
    auto Qcmp = [&](const intg &lhs, const intg &rhs) {
        //周围邻居少的先放置
        //权重大的先放置
        //资源少放置
        //候选的少放置 
        auto lc=circuit.get_vertex(lhs);
        auto rc=circuit.get_vertex(rhs);
            // return lhs<rhs;
        intg lg = circuit.g_set[lhs].size();
        intg rg = circuit.g_set[rhs].size();

        double lcs =(double)circuit.get_vertex(lhs)->cddt.size()*qp_1 + qp_2*lg;
        double rcs =(double)circuit.get_vertex(rhs)->cddt.size()*qp_1 + qp_2*rg;

        if (lcs == rcs) {
            return lhs<rhs;
        }

        return lcs < rcs;
    };


/*
    开始初始划分
*/
    using RType = std::set<pair<intg, intg>, decltype(Rcmp)>;
    using QType = std::set<intg, decltype(Qcmp)>;

    intg fixed_node_cnt=0;
    //不按照BFS去处理比较好
    QType Q(Qcmp);
    // set<intg> Q;
    vector<RType> R(circuit.num_vertex, RType(Rcmp));
    for (auto &c : circuit.get_all_vertex()) {
        if (c->is_fixed()){
            fixed_node_cnt++;
            continue;
        }
        c->fpga_node = nullptr;
        Q.emplace(c->name);
    }

    Tensor<intg> v_cddt(hop_fpga.num_vertex);
    // int loop=1;while (loop--) {
    while (!Q.empty()) {
        auto node_vj = *(Q.begin());
        Q.erase(Q.begin());
        auto c = circuit.get_vertex(node_vj);
        //是否已经确定了FPGA,已经确定就离开
        if (c->assigned())      
            continue;

        //没有候选节点，这个节点标记上defer后续重新处理
        if (c->cddt.size() == 0) {    
            circuit_defer.emplace(c);  
            c->defer();
            continue;
        }
        //计算node的R,计算对于node来说每个候选FPGA放置后的cut-size增加
        //计算吗每个候选fpga中fpga及其邻居的资源，为R中FPGA的先后顺序 
        R[c->name].clear();
        c->reset_cut_increment();
        c->reset_connect_increment();
        for (auto &ff : c->cddt) {
            calculate_fpga_neighbor_free_space(ff->name,c);
            c->calculate_cut_increment(ff); //鏀规垚璁＄畻hop澶у皬锟�?
            c->calculate_connect_increment(ff);
            R[c->name].emplace(make_pair(c->name, ff->name));
        }
        // assert(node_vj == c->name);
        bool traceback = false;
        FPGANode *f = nullptr;      
        intg fpga_vj;
        do {
            if (R[node_vj].size() <= 0) {
                circuit_defer.emplace(c); 
                c->defer();
                goto trace_end;
            }
            fpga_vj = (*(R[node_vj].begin())).second;
            f = hop_fpga.get_vertex(fpga_vj);
            R[node_vj].erase(R[node_vj].begin());
        } while (valid(f,c)<=0);   

        //给nodec的邻居不能放置的FPGA(非邻居的FPGA)标上1
        v_cddt.clear();
        // the fpga node, which didn't directly connect to f
        for (auto &indirect_f : hop_fpga.get_all_vertex()) {
            //如果fpga不是当前fpga的邻居那就对应的位置上计数
            if (hop_fpga.g_set[f->name].count(indirect_f) <= 0)
                v_cddt.at(indirect_f->name) += 1;
        }
        v_cddt.at(fpga_vj) = 0;

        //更新未处理的不需要refer的邻居的Cddt
        for (auto &neighbor : circuit.g[c->name]) {
            
            if (neighbor->assigned() || neighbor->should_defer)
                continue;

            neighbor->tsr_cddt -= v_cddt;

            if (neighbor->tsr_cddt.all_zero()) {
                traceback = true;
            }
        }
        //trackback说明当前放置不合适，需要撤回当前的放置
        if (traceback) {
            c->cddt.erase(f);
            //鐢╟ddt鏇存柊tsr
            c->flush_cddt_to_tsr();
            c->fpga_node = nullptr;
            Q.emplace(c->name);
            //鎭拷?閭诲眳鐨刢ddt
            for (auto &neighbor : circuit.g[c->name]) {
                if (neighbor->assigned())
                    continue;

                neighbor->tsr_cddt += v_cddt;
            }
            R[c->name].clear();
        } else {
            f->add_circuit(c);
            c->add_fpga(f);
            hop_after_parition+=c->cut_increment_map[f->name];
            for (auto &neighbor : circuit.g[c->name]) {
                if (neighbor->assigned() || neighbor->should_defer)
                    continue;
                Q.erase(neighbor->name);
                neighbor->flush_tsr_to_cddt(hop_fpga.v);
                Q.emplace(neighbor->name);
            }
        }
    trace_end:
        if (traceback)
            log("Traceback");
        else
            log("Don't traceback");
    }
    // set<intg> Q_before(Q.begin(),Q.end());
    // partition_after(Q_before);
    //释放固定的CV
    free_fixed_CV();
    #ifdef LOG
    cout<<"fixed_node_cnt:"<<fixed_node_cnt<<endl;
    cout<<"CV node's num:"<<circuit_defer.size()<<endl;
    if(circuit_defer.size()==0){
        debug_resouce();
        cout<<"cur hop:"<<hop_after_parition<<endl;
    }
    #endif
}

intg DB::try_move_CV_space_increment(CircuitNode *c,FPGANode *f){
    //求移动后的收益
    double dec = 0;
    double max_capa=0;
    for(int i=0;i<8;i++){
        if(c->node_capacity[i]==0||f->fpga_capacity[i]==0)  continue;
        max_capa=std::max(max_capa,(double)f->fpga_capacity[i]);
        dec+=abs((double)f->usage[i]-ave_source[i])/(double)f->fpga_capacity[i];   
        dec-=abs((double)f->usage[i]+c->node_capacity[i]-ave_source[i])/(double)f->fpga_capacity[i];   
    }
    return max_capa*dec;
}
intg DB::move_originf_connect_increment(CircuitNode *c){
    //输出的结果是减少的hop
    intg dec=0;
    auto origin_fpga=c->fpga_node;
    for(auto net:c->src_nets){
        if(net->used_fpga_node.size()<=1) continue;
        int temp_dec=0;
        for(auto fpga_:net->used_fpga_node)
            temp_dec+=fpga_weight[fpga_->name][origin_fpga->name];
        dec+=temp_dec*net->weight;
    }
    for(auto net:c->drain_nets){
        auto src_fpga=net->source->fpga_node;
        if(src_fpga==nullptr||origin_fpga==src_fpga) continue;
        if(net->used_fpga_node_count[origin_fpga->name]==1)
            dec+=net->weight*fpga_weight[src_fpga->name][origin_fpga->name];
    }
    return dec;
}
void DB::force_CV(double argu_conn,vector<CircuitNode*> &re_q){
        timer.output_time("start force_CV");
    vector<CircuitNode *> fpga_req[hop_fpga.num_vertex];
    vector<int> fpga_req_index(hop_fpga.num_vertex,0);
    intg force_repair_cnt=0;
    auto temp_defer_set=circuit_defer;
    for(auto &defer_node:temp_defer_set){
        if(defer_node->fpga_node!=nullptr)continue;
        Tensor<intg> cddt_move(hop_fpga.num_vertex,1);
        if (try_move(defer_node,cddt_move))
            continue;
        //给要去的FPGA排序
        auto place_fpga_cmp=[&](FPGANode * lhs,FPGANode* rhs){
            //寻找剩余空间多的
            return free_space_function(lhs,defer_node)>free_space_function(rhs,defer_node);
        };
        set<FPGANode *,decltype(place_fpga_cmp)> cddt_fpga_set(place_fpga_cmp); 
        for(int i=0;i<hop_fpga.num_vertex;++i){
            if(cddt_move.at(i)==1){
                cddt_fpga_set.emplace(hop_fpga.v[i]);
            }
        }
        for(auto cddt_fpga:cddt_fpga_set){
            int i=cddt_fpga->name;
            if(cddt_move.at(i)==1){
                int need_copy=valid(hop_fpga.v[i],defer_node);
                if(need_copy>0){
                    defer_node->calculate_cut_increment(hop_fpga.v[i]);
                    hop_after_parition+=defer_node->cut_increment_map[i];

                    hop_fpga.v[i]->add_circuit(defer_node);
                    defer_node->add_fpga(hop_fpga.v[i]);
                    circuit_defer.erase(defer_node);
                    force_repair_cnt++;
                    break;
                }
                else{
                    auto &cddt_move_nodes=hop_fpga.v[i]->placed_circuit;
                    //如果已经有排序的序列那就使用，否则不使用
                    if(fpga_req[i].size()==0){
                        fpga_req[i].assign(cddt_move_nodes.begin(),cddt_move_nodes.end());
                        sort(fpga_req[i].begin(), fpga_req[i].end(),
                            [&](CircuitNode *lhs, CircuitNode *rhs) {
                                //move net's num is less,capacity is big
                                const auto &lns = lhs->nets.size();
                                const auto &rns = rhs->nets.size();
                                const auto &lgs = circuit.g_set[lhs->name].size();
                                const auto &rgs = circuit.g_set[rhs->name].size();
                                // const auto &l_src=get_max_src(lhs);
                                // const auto &r_src=get_max_src(rhs);

                                if (lns  == rns ) {
                                    return lgs < rgs;
                                }
                                return lns-get_node_powerful_capacity(lhs) < rns-get_node_powerful_capacity(rhs) ;
                            });
                    }
                    // for(auto c:fpga_req[i]){
                    int k=fpga_req_index[i];
                    for(;k<fpga_req[i].size();k++){
                        auto c=fpga_req[i][k];
                        if(circuit.g_set[defer_node->name].count(c))continue;
                        auto origin_fpga=c->fpga_node;
                        if(origin_fpga!=hop_fpga.v[i])  continue;
                        Tensor<intg> cddt_move_(hop_fpga.num_vertex,1);
                        cddt_move_.at_is(origin_fpga->name,0);
                        if (try_move(c,cddt_move_))
                            continue;
                        intg tmp_dec = std::numeric_limits<intg>::min();
                        // intg tmp_dec = 0;
                        FPGANode *refine_f = nullptr;
                        intg move_increment=move_originf_connect_increment(c);
                        intg move_refine_f_increment=0;
                        origin_fpga->remove_circuit(c);
                        c->remove_fpga(origin_fpga);
                        for(int j=0;j<hop_fpga.num_vertex;++j){
                            if(cddt_move_.at(j)==1){
                                auto nf=hop_fpga.v[j];
                                if (nf == origin_fpga)
                                    continue;
                                intg td=0;
                                c->calculate_cut_increment(nf);
                                intg temp_move_increment=-c->cut_increment_map[nf->name]+move_increment;
                                td+=(intg)argu_conn*temp_move_increment;
                                td+=free_space_function(nf,c);
                                // td+=try_move_CV_space_increment(c, nf);
                                // //td求得是c移动到nf后的收益     
                                if (tmp_dec < td ) {
                                    if(valid(nf,c)>0){ 
                                        tmp_dec = td;
                                        refine_f = nf;
                                        move_refine_f_increment=temp_move_increment;
                                    }
                                } 
                                else if (tmp_dec == td && refine_f != nullptr ) {
                                    if(valid(nf,c)>0){ 
                                        if (hop_fpga.g_set[nf->name].size() <
                                            hop_fpga.g_set[refine_f->name].size()) {
                                            tmp_dec = td;
                                            refine_f = nf;
                                            move_refine_f_increment=temp_move_increment;
                                        }
                                    }
                                }
                            }
                        }//for(int j=0;j<hop_fpga.num_vertex;++j)
                        if(refine_f==nullptr){
                            origin_fpga->add_circuit(c);
                            c->add_fpga(origin_fpga);
                            continue;
                        }
                        refine_f->add_circuit(c);
                        c->add_fpga(refine_f);
                        hop_after_parition-=move_refine_f_increment;
                        if(valid(origin_fpga,defer_node)>0){
                            defer_node->calculate_cut_increment(origin_fpga);
                            hop_after_parition+=defer_node->cut_increment_map[origin_fpga->name];

                            origin_fpga->add_circuit(defer_node);
                            defer_node->add_fpga(origin_fpga);
                            break;
                        }
                    }//for(auto c:fpga_req[hop_fpga.v[i]])   
                    fpga_req_index[i]=k;              
                    //移动成功  
                    if(defer_node->fpga_node!=nullptr){
                        circuit_defer.erase(defer_node);
                        re_q.push_back(defer_node);
                        force_repair_cnt++;
                        break;
                    }
                }//if(nedd>0)
            }//if(cddt_move.at(i)==1)
        }//for(int i=0;i<hop_fpga.num_vertex;++i)
    }
    #ifdef LOG
    cout<<"preprocess successfully:"<<force_repair_cnt<<endl;
    #endif
        timer.output_time("finish force_CV");
}
double DB::get_node_powerful_capacity(CircuitNode *c){
    double powerful_capacity=0;
    for(int i=0;i<8;i++){
        powerful_capacity+=c->node_capacity[i]*rate_source[i];
    }  
    return powerful_capacity;
}
void DB::refine_CV() {
    #ifdef LOG
    cout << "Start CV refine." << endl;
    #endif
    if(circuit_defer.size()==0) return;
    // Refine with move-based. We try to move the boundary node, which is the
    // node that in the FPGA node that different with other node in the same
    // net.
    intg repair_cnt_therehold=circuit_defer.size()/10;
    int refine_CV_cnt=0;
    //至少为2，设置为1过不了case1
    int block_flag=1;
    intg direct_repair_cnt=0;
    intg refine_cnt=0;

    int fpga_capacity_sum=0;
    for(int i=0;i<8;++i)
        fpga_capacity_sum+=fpga.get_vertex(0)->fpga_capacity[i];
    double node_num_digit = log10(circuit.num_vertex) + 1;
    double argument_conn=log10(fpga.get_vertex(0)->wire_capacity) + 1;
    double capacity_digit = log10(fpga_capacity_sum) + 1;
    argument_conn=pow(10, capacity_digit-argument_conn+1);
    double node_capa_neighbor_ratio=1;
    //
    //    //fpga_sensitive的排序?
    vector<CircuitNode *> re_q;
    unordered_map<CircuitNode *,set<FPGANode *>> cddt_refine_f;
    for (auto &c : circuit.get_all_vertex()) {
        if(c->fpga_node==nullptr)   continue;
        re_q.emplace_back(c);
    }
    auto temp_place_cmp=[&](CircuitNode * lhs,CircuitNode* rhs){
        //鎸戦€夐偦灞呭ぇ锛岃祫婧愬ぇ鐨勫厛澶勭悊
        double l_powerful_capacity=get_node_powerful_capacity(lhs);
        double r_powerful_capacity=get_node_powerful_capacity(rhs);
        return node_capa_neighbor_ratio*l_powerful_capacity+(double)circuit.g[lhs->name].size()
                >(double)node_capa_neighbor_ratio*r_powerful_capacity+circuit.g[rhs->name].size();
    };
    vector<set<CircuitNode *,decltype(temp_place_cmp)>> fpga_sensitive_node(hop_fpga.num_vertex,set<CircuitNode *,decltype(temp_place_cmp)>(temp_place_cmp));
    unordered_map<CircuitNode *,unordered_set<FPGANode *>> node_sensitive_fpga;
    for(auto c:circuit.get_all_vertex()){
        if(c->fpga_node==nullptr){
                for(auto sensitive_fpga_:node_sensitive_fpga[c])
                    fpga_sensitive_node[sensitive_fpga_->name].erase(c);
                node_sensitive_fpga[c].clear();
                Tensor<intg> cddt_move(hop_fpga.num_vertex,1);
                if (try_move(c,cddt_move)){
                    continue;
                }
                for(int i=0;i<hop_fpga.num_vertex;++i)
                if(cddt_move.at(i)==1){
                    fpga_sensitive_node[i].emplace(c);
                    node_sensitive_fpga[c].emplace(fpga.v[i]);
                }//for(int i=0;i<hop_fpga.num_vertex;++i)
        }
    }//for(auto neighbor:circuit.g[c->name])

    while(circuit_defer.size()>0) {
        //        //先获得满足交换条件的node对?
        if(circuit_defer.size()==0)
            break;
        #ifdef LOG
        cout<<"------------------------"<<refine_CV_cnt+1<<"--------------------------"<<endl;
        #endif
        timer.output_time("Once refine");
        // Opt
        intg tmp_dec = std::numeric_limits<intg>::max();
        FPGANode *refine_f = nullptr;
        //鑾峰緱鏁忔劅fpga
        vector<intg> fpga_power(hop_fpga.num_vertex,0);
        for(int i=0;i<hop_fpga.num_vertex;++i){
            fpga_power[i]=fpga_sensitive_node[i].size();
        }
        // fpga_sensitive_node.clear();
        // node_sensitive_fpga.clear();
        //获得敏感fpga
        sort(re_q.begin(), re_q.end(),
             [&](CircuitNode *lhs, CircuitNode *rhs) {
                // //鍙備笌鐨刵et灏戠殑鍏堢Щ锟�?锟斤拷鍏讹拷?锟�?锟斤拷灞呭皯鐨勫厛绉诲姩
                // //鎵€鍦‵PGA绌洪棿绱у紶鐨勫厛绉诲姩
                 const auto &lns = lhs->nets.size();
                 const auto &rns = rhs->nets.size();
                 const auto &lhs_power=fpga_power[lhs->fpga_node->name];
                 const auto &rhs_power=fpga_power[rhs->fpga_node->name];
                 const auto &lgs = circuit.g_set[lhs->name].size();
                 const auto &rgs = circuit.g_set[rhs->name].size();
                if(lhs_power==rhs_power){
                    return lns<rns;
                }
                return lhs_power > rhs_power;     
             });

        // set<intg> sets[re_q.size()*hop_fpga.num_vertex];
        intg sets_cnt=0;
        direct_repair_cnt=0;
        refine_cnt=0;
        auto re_q_=re_q;
        for(int i=0;i<re_q_.size();i++){
            auto c=re_q_[i];
        // for (auto &c : re_q) {
            if(c->fpga_node==nullptr){continue;}
            if(circuit_defer.size()==0){
                break;
            }
            Tensor<intg> cddt_move(hop_fpga.num_vertex,1);
            cddt_move.at_is(c->fpga_node->name,0);
            if(try_move(c,cddt_move))   continue;
            cddt_refine_f[c].clear();
            for(int i=0;i<hop_fpga.num_vertex;++i){
                if(cddt_move.at(i)==1)
                    cddt_refine_f[c].emplace(hop_fpga.v[i]);
            }
            if(cddt_refine_f[c].size()<=0)  continue;

            // tmp_dec = std::numeric_limits<intg>::min();
            tmp_dec = 0;
            refine_f = nullptr;
            intg move_increment=move_originf_connect_increment(c);
            int move_refine_f_increment=0;
            auto remove_fpga=c->fpga_node;
            remove_fpga->remove_circuit(c);
            c->remove_fpga(remove_fpga);
            for (auto &nf : cddt_refine_f[c]) {
                if (nf == remove_fpga)
                    continue;

                c->calculate_cut_increment(nf);
                intg temp_move_increment=-c->cut_increment_map[nf->name]+move_increment;
                intg td=0;
                td+=(intg)argument_conn*temp_move_increment;
                td+=free_space_function(nf,c);
                // td+=try_move_CV_space_increment(c, nf);
                if (tmp_dec < td ) {
                    if(valid(nf,c)>0){ 
                        tmp_dec = td;
                        refine_f = nf;
                        move_refine_f_increment=temp_move_increment;
                    }
                } 
                else if (tmp_dec == td && refine_f != nullptr ) {
                    if(valid(nf,c)>0){ 
                    if (hop_fpga.g_set[nf->name].size() <
                        hop_fpga.g_set[refine_f->name].size()) {
                        tmp_dec = td;
                        refine_f = nf;
                        move_refine_f_increment=temp_move_increment;
                    }
                    }
                }
            }//for (auto &nf : cddt_refine_f[c])
            if (refine_f == nullptr){
                remove_fpga->add_circuit(c);
                c->add_fpga(remove_fpga);
                continue;
            }
            intg all_profit=tmp_dec;
            if (all_profit > 0||fpga_sensitive_node[remove_fpga->name].size()>0) {
                    refine_f->add_circuit(c);//47696
                    c->add_fpga(refine_f);
                    hop_after_parition-=move_refine_f_increment;
                    refine_cnt++;
                    if(1){
                        //重新计算敏感fpga
                        for(auto neighbor:circuit.g[c->name]){
                            if(neighbor->fpga_node==nullptr){
                                    //擦除所有的敏感fpga
                                    for(auto sensitive_fpga_:node_sensitive_fpga[neighbor])
                                        fpga_sensitive_node[sensitive_fpga_->name].erase(neighbor);
                                    node_sensitive_fpga[neighbor].clear();
                                    //重新计算
                                    Tensor<intg> cddt_move(hop_fpga.num_vertex,1);
                                    if (try_move(neighbor,cddt_move)){
                                        continue;
                                    }
                                    for(int i=0;i<hop_fpga.num_vertex;++i)
                                    if(cddt_move.at(i)==1){
                                        fpga_sensitive_node[i].emplace(neighbor);
                                        node_sensitive_fpga[neighbor].emplace(fpga.v[i]);
                                    }//for(int i=0;i<hop_fpga.num_vertex;++i)
                            }
                        }//for(auto neighbor:circuit.g[c->name])
                        if(fpga_sensitive_node[remove_fpga->name].size()>0){
                            for(auto CV_node_:fpga_sensitive_node[remove_fpga->name]){
                                if(CV_node_->fpga_node!=nullptr)    continue;
                                int valid_val=valid(remove_fpga,CV_node_);
                                if(valid_val>0){
                                    direct_repair_cnt++;
                                    circuit_defer.erase(CV_node_);
                                    //更新hop，并修复
                                    CV_node_->calculate_cut_increment(remove_fpga);
                                    hop_after_parition+=CV_node_->cut_increment_map[remove_fpga->name];

                                    remove_fpga->add_circuit(CV_node_);
                                    CV_node_->add_fpga(remove_fpga);
                                    re_q.push_back(CV_node_);
                                    //去除对应的敏感fpga
                                    for(auto sensitive_fpga_:node_sensitive_fpga[CV_node_])
                                        fpga_sensitive_node[sensitive_fpga_->name].erase(CV_node_); 
                                    node_sensitive_fpga[CV_node_].clear();
                                    break;
                                }//if(valid(removed_fpga,CV_node_)>0)
                                else if(valid_val==0){
                                    break;
                                }
                            }//for(auto CV_node_:fpga_sensitive_node[removed_fpga])
                        }//if(fpga_sensitive_node.count(remove_fpga)>0&&fpga_sensitive_node[remove_fpga].size()>0)
                    }//if(circuit_defer.size())
            }//if (all_profit > 0)
            else{
                remove_fpga->add_circuit(c);
                c->add_fpga(remove_fpga);
            }
        }//for (auto &c : re_q)

        if(circuit_defer.size()!=0&&direct_repair_cnt<=repair_cnt_therehold){
            //处理circuit_defer还在但是修复数目为0的情况
            #ifdef LOG
                cout<<"defer is live, direct repair cnt="<<direct_repair_cnt<<endl;
                //debug_resouce();
            #endif
            block_flag=0;
        }
        if(block_flag<=0){
            #ifdef LOG
            cout<<"block flag:"<<block_flag<<endl;
            #endif
            if(circuit_defer.size()){
                force_CV(argument_conn,re_q);   
                if(circuit_defer.size()==0) break; 
            }
            block_flag=1;
        }

        #ifdef LOG
        cout<<"direct repair cnt:"<<direct_repair_cnt<<endl;
        cout<<"refine cnt:"<<refine_cnt<<endl;
        if(circuit_defer.size()==0){
            // debug_resouce();
            cout<<"cur_hop_length:"<<hop_after_parition<<endl;
        }
        cout<<"--------------------------------------------------\n";
        #endif
        refine_CV_cnt++;
    }//for (int i = 0; i < 5; ++i)
        
    // if (timer.timeout())
    //     return;

#ifdef LOG_DB
    output_loss();
#endif
}

bool DB::topo_valid_undercopy(FPGANode *f,CircuitNode *c){
    /*这个函数是用来判断c放在f会不会发生topo违例的*/
    //对于驱动c的点src，保证f在src的copy_fpga或者fpga上
    //否则不要和fpga发生违例
    //对于c驱动的点drain_，保证drain_以及copy_fpga在c的copy_fpga或fpga上
    for(auto &n:c->src_nets){   
        //n作为src_net要驱动所有被驱动节点的fpga_node和copy_fpga
        //这些点要么在if_f的范围内，要么在c的copy_position   
        intg net_weight=n->weight;
        for(auto &fpga_:n->used_all_fpga){
            if(c->copy_position.count(fpga_))continue;
            if(fpga_weight[f->name][fpga_->name]==-1)
                return 0;
        } 
    }
    for(auto &n:c->drain_nets){    
        auto source=n->source;       
        auto copyf_node=source->copy_position;
        //如果c移动到if_f首先要考虑if_f在所有source_node的copy_position中或者与sourcefpga_node不违例
        if(n->used_all_fpga.count(f)||source->fpga_node==nullptr) continue;
        if(!source->copy_position.count(f)){
                if(fpga_weight[f->name][source->fpga_node->name]==-1)
                    return 0;
        }
        //求驱动copy_position的hop增加
        for(auto fpga_:c->copy_position){
            if(!source->copy_position.count(fpga_)){
                #ifdef LOG
                assert(fpga_weight[source->fpga_node->name][fpga_->name]!=-1);
                #endif
            }
        }// for(auto fpga_:copy_position)
    }
    return 1;
}

intg DB::valid_undercopy(FPGANode *f,CircuitNode *c){
    //c没有fpga
    //topo_valid 是-2是最严重的问题
    if(c->copy_position.count(f))   return 0;
    if(topo_valid_undercopy(f,c)==0){
        return -2;
    }
    //return (usage < fpga_capacity)
    assert(c->fpga_node==nullptr);
    assert(f!=nullptr);
    for(int i=0;i<8;i++)
        if(f->usage[i]+c->node_capacity[i]>f->fpga_capacity[i])
            return 0;
    vector<intg> used_wire(hop_fpga.num_vertex,0);

    for(auto &net:c->src_nets){   
        if(net->used_all_fpga.size()==0||net->used_all_fpga.count(f)&&net->used_all_fpga.size()==1)
            continue;
        //挑选出不被copy所覆盖，也跟当前的FPGA无关的FPGA
        auto c_copy_fpga=c->copy_position;
        int result =0;
        for(auto fpga_:net->used_all_fpga){
            if(c_copy_fpga.count(fpga_)||fpga_==f)
                continue;
            else{
                if(fpga_->usage_wire+net->weight+used_wire[fpga_->name]>fpga_->wire_capacity)
                    return -1;
                else
                    used_wire[fpga_->name]+=net->weight; 
                result=1;
            }
        }
        if(result==1)
            if(f->usage_wire+net->weight+used_wire[f->name]>f->wire_capacity)
                return -1;
            else
                used_wire[f->name]+=net->weight;
    }
    for(auto &net:c->drain_nets){ 
        auto source_c=net->source;
        auto source_f=source_c->fpga_node;
        auto source_copyf=source_c->copy_position;
        if(source_f==nullptr||net->used_all_fpga.count(f))
            continue;
        //如果src放置那还新增了fpga那自身肯定要处理
        if(f->usage_wire+net->weight+used_wire[f->name]>f->wire_capacity)
            return -1;
        else
            used_wire[f->name]+=net->weight;
        
        //如果source_f的fpga都是src-fpga或者src-copy的，那就会产生新的割
        #ifdef LOG
        assert(net->used_all_fpga.size()>=(1+source_copyf.size()));
        #endif
        if(net->used_all_fpga.size()>1+source_copyf.size())continue;
        if(source_f->usage_wire+net->weight+used_wire[source_f->name]>source_f->wire_capacity)
            return -1;
        else    
            used_wire[source_f->name]+=net->weight;
    }
    return 1;
}

intg DB::try_move_undercopy(CircuitNode* c,Tensor<intg> &cddt_move){
    //c移动去只要保证去src点的copy-circuit或者不与src_fpga发生违例
    //同时drain的copy_circuit和drain_fpga要么在c的copy_circuit上
    //要么在c的fpga上
    auto origin_fpga=c->fpga_node;
    c->remove_fpga(origin_fpga);
    auto fpga_list=hop_fpga.v;
    for(int i=0;i<hop_fpga.num_vertex;i++){
        if(fpga_list[i]==origin_fpga) continue;
        if(topo_valid_undercopy(fpga_list[i],c)>0)
            cddt_move.at(i)=1;
    }
    c->add_fpga(origin_fpga);
    return cddt_move.all_zero();    
}   


void DB::add_source(CircuitNode *c){
    all_copy_cnt++;
    for(int i=0;i<8;i++){
        sum_source[i]+=c->node_capacity[i];
        ave_source[i]=(double)sum_source[i]/hop_fpga.num_vertex;
    }
}
void DB::remove_source(CircuitNode *c){
    all_copy_cnt--;
    for(int i=0;i<8;i++){
        sum_source[i]-=c->node_capacity[i];
        ave_source[i]=(double)sum_source[i]/hop_fpga.num_vertex;
    }
}
bool DB::copy_valid(FPGANode *f,CircuitNode *c){
    //保证复制节点和原始节点不会出现在一起
    auto origin_fpga=c->fpga_node;
    auto c_copyf=c->copy_position;
    if(c->fpga_node==f||c_copyf.count(f))
        return false;
    for(auto net_:c->drain_nets){
        auto src_node=net_->source;
        auto src_fpga=src_node->fpga_node;
        //如果f跟src点会产生违例
        if(hop_fpga.g_set[src_fpga->name].count(f)==0)
            return false;
    }
    //资源上
    for(int i=0;i<8;i++){
        if(f->usage[i]+c->node_capacity[i]>f->fpga_capacity[i])
            return false;
    }
    //对外互联数目上
    vector<intg> used_wire(hop_fpga.num_vertex,0);
        //复制Node到当前FPGA，要考虑Node的src到Node产生的对外互联数目的增加 
    for(auto drain_:c->drain_nets){
        auto src_node=drain_->source;
        auto src_fpga=src_node->fpga_node;
        auto src_copy_fpga=src_node->copy_position;
        //这个情况下不会产生割
        if(src_fpga==f||src_copy_fpga.count(f))
            continue;
        //本来已经产生过割了
        if(drain_->used_copy_fpga.count(f)||drain_->used_fpga_node.count(f))
            continue;
        if(used_wire[f->name]+f->usage_wire+drain_->weight>f->wire_capacity)
            return false;
        used_wire[f->name]+=drain_->weight;
        //
        //result=0说明原来没有过割，要产生新的割
        bool result=0;
        for(auto fpga_:drain_->used_copy_fpga)
        if(fpga_!=src_fpga&&!src_copy_fpga.count(fpga_)){
            result=1;
            break;
        }
        if(result==0)
            for(auto fpga_:drain_->used_fpga_node)
            if(fpga_!=src_fpga&&!src_copy_fpga.count(fpga_)){
                result=1;
                break;
            }
        if(result==0){
            if(used_wire[src_fpga->name]+src_fpga->usage_wire+drain_->weight>src_fpga->wire_capacity)
                return false;
            used_wire[src_fpga->name]+=drain_->weight;
        }
    }//for(auto drain_:c->drain_nets)
    
    //复制c点，对this之外的被驱动点是没有影响的
    for(auto src_:c->src_nets){
        auto src_fpga=c->fpga_node;
        auto src_copy_fpga=c->copy_position;
        if(src_->used_copy_fpga.count(f)==0&&src_->used_fpga_node.count(f)==0)
            continue;
        if(src_copy_fpga.count(f)||src_fpga==f)     
            continue;    
        if(used_wire[f->name]+f->usage_wire+src_->weight>f->wire_capacity)
            return false;
        used_wire[f->name]+=src_->weight;
        //result=0说明没有没有和其他的FPGA有割，那就可以去除掉
        bool result=0;
        for(auto fpga_:src_->used_copy_fpga)
        if(fpga_!=src_fpga&&!src_copy_fpga.count(fpga_)&&fpga_!=f){
            result=1;
            break;
        }
        if(result==0)
            for(auto fpga_:src_->used_fpga_node)
            if(fpga_!=src_fpga&&!src_copy_fpga.count(fpga_)&&fpga_!=f){
                result=1;
                break;
            }
        if(result==0){
            if(used_wire[src_fpga->name]+src_fpga->usage_wire+src_->weight>src_fpga->wire_capacity)
                return false;
            used_wire[src_fpga->name]+=src_->weight;
        }
    }
    return true;
}
int DB::remove_copy_valid(CircuitNode *c,FPGANode *f){
    #ifdef LOG
    assert(c->copy_position.count(f));
    #endif
    for(auto neighbor_:dir_circuit.g[c->name]){
        if(neighbor_->fpga_node==f){
        //test if neighbor_ topo_invalid 
        for(auto neighbor__:circuit.g[neighbor_->name]){
            if(fpga_weight[neighbor__->fpga_node->name][neighbor_->fpga_node->name]==-1)
                return 0;
        }    
        }
    }
    return 1;
}
bool DB::try_copy(CircuitNode *c,Tensor<intg> &cddt_copy){
    //c的可以复制位置要保证在src的copy_position或者不和src_fpga违例
    intg max_cnt=0;
    for(auto net_:c->drain_nets){
        auto src_node=net_->source;
        auto src_fpga=src_node->fpga_node;
        auto src_copyf=src_node->copy_position;
        if(src_fpga==nullptr) continue;        
        if(!src_node->copy_position.count(src_fpga))
            cddt_copy.at(src_fpga->name)+=1;
        for(auto copy_f:src_node->copy_position){
            cddt_copy.at(copy_f->name)+=1;
        }
        for(auto neighbor_f:hop_fpga.g[src_fpga->name]){
            if(neighbor_f==src_fpga||src_node->copy_position.count(neighbor_f)) continue;
            cddt_copy.at(neighbor_f->name)+=1;
        }
        ++max_cnt;
    }
    cddt_copy.at(c->fpga_node->name)=0;
    int legal_f=0;
    for(int i=0;i<hop_fpga.num_vertex;i++){
        if(cddt_copy.at(i)!=max_cnt)
            cddt_copy.at(i)=0;
        else
            cddt_copy.at(i)=1,legal_f++;
    }
    if(legal_f==0)
        return false;
    else
        return true;
}


intg DB::try_move_increment_undercopy(CircuitNode *c,FPGANode *f){
    //记录没有fpga的c放置在f上的hop增长
    c->calculate_cut_increment_undercopy(f);
    return (c->cut_increment_undercopy_map[f->name]);
}
intg DB::free_space_function(FPGANode *f,CircuitNode *c) {
    //f的space与ave的距离
    //求移动后的收益
    // double dec = 0;
    // double max_capa=0;
    // for(int i=0;i<8;i++){
    //     if(c->node_capacity[i]==0||f->fpga_capacity[i]==0)  continue;
    //     max_capa=std::max(max_capa,(double)f->fpga_capacity[i]);
    //     dec+=(ave_source[i]-f->usage[i])/(double)c->node_capacity[i];
    // }
    // return max_capa*dec;
    double dec = 0;
    for(int i=0;i<8;i++){
        if(c->node_capacity[i]==0||f->fpga_capacity[i]==0)  continue;
        // dec+=(ave_source[i]-f->usage[i])*((double)ave_source[i]/(double)f->fpga_capacity[i]);
        dec+=1.0f*(ave_source[i]-f->usage[i])*rate_source[i];
    }
    return dec;
}
double DB::get_score(int total_hop){
    return (int)total_hop*(1.0f+0.2f*((double)timer.get_time()/3600));
}
bool DB::is_src_or_drain(CircuitNode *c,CircuitNode *verified){
    //检测verified是c的src还是drain，src输出0，drain输出1
    if(circuit.g_set[c->name].count(verified)!=0){
        if(dir_circuit.g_set[c->name].count(verified)!=0){
            //verified是c->name的drain
            return 1;
        }
        return 0;
    }
    return -1;
}
intg DB::try_move_copy_node_increment(CircuitNode *c,FPGANode *f){
    intg copy_increment=0;
    #ifdef LOG
    assert(c->copy_position.count(f));
    #endif
    for(auto drain_:c->drain_nets){
        auto src_node=drain_->source;
        auto src_fpga=src_node->fpga_node;
        auto src_copy_fpga=src_node->copy_position;
        if(src_fpga==nullptr)   continue;
    
        if(src_fpga==f||src_copy_fpga.count(f))
            continue;

        if(drain_->used_fpga_node.count(f)==0&&drain_->used_copy_fpga_count[f->name]==1)
            copy_increment+=drain_->weight*fpga_weight[f->name][src_fpga->name];
    }//for(auto drain_:c->drain_nets)

    for(auto src_:c->src_nets){
        auto src_fpga=c->fpga_node;
        auto src_copy_fpga=c->copy_position; 
        if(fpga_weight[f->name][src_fpga->name]<0) continue;;
        if(src_->used_fpga_node_count[f->name]+src_->used_copy_fpga_count[f->name]>1){
            copy_increment-=src_->weight*fpga_weight[f->name][src_fpga->name];
        }
    }
    return copy_increment;
}
void DB::mixed_refine(int loop_time,int copy_button){
    #ifdef LOG
    // last_hop=output_loss();
    cout << "Start refine." << endl;
    #endif
    intg last_hop=hop_after_parition;
    intg copy_cnt=0;
    intg debug_increment=0;
    int refine_button=1;
    int copy_start_time=2;
    double space_augu=0.1;
    if (timer.timeout())
        return ;

    vector<int> node_weight(circuit.num_vertex,0);
    for(auto c:circuit.get_all_vertex()){
        intg weight=0;
        intg rela_neighbor=0;
        for(auto net_:c->nets)
            weight+=net_->weight;
        weight*=circuit.g_set[c->name].size(); 
        node_weight[c->name]=weight;
    }
    vector<CircuitNode *> re_q=circuit.get_all_vertex();
    // unordered_map<CircuitNode *,set<FPGANode *>> cddt_refine_f;
    sort(re_q.begin(), re_q.end(),
            [&](const CircuitNode *lhs, const CircuitNode *rhs) {
            //参与的net少的先移动，其次是邻居少的先移动
            const auto &lns = lhs->nets.size();
            const auto &rns = rhs->nets.size();
            const auto &lgs = circuit.g_set[lhs->name].size();
            const auto &rgs = circuit.g_set[rhs->name].size();
            // intg l_weight=0;
            // intg l_rela_neighbor=0;
            // for(auto net_:lhs->nets)
            //     l_weight+=net_->weight;
            // l_weight*=lgs;
            // intg r_weight=0;
            // intg r_rela_neighbor=0;
            // for(auto net_:rhs->nets)
            //     r_weight+=net_->weight;
            // r_weight*=rgs;
            int l_weight=node_weight[lhs->name];
            int r_weight=node_weight[rhs->name];
            if (l_weight == r_weight) {
                return lgs < rgs;
            }
            return l_weight > r_weight;
            });

    for (int i = 0; i < loop_time; ++i) {
        //先获得满足交换条件的node对
        #ifdef LOG
        cout<<"------------------------"<<i+1<<"--------------------------\n";
        #endif
        if (timer.timeout())
            break;
    ///////////////////////////////
    /////refine              //////
    //////////////////////////////
        double start_score=get_score(last_hop);
        // Opt
        set<FPGANode *>cddt_refine_f;
        intg tmp_dec = std::numeric_limits<intg>::max();
        FPGANode *refine_f = nullptr;
        // FPGANode *perfect_f=nullptr;
        intg refine_cnt=0;
        copy_cnt=0;
        debug_increment=0;
        if(refine_button){
        for (auto &c : re_q) {
            //判断是否可以继续处理
            //if(c->fpga_node==nullptr){continue;}
            Tensor<intg> cddt_move(hop_fpga.num_vertex,0);
            if(try_move_undercopy(c,cddt_move))   continue;
            cddt_refine_f.clear();
            for(int i=0;i<hop_fpga.num_vertex;++i){
                if(cddt_move.at(i)==1)
                    cddt_refine_f.emplace(hop_fpga.v[i]);
            }
            if(cddt_refine_f.size()<=0)  continue;

            ////忽略资源和线的限制寻找收益更大的FPGA
            tmp_dec = std::numeric_limits<intg>::min();
            refine_f = nullptr;
            auto removed_fpga=c->fpga_node;
            removed_fpga->remove_circuit_undercopy(c);
            c->remove_fpga(removed_fpga);
            intg remove_origin_incremet=try_move_increment_undercopy(c,removed_fpga);

            for (auto &nf : cddt_refine_f) {
                if (nf == removed_fpga)
                    continue;
                // //td求得是c移动到nf后的收益 
                intg td =remove_origin_incremet;
                td-=try_move_increment_undercopy(c, nf);
                // intg td = c->try_move(nf);
                if (tmp_dec < td) {
                    if(valid_undercopy(nf,c)>0){
                        tmp_dec = td;
                        refine_f = nf;
                    }
                } 
                else if (tmp_dec == td && refine_f != nullptr) {
                    //如果收益等于最小收益，那就比其他东西
                    if (valid_undercopy(nf,c)>0&&hop_fpga.g_set[nf->name].size() >
                        hop_fpga.g_set[refine_f->name].size()) {
                        tmp_dec = td;
                        refine_f = nf;
                    }
                }
            }//for (auto &nf : cddt_refine_f[c])
            intg all_profit=tmp_dec;
            if (all_profit > 0&&refine_f!=nullptr) {
                    refine_f->add_circuit_undercopy(c);
                    c->add_fpga(refine_f);
                    refine_cnt++;
                    debug_increment+=tmp_dec;
                    continue;
            }//if (all_profit > 0&&refine_f!=nullptr)
            else if(all_profit==0&&refine_f!=nullptr){
                    //如果refine_f比remove_f的空闲空间多，可以移动过去
                    int move_space_incre=free_space_function(refine_f,c);
                    move_space_incre-=free_space_function(removed_fpga,c);

                    if(move_space_incre>0){
                        refine_f->add_circuit_undercopy(c);
                        c->add_fpga(refine_f);
                        refine_cnt++;
                        continue;
                    }else{
                        removed_fpga->add_circuit_undercopy(c);
                        c->add_fpga(removed_fpga);
                    }
            }//else if(all_profit==0&&refine_f!=nullptr)
            else{

                removed_fpga->add_circuit_undercopy(c);
                c->add_fpga(removed_fpga);
            }//else

        }//for (auto &c : re_q)
        }//if(refine_button)

        // if(1.0f*debug_increment/last_hop<=0.01){
        //     //In the end put all the copy node
        //     refine_button=0;
        // }

        //copy start previous
        if(i<copy_start_time&&1.0f*debug_increment/last_hop<=0.05){
            copy_start_time=i;
        }
        int remove_copy_cnt=0;
        if(i>=copy_start_time){
            if (timer.timeout())
                return ;  
    ///////////////////////////////
    ///// remove   copy     //////
    //////////////////////////////
            if(copy_button&&0){//remove copy node
            for(auto &c_:re_q){
                if(c_->copy_position.size()==0) continue;
                vector<FPGANode *> copy_f(c_->copy_position.begin(),c_->copy_position.end());
                sort(copy_f.begin(),copy_f.end(),[&](FPGANode *a,FPGANode *b){
                    return free_space_function(a,c_)<free_space_function(b,c_);
                });
                for(auto copy_f_:copy_f){
                    if(remove_copy_valid(c_,copy_f_)==0) continue;
                    int remove_copy_profit=try_move_copy_node_increment(c_,copy_f_);
                    copy_f_->remove_copy_circuit(c_);
                    c_->remove_copy_fpga(copy_f_);
                    debug_increment+=remove_copy_profit;
                    remove_source(c_);
                    // copy_f_->add_copy_circuit(c_);
                    // c_->add_copy_fpga(copy_f_);
                    // debug_increment-=remove_copy_profit;
                    // add_source(c_);
                    intg best_profit=std::numeric_limits<intg>::min();
                    CircuitNode* best_refine_c=nullptr;
                    FPGANode* best_refine_f=nullptr;
                    vector<CircuitNode *> refine_c_vec(circuit.g[c_->name]);
                    // refine_c_vec.push_back(c_);
                    for(auto c:refine_c_vec){
                        if(dir_circuit.g_set[c_->name].count(c))
                            continue;
                        Tensor<intg> cddt_move(hop_fpga.num_vertex,0);
                        if(try_move_undercopy(c,cddt_move))   continue;
                        cddt_refine_f.clear();
                        for(int i=0;i<hop_fpga.num_vertex;++i){
                            if(cddt_move.at(i)==1)
                                cddt_refine_f.emplace(hop_fpga.v[i]);
                        }
                        if(cddt_refine_f.size()<=0)  continue;

                        ////忽略资源和线的限制寻找收益更大的FPGA
                        tmp_dec = std::numeric_limits<intg>::min();
                        refine_f = nullptr;
                        auto removed_fpga=c->fpga_node;
                        removed_fpga->remove_circuit_undercopy(c);
                        c->remove_fpga(removed_fpga);
                        intg remove_origin_incremet=try_move_increment_undercopy(c,removed_fpga);

                        for (auto &nf : cddt_refine_f) {
                            if (nf == removed_fpga)
                                continue;
                            intg td =remove_origin_incremet;
                            td-=try_move_increment_undercopy(c, nf);
                            // intg td = c->try_move(nf);
                            if (tmp_dec < td) {
                                if(valid_undercopy(nf,c)>0){
                                    tmp_dec = td;
                                    refine_f = nf;
                                }
                            } 
                            else if (tmp_dec == td && refine_f != nullptr) {
                                if (valid_undercopy(nf,c)>0&&hop_fpga.g_set[nf->name].size() <
                                    hop_fpga.g_set[refine_f->name].size()) {
                                    tmp_dec = td;
                                    refine_f = nf;
                                }
                            }
                        }//for (auto &nf : cddt_refine_f[c])
                        intg all_profit=tmp_dec;
                        if (all_profit > best_profit&&refine_f!=nullptr) {
 
                                best_profit=all_profit;
                                best_refine_c=c;
                                best_refine_f=refine_f;
                                // continue;
                        }//if (all_profit > 0&&refine_f!=nullptr)
                        else if(all_profit==best_profit&&refine_f!=nullptr){
                                //濡傛灉refine_f姣攔emove_f鐨勭┖闂茬┖闂达拷?锛屽彲浠ョЩ鍔ㄨ繃锟�?
                                int move_space_incre=free_space_function(refine_f,c);
                                move_space_incre-=free_space_function(removed_fpga,c);

                                if(move_space_incre>0){
                                    best_profit=all_profit;
                                    best_refine_c=c;
                                    best_refine_f=refine_f;
                                }
                        }
                        removed_fpga->add_circuit_undercopy(c);
                        c->add_fpga(removed_fpga);
                    }//for(auto c:circuit.g[c_->name])

                    // //calculate selfnode can move to copy_f
                    // auto removed_fpga=c_->fpga_node;
                    // removed_fpga->remove_circuit_undercopy(c_);
                    // c_->remove_fpga(removed_fpga);
                    // if(valid_undercopy(copy_f_,c_)>0&&0){
                    // int self_profit=try_move_increment_undercopy(c_,removed_fpga)-try_move_increment_undercopy(c_, copy_f_);
                    //     if (self_profit > best_profit) {
                    //             best_profit=self_profit;
                    //             best_refine_c=c_;
                    //             best_refine_f=copy_f_;
                    //     }//if (all_profit > 0&&refine_f!=nullptr)
                    //     else if(self_profit ==best_profit){
                    //         int move_space_incre=free_space_function(copy_f_,c_);
                    //         move_space_incre-=free_space_function(removed_fpga,c_);
                    //         if(move_space_incre>0){
                    //             best_profit=self_profit;
                    //             best_refine_c=c_;
                    //             best_refine_f=copy_f_;
                    //         }
                    //     }
                    // }
                    // removed_fpga->add_circuit_undercopy(c_);
                    // c_->add_fpga(removed_fpga);


                    //refine
                    if(best_profit>=remove_copy_profit&&best_refine_f!=nullptr){
                        auto origin_fpga=best_refine_c->fpga_node;
                        origin_fpga->remove_circuit_undercopy(best_refine_c);
                        best_refine_c->remove_fpga(origin_fpga);
                        best_refine_f->add_circuit_undercopy(best_refine_c);
                        best_refine_c->add_fpga(best_refine_f);
                        debug_increment+=best_profit;
                        refine_cnt++;
                        remove_copy_cnt++;
                    }else{
                        copy_f_->add_copy_circuit(c_);
                        c_->add_copy_fpga(copy_f_);
                        debug_increment-=remove_copy_profit;
                        add_source(c_);
                    }
                }//for(auto copy_f_:c_->copy_position)
            }
            }//if(copy_button)

    ///////////////////////////////
    /////copy                ////////
    //////////////////////////////
            if(copy_button){//copy 
            for(auto &c:re_q){
                auto copy_c=c;
                    Tensor<intg> cddt_copy(hop_fpga.num_vertex,0);
                    if(!try_copy(copy_c,cddt_copy))   continue;
                    vector<FPGANode*> cddt_copy_f;
                    for(int i=0;i<hop_fpga.num_vertex;++i){
                        // if(cddt_copy.at(i)==1&&advise_copy(hop_fpga.v[i],c))
                        if(cddt_copy.at(i)==1)
                            cddt_copy_f.emplace_back(hop_fpga.v[i]);
                    }
                    if(cddt_copy_f.size()<=0){  
                        continue;
                    }

                    intg tmp_dec_copy = std::numeric_limits<intg>::min();
                    FPGANode* perfect_copy_f=nullptr;
                    int perfect_copy_profit=0;
                    // tmp_dec = std::numeric_limits<intg>::max();
                    FPGANode* copy_f = nullptr;
                    for (auto &nf : cddt_copy_f) {
                        if (nf == copy_c->fpga_node)
                            continue;
                        
                        copy_c->calculate_copy_increment(nf);
                        intg td_copy = -copy_c->copy_increment_map[nf->name];
                        // intg td = copy_c->try_move(nf);
                        // //td求得是c复制到nf后的收益 
                        if (tmp_dec_copy < td_copy) {
                            if(copy_valid(nf,copy_c)>0){
                            tmp_dec_copy = td_copy;
                                copy_f = nf;
                            }
                        } 
                        else if (tmp_dec_copy == td_copy && copy_f != nullptr) {
                            //如果收益等于最小收益，那就比其他东西
                            if (copy_valid(nf,copy_c)>0&&hop_fpga.g_set[nf->name].size() >
                                hop_fpga.g_set[copy_f->name].size()) {
                                tmp_dec_copy = td_copy;
                                copy_f = nf;
                            }
                        }

                        if(perfect_copy_profit<td_copy){
                            perfect_copy_profit=td_copy;
                            perfect_copy_f=nf; 
                        }else if(perfect_copy_profit==td_copy&&perfect_copy_f!=nullptr){
                            if (hop_fpga.g_set[nf->name].size() >
                                hop_fpga.g_set[perfect_copy_f->name].size()) {
                                perfect_copy_profit=td_copy;
                                perfect_copy_f=nf;
                            }
                        }

                    }//for (auto &nf : cddt_refine_f[c])
                    if (copy_f == nullptr)
                        continue;
                    intg copy_profit=tmp_dec_copy;
                    if (copy_profit > 0) {
                        copy_f->add_copy_circuit(copy_c);
                        copy_c->add_copy_fpga(copy_f);
                        add_source(copy_c); 
                        copy_cnt++;
                        debug_increment+=copy_profit;
                        continue;
                    }//if (all_profit > 0)  
                    else if (0){
                        //try to remove some copy node to copy
                        // if(perfect_copy_f==nullptr) continue;
                        // vector<CircuitNode*> all_copy_c(perfect_copy_f->copy_circuit.begin(),perfect_copy_f->copy_circuit.end());
                        // sort(all_copy_c.begin(),all_copy_c.end(),[&](CircuitNode* a,CircuitNode* b){
                        //     return try_move_copy_node_increment(a,perfect_copy_f)>try_move_copy_node_increment(b,perfect_copy_f);
                        // });
                        // for(auto cddt_remove_c:all_copy_c){
                        //     if(remove_copy_valid(cddt_remove_c,perfect_copy_f)<=0) continue;
                        //     int remove_copy_profit=try_move_copy_node_increment(cddt_remove_c,perfect_copy_f);
                        //     if(remove_copy_profit+perfect_copy_profit<=0) continue;
                        //     perfect_copy_f->remove_copy_circuit(cddt_remove_c);
                        //     cddt_remove_c->remove_copy_fpga(perfect_copy_f);
                        //     if(copy_valid(perfect_copy_f,copy_c)>0){
                        //         perfect_copy_f->add_copy_circuit(copy_c);
                        //         copy_c->add_copy_fpga(perfect_copy_f);
                        //         remove_source(cddt_remove_c);
                        //         add_source(copy_c);
                        //         remove_copy_cnt++;
                        //         debug_increment+=remove_copy_profit+perfect_copy_profit;
                        //     }else{
                        //         perfect_copy_f->add_copy_circuit(cddt_remove_c);
                        //         cddt_remove_c->add_copy_fpga(perfect_copy_f);
                        //     }
                        // }
                    }
            }//for(auto &c:re_q)
            }//if(copy_button)
        }//if(i>=copy_start_time)

        intg cur_hop=last_hop-debug_increment;
        hop_after_parition=cur_hop;
        #ifdef LOG
        cout<<"refine cnt:"<<refine_cnt<<endl;
        cout<<"increment:"<<debug_increment<<endl;
        cout<<"real_incremetn:"<<last_hop-cur_hop<<endl;
        cout<<"copy_cnt:"<<copy_cnt<<endl;
        cout<<"remove_copy_cnt:"<<remove_copy_cnt<<endl;
        if(circuit_defer.size()==0){
            cout<<"cur_hop_length:"<<cur_hop<<endl;
        }
        cout<<"cur_score:"<<(int)get_score(cur_hop)<<endl;
        cout<<"--------------------------------------------------\n";
        #endif
        last_hop=cur_hop;
        double end_score=get_score(cur_hop);
        if((start_score<end_score)&&i>=copy_start_time+1){
            #ifdef LOG
                cout<<"hop break"<<endl;
            #endif
            break;
        }
    }//for (int i = 0; i < 5; ++i)


    // if(remove_button){
    //     remove_some_node(re_q);
    // }

#ifdef LOG_DB
    output_loss();
#endif
}
void DB::remove_some_node(vector<CircuitNode *> &re_q){
    // for(auto &c:re_q){
    //     for(auto copy_c:c->copy_position){
            
    //     }
    // }
    int node_num=circuit.num_vertex/5;
    int remove_cnt=0;
    for(int i=re_q.size()-1;i>=re_q.size()-node_num;i--){
        remove_cnt++;
        auto c=re_q[i];
        auto origin_f=c->fpga_node;
        origin_f->remove_circuit_undercopy(c);
        c->remove_fpga(origin_f);
        circuit_defer.emplace(c);
        hop_after_parition-=try_move_increment_undercopy(c,origin_f);
    }
    #ifdef LOG
    cout<<"remove_cnt:"<<remove_cnt<<endl;
    #endif
}
// void DB::build_merge_graph() {
// //    Graph() {} void init(vector<T *> &vv, vector<vector<T *>> &gg,vector<vector<intg>> &e_weight);
// }

// void DB::merge_vertex(MergeNode *c,MergeNode *new_node,vector<MergeNode *> &g,vector<vector<intg>> &e_weight){
//     new_node->merge(c);
//     for(auto neighbor:merge_graph.g[c->name]){
//         if(neighbor->fpga_node==c->fpga_node){
//             new_node->merge(neighbor);
//         }
//     }
// }

intg DB::try_multimove_undercopy(unordered_set<CircuitNode *> &union_c,Tensor<intg> &cddt_move){
    for(auto c_:union_c){
        Tensor<intg> tmp_cddt(hop_fpga.num_vertex,0);
        if(try_move_undercopy(c_,tmp_cddt))
            return 0;
        cddt_move+=tmp_cddt;
    }
    return 1;
}
intg DB::try_multimove_increment_undercopy(MergeNode &c, FPGANode *if_f){
    intg cut_increment = 0;
    // for(auto c_:c.belong_nodes){
    //     c_->calculate_cut_increment_undercopy(if_f);
    //     cut_increment+=c_->cut_increment_undercopy_map[if_f->name];
    // }
    // return cut_increment;
    for(auto &n:c.src_nets){   
        //this放在if_f上，那就要计算被驱动点的fpga_node和copy_node和当前点的if_f以及copy_position的关系
        // if(c.common_nets.count(n)) continue;
        intg net_weight=n->weight;
        auto src_node=n->source;
        for(auto &fpga_:n->used_all_fpga){
            if(src_node->copy_position.count(fpga_)||fpga_==if_f)continue;
            #ifdef LOG
                assert(fpga_weight[if_f->name][fpga_->name]!=-1);
            #endif
            intg cur_weight = fpga_weight[if_f->name][fpga_->name];
            cut_increment+=cur_weight*net_weight;
        } 
    }

    for(auto &n:c.drain_nets){    
        auto source_=n->source;       
        auto copyf_node=source_->copy_position;
        //如果c移动到if_f首先要考虑if_f在所有source_node的copy_position中或者与sourcefpga_node不违例
        if(n->used_all_fpga.count(if_f)||source_->fpga_node==nullptr) continue;
        if(!source_->copy_position.count(if_f)&&if_f!=source_->fpga_node){
            #ifdef LOG
                assert(fpga_weight[if_f->name][source_->fpga_node->name]!=-1);
            #endif
            cut_increment+=n->weight*fpga_weight[if_f->name][source_->fpga_node->name];
        }
    }

    // for(auto &n:c.common_nets){    

    // }
    return cut_increment;
}
int DB::free_space_function(MergeNode &c,FPGANode *if_f){
    double dec = 0;
    for(int i=0;i<8;i++){
        if(c.node_capacity[i]==0||if_f->fpga_capacity[i]==0)  continue;
        // dec+=(ave_source[i]-f->usage[i])*((double)ave_source[i]/(double)f->fpga_capacity[i]);
        dec+=1.0f*(ave_source[i]-if_f->usage[i])*rate_source[i];
    }
    return dec;
}
int DB::topo_multi_valid_undercopy(FPGANode* if_f,MergeNode &c){
    for(auto node_:c.belong_nodes){
        int topo_val=topo_valid(if_f,node_);
        if(topo_val<=0)
            return topo_val;
    }
    return 1;
}
int DB::multi_valid_undercopy(FPGANode *if_f,MergeNode &c){
    if(topo_multi_valid_undercopy(if_f,c)==0){
        return -2;
    }
    //return (usage < fpga_capacity)
    #ifdef LOG
        assert(c.fpga_node==nullptr);
        assert(if_f!=nullptr);
    #endif
    for(int i=0;i<8;i++)
        if(if_f->usage[i]+c.node_capacity[i]>if_f->fpga_capacity[i])
            return 0;
    vector<intg> used_wire(hop_fpga.num_vertex,0);

    for(auto &net:c.src_nets){   
        // if(c.common_nets.count(net)) continue;
        auto src_node=net->source;
        if(net->used_all_fpga.size()==0||net->used_all_fpga.count(if_f)&&net->used_all_fpga.size()==1)
            continue;

        auto c_copy_fpga=src_node->copy_position;
        int result =0;
        for(auto fpga_:net->used_all_fpga){
            if(c_copy_fpga.count(fpga_)||fpga_==if_f)
                continue;
            else{
                if(fpga_->usage_wire+net->weight+used_wire[fpga_->name]>fpga_->wire_capacity)
                    return -1;
                else
                    used_wire[fpga_->name]+=net->weight; 
                result=1;
            }
        }
        if(result==1)
            if(if_f->usage_wire+net->weight+used_wire[if_f->name]>if_f->wire_capacity)
                return -1;
            else
                used_wire[if_f->name]+=net->weight;
    }
    for(auto &net:c.drain_nets){ 
        auto source_c=net->source;
        auto source_f=source_c->fpga_node;
        auto source_copyf=source_c->copy_position;
        if(source_f==nullptr||net->used_all_fpga.count(if_f))
            continue;
        //濡傛灉src鏀剧疆閭ｈ繕鏂帮拷?浜唂pga閭ｈ嚜锟�?锟斤拷瀹氾拷?澶勭悊
        if(if_f->usage_wire+net->weight+used_wire[if_f->name]>if_f->wire_capacity)
            return -1;
        else
            used_wire[if_f->name]+=net->weight;
        
        //濡傛灉source_f鐨刦pga閮芥槸src-fpga鎴栬€卻rc-copy鐨勶紝閭ｅ氨浼氫骇鐢熸柊鐨勫壊
        #ifdef LOG
        assert(net->used_all_fpga.size()>=(1+source_copyf.size()));
        #endif
        if(net->used_all_fpga.size()>1+source_copyf.size())continue;
        if(source_f->usage_wire+net->weight+used_wire[source_f->name]>source_f->wire_capacity)
            return -1;
        else    
            used_wire[source_f->name]+=net->weight;
    }

    // for(auto net:c.common_nets){
    //     auto src_node=net->source;
    //     #ifdef LOG
    //         assert(c.belong_nodes.count(src_node)); 
    //     #endif

    // }
    return 1;    
}
void DB::multi_refine(int loop_time){
    #ifdef LOG
        cout<<"multi_refine"<<endl;
    #endif

    intg last_hop=hop_after_parition;
    intg copy_cnt=0;
    intg debug_increment=0;
    int copy_start_time=2;
    double space_augu=0.1;
    if (timer.timeout())
        return ;

    vector<int> node_weight(circuit.num_vertex,0);
    for(auto c:circuit.get_all_vertex()){
        #ifdef LOG
            assert(c->fpga_node!=nullptr);
        #endif
        intg weight=0;
        intg rela_neighbor=0;
        for(auto net_:c->nets)
            weight+=net_->weight;
        weight*=circuit.g_set[c->name].size(); 
        node_weight[c->name]=weight;
    }
    vector<CircuitNode *> re_q=circuit.get_all_vertex();
    // unordered_map<CircuitNode *,set<FPGANode *>> cddt_refine_f;
    sort(re_q.begin(), re_q.end(),
            [&](const CircuitNode *lhs, const CircuitNode *rhs) {
            //鍙備笌鐨刵et灏戠殑鍏堢Щ锟�?锟斤拷鍏讹拷?锟�?锟斤拷灞呭皯鐨勫厛绉诲姩
            const auto &lns = lhs->nets.size();
            const auto &rns = rhs->nets.size();
            const auto &lgs = circuit.g_set[lhs->name].size();
            const auto &rgs = circuit.g_set[rhs->name].size();
            int l_weight=node_weight[lhs->name];
            int r_weight=node_weight[rhs->name];
            if (l_weight == r_weight) {
                return lgs < rgs;
            }
            return l_weight > r_weight;
            });

    for (int i = 0; i < loop_time; ++i) {
        //鍏堣幏寰楁弧瓒充氦鎹㈡潯浠剁殑node锟�?
        #ifdef LOG
        cout<<"------------------------"<<i+1<<"--------------------------\n";
        #endif
        if (timer.timeout())
            break;

        double start_score=get_score(last_hop);
        // Opt
        intg tmp_dec = std::numeric_limits<intg>::max();
        FPGANode *refine_f = nullptr;
        // FPGANode *perfect_f=nullptr;
        intg refine_cnt=0;
        copy_cnt=0;
        debug_increment=0;
        for (auto &c : re_q) {
            #ifdef LOG
                assert(c->fpga_node!=nullptr);
            #endif
            unordered_set<CircuitNode *> union_c;
            set<FPGANode *> cddt_refine_f;
            MergeNode best_profit_mergenode(0,c,union_c,nullptr);
            int best_profit=0;
            FPGANode *best_refine_f=nullptr;
            for(auto net_:c->nets){
                // //choose one net_ is best to move
                // if(net_==nullptr||c_net->weight*c_net->used_all_fpga.size()>net_->weight*net_->used_all_fpga.size()&&
                // c_net->used_all_fpga.size()!=c_net->source->copy_position.size()+1){
                //     #ifdef LOG
                //         assert(!c_net->source->copy_position.count(c_net->source->fpga_node));
                //     #endif
                //     net_=c_net;
                // }
                if(net_->used_all_fpga.size()!=net_->source->copy_position.size()+1) continue;
                union_c.clear();
                cddt_refine_f.clear();
                auto removed_fpga=c->fpga_node;
                //记录下对于当前net可以移动的节点集合
                for(auto neighbor_:net_->net_cell){
                    if(neighbor_->fpga_node==c->fpga_node)
                        union_c.emplace(neighbor_);
                    else
                        continue;
                }
                Tensor<intg> cddt_move(hop_fpga.num_vertex,0);
                if(try_multimove_undercopy(union_c,cddt_move))   continue;
                cddt_refine_f.clear();
                cddt_move.at(removed_fpga->name)=0;
                for(int i=0;i<hop_fpga.num_vertex;++i){
                    if(cddt_move.at(i)==union_c.size())
                        cddt_refine_f.emplace(hop_fpga.v[i]);
                }
                if(cddt_refine_f.size()<=0)  continue;

                MergeNode merge_node(0,c,union_c,removed_fpga);
                tmp_dec = std::numeric_limits<intg>::min();
                refine_f = nullptr;

                merge_node.remove_fpga_circuit(removed_fpga);    
                intg remove_origin_incremet=try_multimove_increment_undercopy(merge_node,removed_fpga);

                for (auto &nf : cddt_refine_f) {
                    intg td =remove_origin_incremet;
                    td-=try_multimove_increment_undercopy(merge_node, nf);
                    if (tmp_dec < td) {
                        if(multi_valid_undercopy(nf,merge_node)>0){
                            tmp_dec = td;
                            refine_f = nf;
                        }
                    } 
                    else if (tmp_dec == td && refine_f != nullptr) {
                        //濡傛灉鏀剁泭绛変簬鏈€灏忔敹鐩婏紝閭ｅ氨姣斿叾浠栦笢锟�?
                        if (multi_valid_undercopy(nf,merge_node)>0&&hop_fpga.g_set[nf->name].size() <
                            hop_fpga.g_set[refine_f->name].size()) {
                            tmp_dec = td;
                            refine_f = nf;
                        }
                    }
                }
                intg all_profit=tmp_dec;
                if (all_profit > best_profit&&refine_f!=nullptr) {
                        // refine_f->add_circuit_undercopy(c);
                        best_profit_mergenode = merge_node;
                        best_profit = all_profit;
                        best_refine_f = refine_f;
                        merge_node.add_fpga_circuit(removed_fpga);
                        continue;
                }//if (all_profit > 0&&refine_f!=nullptr)
                else if(all_profit==best_profit&&refine_f!=nullptr){
                        //濡傛灉refine_f姣攔emove_f鐨勭┖闂茬┖闂达拷?锛屽彲浠ョЩ鍔ㄨ繃锟�?
                        int move_space_incre=free_space_function(merge_node,refine_f);
                        move_space_incre-=free_space_function(merge_node,removed_fpga);

                        if(move_space_incre>0){
                            // refine_f->add_circuit_undercopy(c);
                            best_profit_mergenode = merge_node;
                            best_profit = all_profit;
                            best_refine_f = refine_f;
                            merge_node.add_fpga_circuit(removed_fpga);
                            continue;
                        }else{
                            // removed_fpga->add_circuit_undercopy(c);
                            merge_node.add_fpga_circuit(removed_fpga);
                        }
                }//else if(all_profit==0&&refine_f!=nullptr)
                else{
                    //refine_f涓嶄负nullptr锛宎llprofit灏忎簬0锛岋拷?鏄庯拷?鍘讳竴锟�?锟斤拷鏃风殑绌洪棿锛屽彲浠ヨ€冭檻鍏ㄤ綋绉绘皯璇曡瘯
                    merge_node.add_fpga_circuit(removed_fpga);
                }//else
            }//for(auto net_:c->nets)

            auto removed_fpga=best_profit_mergenode.fpga_node;
            if(removed_fpga==nullptr) continue;
            best_profit_mergenode.remove_fpga_circuit(removed_fpga);
            if (best_profit > 0&&best_refine_f!=nullptr) {
                    // refine_f->add_circuit_undercopy(c);
                    best_profit_mergenode.add_fpga_circuit(best_refine_f);
                    refine_cnt+=best_profit_mergenode.belong_nodes.size();
                    debug_increment+=best_profit;
                    continue;
            }//if (all_profit > 0&&refine_f!=nullptr)
            else if(best_profit==0&&best_refine_f!=nullptr){

                    int move_space_incre=free_space_function(best_profit_mergenode,best_refine_f);
                    move_space_incre-=free_space_function(best_profit_mergenode,removed_fpga);

                    if(move_space_incre>0){
                        // refine_f->add_circuit_undercopy(c);
                        best_profit_mergenode.add_fpga_circuit(best_refine_f);
                        refine_cnt+=union_c.size();
                        continue;
                    }else{
                        // removed_fpga->add_circuit_undercopy(c);
                        best_profit_mergenode.add_fpga_circuit(removed_fpga);
                    }
            }//else if(all_profit==0&&refine_f!=nullptr)
            else{
                //refine_f涓嶄负nullptr锛宎llprofit灏忎簬0锛岋拷?鏄庯拷?鍘讳竴锟�?锟斤拷鏃风殑绌洪棿锛屽彲浠ヨ€冭檻鍏ㄤ綋绉绘皯璇曡瘯
                best_profit_mergenode.add_fpga_circuit(removed_fpga);
            }//else
        }//for (auto &c : re_q)

        //Start copy
        if(i<copy_start_time&&1.0f*debug_increment/last_hop<=0.05){
            copy_start_time=i;
        }
        //璇曢獙锟�?锟斤拷鏃跺紕姣旇緝濂斤紝杩樻槸涔嬪悗寮勬瘮杈冨ソ
        if(i>=copy_start_time){
            if (timer.timeout())
                return ;  
            if(1){//澶嶅埗寮€锟�
            for(auto &c:re_q){
            // for(auto &c:wait_for_copy){
            // for(auto &c:vector_copy){
                auto copy_c=c;
                    Tensor<intg> cddt_copy(hop_fpga.num_vertex,0);
                    if(!try_copy(copy_c,cddt_copy))   continue;
                    vector<FPGANode*> cddt_copy_f;
                    for(int i=0;i<hop_fpga.num_vertex;++i){
                        // if(cddt_copy.at(i)==1&&advise_copy(hop_fpga.v[i],c))
                        if(cddt_copy.at(i)==1)
                            cddt_copy_f.emplace_back(hop_fpga.v[i]);
                    }
                    if(cddt_copy_f.size()<=0){  
                        continue;
                    }
                    ////蹇界暐璧勬簮鍜岀嚎鐨勯檺鍒讹拷?鎵炬敹鐩婃洿澶х殑FPGA
                    intg tmp_dec_copy = std::numeric_limits<intg>::min();
                    // tmp_dec = std::numeric_limits<intg>::max();
                    FPGANode* copy_f = nullptr;
                    for (auto &nf : cddt_copy_f) {
                        if (nf == copy_c->fpga_node)
                            continue;
                        // //td姹傚緱鏄痗澶嶅埗鍒皀f鍚庣殑鏀剁泭 
                        copy_c->calculate_copy_increment(nf);
                        intg td_copy = -copy_c->copy_increment_map[nf->name];
                        // intg td = copy_c->try_move(nf);
                        if (tmp_dec_copy < td_copy) {
                            if(copy_valid(nf,copy_c)>0){
                            tmp_dec_copy = td_copy;
                                copy_f = nf;
                            }
                        } 
                        else if (tmp_dec_copy == td_copy && copy_f != nullptr) {
                            //濡傛灉鏀剁泭绛変簬鏈€灏忔敹鐩婏紝閭ｅ氨姣斿叾浠栦笢锟�?
                            if (copy_valid(nf,copy_c)>0&&hop_fpga.g_set[nf->name].size() <
                                hop_fpga.g_set[copy_f->name].size()) {
                                tmp_dec_copy = td_copy;
                                copy_f = nf;
                            }
                        }
                    }//for (auto &nf : cddt_refine_f[c])
                    if (copy_f == nullptr)
                        continue;
                    intg copy_profit=tmp_dec_copy;
                    if (copy_profit > 0) {
                        copy_f->add_copy_circuit(copy_c);
                        copy_c->add_copy_fpga(copy_f);
                        add_source(copy_c); 
                        copy_cnt++;
                        debug_increment+=copy_profit;
                        continue;
                    }//if (all_profit > 0)  
                    else if (copy_profit==0&&0){
                        copy_f->add_copy_circuit(copy_c);
                        copy_c->add_copy_fpga(copy_f);
                        add_source(copy_c); 
                        copy_cnt++;
                        debug_increment+=copy_profit;
                        continue;
                    }//else if (copy_profit==0&&free_space_function(copy_f,c)>=0)
            }//for(auto &c:re_q)
            }//if(1)
        }//if(i>=copy_start_time)


        intg cur_hop=last_hop-debug_increment;
        hop_after_parition=cur_hop;
        #ifdef LOG
        cout<<"refine cnt:"<<refine_cnt<<endl;
        cout<<"increment:"<<debug_increment<<endl;
        cout<<"real_incremetn:"<<last_hop-cur_hop<<endl;
        cout<<"copy_cnt:"<<copy_cnt<<endl;
        if(circuit_defer.size()==0){
            cout<<"cur_hop_length:"<<cur_hop<<endl;
        }
        cout<<"cur_score:"<<(int)get_score(cur_hop)<<endl;
        cout<<"--------------------------------------------------\n";
        #endif
        last_hop=cur_hop;
        double end_score=get_score(cur_hop);
        if((start_score<end_score)&&i>=copy_start_time+1){
            #ifdef LOG
                cout<<"hop break"<<endl;
            #endif
            break;
        }
    }//for (int i = 0; i < 5; ++i)
        
    // if (timer.timeout())
    //     return;

#ifdef LOG_DB
    output_loss();
#endif
}

intg DB::output_loss() {
    intg total_hop=0,fixed_hop=0;
    for(auto net:nets){     
        FPGANode *src_fpga=net->source->fpga_node;
        if(src_fpga==nullptr) continue;
        for(auto drain_fpga:net->used_all_fpga){
            if(drain_fpga==src_fpga||net->source->copy_position.count(drain_fpga))  continue;
                total_hop+=net->weight*fpga_weight[src_fpga->name][drain_fpga->name];
        }
    } 
    return total_hop;  
    //cout<<"hop length: "<<total_hop<<endl;
}
void DB::force_place(intg c,intg f){
    c--;f--;
    fpga.v[f]->add_circuit(circuit.v[c]);
    circuit.v[c]->add_fpga(fpga.v[f]);
    cout<<"force place "<<c+1<<" to "<<f+1<<endl;
}
void DB::force_remove(intg c,intg f){
    c--;f--;
    fpga.v[f]->remove_circuit(circuit.v[c]);
    circuit.v[c]->remove_fpga(fpga.v[f]);
    cout<<"force remove "<<c+1<<" to "<<f+1<<endl;
}
void DB::force_exchange(intg a,intg b){
    a--;b--;
    intg a_f=circuit.v[a]->fpga_node->name+1;
     intg b_f=circuit.v[b]->fpga_node->name+1;
    force_remove(a+1,a_f);
    force_remove(b+1,b_f);
    force_place(a+1,b_f);
    force_place(b+1,a_f);  
}
void DB::output(fstream &out) {
    // force_exchange(11,2);
    // force_place(4,3);
    // force_place(4,1);
    // force_place(3,2);
    // auto fpga_1=hop_fpga.v[0];
    // auto node_copy_6=circuit.v[5];
    // auto node_copy_12=circuit.v[11];
    // fpga_1->remove_copy_circuit(node_copy_6);
    // node_copy_6->remove_copy_fpga(fpga_1);
    // fpga_1->add_copy_circuit(node_copy_12);
    // node_copy_12->add_copy_fpga(fpga_1);

    stringstream ss;
    std::set<intg> all_set;
    std::multiset<intg> copy_set;
    for(auto cur_fpga:fpga.get_all_vertex()){
        ss<<"FPGA"<<cur_fpga->name+1<<": ";
        for(auto node:cur_fpga->placed_circuit){
            ss<<"g"<<node->name+1<<" ";
            #ifdef LOG
            all_set.insert(node->name);
            #endif
        }
        for(auto node:cur_fpga->copy_circuit){
            ss<<"g"<<node->name+1<<"* ";
            #ifdef LOG
            copy_set.insert(node->name);
            #endif
        }
        ss<<endl;
    }
    out<<ss.rdbuf();

    #ifdef LOG
    debug_resouce();
    if (all_set.size() == circuit.num_vertex){
        intg cur_hop=output_loss();
        cout<<"Total Hop Length = "<<cur_hop<<endl;
        cout<<"score = "<<(int)(cur_hop*(1.0f+0.2f*((double)timer.get_time()/3600)))<<endl;
        cout <<"copy_num:"<<copy_set.size()<<endl;
        return;
    }
    cout << "=========================================================" << endl;
    cout << "[GY] Total Not Set: " << all_set.size() << "/" << circuit.num_vertex
         <<"noassigned number "<<circuit.num_vertex-all_set.size()<<" copy_num:"<<copy_set.size()<<endl;
        for(auto net:nets){    
            FPGANode *src_fpga=net->source->fpga_node;
            if(src_fpga==nullptr)continue;
            for(auto drain_fpga:net->used_fpga_node){
                if(src_fpga!=drain_fpga){
                    if(fpga_weight[src_fpga->name][drain_fpga->name]==-1){
                       
                        for(auto cell:net->net_cell)
                            if(cell->fpga_node==drain_fpga)
                                cout<<"\nnet's src:"<<net->source->name+1<<" "<<cell->name+1<<"topo Wrong"<<endl;
                    }
                }
            }
            for(auto copy_circuit_:net->copy_circuit){
                for(auto copy_fpga_:copy_circuit_->copy_position){
                    if(fpga_weight[copy_fpga_->name][src_fpga->name]==-1)
                        cout<<"\nnet's src:"<<net->source->name+1<<" copy node:"<<copy_circuit_->name+1<<" topo Wrong"<<endl;
                }
            }
        } 
    cout << "\n=========================================================" << endl;
    #endif
}
void DB::output() {
    // force_exchange(11,2);
    // force_place(4,3);
    // force_place(4,1);
    std::set<intg> all_set;
    std::multiset<intg> copy_set;
    for(auto cur_fpga:fpga.get_all_vertex()){
        for(auto node:cur_fpga->placed_circuit){
            all_set.insert(node->name);
        }
        for(auto node:cur_fpga->copy_circuit){
            copy_set.insert(node->name);
        }
    }


    debug_resouce();
    if (all_set.size() == circuit.num_vertex){
        cout<<"Total Hop Length = :"<<output_loss()<<endl;
        return;
    }
    cout << "=========================================================" << endl;
    cout << "[GY] Total Not Set: " << all_set.size() << "/" << circuit.num_vertex
         <<"noassigned number "<<circuit.num_vertex-all_set.size()<<" copy_num:"<<copy_set.size()<<endl;
        for(auto net:nets){     
            FPGANode *src_fpga=net->source->fpga_node;
            if(src_fpga==nullptr)continue;
            for(auto drain_fpga:net->used_fpga_node){
                if(src_fpga!=drain_fpga){
                    if(fpga_weight[src_fpga->name][drain_fpga->name]==-1){
                       
                        for(auto cell:net->net_cell)
                            if(cell->fpga_node==drain_fpga)
                                cout<<"\nnet's src:"<<net->source->name+1<<" "<<cell->name+1<<"topo Wrong"<<endl;
                    }
                }
            }
            for(auto copy_circuit_:net->copy_circuit){
                for(auto copy_fpga_:copy_circuit_->copy_position){
                    if(fpga_weight[copy_fpga_->name][src_fpga->name]==-1)
                        cout<<"\nnet's src:"<<net->source->name+1<<" copy node:"<<copy_circuit_->name+1<<" topo Wrong"<<endl;
                }
            }
        } 
    cout << "\n=========================================================" << endl;
}
void DB::debug_resouce(){
    cout << "====================resource==============================" << endl;
    for(auto fpga_:fpga.v){
        cout<<"Resources of FPGA"<<fpga_->name+1<<": [ ";
        for(int i=0;i<8;i++)
            cout<<fpga_->usage[i]<<" ";
        cout<<"]\n";
    }
    for(auto fpga_:fpga.v){
        cout<<"Total extern cut of FPGA"<<fpga_->name+1<<": "<<fpga_->usage_wire;
        if(fpga_->usage_wire>fpga_->wire_capacity)
            cout<<" wire overflow\n";
        else
            cout<<endl;
    }
    cout << "=========================================================" << endl;

}
//杈撳嚭鍥綠'鐨勪俊锟�?
template<class T,class I>
void print_G(
            unordered_map<T *, intg> &v_map,
            vector<T *> &vv,
            vector<vector<T *>> &gg,
            vector<vector<I>> &wweight,
            std::string s){
    #ifdef LOG_DB
    cout << "=============================================\n" ;
    cout << "The weigth of "<<s<<endl ;
    for(int i=0;i<gg.size();i++){
        cout<<"node :"<<i+1<<" ";
        for(int j=0;j<gg[i].size();j++)
        {
            cout<<wweight[i][j]<<" ";
        }
        cout<<endl;
    }
    cout << "=============================================\n" ;
    cout << "The neibor of "<<s<<endl ;
    for(int i=0;i<gg.size();i++){
        cout<<"node :"<<i+1<<" ";
        for(int j=0;j<gg[i].size();j++)
        {
            cout<<v_map[gg[i][j]]+1<<" ";
        }
        cout<<endl;
    }
    cout << "=============================================\n" ;
    #endif
}

// SPFA绠楁硶鎵惧埌浠巖oot_index鍒版墍鏈夐《鐐圭殑鏈€闀胯窛锟�?
template<class T,class I,class D>
void SPFA(
        const unordered_map<T *, I> &vv_map,
        const vector<T*> &vv,
          std::vector<D> &dist,
          std::vector<intg> &parents,
        const std::vector<std::vector<T*>> &gg,
          const std::vector<std::vector<D>> &wweight,
          intg root_index) {
    D eps=(D)1e-6;
    std::queue<intg> q;
    std::vector<bool> inQueue(gg.size(), false);      //姣忎釜鍏冪礌锟�?锟斤拷鍦ㄩ槦鍒椾腑
    dist[root_index] = (D)0;
    q.push(root_index);
    inQueue[root_index] = true;

    while (!q.empty()) {
        intg u = q.front();
        q.pop();
        inQueue[u] = false;

        for (size_t i = 0; i < gg[u].size(); ++i) {
            T *v = gg[u][i];
            intg v_index=vv_map.at(v);
            if (v_index != u) { // 閬垮厤锟�?锟斤拷
                D weight = wweight[u][i];
                if (dist[u] + weight-dist[v_index] < eps) {
                    dist[v_index] = dist[u] + weight;
                    parents[v_index] = u;
                    if (!inQueue[v_index]) {
                        q.push(v_index);
                        inQueue[v_index] = true;
                    }
                }
            }
        }
    }
}
void DB::debugnet(){
    for(auto net:nets){
        auto fpga_list=fpga.v;
        cout<<"src "<<net->source->name+1<<" weight "<<net->weight
        <<" netsize "<<net->net_cell.size()<<endl;
        cout<<"the situation used fpga:\n";
        for(auto fpga_:net->used_fpga_node_count){
            cout<<fpga_list[fpga_.first]->name+1<<" "<<fpga_.second<<endl;
        }
    }
}


void DB::refine_under_copy(){
    #ifdef LOG
    cout << "Start refine." << endl;
    #endif
    // Refine with move-based. We try to move the boundary node, which is the
    // node that in the FPGA node that different with other node in the same
    // net.
    intg total_dec = 0;
    intg total_topo_vio = 0;
    intg debug_increment=0;
    intg best_hop=std::numeric_limits<intg>::max();
    intg last_refine_cnt=1;
    for (int i = 0; i < 1; ++i) {
        //鍏堣幏寰楁弧瓒充氦鎹㈡潯浠剁殑node锟�?
        #ifdef LOG
        cout<<"------------------------"<<i+1<<"--------------------------\n";
        #endif
        if (timer.timeout())
            break;
        // Opt
        intg tmp_dec = std::numeric_limits<intg>::max();
        FPGANode *refine_f = nullptr;

        //瀛樺偍锟�?锟斤拷绉诲姩鐨勭偣c
        vector<CircuitNode *> re_q;
        unordered_map<CircuitNode *,set<FPGANode *>> cddt_refine_f;
        for (auto &c : circuit.get_all_vertex()) {
            if(c->fpga_node==nullptr)   continue;
            Tensor<intg> cddt_move(hop_fpga.num_vertex,0);
            if (try_move_undercopy(c,cddt_move))
                continue;
            re_q.emplace_back(c);
            for(int i=0;i<hop_fpga.num_vertex;++i){
                if(cddt_move.at(i)==1)
                    cddt_refine_f[c].emplace(hop_fpga.v[i]);
            }
        }
        sort(re_q.begin(), re_q.end(),
             [&](const CircuitNode *lhs, const CircuitNode *rhs) {
                //鍙備笌鐨刵et灏戠殑鍏堢Щ锟�?锟斤拷鍏讹拷?锟�?锟斤拷灞呭皯鐨勫厛绉诲姩
                 const auto &lns = lhs->nets.size();
                 const auto &rns = rhs->nets.size();
                 const auto &lgs = circuit.g_set[lhs->name].size();
                 const auto &rgs = circuit.g_set[rhs->name].size();
                 if (lns == rns) {
                     return lgs < rgs;
                 }
                 return lns < rns;
             });
        unordered_map<intg,pair<intg,intg>> query_profit_best;   //circuit if_f_fpga profit
        auto temp_place_cmp=[&](intg lhs,intg rhs){
            //璋佺殑鏀剁泭楂樼敤锟�?
            auto l_ans=query_profit_best[lhs].second;
            auto r_ans=query_profit_best[rhs].second;
            
            if(l_ans==r_ans){
                //鍒ゆ柇鍓╀綑绌洪棿
                if(fpga_neighbor_free_space[lhs]==fpga_neighbor_free_space[rhs])
                    return lhs<rhs;
                return fpga_neighbor_free_space[lhs]>fpga_neighbor_free_space[rhs];
            }
            return l_ans>r_ans;
        };
        vector<set<intg,decltype(temp_place_cmp)>> temp_place_sets(hop_fpga.num_vertex,set<intg,decltype(temp_place_cmp)>(temp_place_cmp));
        // set<intg> sets[re_q.size()*hop_fpga.num_vertex];
        intg sets_cnt=0;
        intg refine_cnt=0;
        intg wait_place_cnt=0;
        intg exchange_cnt=0;
        debug_increment=0;
        for (auto &c : re_q) {
            //鍒ゆ柇锟�?锟斤拷锟�?锟斤拷缁х画澶勭悊
            //if(c->fpga_node==nullptr){continue;}
            Tensor<intg> cddt_move(hop_fpga.num_vertex,0);
            if(try_move_undercopy(c,cddt_move))   continue;
            cddt_refine_f[c].clear();
            for(int i=0;i<hop_fpga.num_vertex;++i){
                if(cddt_move.at(i)==1)
                    cddt_refine_f[c].emplace(hop_fpga.v[i]);
            }
            if(cddt_refine_f[c].size()<=0)  continue;
            ////蹇界暐璧勬簮鍜岀嚎鐨勯檺鍒讹拷?鎵炬敹鐩婃洿澶х殑FPGA
            tmp_dec = std::numeric_limits<intg>::min();
            refine_f = nullptr;
            auto removed_fpga=c->fpga_node;
            removed_fpga->remove_circuit_undercopy(c);
            c->remove_fpga(removed_fpga);
            intg remove_origin_incremet=try_move_increment_undercopy(c,removed_fpga);
            for (auto &nf : cddt_refine_f[c]) {
                if (nf == removed_fpga)
                    continue;
                // //td姹傚緱鏄痗绉诲姩鍒皀f鍚庣殑鏀剁泭 
                intg td =remove_origin_incremet-try_move_increment_undercopy(c, nf);
                // intg td = c->try_move(nf);
                if (tmp_dec < td) {
                    // if(valid_undercopy(nf,c)>0){
                    //     tmp_dec = td;
                    //     refine_f = nf;
                    // }
                        tmp_dec = td;
                        refine_f = nf;
                } 
                else if (tmp_dec == td && refine_f != nullptr) {
                    //濡傛灉鏀剁泭绛変簬鏈€灏忔敹鐩婏紝閭ｅ氨姣斿叾浠栦笢锟�?
                    // if (valid_undercopy(nf,c)&&hop_fpga.g_set[nf->name].size() <
                    //     hop_fpga.g_set[refine_f->name].size()) {
                    //     tmp_dec = td;
                    //     refine_f = nf;
                    // }
                    if (hop_fpga.g_set[nf->name].size() <
                        hop_fpga.g_set[refine_f->name].size()) {
                        tmp_dec = td;
                        refine_f = nf;
                    }
                }
            }//for (auto &nf : cddt_refine_f[c])
            if (refine_f == nullptr){
                removed_fpga->add_circuit_undercopy(c);
                c->add_fpga(removed_fpga);
                continue;
            }
            intg all_profit=tmp_dec;
            if (all_profit > 0) {
                    if(valid_undercopy(refine_f,c)){
                        refine_f->add_circuit_undercopy(c);
                        c->add_fpga(refine_f);
                        refine_cnt++;
                        debug_increment+=all_profit;
                        continue;
                    }else{
                    //鏃犳晥锟�?锟斤拷瀛樿捣锟�?,绛夊緟浜ゆ崲
                        if(query_profit_best.count(c->name)<=0){
                            temp_place_sets[refine_f->name].emplace(c->name);
                            query_profit_best[c->name]={refine_f->name,tmp_dec};     
                            calculate_fpga_neighbor_free_space(refine_f->name,c);
                        }
                        else if(query_profit_best[c->name].second<tmp_dec){
                            auto origin_refine_f=query_profit_best[c->name].first;
                            temp_place_sets[origin_refine_f].erase(c->name);
                            query_profit_best[c->name]={refine_f->name,tmp_dec};
                            calculate_fpga_neighbor_free_space(refine_f->name,c);
                        }
                        removed_fpga->add_circuit_undercopy(c);
                        c->add_fpga(removed_fpga);
                    }
            }//if (all_profit > 0)
            else{
                //鏀剁泭浣庯紝浣嗘槸浜ゆ崲鏀剁泭锟�?
                //瀵绘壘涓€锟�?锟斤拷锟�?emoved_fpga鐨勮妭鐐筺ode_want_come锛岋拷?褰撳墠鑺傜偣c鍘籸efine_c
                // if(valid_undercopy(refine_f,c)){
                //     for(auto want_come_here:temp_place_sets[removed_fpga->name]){
                //         // auto want_come_here=*temp_place_sets[removed_fpga->name].begin();
                //         auto node_want_come=circuit.v[want_come_here];
                //         auto node_want_come_fpga=node_want_come->fpga_node;
                //         if(query_profit_best[want_come_here].second+all_profit>0){
                //             node_want_come_fpga->remove_circuit_undercopy(node_want_come);
                //             node_want_come->remove_fpga(node_want_come_fpga);
                //             if(valid_undercopy(removed_fpga,node_want_come)){
                //                 removed_fpga->add_circuit_undercopy(node_want_come);
                //                 node_want_come->add_fpga(removed_fpga);
                //                 temp_place_sets[removed_fpga->name].erase(want_come_here);
                //                 refine_f->add_circuit_undercopy(c);
                //                 c->add_fpga(refine_f);
                //                 exchange_cnt++;
                //                 break;
                //             }else{
                //                 node_want_come_fpga->add_circuit_undercopy(node_want_come);
                //                 node_want_come->add_fpga(node_want_come_fpga);  
                //             }//if(valid_undercopy(removed_fpga,node_want_come))
                //         }else{
                //             //鍚庣画鐨勭偣閮戒笉锟�?锟斤拷鏀剁泭鏇村ぇ锟�?
                //             break;
                //         }//if(query_profit_best[want_come_here].second+all_profit>0)
                //     }
                // }
                if(c->fpga_node==nullptr){
                    removed_fpga->add_circuit_undercopy(c);
                    c->add_fpga(removed_fpga);
                }
            }
        }//for (auto &c : re_q)
        #ifdef LOG
        cout<<"refine cnt:"<<refine_cnt<<endl;
        cout<<"increment:"<<debug_increment<<endl;
        cout<<"exchange_cnt:"<<exchange_cnt<<endl;
        int cur_hop;
        if(circuit_defer.size()==0){
            cur_hop=output_loss();
            cout<<"cur_hop_length:"<<cur_hop<<endl;
        }
        cout<<"--------------------------------------------------\n";
        #endif
        if((double)refine_cnt/last_refine_cnt<0.1)
            break;
        else
            last_refine_cnt=refine_cnt;
    }//for (int i = 0; i < 5; ++i)
        
    // if (timer.timeout())
    //     return;

#ifdef LOG_DB
    output_loss();
#endif
}


intg DB::try_move_increment(CircuitNode *c,FPGANode *f){
    intg dec = 0;
    auto origin_fpga=c->fpga_node;
    origin_fpga->remove_circuit(c);
    c->remove_fpga(origin_fpga);
    c->calculate_cut_increment(f);
    dec-=(c->cut_increment_map[f->name]);
    c->calculate_cut_increment(origin_fpga);
    dec+=(c->cut_increment_map[origin_fpga->name]);
    origin_fpga->add_circuit(c);
    c->add_fpga(origin_fpga);
    return dec;
}
void DB::refine() {
    #ifdef LOG
    cout << "Start refine." << endl;
    #endif  
    // Refine with move-based. We try to move the boundary node, which is the
    // node that in the FPGA node that different with other node in the same
    // net.
    intg total_dec = 0;
    intg total_topo_vio = 0;
    intg best_hop=std::numeric_limits<intg>::max();
    intg last_refine_cnt=1;
    for (int i = 0; i < 1; ++i) {
        //鍏堣幏寰楁弧瓒充氦鎹㈡潯浠剁殑node锟�?
        #ifdef LOG
        cout<<"------------------------"<<i+1<<"--------------------------\n";
        #endif
        if (timer.timeout())
            break;
        // Opt
        intg tmp_dec = std::numeric_limits<intg>::max();
        FPGANode *refine_f = nullptr;

        //瀛樺偍锟�?锟斤拷绉诲姩鐨勭偣c
        vector<CircuitNode *> re_q;
        unordered_map<CircuitNode *,set<FPGANode *>> cddt_refine_f;
        for (auto &c : circuit.get_all_vertex()) {
            if(c->fpga_node==nullptr)   continue;
            Tensor<intg> cddt_move(hop_fpga.num_vertex,1);
            cddt_move.at_is(c->fpga_node->name,0);
            if (try_move(c,cddt_move))
                continue;
            re_q.emplace_back(c);
            for(int i=0;i<hop_fpga.num_vertex;++i){
                if(cddt_move.at(i)==1)
                    cddt_refine_f[c].emplace(hop_fpga.v[i]);
            }
        }
        sort(re_q.begin(), re_q.end(),
             [&](const CircuitNode *lhs, const CircuitNode *rhs) {
                //鍙備笌鐨刵et灏戠殑鍏堢Щ锟�?锟斤拷鍏讹拷?锟�?锟斤拷灞呭皯鐨勫厛绉诲姩
                 const auto &lns = lhs->nets.size();
                 const auto &rns = rhs->nets.size();
                 const auto &lgs = circuit.g_set[lhs->name].size();
                 const auto &rgs = circuit.g_set[rhs->name].size();
                 if (lns == rns) {
                     return lgs < rgs;
                 }
                 return lns < rns;
             });
        unordered_map<intg,pair<intg,intg>> query_profit_best;   //circuit if_f_fpga profit
        auto temp_place_cmp=[&](intg lhs,intg rhs){
            //璋佺殑鏀剁泭楂樼敤锟�?
            auto l_ans=query_profit_best[lhs].second;
            auto r_ans=query_profit_best[rhs].second;
            
            if(l_ans==r_ans){
                if(fpga_neighbor_free_space[lhs]==fpga_neighbor_free_space[rhs])
                    return lhs<rhs;
                return fpga_neighbor_free_space[lhs]>fpga_neighbor_free_space[rhs];
            }
            return l_ans>r_ans;
        };
        vector<set<intg,decltype(temp_place_cmp)>> temp_place_sets(hop_fpga.num_vertex,set<intg,decltype(temp_place_cmp)>(temp_place_cmp));
        // set<intg> sets[re_q.size()*hop_fpga.num_vertex];
        intg sets_cnt=0;
        intg refine_cnt=0;
        intg wait_place_cnt=0;
        for (auto &c : re_q) {
            //鍒ゆ柇锟�?锟斤拷锟�?锟斤拷缁х画澶勭悊
            //if(c->fpga_node==nullptr){continue;}
            Tensor<intg> cddt_move(hop_fpga.num_vertex,1);
            cddt_move.at_is(c->fpga_node->name,0);
            if(try_move(c,cddt_move))   continue;
            cddt_refine_f[c].clear();
            for(int i=0;i<hop_fpga.num_vertex;++i){
                if(cddt_move.at(i)==1)
                    cddt_refine_f[c].emplace(hop_fpga.v[i]);
            }
            if(cddt_refine_f[c].size()<=0)  continue;
            ////蹇界暐璧勬簮鍜岀嚎鐨勯檺鍒讹拷?鎵炬敹鐩婃洿澶х殑FPGA
            tmp_dec = std::numeric_limits<intg>::min();
            // tmp_dec = std::numeric_limits<intg>::max();
            refine_f = nullptr;
            for (auto &nf : cddt_refine_f[c]) {
                if (nf == c->fpga_node)
                    continue;
                // //td姹傚緱鏄痗绉诲姩鍒皀f鍚庣殑鏀剁泭 
                intg td = try_move_increment(c, nf);
                // intg td = c->try_move(nf);
                if (tmp_dec < td) {
                    tmp_dec = td;
                    refine_f = nf;
                } 
                else if (tmp_dec == td && refine_f != nullptr) {
                    //濡傛灉鏀剁泭绛変簬鏈€灏忔敹鐩婏紝閭ｅ氨姣斿叾浠栦笢锟�?
                    if (hop_fpga.g_set[nf->name].size() <
                        hop_fpga.g_set[refine_f->name].size()) {
                        tmp_dec = td;
                        refine_f = nf;
                    }
                }
            }//for (auto &nf : cddt_refine_f[c])
            if (refine_f == nullptr)
                continue;
            intg all_profit=tmp_dec;
            auto removed_fpga=c->fpga_node;
            if (all_profit > 0) {
                removed_fpga->remove_circuit(c);
                c->remove_fpga(removed_fpga);
                if(valid(refine_f,c)>0){
                    refine_f->add_circuit(c);//47696
                    c->add_fpga(refine_f);
                    refine_cnt++;
                    // //娑堥櫎瀹孋V锛岄偅灏卞皾璇曟斁锟�?锟斤拷鍓嶆病鏀剧疆杩囩殑
                    // if(temp_place_sets[removed_fpga->name].size()>0){
                    //     auto index_wait_node=*(temp_place_sets[removed_fpga->name].begin());
                    //     temp_place_sets[removed_fpga->name].erase(index_wait_node);
                    //     auto c_wait_node=circuit.v[index_wait_node];
                    //     auto fpga_wait_node=c_wait_node->fpga_node;
                    //     fpga_wait_node->remove_circuit(c_wait_node);
                    //     c_wait_node->remove_fpga(fpga_wait_node);
                    //     if(valid(removed_fpga,c_wait_node)>0){
                    //         //绉诲姩锟�?
                    //         removed_fpga->add_circuit(c_wait_node);
                    //         c_wait_node->add_fpga(removed_fpga);
                    //         wait_place_cnt++;
                            
                    //     }else{
                    //         fpga_wait_node->add_circuit(c_wait_node);
                    //         c_wait_node->add_fpga(fpga_wait_node);
                    //         temp_place_sets[removed_fpga->name].emplace(index_wait_node);
                    //     }
                    // }
                    continue;
                }//if(valid(refine_f,c)>0)
                else{
                    // if(query_profit_best.count(c->name)<=0){
                    //     temp_place_sets[refine_f->name].emplace(c->name);
                    //     query_profit_best[c->name]={refine_f->name,tmp_dec};     
                    //     calculate_fpga_neighbor_free_space(refine_f->name,c);
                    // }
                    // else if(query_profit_best[c->name].second<tmp_dec){
                    //     auto origin_refine_f=query_profit_best[c->name].first;
                    //     temp_place_sets[origin_refine_f].erase(c->name);
                    //     query_profit_best[c->name]={refine_f->name,tmp_dec};
                    //     calculate_fpga_neighbor_free_space(refine_f->name,c);
                    // }
                    removed_fpga->add_circuit(c);
                    c->add_fpga(removed_fpga);
                }
            }//if (all_profit > 0)
        }//for (auto &c : re_q)
        #ifdef LOG
        cout<<"refine cnt:"<<refine_cnt<<endl;
        int cur_hop;
        if(circuit_defer.size()==0){
            cur_hop=output_loss();
            cout<<"cur_hop_length:"<<cur_hop<<endl;
        }
        cout<<"--------------------------------------------------\n";
        #endif
        if((double)refine_cnt/last_refine_cnt<0.1)
            break;
        else
            last_refine_cnt=refine_cnt;
    }//for (int i = 0; i < 5; ++i)
        
    // if (timer.timeout())
    //     return;

#ifdef LOG_DB
    output_loss();
#endif
}

intg DB::exchange_valid(CircuitNode *origin,FPGANode *origin_fpga,CircuitNode *new_come,FPGANode *new_fpga){
    /*
        锟�?锟斤拷锟�?锟斤拷锟�?
        姣旓拷?瑕佹妸5鍙风Щ鍔ㄥ埌fpga2锛孋V锟�?2鍙锋斁锟�?锟斤拷fpga3
        濡傛灉锟�?5鍙锋病鏀剧疆鐨勬椂鍊欙拷?锟�?2鍙锋槸鍚︽湁鏁堝氨浼氬嚭锟�?
        搴旓拷?锟�?5鍙锋斁锟�?锟斤拷鍐嶏拷?锟�?2鍙锋槸鍚︽湁锟�?
    */
   //閽堬拷?CV鏀剧疆鍦烘櫙
    if(new_come->fpga_node==nullptr){
        // std::cout<<"1\n";
        origin_fpga->remove_circuit(origin);
        origin->remove_fpga(origin_fpga);
        
        //灏濊瘯鏀剧疆new_come鍒皁rigin_fpga 
        intg result_unstrict=valid(origin_fpga,new_come);
        if(result_unstrict<=0){
            //涓嶄弗鏍兼斁锟�?锟斤拷杩囦笉锟�?
            origin_fpga->add_circuit(origin);
            origin->add_fpga(origin_fpga);
            return 0;
        }
        //濡傛灉锟�?锟斤拷鏀剧疆鐨勮瘽
        origin_fpga->add_circuit(new_come);
        new_come->add_fpga(origin_fpga);
        result_unstrict=valid(new_fpga,origin);
        if(result_unstrict<=0){
            //鏀剧疆new_come鍒皀ew_fpga涔嬪悗鏃犳晥
            origin_fpga->remove_circuit(new_come);
            new_come->remove_fpga(origin_fpga);
            origin_fpga->add_circuit(origin);
            origin->add_fpga(origin_fpga);
            return 0;
        }      
        origin_fpga->remove_circuit(new_come);
        new_come->remove_fpga(origin_fpga);
        origin_fpga->add_circuit(origin);
        origin->add_fpga(origin_fpga);    
    }
    else{
        origin_fpga->remove_circuit(origin);
        origin->remove_fpga(origin_fpga);
        new_fpga->remove_circuit(new_come);
        new_come->remove_fpga(new_fpga);
        //灏濊瘯鏀剧疆new_come鍒皁rigin_fpga 
        intg result_unstrict=valid(origin_fpga,new_come);
        if(result_unstrict<=0){
            //涓嶄弗鏍兼斁锟�?锟斤拷杩囦笉锟�?,鎭拷?
            origin_fpga->add_circuit(origin);
            origin->add_fpga(origin_fpga);
            new_fpga->add_circuit(new_come);
            new_come->add_fpga(new_fpga);
            return 0;
        }
        //濡傛灉锟�?锟斤拷鏀剧疆鐨勮瘽
        origin_fpga->add_circuit(new_come);
        new_come->add_fpga(origin_fpga);
        result_unstrict=valid(new_fpga,origin);
        if(result_unstrict<=0){
            //鏀剧疆new_come鍒皀ew_fpga涔嬪悗鏃犳晥锛屾仮锟�?
            origin_fpga->remove_circuit(new_come);
            new_come->remove_fpga(origin_fpga);
            origin_fpga->add_circuit(origin);
            origin->add_fpga(origin_fpga);
            new_fpga->add_circuit(new_come);
            new_come->add_fpga(new_fpga);
            return 0;
        }      
        //鎭拷?
        origin_fpga->remove_circuit(new_come);
        new_come->remove_fpga(origin_fpga);
        origin_fpga->add_circuit(origin);
        origin->add_fpga(origin_fpga);
        new_fpga->add_circuit(new_come);
        new_come->add_fpga(new_fpga);     
    }
    return 1;
}

//杩樻病寰堝ソ鐨勶拷?鐞嗭拷?锟�?锟斤拷鐐规湭鏀剧疆鐨勬儏鍐垫€庝箞澶勭悊锟�?
intg DB::topo_vio(FPGANode *f, intg c) {
    //缁燂拷?c鏀惧湪f鍛ㄥ洿鐨勫悎娉曟€э紝-1闈炴硶
    //鏁板瓧浠ｈ〃c鍛ㄥ洿鏈夛拷?灏戞湭鏀剧疆鐨勮妭锟�?
    //濡傛灉c鍛ㄥ洿锟�?锟斤拷锟�?锟斤拷鑺傜偣寰堬拷?
    intg vio = 0;
    for (auto &neighbor : circuit.g_set[c]) {
        if (neighbor->fpga_node == nullptr) {
            //濡傛灉c鐨勯偦灞呬篃鏈変竴涓猣pga娌℃斁锟�?
            vio += 1;
            continue;
        }
        //濡傛灉c鐨勯偦灞呯殑fpga涓嶆槸f鎴栬€呬笉鏄痜鐨勯偦灞呴偅涔堣偗瀹氭槸vio+1
        if (hop_fpga.g_set[f->name].count(
                hop_fpga.get_vertex(neighbor->fpga_node->name)) <= 0 &&
            f != neighbor->fpga_node)
            return -1;
    }
    return vio;
}

void DB::get_bfs_order(){
    vector<bool> done;
    done.resize(circuit.num_vertex, false);
    bfs_order.resize(circuit.num_vertex,-1);
    // node id, depth
    queue<pair<intg, intg>> q;
    for(int i=0;i<fixed_nodes.size();i++){
        q.push(make_pair(fixed_nodes[i].first, 0));
        done[fixed_nodes[i].first]=true;
    }
    while (!q.empty()) {
        auto top = q.front();
        auto &id = top.first;
        auto &depth = top.second;
        q.pop();
        bfs_order[id]=depth;
        for (auto &neighbor : circuit.g[id]) {
            auto neighbor_id = neighbor->name;
            if (done[neighbor_id] == false) {
                q.push(make_pair(neighbor_id, depth + 1));
                done[neighbor_id] = true;
            }
        }
    }
}

void DB::refine_copy(){
    #ifdef LOG
    cout << "Start copy refine." << endl;
    #endif
    // Refine with move-based. We try to move the boundary node, which is the
    // node that in the FPGA node that different with other node in the same
    // net.
    unordered_set<CircuitNode *> copy_node;
    intg total_dec = 0;
    intg total_topo_vio = 0;
    intg best_hop=std::numeric_limits<intg>::max();
    intg last_copy_cnt=1;
    intg copy_hop_increment;
    for (int i = 0; i < 1; ++i) {
        //鍏堣幏寰楁弧瓒充氦鎹㈡潯浠剁殑node锟�?
        #ifdef LOG
        cout<<"------------------------"<<i+1<<"--------------------------\n";
        #endif
        if (timer.timeout())
            break;
        // Opt
        intg tmp_dec = std::numeric_limits<intg>::max();
        FPGANode *copy_f = nullptr;
        copy_hop_increment=0;
        //瀛樺偍锟�?锟斤拷绉诲姩鐨勭偣c
        vector<CircuitNode *> re_q;
        unordered_map<CircuitNode *,set<FPGANode *>> cddt_copy_f;
        for (auto &c : circuit.get_all_vertex()) {
            if(c->fpga_node==nullptr||no_indegree_nodes.count(c))   continue;
            Tensor<intg> cddt_move(hop_fpga.num_vertex,0);
            if(!try_copy(c,cddt_move))
                continue;
            re_q.emplace_back(c);
            for(int i=0;i<hop_fpga.num_vertex;++i){
                if(cddt_move.at(i)==1)
                    cddt_copy_f[c].emplace(hop_fpga.v[i]);
            }
        }
        sort(re_q.begin(), re_q.end(),
             [&](const CircuitNode *lhs, const CircuitNode *rhs) {
                //鍙備笌鐨刵et灏戠殑鍏堢Щ锟�?锟斤拷鍏讹拷?锟�?锟斤拷灞呭皯鐨勫厛绉诲姩
                 const auto &lns = lhs->nets.size();
                 const auto &rns = rhs->nets.size();
                 const auto &lgs = circuit.g_set[lhs->name].size();
                 const auto &rgs = circuit.g_set[rhs->name].size();
                 if (lns == rns) {
                     return lgs < rgs;
                 }
                 return lns < rns;
             });
        unordered_map<intg,pair<intg,intg>> query_profit_best;   //circuit if_f_fpga profit
        auto temp_place_cmp=[&](intg lhs,intg rhs){
            //璋佺殑鏀剁泭楂樼敤锟�?
            auto l_ans=query_profit_best[lhs].second;
            auto r_ans=query_profit_best[rhs].second;
            
            if(l_ans==r_ans){
                if(fpga_neighbor_free_space[lhs]==fpga_neighbor_free_space[rhs])
                    return lhs<rhs;
                return fpga_neighbor_free_space[lhs]>fpga_neighbor_free_space[rhs];
            }
            return l_ans>r_ans;
        };
        vector<set<intg,decltype(temp_place_cmp)>> temp_place_sets(hop_fpga.num_vertex,set<intg,decltype(temp_place_cmp)>(temp_place_cmp));
        // set<intg> sets[re_q.size()*hop_fpga.num_vertex];
        intg sets_cnt=0;
        intg copy_cnt=0;
        intg wait_place_cnt=0;
        for (auto &c : re_q) {
            //鍒ゆ柇锟�?锟斤拷锟�?锟斤拷缁х画澶勭悊
            //if(c->fpga_node==nullptr){continue;}
            Tensor<intg> cddt_move(hop_fpga.num_vertex,0);
            if(!try_copy(c,cddt_move))   continue;
            cddt_copy_f[c].clear();
            for(int i=0;i<hop_fpga.num_vertex;++i){
                if(cddt_move.at(i)==1)
                    cddt_copy_f[c].emplace(hop_fpga.v[i]);
            }
            if(cddt_copy_f[c].size()<=0)  continue;
            ////蹇界暐璧勬簮鍜岀嚎鐨勯檺鍒讹拷?鎵炬敹鐩婃洿澶х殑FPGA
            tmp_dec = std::numeric_limits<intg>::min();
            // tmp_dec = std::numeric_limits<intg>::max();
            copy_f = nullptr;
            for (auto &nf : cddt_copy_f[c]) {
                if (nf == c->fpga_node)
                    continue;
                // //td姹傚緱鏄痗澶嶅埗鍒皀f鍚庣殑鏀剁泭 
                c->calculate_copy_increment(nf);
                intg td = -c->copy_increment_map[nf->name];
                // intg td = c->try_move(nf);
                if (tmp_dec < td) {
                    tmp_dec = td;
                    copy_f = nf;
                } 
                else if (tmp_dec == td && copy_f != nullptr) {
                    //濡傛灉鏀剁泭绛変簬鏈€灏忔敹鐩婏紝閭ｅ氨姣斿叾浠栦笢锟�?
                    if (hop_fpga.g_set[nf->name].size() <
                        hop_fpga.g_set[copy_f->name].size()) {
                        tmp_dec = td;
                        copy_f = nf;
                    }
                }
            }//for (auto &nf : cddt_refine_f[c])
            if (copy_f == nullptr)
                continue;
            intg all_profit=tmp_dec;
            if (all_profit > 0) {
                if(copy_valid(copy_f,c)>0){
                    copy_f->add_copy_circuit(c);
                    c->add_copy_fpga(copy_f);
                    add_source(c); 
                    copy_cnt++;
                    copy_node.emplace(c);
                    copy_hop_increment+=all_profit;
                    continue;
                }//if(valid(refine_f,c)>0)
                else{
                    // if(query_profit_best.count(c->name)<=0){
                    //     temp_place_sets[refine_f->name].emplace(c->name);
                    //     query_profit_best[c->name]={refine_f->name,tmp_dec};     
                    //     calculate_fpga_neighbor_free_space(refine_f->name,c);
                    // }
                    // else if(query_profit_best[c->name].second<tmp_dec){
                    //     auto origin_refine_f=query_profit_best[c->name].first;
                    //     temp_place_sets[origin_refine_f].erase(c->name);
                    //     query_profit_best[c->name]={refine_f->name,tmp_dec};
                    //     calculate_fpga_neighbor_free_space(refine_f->name,c);
                    // }
                }
            }//if (all_profit > 0)
        }//for (auto &c : re_q)
        #ifdef LOG
        cout<<"copy cnt:"<<copy_cnt<<endl;
        cout<<"copy increment:"<<copy_hop_increment<<endl;
        int cur_hop;
        if(circuit_defer.size()==0){
            cur_hop=output_loss();
            cout<<"cur_hop_length:"<<cur_hop<<endl;
        }
        cout<<"--------------------------------------------------\n";
        #endif
        if((double)copy_cnt/last_copy_cnt<0.1)
            break;
        else
            last_copy_cnt=copy_cnt;
    }//for (int i = 0; i < 5; ++i)
        
#ifdef LOG_DB
    output_loss();
#endif
}

}  // namespace topart
