#include "node.h"

namespace topart {

//用于FPGA之间判断各自的关系,是个矩阵，在bd.cpp的build_hop_fpga中初始化
unordered_map<FPGANode *, intg> fpga_to_index;
vector<FPGANode*> index_to_fpga;       
vector<vector<intg>> fpga_weight;      //合法的有距离，不合法只有-1
vector<vector<intg>> fpga_parents;     //fpga[i]是以i为起始点,所有点的父亲

void MergeNode::remove_fpga_circuit(FPGANode *fn){
    for(auto c_:belong_nodes){
        fn->remove_circuit_undercopy(c_);
        c_->remove_fpga(fn);
    }
}
void MergeNode::add_fpga_circuit(FPGANode *fn){
    for(auto c_:belong_nodes){
        fn->add_circuit_undercopy(c_);
        c_->add_fpga(fn);
    }
}

void MergeNode::calculate_copy_increment(FPGANode *if_f){

}

void MergeNode::calculate_cut_increment_undercopy(FPGANode *if_f){
    
}

void CircuitNode::set_fixed(FPGANode *fn) {
    fixed = true;
    fpga_node = fn;

    // update Net
    for (auto &n : nets) {
        n->add_fpga(fn);
    }

}

void CircuitNode::add_fpga(FPGANode *fn) {
    fpga_node = fn;

    // update Net
    for (auto &n : nets) {
        n->add_fpga(fn);
    }

}

void CircuitNode::remove_fpga(FPGANode *fn) {

    //this->fpga_node=nullptr;
    for (auto &n : nets) {
        n->remove_fpga(this->fpga_node);
    }
    this->fpga_node=nullptr;

}
void CircuitNode::add_copy_fpga(FPGANode *fn){
    copy_position.emplace(fn);
    for (auto &n : nets) {
        n->add_copy_circuit(this);
        n->add_copy_fpga(fn);
    }
}
void CircuitNode::remove_copy_fpga(FPGANode *fn){
    copy_position.erase(fn);
    for (auto &n : nets) {
        n->remove_copy_circuit(this);
        n->remove_copy_fpga(fn);
    }
}

void CircuitNode::flush_tsr_to_cddt(vector<FPGANode *> &mapping) {
    cddt.clear();
    for (intg i = 0; i < tsr_cddt.v.size(); ++i) {
        if (tsr_cddt.at(i) > 0)
            cddt.emplace(mapping[i]);
        else {
            tsr_cddt.at(i) = 0;
        }
    }
}

void CircuitNode::flush_cddt_to_tsr() {
    tsr_cddt.clear();
    for (auto &f : cddt) {
        tsr_cddt.at(f->name) += 1;
    }
}
void CircuitNode::calculate_cut_increment(              //处理hop增长数
    FPGANode *if_f) {
    intg cut_increment = 0;

    for(auto &n:src_nets){      //对于当前节点作为src节点，对于每一个被驱动点都更新hop
        intg net_weight=n->weight;
        for(auto &fpga:n->used_fpga_node){
            intg cur_weight = fpga_weight[if_f->name][fpga->name];
            cut_increment+=cur_weight*net_weight;
        }
    }
    for(auto &n:drain_nets){    //node是被驱动节点，src未放置或if-f已经出现过就不再产生hop
        auto source=n->source;       
        if(source->fpga_node==nullptr||n->used_fpga_node.count(if_f)>0) continue;
        cut_increment+=n->weight*fpga_weight[if_f->name][source->fpga_node->name];
    }
    cut_increment_map[if_f->name] = cut_increment;
}

void CircuitNode::calculate_connect_increment(          //专门处理对外互联数目
    FPGANode *if_f) {
    //仅限于c点没有放置过
    std::array<intg,8> resource_incirement;
    unordered_map<FPGANode*,intg> wire_demand;
    if_f->get_demand(this,resource_incirement,wire_demand);

    intg total_connect_increment = 0;
    for(auto wire_:wire_demand){
        auto f_connect=wire_.first;
        auto num_connect=wire_.second;
        total_connect_increment+=num_connect;
    }
    connect_increment_map[if_f->name] = total_connect_increment;
}
void CircuitNode::calculate_copy_increment(FPGANode *f){//复制的hop
    intg copy_increment=0;
    //如果原始点和复制点在一起不影响对外互联
    if(this->fpga_node==f){
        copy_increment_map[f->name]=0;
        return ;
    }

    //复制Node到当前FPGA，要考虑Node的src到Node产生的对外互联数目的增加 
    for(auto drain_:this->drain_nets){
        auto src_node=drain_->source;
        auto src_fpga=src_node->fpga_node;
        auto src_copy_fpga=src_node->copy_position;
        if(src_fpga==nullptr)   continue;
        //这个情况下不会产生割
        if(src_fpga==f||src_copy_fpga.count(f))
            continue;
        //本来已经产生过割了
        if(drain_->used_fpga_node.count(f)||drain_->used_copy_fpga.count(f))
            continue;
        copy_increment+=drain_->weight*fpga_weight[f->name][src_fpga->name];
    }//for(auto drain_:c->drain_nets)
    
    //复制c点，对this之外的被驱动点是没有影响的
    for(auto src_:this->src_nets){
        auto src_fpga=this->fpga_node;
        auto src_copy_fpga=this->copy_position;
        if(!src_->used_all_fpga.count(f))
            continue;
        if(src_copy_fpga.count(f)||src_fpga==f)     
            continue;    
        copy_increment-=src_->weight*fpga_weight[f->name][src_fpga->name];
    }
    copy_increment_map[f->name]=copy_increment;
}
void CircuitNode::calculate_cut_increment_undercopy(             
    FPGANode *if_f) {
    //建立在当前点没有fpga_node的情况下
    intg cut_increment = 0;

    for(auto &n:src_nets){   
        //this放在if_f上，那就要计算被驱动点的fpga_node和copy_node和当前点的if_f以及copy_position的关系
        intg net_weight=n->weight;
        for(auto &fpga_:n->used_all_fpga){
            if(copy_position.count(fpga_)||fpga_==if_f)continue;
            #ifdef LOG
                assert(fpga_weight[if_f->name][fpga_->name]!=-1);
            #endif
            intg cur_weight = fpga_weight[if_f->name][fpga_->name];
            cut_increment+=cur_weight*net_weight;
        } 
    }

    for(auto &n:drain_nets){    
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
        #ifdef LOG
            assert(!copy_position.count(if_f));
        #endif
        //为啥加这一段会出问题呢？
        //求驱动copy_position的hop增加
        // for(auto fpga_:copy_position){
        //     if(!source_->copy_position.count(fpga_)&&fpga_!=source_->fpga_node){
        //         #ifdef LOG
        //         assert(fpga_!=if_f);
        //         assert(fpga_weight[source_->fpga_node->name][fpga_->name]!=-1);
        //         #endif
        //         cut_increment+=n->weight*fpga_weight[source_->fpga_node->name][fpga_->name];
        //     }
        // }// for(auto fpga_:copy_position)
    }
    cut_increment_undercopy_map[if_f->name] = cut_increment;
}
intg CircuitNode::try_move(FPGANode *f) {
    intg dec = 0;
        for (auto &n : nets) {
            intg cur_dec = n->estimate_increase_cut_size_refine(f,this->fpga_node,(n->source==this));
            if(cur_dec==std::numeric_limits<intg>::min())
                return std::numeric_limits<intg>::min();
            dec+=cur_dec;
        }
//     f->remove_circuit(this);
//     this->remove_fpga(f);
//     this->calculate_cut_increment(f);
//     dec=cut_increment_map[f->name];
//     f->add_circuit(this);
//     this->add_fpga(f);
    return dec;
}

//一定是先使用这个函数再使用add_fpga,remove_fpga
void FPGANode::add_circuit(CircuitNode* c) { 
    for (size_t i = 0; i < 8; ++i) {
        usage[i] = usage[i] + c->node_capacity[i];
    }
    //检验所有有关的超边负责处理所有节点的wire
    for(auto &net:c->src_nets){     
        //net里面还没有fpga，或者只有自身那肯定不产生割
        if(net->used_fpga_node.size()==0||net->used_fpga_node.size()==1&&net->used_fpga_node.count(this))
            continue;
        //剩下的情况下都会产生割
        usage_wire+=net->weight;
        for(auto fpga_:net->used_fpga_node){
            if(fpga_==this) continue;
            fpga_->usage_wire+=net->weight;
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
        if(source_f==nullptr||net->used_fpga_node.count(this))
            continue;
        //如果src放置那还新增了fpga那自身肯定要处理
        this->usage_wire+=net->weight;
        //如果只有src_fpga一个fpga，那src_fpga也处理
        if(net->used_fpga_node.size()==1)
            source_f->usage_wire+=net->weight;
    }   
    placed_circuit.emplace(c);
}
void FPGANode::add_circuit_undercopy(CircuitNode* c) {
    for (size_t i = 0; i < 8; ++i) {
        usage[i] = usage[i] + c->node_capacity[i];
    }
    //检验所有有关的超边负责处理所有节点的wire
    for(auto &net:c->src_nets){     
        //net里面还没有fpga，或者只有自身那就没有割
        if(net->used_all_fpga.size()==0||net->used_all_fpga.count(this)&&net->used_all_fpga.size()==1)
            continue;
        //挑选出不被copy所覆盖，也跟当前的FPGA无关的FPGA
        auto c_copy_fpga=c->copy_position;
        int result =0;
        for(auto fpga_:net->used_all_fpga){
            if(c_copy_fpga.count(fpga_)||fpga_==this)
                continue;
            else{
                fpga_->usage_wire+=net->weight;   
                result=1;
            }
        }
        if(result==1)
            usage_wire+=net->weight;
    }
    for(auto &net:c->drain_nets){ 
        auto source_c=net->source;
        auto source_f=source_c->fpga_node;
        auto source_copyf=source_c->copy_position;
        //如果在this的地方
        if(source_f==nullptr||net->used_all_fpga.count(this))
            continue;
        // if(source_f==this||source_copyf.count(this))    continue;
        this->usage_wire+=net->weight;
        
        //如果source_f的fpga都是src-fpga或者src-copy的，那就会产生新的割
        #ifdef LOG
        assert(net->used_all_fpga.size()>=(1+source_copyf.size()));
        #endif
        if(net->used_all_fpga.size()>1+source_copyf.size())continue;
        source_f->usage_wire+=net->weight;
    }   
    placed_circuit.emplace(c);  
}
void FPGANode::remove_circuit(CircuitNode* c) { 
    for (size_t i = 0; i < 8; ++i) {
        usage[i] = usage[i] - c->node_capacity[i];
    }
    //检验所有有关的超边负责处理所有节点的wire
    // 默认used_fpga_node至少有一个this
    for(auto &net:c->src_nets){     
        //如果移除的点作为net的src点
        //如果net只有一个fpga，那肯定是当前fpga，自然没有对外互联数
        //如果net有多个fpga，那所有fpga的对外互联数目都被复原
        if(net->used_fpga_node.size()==1)
            continue;
        usage_wire-=net->weight;
        for(auto fpga_:net->used_fpga_node)
            if(fpga_!=this)
                fpga_->usage_wire-=net->weight;
    }
    for(auto &net:c->drain_nets){ 
        auto source_c=net->source;
        auto source_f=source_c->fpga_node;
        //如果net只有一个fpga那就不用处理
        //如果this占了两次及以上，那就不用处理自己
        if(source_f==nullptr||net->used_fpga_node.size()==1)
            continue;
        //如果删去自己后导致fpga数目变化，那自身肯定要处理
        if(net->used_fpga_node_count[this->name]<=1)
            this->usage_wire-=net->weight;
        //如果删去自己后，刚好让割消失，那src也要变化
        if(net->used_fpga_node.size()==2&&net->used_fpga_node_count[this->name]==1)
            source_f->usage_wire-=net->weight;
    }   
    placed_circuit.erase(c);
}
void FPGANode::remove_circuit_undercopy(CircuitNode* c) {
    for (size_t i = 0; i < 8; ++i) {
        usage[i] = usage[i] - c->node_capacity[i];
    }
    //检验所有有关的超边负责处理所有节点的wire
    // 默认used_fpga_node至少有一个this
    for(auto &net:c->src_nets){     
        auto source_f=c->fpga_node;
        auto source_copyf=c->copy_position;
        auto c_copyf=c->copy_position;
        //如果没有产生过割就跳过,即result=0
        #ifdef LOG
            assert(source_f==this);
            assert(net->used_all_fpga.size()>=(1+source_copyf.size()));
        #endif
        if(net->used_all_fpga.size()==1+source_copyf.size())
            continue; 
        //如果产生过割，那么所有产生过割的fpga都要除去割
        bool result=0;
        for(auto fpga_:net->used_all_fpga){
            if(fpga_!=source_f&&!c_copyf.count(fpga_)){
                fpga_->usage_wire-=net->weight;  
                result=1;
            }
        }
        if(result==1)
            usage_wire-=net->weight;
    }
    for(auto &net:c->drain_nets){ 
        auto source_c=net->source;
        auto source_f=source_c->fpga_node;
        auto source_copyf=source_c->copy_position;
        //如果存在过就不处理
        if(source_f==nullptr||this==source_f||source_copyf.count(this))
            continue;
        //如果删去自己后导致fpga数目变化，那自身肯定要处理
        if(!net->used_copy_fpga.count(this)&&net->used_fpga_node_count[this->name]<=1){
            this->usage_wire-=net->weight;
            //如果除了source-f和source_copyf之外只有this，那src的割也要去除
            #ifdef LOG
            assert(!source_copyf.count(source_f));
            //assert(net->used_all_fpga.size()>=2+source_copyf.size());
            #endif
            if(net->used_all_fpga.size()>2+source_copyf.size())
                continue; 
            source_f->usage_wire-=net->weight;
        }
    }   
    placed_circuit.erase(c);
}
//先找个，再处理circuit
void FPGANode::add_copy_circuit(CircuitNode* c){
    //复制产生后，判断net被多少FPGA割的条件变了
    //是否产生割也变了
    //复制后对于drain_node，是否产生hop也要看复制的src_node是不是在一起了
    //对于src_node，产生的hop也要看对于复制的drain_node也会产生新的割
    for (size_t i = 0; i < 8; ++i) {
        usage[i] = usage[i] + c->node_capacity[i];
    }

    //如果原始点和复制点在一起不影响对外互联
    if(c->fpga_node==this)  return;

    //复制Node到当前FPGA，要考虑Node的src到Node产生的对外互联数目的增加 
    for(auto drain_:c->drain_nets){
        auto src_node=drain_->source;
        auto src_fpga=src_node->fpga_node;
        auto src_copy_fpga=src_node->copy_position;
        //如果放置在源点或者复制点所在的FPGA则跳过
        if(src_fpga==nullptr) continue;
        if(src_fpga==this||src_copy_fpga.count(this))
            continue;
        //本来已经产生过割了
        if(drain_->used_all_fpga.count(this))
            continue;
        // if(drain_->used_copy_fpga.count(this)||drain_->used_fpga_node.count(this))
        //     continue;
        usage_wire+=drain_->weight;
        //如果没产生过割(used_copy_fpga里的和used_fpga_node里面的都是src_fpga和src_copy_fpga内的)，那就给src加上割
        #ifdef LOG
        assert(drain_->used_all_fpga.size()>=src_copy_fpga.size()+1);
        #endif
        if(drain_->used_all_fpga.size()>src_copy_fpga.size()+1) continue;
        src_fpga->usage_wire+=drain_->weight;
    }//for(auto drain_:c->drain_nets)
    
    //复制c点，对this之外的被驱动点是没有影响的
    for(auto src_:c->src_nets){
        auto src_fpga=c->fpga_node;
        auto src_copy_fpga=c->copy_position;
        //如果当前fpga没有产生割那就跳过
        if(!src_->used_all_fpga.count(this))
            continue;
        // if(!src_->used_copy_fpga.count(this)&&!src_->used_fpga_node.count(this))
        //     continue;
        if(src_copy_fpga.count(this)||src_fpga==this)     
            continue;    
        //跟当前节点产生过割
        usage_wire-=src_->weight;
        //result=0说明没有没有和其他的FPGA有割，那就可以去除掉
        #ifdef LOG
        assert(src_->used_all_fpga.size()>=src_copy_fpga.size()+2);
        #endif
        if(src_->used_all_fpga.size()>src_copy_fpga.size()+2) continue;
        src_fpga->usage_wire-=src_->weight;
    }
    copy_circuit.emplace(c);
}
void FPGANode::remove_copy_circuit(CircuitNode* c){
    for (size_t i = 0; i < 8; ++i) {
        usage[i] = usage[i] - c->node_capacity[i];
    }
    for(auto drain_:c->drain_nets){
        auto src_node=drain_->source;
        auto src_fpga=src_node->fpga_node;
        auto &src_copy_fpga=src_node->copy_position;
        //判断this是否在当前net只有一个，如果有多个就跳过
        if(src_fpga==nullptr)   continue;
        // if(drain_->used_all_fpga.count(this))
        //     continue;
        if(drain_->used_fpga_node.count(this)||drain_->used_copy_fpga_count[this->name]>1)
            continue;
        if(src_fpga==this||src_copy_fpga.count(this))
            continue;
        usage_wire-=drain_->weight;  
        //如果只跟this有割那就减去src的割
        #ifdef LOG
        assert(drain_->used_all_fpga.size()>=src_copy_fpga.size()+2);
        #endif
        if(drain_->used_all_fpga.size()>src_copy_fpga.size()+2) continue;
        src_fpga->usage_wire-=drain_->weight;
    }//for(auto drain_:c->drain_nets)
    
    //对于当前节点驱动的节点,如果复制所在的this里存在被驱动点或者其复制点，那就加上割
    //如果这导致src_fpga里面仅仅只有this，那就加上src_fpga的割
    for(auto src_:c->src_nets){
        auto src_fpga=c->fpga_node;
        auto &src_copy_fpga=c->copy_position;
        //如果没有c的copy位置没有其他节点存在，那就跳过
        if(src_->used_fpga_node.count(this)==0&&src_->used_copy_fpga_count[this->name]==1)
            continue;
        usage_wire+=src_->weight;
        //如果src原来没有产生过割那就加上割
        //src原来没有产生过割指的是src_->used_all_fpga.size()==src_copy_fpga.size()+1
        #ifdef LOG
        assert(!src_copy_fpga.count(src_fpga));
        #endif
        if(src_->used_all_fpga.size()>src_copy_fpga.size()+1) continue;
        src_fpga->usage_wire+=src_->weight;
    }
    copy_circuit.erase(c);
}
// intg FPGANode::valid(CircuitNode *c){  
//     //return (usage < fpga_capacity)
//     for(int i=0;i<8;i++)
//         if(usage[i]+c->node_capacity[i]>fpga_capacity[i])
//             return 0;

//     vector<intg> used_wire(index_to_fpga.size(),0);
//     for(auto &net:c->src_nets){   
//         //对作为src的点来说，一旦放置后会产生割，那所有fpga都要检验   

//         //net里面还没有fpga，或者只有自身那肯定不产生割
//         if(net->used_fpga_node.size()==0||net->used_fpga_node.size()==1&&net->used_fpga_node.count(this))
//             continue;

//         //剩下的情况下都会产生割
//         if(usage_wire+net->weight+used_wire[this->name]>wire_capacity)
//             return -1;
//         else
//             used_wire[this->name]+=net->weight;
//         for(auto fpga_:net->used_fpga_node){
//             if(fpga_==this) continue;
//             if(fpga_->usage_wire+net->weight+used_wire[fpga_->name]>fpga_->wire_capacity)
//                 return -1;
//             else
//                 used_wire[fpga_->name]+=net->weight;
//         }
//     }
//     for(auto &net:c->drain_nets){ 
//         auto source_c=net->source;
//         auto source_f=source_c->fpga_node;
//         //对于drain点来说
//         //如果src已经放置，当前点放置后产生了第二个fpga，那要给自身和src都添加割
//         //如果src已经放置，当前点放置后没有产生新的fpga，那就不影响
//         //如果src已经放置，当前点放置后产生了新的fpga那就只处理自己
//         //如果source点未放放置，也要保证所有drain_fpga之间在hop_fpga上形成一个菊花图,偷个懒，写成相邻
//         if(net->used_fpga_node.count(this))
//             continue;
//         if(source_f==nullptr){
//             for(auto fpga_:net->used_fpga_node)
//                 if(fpga_weight[fpga_->name][this->name]==-1)
//                     return 0;
//         }
//         //如果src放置那还新增了fpga那自身肯定要处理
//         if(!net->used_fpga_node.count(this)){
//             if(usage_wire+net->weight+used_wire[this->name]>wire_capacity)
//                 return -1;
//             else
//                 used_wire[this->name]+=net->weight;
//             //如果只有src_fpga一个fpga，那src_fpga也处理
//             if(net->used_fpga_node.size()==1)
//                 if(source_f->usage_wire+net->weight+used_wire[source_f->name]>source_f->wire_capacity)
//                     return -1;
//                 else    
//                     used_wire[source_f->name]+=net->weight;
//         }
//     }
//     return 1;
// }
void FPGANode::get_demand(CircuitNode* c,std::array<intg,8> &resource_demand,unordered_map<FPGANode*,intg>& wire_demand){
    //统计c想要放在当前FPGA需要"空余出"资源数量各自至少是多少
    //负数就是盈余
    for(int i=0;i<8;i++)
        if(usage[i]+c->node_capacity[i]>fpga_capacity[i])
            resource_demand[i]=usage[i]+c->node_capacity[i]-fpga_capacity[i];
    
    //统计出添加c点后产生的used_wire
    vector<intg> used_wire(index_to_fpga.size(),0);
    for(auto &net:c->src_nets){   
        //对作为src的点来说，一旦放置后会产生割，那所有fpga都要检验   

        //net里面还没有fpga，或者只有自身那肯定不产生割
        if(net->used_fpga_node.size()==0||net->used_fpga_node.size()==1&&net->used_fpga_node.count(this))
            continue;

        //剩下的情况下都会产生割,
        used_wire[this->name]+=net->weight;
        for(auto fpga_:net->used_fpga_node){
            if(fpga_==this) continue;
            used_wire[fpga_->name]+=net->weight;
        }
    }
    for(auto &net:c->drain_nets){ 
        auto source_c=net->source;
        auto source_f=source_c->fpga_node;
        if(source_f==nullptr||net->used_fpga_node.count(this))
            continue;
        //如果src放置那还新增了fpga那自身肯定要处理
        if(!net->used_fpga_node.count(this)){
            used_wire[this->name]+=net->weight;
            //如果只有src_fpga一个fpga，那src_fpga也处理
            if(net->used_fpga_node.size()==1)
                used_wire[source_f->name]+=net->weight;
        }
    }  
    //根据used_wire多的记录下需求
    for(int i=0;i<used_wire.size();i++){
        auto cur_fpga = index_to_fpga[i];
        if(used_wire[i]>0)
            if(wire_demand.count(cur_fpga))
                wire_demand[cur_fpga]+=used_wire[i];
            else
                wire_demand[cur_fpga]=used_wire[i];
    }
}
std::pair<double,double> FPGANode::free_space(CircuitNode *base) {    
                                                //返回剩余资源是base的倍数，和剩余总资源
    //记录base的资源是多少倍
    double max_time=0;   
    double all_src=0;
    for(int i=0;i<8;i++){
        if(fpga_capacity[i]==0)
            continue;
        if(fpga_capacity[i]-usage[i]==0){
            max_time=1000000000;continue;}
        if(base->node_capacity[i]!=0){
            max_time=std::max(max_time,(double)base->node_capacity[i]/(fpga_capacity[i]-usage[i]));
            all_src+=(double)base->node_capacity[i]/(fpga_capacity[i]-usage[i]);
        }
    }    
    return std::make_pair(max_time,all_src);
}
intg FPGANode::free_wire() {   
    return wire_capacity-usage_wire;
}
// intg FPGANode::exchange_valid(CircuitNode *origin,CircuitNode *new_come,FPGANode *new_fpga){
//     /*
//         犯了个错误
//         比如要把5号移动到fpga2，CV的2号放置到fpga3
//         如果是5号没放置的时候检测2号是否有效就会出错
//         应该等5号放置后再检测2号是否有效
//     */
//    //针对CV放置场景
//     if(new_come->fpga_node==nullptr){
//         // std::cout<<"1\n";
//         this->remove_circuit(origin);
//         origin->remove_fpga(this);
        
//         //尝试放置new_come到this 
//         intg result_unstrict=this->valid(new_come);
//         if(result_unstrict<=0){
//             //不严格放置都过不去
//             this->add_circuit(origin);
//             origin->add_fpga(this);
//             return result_unstrict;
//         }
//         //如果可以放置的话
//         this->add_circuit(new_come);
//         new_come->add_fpga(this);
//         result_unstrict=new_fpga->valid(origin);
//         if(result_unstrict<=0){
//             //放置new_come到new_fpga之后无效
//             this->remove_circuit(new_come);
//             new_come->remove_fpga(this);
//             this->add_circuit(origin);
//             origin->add_fpga(this);
//             return result_unstrict;
//         }      
//         this->remove_circuit(new_come);
//         new_come->remove_fpga(this);
//         this->add_circuit(origin);
//         origin->add_fpga(this);    
//     }
//     else{
//         this->remove_circuit(origin);
//         origin->remove_fpga(this);
//         new_fpga->remove_circuit(new_come);
//         new_come->remove_fpga(new_fpga);
//         //尝试放置new_come到this 
//         intg result_unstrict=this->valid(new_come);
//         if(result_unstrict<=0){
//             //不严格放置都过不去,恢复
//             this->add_circuit(origin);
//             origin->add_fpga(this);
//             new_fpga->add_circuit(new_come);
//             new_come->add_fpga(new_fpga);
//             return result_unstrict;
//         }
//         //如果可以放置的话
//         this->add_circuit(new_come);
//         new_come->add_fpga(this);
//         result_unstrict=new_fpga->valid(origin);
//         if(result_unstrict<=0){
//             //放置new_come到new_fpga之后无效，恢复
//             this->remove_circuit(new_come);
//             new_come->remove_fpga(this);
//             this->add_circuit(origin);
//             origin->add_fpga(this);
//             new_fpga->add_circuit(new_come);
//             new_come->add_fpga(new_fpga);
//             return result_unstrict;
//         }      
//         //恢复
//         this->remove_circuit(new_come);
//         new_come->remove_fpga(this);
//         this->add_circuit(origin);
//         origin->add_fpga(this);
//         new_fpga->add_circuit(new_come);
//         new_come->add_fpga(new_fpga);     
//     }
//     return 1;
// }
intg FPGANode::exchange_valid_test(CircuitNode *origin,CircuitNode *new_come){   
    //new_come可能是浮空的，也可能不是浮空的

    //线的比较，只需要比较交叉fpga的线容量即可
    for(int i=0;i<8;i++)
        if(usage[i]+new_come->node_capacity[i]-origin->node_capacity[i]
                                                >fpga_capacity[i])
            return 0;
    //获取相交的Net单独处理
    unordered_set<Net*> intersect_nets;
    //先获得移除origin之后空余出来的线
    vector<intg> used_wire(index_to_fpga.size(),0);
    for(auto &net:origin->src_nets){     
        //如果移除的点作为net的src点
        //如果net只有一个fpga，那肯定是当前fpga，自然没有对外互联数
        //如果net有多个fpga，那所有fpga的对外互联数目都被复原
        if(net->used_fpga_node.size()==1)
            continue;
        used_wire[this->name]-=net->weight;
        for(auto fpga_:net->used_fpga_node)
            if(fpga_!=this)
                used_wire[fpga_->name]-=net->weight;
    }
    for(auto &net:origin->drain_nets){ 
        auto source_c=net->source;
        auto source_f=source_c->fpga_node;
        //如果net只有一个fpga那就不用处理，如果只有两个那就都要处理
        //如果有3个及以上的fpga那就只处理自己就好
        if(net->used_fpga_node.size()==0)
            continue;
        if(net->used_fpga_node_count[this->name]>1)
            used_wire[this->name]-=net->weight;
        if(net->used_fpga_node.size()==2&&net->used_fpga_node_count[this->name]==1)
            used_wire[source_f->name]-=net->weight;
    }   

    //判断添加new_come之后是否会超出线容量
    for(auto &net:new_come->src_nets){   
        if(origin->nets.count(net)){
            intersect_nets.emplace(net);
            continue;
        }
        //对作为src的点来说，一旦放置后会产生割，那所有fpga都要检验   
        //net里面还没有fpga，或者只有自身那肯定不产生割
        if(net->used_fpga_node.size()==0||net->used_fpga_node.size()==1&&net->used_fpga_node.count(this))
            continue;
        //剩下的情况下都会产生割
        if(usage_wire+net->weight+used_wire[this->name]>wire_capacity){
            return -1;
        }
        else
            used_wire[this->name]+=net->weight;
        for(auto fpga_:net->used_fpga_node){
            if(fpga_==this) continue;
            if(fpga_->usage_wire+net->weight+used_wire[fpga_->name]>fpga_->wire_capacity)
                return -1;
            else
                used_wire[fpga_->name]+=net->weight;
        }
    }
    for(auto &net:new_come->drain_nets){ 
        if(origin->nets.count(net)){    
            intersect_nets.emplace(net);
            continue;
        }
        auto source_c=net->source;
        auto source_f=source_c->fpga_node;
        //对于drain点来说
        //如果source未放置那无所谓
        //如果src就是origin，那么自身的割已经被去除了，因此不用处理 
        //如果src不是origin并且this的数目大于1，那也不影响
        //如果src已经放置，当前点放置后产生了第二个fpga，那要给自身和src都添加割
        //如果src已经放置，当前点放置后没有产生新的fpga，那就不影响
        //如果src已经放置，当前点放置后产生了新的fpga那就只处理自己
        if(source_f==nullptr||source_c==origin||source_c!=origin&&net->used_fpga_node_count[this->name]>1)
            continue;
        //如果src放置那还新增了fpga那自身肯定要处理
        if(!net->used_fpga_node.count(this)){
            if(usage_wire+net->weight+used_wire[this->name]>wire_capacity)
                return -1;
            else
                used_wire[this->name]+=net->weight;
            //如果只有src_fpga一个fpga，那src_fpga也处理
            if(net->used_fpga_node.size()==1)
                if(source_f->usage_wire+net->weight+used_wire[source_f->name]>source_f->wire_capacity)
                    return -1;
                else    
                    used_wire[source_f->name]+=net->weight;
        }
    }

    //如果origin和new_come在一个超边 
    for(auto net:intersect_nets){
        if(new_come->src_nets.count(net)){
            if(net->used_fpga_node.size()==1)
                continue;
            //剩下的情况下都会产生割
            if(usage_wire+net->weight+used_wire[this->name]>wire_capacity)
                return -1;
            else
                used_wire[this->name]+=net->weight;
            for(auto fpga_:net->used_fpga_node){
                if(fpga_==this) continue;
                if(fpga_->usage_wire+net->weight+used_wire[fpga_->name]>fpga_->wire_capacity)
                    return -1;
                else
                    used_wire[fpga_->name]+=net->weight;
            }      
        }
        else{
                auto source_c=net->source;
                auto source_f=source_c->fpga_node;
                //origin驱动了new-come
                if(net->source==origin){

                }
                //new_come和origin都被驱动
                else{
                    //如果source_c是origin或source_f=nullptr那就直接跳过
                    if(source_f==nullptr||source_c==origin)
                        continue;
                    //source_f一定不是this，因此放置之后一定会产生割
                    if(!net->used_fpga_node.count(this)){
                        if(usage_wire+net->weight+used_wire[this->name]>wire_capacity)
                            return -1;
                        else
                            used_wire[this->name]+=net->weight;
                        //如果有两个fpga，意味着移除后就剩下source_f了，重新加回来之后都有source_f存在
                        if(net->used_fpga_node.size()==2)
                            if(source_f->usage_wire+net->weight+used_wire[source_f->name]>source_f->wire_capacity)
                                return -1;
                            else    
                                used_wire[source_f->name]+=net->weight;
                    }
                }           
        }
    }
    return 1;
}
}  // namespace topart

