#include "node.h"

namespace topart {

//����FPGA֮���жϸ��ԵĹ�ϵ,�Ǹ�������bd.cpp��build_hop_fpga�г�ʼ��
unordered_map<FPGANode *, intg> fpga_to_index;
vector<FPGANode*> index_to_fpga;       
vector<vector<intg>> fpga_weight;      //�Ϸ����о��룬���Ϸ�ֻ��-1
vector<vector<intg>> fpga_parents;     //fpga[i]����iΪ��ʼ��,���е�ĸ���

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
void CircuitNode::calculate_cut_increment(              //����hop������
    FPGANode *if_f) {
    intg cut_increment = 0;

    for(auto &n:src_nets){      //���ڵ�ǰ�ڵ���Ϊsrc�ڵ㣬����ÿһ���������㶼����hop
        intg net_weight=n->weight;
        for(auto &fpga:n->used_fpga_node){
            intg cur_weight = fpga_weight[if_f->name][fpga->name];
            cut_increment+=cur_weight*net_weight;
        }
    }
    for(auto &n:drain_nets){    //node�Ǳ������ڵ㣬srcδ���û�if-f�Ѿ����ֹ��Ͳ��ٲ���hop
        auto source=n->source;       
        if(source->fpga_node==nullptr||n->used_fpga_node.count(if_f)>0) continue;
        cut_increment+=n->weight*fpga_weight[if_f->name][source->fpga_node->name];
    }
    cut_increment_map[if_f->name] = cut_increment;
}

void CircuitNode::calculate_connect_increment(          //ר�Ŵ�����⻥����Ŀ
    FPGANode *if_f) {
    //������c��û�з��ù�
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
void CircuitNode::calculate_copy_increment(FPGANode *f){//���Ƶ�hop
    intg copy_increment=0;
    //���ԭʼ��͸��Ƶ���һ��Ӱ����⻥��
    if(this->fpga_node==f){
        copy_increment_map[f->name]=0;
        return ;
    }

    //����Node����ǰFPGA��Ҫ����Node��src��Node�����Ķ��⻥����Ŀ������ 
    for(auto drain_:this->drain_nets){
        auto src_node=drain_->source;
        auto src_fpga=src_node->fpga_node;
        auto src_copy_fpga=src_node->copy_position;
        if(src_fpga==nullptr)   continue;
        //�������²��������
        if(src_fpga==f||src_copy_fpga.count(f))
            continue;
        //�����Ѿ�����������
        if(drain_->used_fpga_node.count(f)||drain_->used_copy_fpga.count(f))
            continue;
        copy_increment+=drain_->weight*fpga_weight[f->name][src_fpga->name];
    }//for(auto drain_:c->drain_nets)
    
    //����c�㣬��this֮��ı���������û��Ӱ���
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
    //�����ڵ�ǰ��û��fpga_node�������
    intg cut_increment = 0;

    for(auto &n:src_nets){   
        //this����if_f�ϣ��Ǿ�Ҫ���㱻�������fpga_node��copy_node�͵�ǰ���if_f�Լ�copy_position�Ĺ�ϵ
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
        //���c�ƶ���if_f����Ҫ����if_f������source_node��copy_position�л�����sourcefpga_node��Υ��
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
        //Ϊɶ����һ�λ�������أ�
        //������copy_position��hop����
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

//һ������ʹ�����������ʹ��add_fpga,remove_fpga
void FPGANode::add_circuit(CircuitNode* c) { 
    for (size_t i = 0; i < 8; ++i) {
        usage[i] = usage[i] + c->node_capacity[i];
    }
    //���������йصĳ��߸��������нڵ��wire
    for(auto &net:c->src_nets){     
        //net���滹û��fpga������ֻ�������ǿ϶���������
        if(net->used_fpga_node.size()==0||net->used_fpga_node.size()==1&&net->used_fpga_node.count(this))
            continue;
        //ʣ�µ�����¶��������
        usage_wire+=net->weight;
        for(auto fpga_:net->used_fpga_node){
            if(fpga_==this) continue;
            fpga_->usage_wire+=net->weight;
        }
    }
    for(auto &net:c->drain_nets){ 
        auto source_c=net->source;
        auto source_f=source_c->fpga_node;
        //����drain����˵
        //���sourceδ����������ν
        //���src�Ѿ����ã���ǰ����ú�����˵ڶ���fpga����Ҫ�������src����Ӹ�
        //���src�Ѿ����ã���ǰ����ú�û�в����µ�fpga���ǾͲ�Ӱ��
        //���src�Ѿ����ã���ǰ����ú�������µ�fpga�Ǿ�ֻ�����Լ�
        if(source_f==nullptr||net->used_fpga_node.count(this))
            continue;
        //���src�����ǻ�������fpga������϶�Ҫ����
        this->usage_wire+=net->weight;
        //���ֻ��src_fpgaһ��fpga����src_fpgaҲ����
        if(net->used_fpga_node.size()==1)
            source_f->usage_wire+=net->weight;
    }   
    placed_circuit.emplace(c);
}
void FPGANode::add_circuit_undercopy(CircuitNode* c) {
    for (size_t i = 0; i < 8; ++i) {
        usage[i] = usage[i] + c->node_capacity[i];
    }
    //���������йصĳ��߸��������нڵ��wire
    for(auto &net:c->src_nets){     
        //net���滹û��fpga������ֻ�������Ǿ�û�и�
        if(net->used_all_fpga.size()==0||net->used_all_fpga.count(this)&&net->used_all_fpga.size()==1)
            continue;
        //��ѡ������copy�����ǣ�Ҳ����ǰ��FPGA�޹ص�FPGA
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
        //�����this�ĵط�
        if(source_f==nullptr||net->used_all_fpga.count(this))
            continue;
        // if(source_f==this||source_copyf.count(this))    continue;
        this->usage_wire+=net->weight;
        
        //���source_f��fpga����src-fpga����src-copy�ģ��Ǿͻ�����µĸ�
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
    //���������йصĳ��߸��������нڵ��wire
    // Ĭ��used_fpga_node������һ��this
    for(auto &net:c->src_nets){     
        //����Ƴ��ĵ���Ϊnet��src��
        //���netֻ��һ��fpga���ǿ϶��ǵ�ǰfpga����Ȼû�ж��⻥����
        //���net�ж��fpga��������fpga�Ķ��⻥����Ŀ������ԭ
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
        //���netֻ��һ��fpga�ǾͲ��ô���
        //���thisռ�����μ����ϣ��ǾͲ��ô����Լ�
        if(source_f==nullptr||net->used_fpga_node.size()==1)
            continue;
        //���ɾȥ�Լ�����fpga��Ŀ�仯��������϶�Ҫ����
        if(net->used_fpga_node_count[this->name]<=1)
            this->usage_wire-=net->weight;
        //���ɾȥ�Լ��󣬸պ��ø���ʧ����srcҲҪ�仯
        if(net->used_fpga_node.size()==2&&net->used_fpga_node_count[this->name]==1)
            source_f->usage_wire-=net->weight;
    }   
    placed_circuit.erase(c);
}
void FPGANode::remove_circuit_undercopy(CircuitNode* c) {
    for (size_t i = 0; i < 8; ++i) {
        usage[i] = usage[i] - c->node_capacity[i];
    }
    //���������йصĳ��߸��������нڵ��wire
    // Ĭ��used_fpga_node������һ��this
    for(auto &net:c->src_nets){     
        auto source_f=c->fpga_node;
        auto source_copyf=c->copy_position;
        auto c_copyf=c->copy_position;
        //���û�в������������,��result=0
        #ifdef LOG
            assert(source_f==this);
            assert(net->used_all_fpga.size()>=(1+source_copyf.size()));
        #endif
        if(net->used_all_fpga.size()==1+source_copyf.size())
            continue; 
        //������������ô���в��������fpga��Ҫ��ȥ��
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
        //������ڹ��Ͳ�����
        if(source_f==nullptr||this==source_f||source_copyf.count(this))
            continue;
        //���ɾȥ�Լ�����fpga��Ŀ�仯��������϶�Ҫ����
        if(!net->used_copy_fpga.count(this)&&net->used_fpga_node_count[this->name]<=1){
            this->usage_wire-=net->weight;
            //�������source-f��source_copyf֮��ֻ��this����src�ĸ�ҲҪȥ��
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
//���Ҹ����ٴ���circuit
void FPGANode::add_copy_circuit(CircuitNode* c){
    //���Ʋ������ж�net������FPGA�����������
    //�Ƿ������Ҳ����
    //���ƺ����drain_node���Ƿ����hopҲҪ�����Ƶ�src_node�ǲ�����һ����
    //����src_node��������hopҲҪ�����ڸ��Ƶ�drain_nodeҲ������µĸ�
    for (size_t i = 0; i < 8; ++i) {
        usage[i] = usage[i] + c->node_capacity[i];
    }

    //���ԭʼ��͸��Ƶ���һ��Ӱ����⻥��
    if(c->fpga_node==this)  return;

    //����Node����ǰFPGA��Ҫ����Node��src��Node�����Ķ��⻥����Ŀ������ 
    for(auto drain_:c->drain_nets){
        auto src_node=drain_->source;
        auto src_fpga=src_node->fpga_node;
        auto src_copy_fpga=src_node->copy_position;
        //���������Դ����߸��Ƶ����ڵ�FPGA������
        if(src_fpga==nullptr) continue;
        if(src_fpga==this||src_copy_fpga.count(this))
            continue;
        //�����Ѿ�����������
        if(drain_->used_all_fpga.count(this))
            continue;
        // if(drain_->used_copy_fpga.count(this)||drain_->used_fpga_node.count(this))
        //     continue;
        usage_wire+=drain_->weight;
        //���û��������(used_copy_fpga��ĺ�used_fpga_node����Ķ���src_fpga��src_copy_fpga�ڵ�)���Ǿ͸�src���ϸ�
        #ifdef LOG
        assert(drain_->used_all_fpga.size()>=src_copy_fpga.size()+1);
        #endif
        if(drain_->used_all_fpga.size()>src_copy_fpga.size()+1) continue;
        src_fpga->usage_wire+=drain_->weight;
    }//for(auto drain_:c->drain_nets)
    
    //����c�㣬��this֮��ı���������û��Ӱ���
    for(auto src_:c->src_nets){
        auto src_fpga=c->fpga_node;
        auto src_copy_fpga=c->copy_position;
        //�����ǰfpgaû�в������Ǿ�����
        if(!src_->used_all_fpga.count(this))
            continue;
        // if(!src_->used_copy_fpga.count(this)&&!src_->used_fpga_node.count(this))
        //     continue;
        if(src_copy_fpga.count(this)||src_fpga==this)     
            continue;    
        //����ǰ�ڵ��������
        usage_wire-=src_->weight;
        //result=0˵��û��û�к�������FPGA�и�ǾͿ���ȥ����
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
        //�ж�this�Ƿ��ڵ�ǰnetֻ��һ��������ж��������
        if(src_fpga==nullptr)   continue;
        // if(drain_->used_all_fpga.count(this))
        //     continue;
        if(drain_->used_fpga_node.count(this)||drain_->used_copy_fpga_count[this->name]>1)
            continue;
        if(src_fpga==this||src_copy_fpga.count(this))
            continue;
        usage_wire-=drain_->weight;  
        //���ֻ��this�и��Ǿͼ�ȥsrc�ĸ�
        #ifdef LOG
        assert(drain_->used_all_fpga.size()>=src_copy_fpga.size()+2);
        #endif
        if(drain_->used_all_fpga.size()>src_copy_fpga.size()+2) continue;
        src_fpga->usage_wire-=drain_->weight;
    }//for(auto drain_:c->drain_nets)
    
    //���ڵ�ǰ�ڵ������Ľڵ�,����������ڵ�this����ڱ�����������临�Ƶ㣬�Ǿͼ��ϸ�
    //����⵼��src_fpga�������ֻ��this���Ǿͼ���src_fpga�ĸ�
    for(auto src_:c->src_nets){
        auto src_fpga=c->fpga_node;
        auto &src_copy_fpga=c->copy_position;
        //���û��c��copyλ��û�������ڵ���ڣ��Ǿ�����
        if(src_->used_fpga_node.count(this)==0&&src_->used_copy_fpga_count[this->name]==1)
            continue;
        usage_wire+=src_->weight;
        //���srcԭ��û�в��������Ǿͼ��ϸ�
        //srcԭ��û�в�������ָ����src_->used_all_fpga.size()==src_copy_fpga.size()+1
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
//         //����Ϊsrc�ĵ���˵��һ�����ú������������fpga��Ҫ����   

//         //net���滹û��fpga������ֻ�������ǿ϶���������
//         if(net->used_fpga_node.size()==0||net->used_fpga_node.size()==1&&net->used_fpga_node.count(this))
//             continue;

//         //ʣ�µ�����¶��������
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
//         //����drain����˵
//         //���src�Ѿ����ã���ǰ����ú�����˵ڶ���fpga����Ҫ�������src����Ӹ�
//         //���src�Ѿ����ã���ǰ����ú�û�в����µ�fpga���ǾͲ�Ӱ��
//         //���src�Ѿ����ã���ǰ����ú�������µ�fpga�Ǿ�ֻ�����Լ�
//         //���source��δ�ŷ��ã�ҲҪ��֤����drain_fpga֮����hop_fpga���γ�һ���ջ�ͼ,͵������д������
//         if(net->used_fpga_node.count(this))
//             continue;
//         if(source_f==nullptr){
//             for(auto fpga_:net->used_fpga_node)
//                 if(fpga_weight[fpga_->name][this->name]==-1)
//                     return 0;
//         }
//         //���src�����ǻ�������fpga������϶�Ҫ����
//         if(!net->used_fpga_node.count(this)){
//             if(usage_wire+net->weight+used_wire[this->name]>wire_capacity)
//                 return -1;
//             else
//                 used_wire[this->name]+=net->weight;
//             //���ֻ��src_fpgaһ��fpga����src_fpgaҲ����
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
    //ͳ��c��Ҫ���ڵ�ǰFPGA��Ҫ"�����"��Դ�������������Ƕ���
    //��������ӯ��
    for(int i=0;i<8;i++)
        if(usage[i]+c->node_capacity[i]>fpga_capacity[i])
            resource_demand[i]=usage[i]+c->node_capacity[i]-fpga_capacity[i];
    
    //ͳ�Ƴ����c��������used_wire
    vector<intg> used_wire(index_to_fpga.size(),0);
    for(auto &net:c->src_nets){   
        //����Ϊsrc�ĵ���˵��һ�����ú������������fpga��Ҫ����   

        //net���滹û��fpga������ֻ�������ǿ϶���������
        if(net->used_fpga_node.size()==0||net->used_fpga_node.size()==1&&net->used_fpga_node.count(this))
            continue;

        //ʣ�µ�����¶��������,
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
        //���src�����ǻ�������fpga������϶�Ҫ����
        if(!net->used_fpga_node.count(this)){
            used_wire[this->name]+=net->weight;
            //���ֻ��src_fpgaһ��fpga����src_fpgaҲ����
            if(net->used_fpga_node.size()==1)
                used_wire[source_f->name]+=net->weight;
        }
    }  
    //����used_wire��ļ�¼������
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
                                                //����ʣ����Դ��base�ı�������ʣ������Դ
    //��¼base����Դ�Ƕ��ٱ�
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
//         ���˸�����
//         ����Ҫ��5���ƶ���fpga2��CV��2�ŷ��õ�fpga3
//         �����5��û���õ�ʱ����2���Ƿ���Ч�ͻ����
//         Ӧ�õ�5�ŷ��ú��ټ��2���Ƿ���Ч
//     */
//    //���CV���ó���
//     if(new_come->fpga_node==nullptr){
//         // std::cout<<"1\n";
//         this->remove_circuit(origin);
//         origin->remove_fpga(this);
        
//         //���Է���new_come��this 
//         intg result_unstrict=this->valid(new_come);
//         if(result_unstrict<=0){
//             //���ϸ���ö�����ȥ
//             this->add_circuit(origin);
//             origin->add_fpga(this);
//             return result_unstrict;
//         }
//         //������Է��õĻ�
//         this->add_circuit(new_come);
//         new_come->add_fpga(this);
//         result_unstrict=new_fpga->valid(origin);
//         if(result_unstrict<=0){
//             //����new_come��new_fpga֮����Ч
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
//         //���Է���new_come��this 
//         intg result_unstrict=this->valid(new_come);
//         if(result_unstrict<=0){
//             //���ϸ���ö�����ȥ,�ָ�
//             this->add_circuit(origin);
//             origin->add_fpga(this);
//             new_fpga->add_circuit(new_come);
//             new_come->add_fpga(new_fpga);
//             return result_unstrict;
//         }
//         //������Է��õĻ�
//         this->add_circuit(new_come);
//         new_come->add_fpga(this);
//         result_unstrict=new_fpga->valid(origin);
//         if(result_unstrict<=0){
//             //����new_come��new_fpga֮����Ч���ָ�
//             this->remove_circuit(new_come);
//             new_come->remove_fpga(this);
//             this->add_circuit(origin);
//             origin->add_fpga(this);
//             new_fpga->add_circuit(new_come);
//             new_come->add_fpga(new_fpga);
//             return result_unstrict;
//         }      
//         //�ָ�
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
    //new_come�����Ǹ��յģ�Ҳ���ܲ��Ǹ��յ�

    //�ߵıȽϣ�ֻ��Ҫ�ȽϽ���fpga������������
    for(int i=0;i<8;i++)
        if(usage[i]+new_come->node_capacity[i]-origin->node_capacity[i]
                                                >fpga_capacity[i])
            return 0;
    //��ȡ�ཻ��Net��������
    unordered_set<Net*> intersect_nets;
    //�Ȼ���Ƴ�origin֮������������
    vector<intg> used_wire(index_to_fpga.size(),0);
    for(auto &net:origin->src_nets){     
        //����Ƴ��ĵ���Ϊnet��src��
        //���netֻ��һ��fpga���ǿ϶��ǵ�ǰfpga����Ȼû�ж��⻥����
        //���net�ж��fpga��������fpga�Ķ��⻥����Ŀ������ԭ
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
        //���netֻ��һ��fpga�ǾͲ��ô������ֻ�������ǾͶ�Ҫ����
        //�����3�������ϵ�fpga�Ǿ�ֻ�����Լ��ͺ�
        if(net->used_fpga_node.size()==0)
            continue;
        if(net->used_fpga_node_count[this->name]>1)
            used_wire[this->name]-=net->weight;
        if(net->used_fpga_node.size()==2&&net->used_fpga_node_count[this->name]==1)
            used_wire[source_f->name]-=net->weight;
    }   

    //�ж����new_come֮���Ƿ�ᳬ��������
    for(auto &net:new_come->src_nets){   
        if(origin->nets.count(net)){
            intersect_nets.emplace(net);
            continue;
        }
        //����Ϊsrc�ĵ���˵��һ�����ú������������fpga��Ҫ����   
        //net���滹û��fpga������ֻ�������ǿ϶���������
        if(net->used_fpga_node.size()==0||net->used_fpga_node.size()==1&&net->used_fpga_node.count(this))
            continue;
        //ʣ�µ�����¶��������
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
        //����drain����˵
        //���sourceδ����������ν
        //���src����origin����ô����ĸ��Ѿ���ȥ���ˣ���˲��ô��� 
        //���src����origin����this����Ŀ����1����Ҳ��Ӱ��
        //���src�Ѿ����ã���ǰ����ú�����˵ڶ���fpga����Ҫ�������src����Ӹ�
        //���src�Ѿ����ã���ǰ����ú�û�в����µ�fpga���ǾͲ�Ӱ��
        //���src�Ѿ����ã���ǰ����ú�������µ�fpga�Ǿ�ֻ�����Լ�
        if(source_f==nullptr||source_c==origin||source_c!=origin&&net->used_fpga_node_count[this->name]>1)
            continue;
        //���src�����ǻ�������fpga������϶�Ҫ����
        if(!net->used_fpga_node.count(this)){
            if(usage_wire+net->weight+used_wire[this->name]>wire_capacity)
                return -1;
            else
                used_wire[this->name]+=net->weight;
            //���ֻ��src_fpgaһ��fpga����src_fpgaҲ����
            if(net->used_fpga_node.size()==1)
                if(source_f->usage_wire+net->weight+used_wire[source_f->name]>source_f->wire_capacity)
                    return -1;
                else    
                    used_wire[source_f->name]+=net->weight;
        }
    }

    //���origin��new_come��һ������ 
    for(auto net:intersect_nets){
        if(new_come->src_nets.count(net)){
            if(net->used_fpga_node.size()==1)
                continue;
            //ʣ�µ�����¶��������
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
                //origin������new-come
                if(net->source==origin){

                }
                //new_come��origin��������
                else{
                    //���source_c��origin��source_f=nullptr�Ǿ�ֱ������
                    if(source_f==nullptr||source_c==origin)
                        continue;
                    //source_fһ������this����˷���֮��һ���������
                    if(!net->used_fpga_node.count(this)){
                        if(usage_wire+net->weight+used_wire[this->name]>wire_capacity)
                            return -1;
                        else
                            used_wire[this->name]+=net->weight;
                        //���������fpga����ζ���Ƴ����ʣ��source_f�ˣ����¼ӻ���֮����source_f����
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

