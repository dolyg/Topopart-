#include "node.h"

namespace topart {


intg Net::estimate_increase_cut_size_refine(FPGANode *if_f,FPGANode* origin,bool src) {
    //增益等于原来origin的hop - if_f贡献的hop，意思是hop减少了多少
    intg result =0;
    //如果这个超边一个都没有，那肯定非法 
    if (used_fpga_node.size() == 0) {
        return std::numeric_limits<intg>::min();
    } 

    //如果是移动源点,那么移动后hop的贡献取决于Net里所有fpga与之的贡献
    if(src){
          //if_f违例检查同时计算if_f的贡献
        for(auto fpga:used_fpga_node){
            if(fpga==origin&&used_fpga_node_count[origin->name]==1) 
                continue;
            if(fpga_weight[fpga->name][if_f->name]==-1){
                //std::cout<<"Net::estimate_increase_cut_size_refine: fpga_weight["<<if_f->name<<"]["<<fpga->name<<"] is -1"<<std::endl;
                return std::numeric_limits<intg>::min();
            }
            result-=fpga_weight[if_f->name][fpga->name];
        }  
        //计算origin的hop贡献
        for(auto fpga:used_fpga_node)
            result+=fpga_weight[origin->name][fpga->name];

    }
    //如果是移动drain
    else{
        //违例检查
        if(fpga_weight[source->fpga_node->name][if_f->name]==-1)
            return std::numeric_limits<intg>::min();
        //如果origin移动后消失那么收益增加
        if(used_fpga_node_count[origin->name]==1)
            result+=fpga_weight[source->fpga_node->name][if_f->name];
        //如果if_f本来不存在那肯定收益减少
        if(used_fpga_node_count[if_f->name]==0)
            result-=fpga_weight[source->fpga_node->name][if_f->name];
    }
    return weight*result;
}

bool Net::try_move(FPGANode *f) {
    //返回Net对应的FPGA的数量
    return used_fpga_node_count[f->name] <= 1;
}

void Net::add_copy_circuit(CircuitNode *c){
    copy_circuit.emplace(c);
}
void Net::remove_copy_circuit(CircuitNode *c){
    copy_circuit.erase(c);
}
void Net::add_copy_fpga(FPGANode *f) {
    used_copy_fpga.emplace(f);
    used_all_fpga.emplace(f);
    if (used_copy_fpga_count.count(f->name) <= 0) {
        used_copy_fpga_count[f->name] = 1;
    } else {
        used_copy_fpga_count[f->name]++;
    }
}
void Net::remove_copy_fpga(FPGANode *f) {
    if (used_copy_fpga_count[f->name] == 1) {
        used_copy_fpga_count.erase(f->name);
        used_copy_fpga.erase(f);
        if(!used_fpga_node.count(f))
            used_all_fpga.erase(f);
    } else {
        used_copy_fpga_count[f->name] -= 1;
    }
}

}  // namespace topart
