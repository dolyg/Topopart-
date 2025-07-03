#include "graph.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <queue>
#include <utility>

#include "log.h"

namespace topart {

using std::cout;
using std::endl;
using std::make_pair;
using std::max;
using std::pair;
using std::queue;

template <class T>
void Graph<T>::init(vector<T *> &vv, vector<vector<T *>> &gg,vector<vector<intg>> &e_weight) {     //根据顶点向量和邻接表建图
    num_vertex = gg.size();
    max_dist = vector<intg>(num_vertex, 0);
    v = vv;
    for (int i = 0; i < vv.size(); ++i) {
        v_map[vv[i]] = i;
    }
    g = gg;
    weight=e_weight;
    g_set.resize(gg.size());
    for (intg i = 0; i < g_set.size(); ++i) {
        g_set[i] = unordered_set<T *>(gg[i].begin(), gg[i].end());
    }
}

template <class T>
void Graph<T>::calculate_max_dist(Func f) {
    for (int i = 0; i < num_vertex; ++i) {
        max_dist[i] = get_max_dist(i, f);
    }
}

template <class T>
T *Graph<T>::get_vertex(intg i) {
    // assert(i < v.size());
    return v[i];
}


template <class T>
intg Graph<T>::get_max_dist(intg start, Func f) {   //f代表根据源点到某个点的距离进行的操作
    // bfs for max dist
    // So we can calculate S hat value for each fpga
    vector<bool> done;
    done.resize(num_vertex, false);

    intg ans = std::numeric_limits<intg>::min();
    // node id, depth
    queue<pair<intg, intg>> q;
    q.push(make_pair(start, 0));

    done[start] = true;
    while (!q.empty()) {
        auto top = q.front();
        auto &id = top.first;
        auto &depth = top.second;
        q.pop();

        if (depth > 0)
            f(start, id, depth);

        ans = max(ans, depth);

        for (auto &neighbor : g[id]) {
            auto neighbor_id = v_map[neighbor];
            if (done[neighbor_id] == false) {
                q.push(make_pair(neighbor_id, depth + 1));
                done[neighbor_id] = true;
            }
        }
    }
    return ans;
}

template <class T>
intg Graph<T>::opt_get_max_dist(intg start, Func f,unordered_set<intg> &object_to_bfs) {   //为了求node的S集合专门优化的
    // bfs for max dist
    // So we can calculate S hat value for each fpga
    vector<bool> done;
    done.resize(num_vertex, false);

    intg ans = std::numeric_limits<intg>::min();
    // node id, depth
    queue<pair<intg, intg>> q;
    q.push(make_pair(start, 0));

    done[start] = true;
    while (!q.empty()) {
        auto top = q.front();
        auto id = top.first;
        auto depth = top.second;
        q.pop();

        if (depth > 0)
            f(start, id, depth);

        ans = max(ans, depth);

        for (auto &neighbor : g[id]) {
            auto neighbor_id = v_map[neighbor];
            if (done[neighbor_id] == false && object_to_bfs.count(neighbor_id)<=0) {
                q.push(make_pair(neighbor_id, depth + 1));
                done[neighbor_id] = true;
            }
        }
    }
    return ans;
}

template <class T>
intg Graph<T>::get_weight(T* src,T* dst){
    return weight[v_map[src]][v_map[dst]];
}

template <class T>
void Graph<T>::get_status() {
#ifdef LOG_DB
    std::cout << "=============================================\n" ;
    log("Max dist list: ");
    for (int i = 0; i < num_vertex; ++i) {
        cout << i << "'s max dist is: " << max_dist[i] << endl;
    }

    std::cout << "=============================================\n" ;
    log("S-hat list: ");
    for (int i = 0; i < num_vertex; ++i) {
        cout << v[i]->status();
    }
    cout << std::flush;
#endif
}


template class Graph<FPGANode>;
template class Graph<CircuitNode>;

}  // namespace topart
