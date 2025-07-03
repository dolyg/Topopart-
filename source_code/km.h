#include "type.h"
#include <vector>
#include <queue>
#include <map>
#include <limits>
#include <numeric>

using std::vector;
using std::queue;
using std::pair;
using std::move;
using std::numeric_limits;

template <typename T>
struct hungarian {  // km
  intg n;
  vector<intg> matchx;
  vector<intg> matchy;
  vector<intg> pre;
  vector<bool> visx;
  vector<bool> visy;
  vector<T> lx;
  vector<T> ly;
  vector<vector<T> > g;
  vector<T> slack;
  T inf;
  T res;
  queue<intg> q;
  intg org_n;
  intg org_m;
	intg max(intg a,intg b){
		return a>b?a:b;
	}
	void swap(intg &a,intg &b){
		intg temp=a;
		a=b;
		b=temp;
	}
  hungarian(intg _n, intg _m) {
    org_n = _n;
    org_m = _m;
    n = max(_n, _m);
    inf = numeric_limits<T>::max();
    res = 0;
    g = vector<vector<T> >(n, vector<T>(n));
    matchx = vector<intg>(n, -1);
    matchy = vector<intg>(n, -1);
    pre = vector<intg>(n);
    visx = vector<bool>(n);
    visy = vector<bool>(n);
    lx = vector<T>(n, -inf);
    ly = vector<T>(n);
    slack = vector<T>(n);
  }

  void addEdge(intg u, intg v, intg w) {
    g[u][v] = max(w, 0);  // 负值还不如不匹配 因此设为0不影响
  }

  bool check(intg v) {
    visy[v] = true;
    if (matchy[v] != -1) {
      q.push(matchy[v]);
      visx[matchy[v]] = true;
      return false;
    }
    while (v != -1) {
      matchy[v] = pre[v];
      swap(v, matchx[pre[v]]);
    }
    return true;
  }

  void bfs(intg i) {
    while (!q.empty()) {
      q.pop();
    }
    q.push(i);
    visx[i] = true;
    while (true) {
      while (!q.empty()) {
        intg u = q.front();
        q.pop();
        for (intg v = 0; v < n; v++) {
          if (!visy[v]) {
            T delta = lx[u] + ly[v] - g[u][v];
            if (slack[v] >= delta) {
              pre[v] = u;
              if (delta) {
                slack[v] = delta;
              } else if (check(v)) {
                return;
              }
            }
          }
        }
      }
      // 没有增广路 修改顶标
      T a = inf;
      for (intg j = 0; j < n; j++) {
        if (!visy[j]) {
          a = std::min(a, slack[j]);
        }
      }
      for (intg j = 0; j < n; j++) {
        if (visx[j]) {  // S
          lx[j] -= a;
        }
        if (visy[j]) {  // T
          ly[j] += a;
        } else {  // T'
          slack[j] -= a;
        }
      }
      for (intg j = 0; j < n; j++) {
        if (!visy[j] && slack[j] == 0 && check(j)) {
          return;
        }
      }
    }
  }

  intg solve(vector<intg> &match) {
    // 初始顶标
    for (intg i = 0; i < n; i++) {
      for (intg j = 0; j < n; j++) {
        lx[i] = max(lx[i], g[i][j]);
      }
    }

    for (intg i = 0; i < n; i++) {
      fill(slack.begin(), slack.end(), inf);
      fill(visx.begin(), visx.end(), false);
      fill(visy.begin(), visy.end(), false);
      bfs(i);
    }

    // custom
    for (intg i = 0; i < n; i++) {
      if (g[i][matchx[i]] > 0) {
        res += g[i][matchx[i]];
      } else {
        matchx[i] = -1;
      }
    }
    // for (intg i = 0; i < org_n; i++) {
    //   cout << matchx[i] + 1 << " ";
    // }
    // cout << "\n";
	match=move(matchx);
	return res;
  }
};
intg bitpart_match(vector<intg> &bi_lv,vector<intg> &bi_rv,
			unordered_map<intg,unordered_map<intg,intg>> &graph,
			vector<pair<intg,intg>> &ans);
intg bitpart_match(vector<intg> &bi_lv,vector<intg> &bi_rv,
			unordered_map<intg,unordered_map<intg,intg>> &graph,
			vector<pair<intg,intg>> &ans) {
	//给的顶点是从0开始的，这里的是1开始的，注意偏置
	intg n=bi_lv.size(),m=bi_rv.size();
	unordered_map<intg,intg> lv_to_index;
	unordered_map<intg,intg> rv_to_index;
	vector<intg> match;
	for(intg i=0;i<n;i++)
		lv_to_index[bi_lv[i]]=i;
	for(intg i=0;i<m;i++)
		rv_to_index[bi_rv[i]]=i;

	hungarian<intg> solver(n, m);
	for(auto &eedge:graph){
		intg u,v,w;
		u=eedge.first;
		for(auto &e:eedge.second){
			v=e.first;
			w=e.second;
			u=lv_to_index[u];
			v=rv_to_index[v];
			solver.addEdge(u,v,w);
		}
	}
  intg max_power=solver.solve(match);
  for(intg i=0;i<match.size();i++){
    if (match[i]<0)
        continue;
	ans.emplace_back(make_pair(bi_lv[i],bi_rv[match[i]]));
  }
  return max_power;
}
