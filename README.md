# Topopart+
It's my implementation topopart+ with 2024EDAcompetition with http://eda.icisc.cn/

The project is modified based on BWbwchen's TopoPart,But I don`t know how to use Github to show its relationship with his project.

The project is a demo,in my project you can see how I implement the topopart+

And the project also show my mind that how I ensure the partition is legally valid,u can see the refine_CV() and Force_CV() to understand my ideas.

And i am a beginner and i don`t know how to use github professionally.

I would appreciate it if you could give more advice.

# technique report 

## Public Case Score Report
| #name  | #Hop   | #time    |
|--------|--------|----------|
| Case01 | 12     | 0.00354529 |
| Case02 | 3127   | 0.181807   |
| Case03 | 14850  | 56.9944    |
| Case04 | 739974 | 651.164    |

| #name  | #Hop    | #time    | #increment |
|--------|---------|----------|------------|
| Case01 | 20      | 0.00354529 | 1.66       |
| Case02 | 3197    | 0.181807   | 1.022      |
| Case03 | 18877   | 56.9944    | 1.27       |
| Case04 | 840357  | 651.164    | 1.14       |

##     Hidden Case score report 

| info_type | runtime(s) | hop length | score  | rank score |
| ------------ | --------- | ---------- | ---------- | ---------- |
|  case1     | 0     | 12      | 12   | 14       |
|  case2     | 0.0        | 3127.0     | 3127.0 | 17.0       |
|  case3     | 66.0       | 14850.0    | 14904.0 | 14.0       |
|  case4     | 610.0      | 739974.0   | 765051.0 | 9.0        |
|  case5     | 0.0        | 1261.0     | 1261.0 | 8.0        |
|  case6     | 1.0        | 1442.0     | 1442.0 | 6.0        |
|  case7     | 28.0       | 29750.0    | 29796.0 | 9.0        |
|  case8     | 50.0       | 26702.0    | 26776.0 | 9.0        |
|  case10    | 740.0      | 672951.0   | 700617.0 | 10.0       |

# How to get all the case
If u want to download all the case(about 10),please see the Case download information.pdf.

# English report
## Algorithm Introduction

### Flowchart:

#### Core Algorithm
This algorithm is a modification based on the work of Mapart.

#### Initial Partitioning
Simultaneously considers the remaining FPGA space and the hop count growth of the entire graph during partitioning.

#### Handling Unplaced Nodes
After initial partitioning, uses `refine_CV` and `force_CV` to attempt placing nodes that failed to be partitioned.

#### Hybrid Optimization
Integrates both node replication and node movement strategies during refinement to ensure consistency, achieving optimization improvements ranging from 3% (Case 02) to 30% (Case 03).

### Innovations

#### Multi - Level Refinement Framework
Addresses limitations in Mapart's refinement phase, which relies on the Fiducia - Mattheyses (FM) algorithm. Experiments showed that refining other nodes while handling violations incrementally could lead to incomplete partitioning, high computational costs, and degraded results. Proposed `refine_CV` and `force_CV` algorithms based on a multi - level framework to efficiently resolve violations for nodes with varying resource requirements, ensuring legality and minimizing negative impacts on final results.

#### Hybrid Node Movement Algorithm
Inspired by the FM algorithm, this approach combines node movement and replication strategies to maintain consistency. Demonstrates significant improvements in certain scenarios, particularly when initial partitioning is less constrained.

#### Space Evaluation Function
Introduces a new `spaceFunction` to assess the quality of remaining FPGA space, considering the spatial overhead of logic replication.

### Experiments
#### Combining Movement and Replication in FM
Directly integrating node replication into the FM algorithm yielded suboptimal results.

# 中文介绍
算法介绍
流程图:

核心算法:
本算法是基于Mapart工作做的修改.
在Partition中，同时考虑FPGA的剩余空间和整个图的hop增长，进行初始划分.

在初始化分后对于未能成功划分的节点使用refine_CV和force_CV去尝试放置.

考虑到逻辑复制和节点移动的一致性，在refine中混合入了节点复制和节点移动两种选择，带来了最高百分之30(case03)，最低百分之3(case02)的优化结果。

创新点:

1.改进了Mapart的工作，Mapart是在refine阶段基于Fiducia Mattheyses算法处理的违例，在实验中发现，如果初始划分后违例节点过去，边细化其他节点，边处理违例节点的做法可能无法达到合法的要求，事件上也消耗巨大也会让细化结果恶化不少。因此提出了基于多级框架的refine_CV,force_CV，该方法针对于违例节点对资源的两种不同需求情况提出，不但可以保证在极快时间内完美处理任何比例的初始化分中的违例节点也可以保证减少对最终结果的坏影响。

2.受到Fiducia Mattheyses算法启发，考虑到移动节点和复制节点的一致性，提出了混合移动节点和复制节点的多级点移动算法，该算法可以在一些情况下提升巨大，在初始划分紧张的情况下提升较小。

3.由于逻辑复制对空间的占用，因此对FPGA的空间预留是有要求的，提出新的spacefunction用来评价FPGA的剩余空间的优劣。

尝试:

1.在FM算法的基础上，尝试让节点选择移动和选择复制，效果不好。

2.尝试假定一个节点已经复制过，去进行菊花图为单位的细化，效果不好

3.尝试假定一个节点的复制点被移除，去进行菊花图为单位复制撤回，效果不好

4.十分希望建立一个流模型，但是考虑不到怎么去处理一个节点移走之后空闲fpga的容量怎么流出来。
