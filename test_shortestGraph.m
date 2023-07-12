% 起始节点序号
s = [1 1 1 2 2 6 6 7 7 3 3 9 9 4 4 11 11 8];
% 目标节点序号
t = [2 3 4 5 6 7 8 5 8 9 10 5 10 11 12 10 12 12];
% 边的权重
weights = [10 10 10 10 10 1 1 1 1 1 1 1 1 1 1 1 1 1];
G = graph(s,t,weights);
plot(G,'EdgeLabel',G.Edges.Weight)

[P,d] = shortestpath(G,3,8)