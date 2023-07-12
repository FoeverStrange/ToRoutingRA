function [Jx, F, Rss] = RA(Xi,Xo,Xc,para)
%     RA 先计算router，再计算资源分配
    Jx = 0;
    F = 0;
    Rss = 0;
% 最短路径计算
sum(Xi(user, server, :))
[P,d] = shortestpath(G,3,8);

end