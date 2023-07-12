function [Jx, F, Rss] = RA(Xi,Xo,Xc,para)

[userNumber, serverNumber, subbandNumber] = size(Xi);
G = para.G;
%     RA 先计算router，再计算资源分配
    Jx = 0;
    F = 0;
    Rss = 0;
% 最短路径计算
Xi_temp = sum(Xi, 3);
Xo_temp = sum(Xo, 3);
% 最短路径，sum(1/H)存储
H_sum = zeros(userNumber, 2);
for user = 1:userNumber
    server_in = find(Xi_temp(user, :) == 1);
    server_out = find(Xo_temp(user, :) == 1);
    server_cp = find(Xc(user, :) == 1);
    
    [~,d] = shortestpath(G,server_in,server_cp);
    H_sum(user,1) = d;
    [~,d] = shortestpath(G,server_cp,server_out);
    H_sum(user,2) = d;
end

[F,res_cra] = cra(x,Fs,Eta_user);        %计算资源分配
[Rss, res_comu] = Rss_second_order_derivative(x,beta_time,beta_enengy,para,sub_bandNumber,serverNumber);

end
function [F,res_cra] = cra(x,Fs,Eta_user)

end
function [Rss, res_comu] = Rss_second_order_derivative(x,beta_time,beta_enengy,para,sub_bandNumber,serverNumber)

end
