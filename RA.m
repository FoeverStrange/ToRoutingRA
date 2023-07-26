function [Jx, F, Rss_i, Rss_o,res_cra,res_comu] = RA(Xi,Xo,Xc,para)

[userNumber, serverNumber, subbandNumber] = size(Xi);
G = para.G;
%     RA 先计算router，再计算资源分配
%     Jx = 0;
%     F = 0;
%     Rss = 0;
% 最短路径计算
Xi_temp = sum(Xi, 3);
Xo_temp = sum(Xo, 3);
% 最短路径，sum(1/H)存储
H_sum = zeros(userNumber, 4);
H_ASL = para.H_ASL;
userP = para.userP;
PANISH = 1*10^-17;
for user = 1:userNumber
    server_in = find(Xi_temp(user, :) == 1);
    server_out = find(Xo_temp(user, :) == 1);
    server_cp = find(Xc(user, :) == 1);
    if isempty(server_in)
        H_sum(user,1) = 2;
        H_sum(user,2) = 1/PANISH;
        H_sum(user,3) = 2;
        H_sum(user,4) = 1/PANISH;
        continue
    end
     if isempty(server_out)
        H_sum(user,1) = 5;
        H_sum(user,2) = 1/PANISH;
        H_sum(user,3) = 5;
        H_sum(user,4) = 1/PANISH;
        continue
    end
    user_out = userP(user);
    temp1 = H_ASL(user,server_in);
    temp2 = H_ASL(user_out,server_out);
    if temp1 == 0
        Hasl_temp = 1/PANISH;
    else
        Hasl_temp = 1/temp1;
    end
    if temp2 == 0
        Has2_temp = 1/PANISH;
    else
%         temp2
%         user_out
%         server_out
        Has2_temp = 1/temp2;
    end
%     Hasl_temp = min(1/temp1,1/PANISH);
%     Has2_temp = min(1/temp2,1/PANISH);
    
    [P,d] = shortestpath(G,server_in,server_cp);
    H_sum(user,1) = length(P);
    H_sum(user,2) = d+Hasl_temp;
    [P,d] = shortestpath(G,server_cp,server_out);
    H_sum(user,3) = length(P);
    H_sum(user,4) = d+Has2_temp;
end

[F,res_cra] = cra(Xc,para);        %计算资源分配
% 通信资源分配
[Rss_i, Rss_o, res_comu] = Rss_second_order_derivative(H_sum,para); 

Jx = res_cra + res_comu;
end
function [F,res_cra] = cra(Xc,para)
    res_cra = 0;
    lamda = para.lamda;
    fs = para.Fs(1);
    beta_time = para.beta_time;
    beta_enengy = para.beta_enengy;
    k = para.k;
    Tu = para.Tu;
    H_ASL = para.H_ASL;
    [userNumber, serverNumber,~] = size(H_ASL);
    F = zeros(userNumber,serverNumber);
    for server = 1:serverNumber
        %     对每个服务器分别优化
        X_server_s = Xc(:,server);
% 找到一个服务器上卸载的所有任务
        positiveElements = X_server_s > 0;
        Xu = find(positiveElements);
        if isempty(Xu)
           continue 
        end
        n = length(Xu);
%         Fs = zeros(userNumber,1);
%         Au = zeros(n,1);
%         Bu = zeros(n,1);
        psi = zeros(n,1);
        for user_p = 1:n
            user = Xu(user_p);
%            Au(user_p) =  beta_time(user)*Tu(user).circle * lamda(user);
%            Bu(user_p) = 2*beta_enengy(user) * k * Tu(user).circle* lamda(user);
           psi(user_p) = Tu(user).circle* lamda(user) *(beta_time(user)* + k*beta_enengy(user)*fs^2);
        end
        psi_sq = sqrt(psi);
        fus_o = psi_sq*fs/sum(psi_sq);
        
        res_cra_temp = psi./fus_o ;
%         res_cra_serverr = sum(res_cra_temp.*lamda(Xu));
        res_cra = res_cra + sum(res_cra_temp);
    end
    
end


function [Rss_i_out, Rss_o_out, res_comu] = Rss_second_order_derivative(H_sum,para)
res_comu = 0;
    lamda = para.lamda;
    [userNumber, ~] = size(H_sum);
    Tu = para.Tu;
    Sigma_square = para.Sigma_square;
    beta_enengy = para.beta_enengy;
    beta_time = para.beta_time;
    RssMax = para.RssMax;
    W = para.W;
    n = 5;
    shrink = para.shrink;
    for user = 1:userNumber
        du = Tu(user).data;
        [phi_i, phie_i,phi_o, phie_o] = getPhie(H_sum, beta_enengy, beta_time, Sigma_square, du,shrink,user);
%         入
        if getPi_s(RssMax, phi_i, phie_i, W) <=0
            Rss_i = RssMax;
        else
            Rss_r = RssMax;
            Rss_l = 0;
            for i = 1:n
                Rss_mid = (Rss_r+Rss_l)/2;
                Pi_s_temp = getPi_s(Rss_mid, phi_i, phie_i, W);
                if Pi_s_temp >0
                    Rss_r = Rss_mid;
                else
                    Rss_l = Rss_mid;
                end
            end
            Rss_i = (Rss_r+Rss_l)/2;
        end
        
        
        if getPi_s(RssMax, phi_o, phie_o, W) <=0
            Rss_o = RssMax;
        else
            Rss_r = RssMax;
            Rss_l = 0;
            for i = 1:n
                Rss_mid = (Rss_r+Rss_l)/2;
                Pi_s_temp = getPi_s(Rss_mid, phi_i, phie_i, W);
                if Pi_s_temp >0
                    Rss_r = Rss_mid;
                else
                    Rss_l = Rss_mid;
                end
            end
            Rss_o = (Rss_r+Rss_l)/2;
        end
        res_comu_temp = (phi_i + phie_i*2^(Rss_i/W))/Rss_i + (phi_o + phie_o*2^(Rss_o/W))/Rss_o;
        Rss_i_out(user) = Rss_i;
        Rss_o_out(user) = Rss_o;
        res_comu = res_comu + res_comu_temp*lamda(user);
    end

end

function [phi_i, phie_i,phi_o, phie_o] = getPhie(H_sum, beta_enengy, beta_time, sigma, du,shrink,user)
    phi_i = beta_enengy(user) * sigma * H_sum(user,2) * du;
    phie_i = H_sum(user,1)*beta_time(user)*du - phi_i;
    du = du/shrink;
    phi_o = beta_enengy(user) * sigma * H_sum(user,4)*du;
    phie_o = H_sum(user,3)*beta_time(user)*du - phi_o;
end
function Pi_s = getPi_s(R, phi, phie, W)
    Pi_s = R * (phie * 2^(R/W)) - phi - phie*2^(R/W);
end