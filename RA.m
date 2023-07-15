function [Jx, F, Rss_i, Rss_o] = RA(Xi,Xo,Xc,para)

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
PANISH = 1*10^-15;
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
    Hasl_temp = min(1/temp1,1/PANISH);
    Has2_temp = min(1/temp2,1/PANISH);
    
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
        Fs = zeros(userNumber,1);
        Au = zeros(n,1);
        Bu = zeros(n,1);
        for user_p = 1:n
            user = Xu(user_p);
           Au(user_p) =  beta_time(user)*Tu(user).circle;
           Bu(user_p) = beta_enengy(user) * k * Tu(user).circle;
        end
        
        fus_o = (Au./Bu).^(1/3);
        if sum(fus_o) <= fs
           Fs(Xu) = fus_o;
           F(:,server) =  Fs;
        else   %解kkt
           % 初始估计
X0 = [ones(n, 1)*fs/n;1];  % 初始估计，可以根据实际情况修改
% 使用fsolve求解
X = fsolve(@(X)KKTfun(X, Au, Bu, fs), X0);
% X就是你的解，其中前n个元素是F，最后一个元素是v。
fus_o = X(1:n);
v = X(n+1);
            Fs(Xu) = fus_o;
            F(:,server) =  Fs;
        end
        res_cra = sum(Au./fus_o + Bu .* fus_o .^2);
        
    end

% 解kkt条件
% 情况1：v = 0,所有fus的和小于fs
% 情况2：v>0,所有fus的和=fs，解一元三次方程组
    
end
function y = KKTfun(X, Au, Bu, fs)
    n = length(X) - 1; % 减1是因为向量X的最后一个元素是v
    F = X(1:n);  % 前n个元素是F
    v = X(n+1);  % 最后一个元素是v
    y = zeros(n + 1, 1); % Initialize y
    y(1:n) = Bu .* (F.^3) + v * (F.^2) - Au;  % Bu.*F.^3+v.*F.^2=Au
    y(n+1) = sum(F) - fs;  % ones(1,n)*F=fs
end

function [Rss_i, Rss_o, res_comu] = Rss_second_order_derivative(H_sum,para)
res_comu = 0;
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
        res_comu = res_comu + (phi_i + phie_i*2^(Rss_i/W))/Rss_i + (phi_o + phie_o*2^(Rss_o/W))/Rss_o;
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