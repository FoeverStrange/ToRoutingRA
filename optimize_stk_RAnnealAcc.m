function [J1, X1i,X1o,X1c, F1, Rss_i, Rss_o] = optimize_stk_HRAcc(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber, ...
    G,shrink,userP ...
    )
    %封装参数
    para.beta_time = beta_time;               %时间/能耗偏好
    para.beta_enengy = beta_enengy;
    para.Tu = Tu;                             %任务数据
    para.W = W;                               %带宽
    para.H_ASL = H_ASL;                       %ASL信道增益矩阵
    para.Ttol = Ttol;                        %时变拓扑约束
    para.H_ISL = H_ISL;                       %ISL信道增益矩阵
    para.Ttol_S = Ttol_S;                       %ISL时变拓扑约束
    para.lamda = lamda;                       %用户优先级参数
    para.Sigma_square = Sigma_square;         %噪声方差
    para.Fs = Fs;                             %服务器运算能力矩阵
    para.RssMax = RssMax;
    para.k = k; 
    para.G = G; 
    para.shrink = shrink;
    para.userP = userP;
    
    [J1, X1i,X1o,X1c, F1,Rss_i, Rss_o] = ta( ...
    userNumber,...              % 用户个数
    serverNumber,...            % 服务器个数
    sub_bandNumber,...          % 子带个数
    para ...                    % 所需参数
    );
end

function [J, Xi,Xo,Xc, F,Rss_i, Rss_o] = ta( ...
    userNumber,...              % 用户个数
    serverNumber,...            % 服务器个数
    sub_bandNumber,...          % 子带个数
    para...                     % 所需参数
)
% 模拟退火参数初始化
T = userNumber/10;  %初始温度
T_min = 10^-1;
alpha = 0.9;
%TA Task allocation,任务分配算法
[Xi,Xo,Xc] = genOriginXH(userNumber, serverNumber,sub_bandNumber,para);    %得到初始X
[J, F,Rss_i, Rss_o] = RA(Xi,Xo,Xc,para);
% 迭代调整Xc，参考H算法
flaj = 0;
iterations = 1;
H_ASL = para.H_ASL;
T = userNumber/10;  %初始温度
T_min = 10^-1;
alpha = 0.84;
while(T>T_min)
    genUs = sum(Xc,1);
    genUs(genUs<2) = 0;
    probabilities = genUs / sum(genUs);
%     [~,maxServerp] = max(genUs);
%     indices = find(genUs(:,maxServerp) > 0);
maxServerp = find(rand <= cumsum(probabilities), 1, 'first');
indices = find(Xc(:,maxServerp) > 0);
% 找到最不重要的user
    user_temp = indices(end);
    Xi_temp = Xi;
    Xc_temp = Xc;
%     把这个user占用的子载波置零
    sub_band_temp = find(Xi_temp(user_temp,maxServerp,:)==1);
    Xi_temp(user_temp,maxServerp,sub_band_temp) = 0;
    Xc_temp(user_temp,maxServerp) = 0;

    user_ASL_vec = H_ASL(user_temp,:);
    sorted_vector = sort(user_ASL_vec, 'descend');
    serverNumTemp = sum(sorted_vector>0);
    Jm = J;
    Fm = F;
    Rss_im = Rss_i;
    Rss_om = Rss_o;
    for n = 1:serverNumTemp
        ASL_temp = sorted_vector(n);
%         找到ASL_temp 对应的server以及空闲subband
        mServer = find(user_ASL_vec == ASL_temp);
        for band = 1:sub_bandNumber
%                 sumResult = sum(Xi(:, randomServer, band));
               if sum(Xi(:, mServer, band)) == 0
                   Xi_temp1 = Xi_temp;
                   Xc_temp1 = Xc_temp;
                   Xi_temp1(user_temp, mServer, band) = 1;
                   Xc_temp1(user_temp, mServer) = 1;
                   [J_t, F_t,Rss_i_t, Rss_o_t,res_cra,res_comu] = RA(Xi_temp1,Xo,Xc_temp1,para);
                   if J_t < Jm*(1-0.00001)
                        Jm = J_t; 
                        Fm = F_t;
                        Rss_im = Rss_i_t;
                        Rss_om = Rss_o_t;
                   end
                   break
               end
        end
    end
    if Jm < J
       J = Jm;
       F = Fm;
       Rss_i = Rss_im;
       Rss_o = Rss_om;
    else
        flaj = flaj + 1;
    end
    picture(iterations,1) = J;
    iterations = iterations +1;
    T = T * alpha;
end
%     figure()
%     hold on
%     plot(picture(:,1));
%     hold off
end
function p = getProbability(delta,t)
    p = exp(delta/t);
end
function [Xi,Xo,Xc] = genOriginXH(userNumber, serverNumber,sub_bandNumber,para)
% 通过路由决定Xi,Xo，并认为Xc为固定
G = para.G;
userP = para.userP;
    Xi = zeros(userNumber, serverNumber,sub_bandNumber);
    Xo = zeros(userNumber, serverNumber,sub_bandNumber);
    Xc = zeros(userNumber, serverNumber);
    Xo_bandtemp = zeros(userNumber, serverNumber,sub_bandNumber);
for user = 1:userNumber
    user_out = userP(user);
    flaj1 = 0;
    while flaj1 ==0
        [P,~] = shortestpath(G,serverNumber+user,serverNumber+user_out);
        randomServer_in = P(2);
        randomServer_out = P(end-1);
        user_in = user;
        band_in_full = 0;
        band_out_full = 0;
        for band_in = 1:sub_bandNumber
%                 sumResult = sum(Xi(:, randomServer, band));
           if sum(Xi(:, randomServer_in, band_in)) == 0 && flaj1 == 0
               for band_out = 1:sub_bandNumber
                   if sum(Xo_bandtemp(:, randomServer_out, band_out)) == 0 && flaj1 == 0
                       Xi(user_in, randomServer_in, band_in) = 1;
                       Xc(user_in, randomServer_in) = 1;
                       Xo_bandtemp(user_out, randomServer_out, band_out) = 1;
                       Xo(user_in, randomServer_out, band_out) = 1;
                       flaj1 = 1;
                       break
                   end
                   if band_out == sub_bandNumber && flaj1 == 0
                      band_out_full = 1; 
                   end
               end
               if flaj1 == 1
                  break 
               end
           end
           
           if band_in == sub_bandNumber && flaj1 == 0
              band_in_full = 1; 
           end
        end
        if flaj1 == 0
            if band_in_full == 1
                G = rmedge(G, serverNumber+user, randomServer_in);
            end
            if band_out_full == 1
                G = rmedge(G, serverNumber+user_out, randomServer_out);
            end
        end
    end
    
end

end