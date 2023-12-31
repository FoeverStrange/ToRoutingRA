clear;
close all;

sub_bandNumber = 4;          %子带个数

% 生成信道增益矩阵以及Ttol矩阵
% 卫星通信信道增益很差，因此考虑任务数据量较小但计算量很大的任务
pathStr = 'D:\sys\Resource_Allocation\STK\STK\Sc_PGSateNet\PGSateNet.sc';
% num_begin = 1;
% num_end = 15;
% stk_aircraft_construct_func(pathStr,num_begin, num_end);

% [H_ASL,Ttol,H_ISL,Ttol_S] = stkIriGenGain(sub_bandNumber,pathStr);

% save UAVsat230715_improve.mat H_ASL Ttol H_ISL Ttol_S;
load("UAVsat230715_improve.mat")

% 图生成
% 获取矩阵的维度
n = size(H_ISL, 1);
% 创建一个空的加权图对象
G = graph();
% 添加节点到图中
G = addnode(G, n);
% 遍历矩阵元素，将大于0的元素添加为边，并设置对应的权重
for i = 2:n
    for j = 1:i-1
        if H_ISL(i, j) > 0
            G = addedge(G, i, j, 1/H_ISL(i, j));
        end
    end
end
[n_user,n_server] = size(H_ASL);
G = addnode(G, n_user);
for i = 1:n_user
    for j = 1:n_server
        if H_ASL(i, j) > 0
            G = addedge(G, n_server+i, j, 1/H_ASL(i, j));
        end
    end
end

% 输出加权图对象
% disp(G);
% plot(G,'EdgeLabel',G.Edges.Weight)

% 设定：用户1和用户2是一对，用户3和用户4是一对……
[userNumber, serverNumber, ~ ] = size(H_ASL);

Fs = 1e8 * ones(serverNumber,1);  %服务器运算能力矩阵

T0.data = [];                      %任务数据大小
T0.circle = [];                    %任务所需时钟周期
Tu = repmat(T0,userNumber,1);
Tu_check = 1:userNumber;
for i = 1:userNumber
    if ismember(i,Tu_check)
        Tu(i).data = 1e3;
        Tu(i).circle = 3000*1e6;
    else
        Tu(i).data = 0;
        Tu(i).circle = 0;
    end
end

% 用户优先级参数，全设定为1
% lamda = ones(userNumber,1);
lamda = rand(userNumber, 1) * 0.8 + 0.2; % 生成在0.2到1之间的随机小数
lamda = sort(lamda, 'descend'); % 降序排列向量
% 优化偏好，默认为0.5，即能耗偏好和时延偏好相同
beta_time = 0.5 * ones(userNumber,1);
beta_enengy = ones(userNumber,1) - beta_time;
% 星间链路速率
RssMax = 10 * 10^3; %10Mb
% RssMin = 0;

Sigma_square = 1e-20;       %噪声方差
W = 100e6;   %系统带宽100MHz
k = 1 * 10^-16;  %芯片能耗系数
shrink = 3;
userP = 1:userNumber;
for user = 1:userNumber
    userP(user) = randi(userNumber);
    if userP(user) == user
        userP(user) = mod(user+1,userNumber)+1;
    end
end


%测试算法质量
data_show = cell(5,3);
data_show(1,1) = {'Algorithm'};
data_show(1,2) = {'Objective'};
data_show(1,3) = {'Computing Time'};

disp('Rand Acc Computing')
tic;
[J1, X1i,X1o,X1c, F1, Rss1_i, Rss1_o] = optimize_stk_randAcc(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber, ...
    G, shrink, userP...
    );
Rand_time = toc;
Rand_objective = J1;
data_show(2,1) = {'RandAcc'};
data_show(2,2) = {Rand_objective};
data_show(2,3) = {Rand_time};

disp('Fix Acc Computing')
tic;
[J2, X2i,X2o,X2c, F2, Rss2_i, Rss2_o] = optimize_stk_fixAcc(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber, ...
    G, shrink, userP...
    );
Fix_time = toc;
Fix_objective = J2;
data_show(3,1) = {'FixAcc'};
data_show(3,2) = {Fix_objective};
data_show(3,3) = {Fix_time};

disp('Greedy Acc Computing')
tic;
[J3, X3i,X3o,X3c, F3, Rss3_i, Rss3_o] = optimize_stk_GreedyAcc(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber, ...
    G, shrink, userP...
    );
Greedy_time = toc;
Greedy_objective = J3;
data_show(4,1) = {'GreedyAcc_最小接入代价'};
data_show(4,2) = {Greedy_objective};
data_show(4,3) = {Greedy_time};

disp('H Acc Computing')
tic;
[J4, X4i,X4o,X4c, F4, Rss4_i, Rss4_o,res_cra,res_comu] = optimize_stk_HAcc(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber, ...
    G, shrink, userP...
    );
H_time = toc;
H_objective = J4;
data_show(5,1) = {'H'};
data_show(5,2) = {H_objective};
data_show(5,3) = {H_time};

disp('HR Acc Computing')
tic;
[J5, X5i,X5o,X5c, F5, Rss5_i, Rss5_o] = optimize_stk_HRAcc(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber, ...
    G, shrink, userP...
    );
HR_time = toc;
HR_objective = J5;
data_show(6,1) = {'RouterH'};
data_show(6,2) = {HR_objective};
data_show(6,3) = {HR_time};

disp('RAnneal Acc Computing')
tic;
[J6, X6i,X6o,X6c, F6, Rss6_i, Rss6_o] = optimize_stk_RAnnealAcc(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber, ...
    G, shrink, userP...
    );
RAnneal_time = toc;
RAnneal_objective = J6;
data_show(7,1) = {'RAnneal'};
data_show(7,2) = {RAnneal_objective};
data_show(7,3) = {RAnneal_time};

disp('HAcc2 Acc Computing')
tic;
[J7, X7i,X7o,X7c, F7, Rss7_i, Rss7_o,res_cra,res_comu] = optimize_stk_HAcc2(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber, ...
    G, shrink, userP...
    );
H2_time = toc;
H2_objective = J7;
data_show(8,1) = {'H2'};
data_show(8,2) = {H2_objective};
data_show(8,3) = {H2_time};

disp('HAcc3 Acc Computing')
tic;
[J8, X8i,X8o,X8c, F8, Rss8_i, Rss8_o,res_cra,res_comu] = optimize_stk_HAcc3(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber, ...
    G, shrink, userP...
    );
H8_time = toc;
H8_objective = J8;
data_show(9,1) = {'H3'};
data_show(9,2) = {H8_objective};
data_show(9,3) = {H8_time};

disp('Mission Compeleted')







