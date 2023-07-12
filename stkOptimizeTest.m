clear;
close all;

sub_bandNumber = 4;          %子带个数

% 生成信道增益矩阵以及Ttol矩阵
% 卫星通信信道增益很差，因此考虑任务数据量较小但计算量很大的任务
pathStr = 'D:\sys\Resource_Allocation\STK\STK\Sc_PGSateNet\PGSateNet.sc';
num_begin = 1;
num_end = 8;
stk_aircraft_construct_func(pathStr,num_begin, num_end);
% [H_ASL,Ttol,H_ISL,Ttol_S] = stkIriGenGain(sub_bandNumber,pathStr);

% save UAVsat230712_improve.mat H_ASL Ttol H_ISL Ttol_S;
load("UAVsat230712_improve.mat")
% 设定：用户1和用户2是一对，用户3和用户4是一对……
[userNumber, serverNumber, ~ ] = size(H);

Fs = 40e10 * ones(serverNumber,1);  %服务器运算能力矩阵

T0.data = [];                      %任务数据大小
T0.circle = [];                    %任务所需时钟周期
Tu = repmat(T0,userNumber,1);
Tu_check = 1:userNumber;
for i = 1:userNumber
    if ismember(i,Tu_check)
        Tu(i).data = 10*8* 1024 * 8;
        Tu(i).circle = 10e9;
    else
        Tu(i).data = 0;
        Tu(i).circle = 0;
    end
end

% 用户优先级参数，全设定为1
lamda = ones(userNumber,1);
% 优化偏好，默认为0.5，即能耗偏好和时延偏好相同
beta_time = 0.5 * ones(userNumber,1);
beta_enengy = ones(userNumber,1) - beta_time;
% 星间链路速率
RssMax = 10 * 10^3;
% RssMin = 0;

Sigma_square = 1e-20;       %噪声方差
W = 100e6;   %系统带宽100MHz
k = 1 * 10^-26;  %芯片能耗系数

%测试算法质量
data_show = cell(5,3);
data_show(1,1) = {'Algorithm'};
data_show(1,2) = {'Objective'};
data_show(1,3) = {'Computing Time'};

disp('Rand Acc Computing')
tic;
[J1, X1, F1, Rss1] = optimize_stk_randAcc(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber ...
    );
Rand_time = toc;
Rand_objective = J1;
data_show(2,1) = {'optimize_stk_hJTORA'};
data_show(2,2) = {Rand_objective};
data_show(2,3) = {Rand_time};







