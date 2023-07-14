clear;
close all;

sub_bandNumber = 4;          %�Ӵ�����

% �����ŵ���������Լ�Ttol����
% ����ͨ���ŵ�����ܲ��˿���������������С���������ܴ������
pathStr = 'D:\sys\Resource_Allocation\STK\STK\Sc_PGSateNet\PGSateNet.sc';
% num_begin = 1;
% num_end = 8;
% stk_aircraft_construct_func(pathStr,num_begin, num_end);

% [H_ASL,Ttol,H_ISL,Ttol_S] = stkIriGenGain(sub_bandNumber,pathStr);

% save UAVsat230712_improve.mat H_ASL Ttol H_ISL Ttol_S;
load("UAVsat230712_improve.mat")

% ͼ����
% ��ȡ�����ά��
n = size(H_ISL, 1);
% ����һ���յļ�Ȩͼ����
G = graph();
% ��ӽڵ㵽ͼ��
G = addnode(G, n);
% ��������Ԫ�أ�������0��Ԫ�����Ϊ�ߣ������ö�Ӧ��Ȩ��
for i = 2:n
    for j = 1:i-1
        if H_ISL(i, j) > 0
            G = addedge(G, i, j, 1/H_ISL(i, j));
        end
    end
end
% �����Ȩͼ����
% disp(G);
% plot(G,'EdgeLabel',G.Edges.Weight)

% �趨���û�1���û�2��һ�ԣ��û�3���û�4��һ�ԡ���
[userNumber, serverNumber, ~ ] = size(H_ASL);

Fs = 4e8 * ones(serverNumber,1);  %������������������

T0.data = [];                      %�������ݴ�С
T0.circle = [];                    %��������ʱ������
Tu = repmat(T0,userNumber,1);
Tu_check = 1:userNumber;
for i = 1:userNumber
    if ismember(i,Tu_check)
        Tu(i).data = 10*8* 1024 * 8;
        Tu(i).circle = 1000*10e9;
    else
        Tu(i).data = 0;
        Tu(i).circle = 0;
    end
end

% �û����ȼ�������ȫ�趨Ϊ1
% lamda = ones(userNumber,1);
lamda = rand(userNumber, 1) * 0.8 + 0.2; % ������0.2��1֮������С��
lamda = sort(lamda, 'descend'); % ������������
% �Ż�ƫ�ã�Ĭ��Ϊ0.5�����ܺ�ƫ�ú�ʱ��ƫ����ͬ
beta_time = 0.5 * ones(userNumber,1);
beta_enengy = ones(userNumber,1) - beta_time;
% �Ǽ���·����
RssMax = 10 * 10^3; %10Mb
% RssMin = 0;

Sigma_square = 1e-20;       %��������
W = 100e6;   %ϵͳ����100MHz
k = 1 * 10^-26;  %оƬ�ܺ�ϵ��
shrink = 10;
userP = 1:userNumber;
for user = 1:userNumber
    userP(user) = randi(userNumber);
end


%�����㷨����
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
data_show(4,1) = {'GreedyAcc_��С�������'};
data_show(4,2) = {Greedy_objective};
data_show(4,3) = {Greedy_time};

disp('H Acc Computing')
tic;
[J4, X4i,X4o,X4c, F4, Rss4_i, Rss4_o] = optimize_stk_HAcc(Fs, Tu, W, RssMax,...
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

disp('Mission Compeleted')







