clear;
close all;

sub_bandNumber = 4;          %�Ӵ�����

% �����ŵ���������Լ�Ttol����
% ����ͨ���ŵ�����ܲ��˿���������������С���������ܴ������
pathStr = 'D:\sys\Resource_Allocation\STK\STK\Sc_PGSateNet\PGSateNet.sc';
num_begin = 1;
num_end = 8;
stk_aircraft_construct_func(pathStr,num_begin, num_end);
% [H_ASL,Ttol,H_ISL,Ttol_S] = stkIriGenGain(sub_bandNumber,pathStr);

% save UAVsat230712_improve.mat H_ASL Ttol H_ISL Ttol_S;
load("UAVsat230712_improve.mat")
% �趨���û�1���û�2��һ�ԣ��û�3���û�4��һ�ԡ���
[userNumber, serverNumber, ~ ] = size(H);

Fs = 40e10 * ones(serverNumber,1);  %������������������

T0.data = [];                      %�������ݴ�С
T0.circle = [];                    %��������ʱ������
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

% �û����ȼ�������ȫ�趨Ϊ1
lamda = ones(userNumber,1);
% �Ż�ƫ�ã�Ĭ��Ϊ0.5�����ܺ�ƫ�ú�ʱ��ƫ����ͬ
beta_time = 0.5 * ones(userNumber,1);
beta_enengy = ones(userNumber,1) - beta_time;
% �Ǽ���·����
RssMax = 10 * 10^3;
% RssMin = 0;

Sigma_square = 1e-20;       %��������
W = 100e6;   %ϵͳ����100MHz
k = 1 * 10^-26;  %оƬ�ܺ�ϵ��

%�����㷨����
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







