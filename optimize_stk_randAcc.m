function [J1, X1i,X1o,X1c, F1, Rss1] = optimize_stk_randAcc(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber ...
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
    
    [J1, X1i,X1o,X1c, F1,Rss1] = ta( ...
    userNumber,...              % 用户个数
    serverNumber,...            % 服务器个数
    sub_bandNumber,...          % 子带个数
    para ...                    % 所需参数
    );
end
function [J, Xi,Xo,Xc, F,Rss] = ta( ...
    userNumber,...              % 用户个数
    serverNumber,...            % 服务器个数
    sub_bandNumber,...          % 子带个数
    para...                     % 所需参数
)
%TA Task allocation,任务分配算法
X = genOriginX(userNumber, serverNumber,sub_bandNumber,para);    %得到初始X
[J, F,Rss] = RA(X,para);

end