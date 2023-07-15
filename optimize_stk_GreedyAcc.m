function [J1, X1i,X1o,X1c, F1, Rss_i, Rss_o] = optimize_stk_GreedyAcc(Fs, Tu, W, RssMax,...
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
%TA Task allocation,任务分配算法
[Xi,Xo,Xc] = genOriginXgreedy(userNumber, serverNumber,sub_bandNumber,para);    %得到初始X
[J, F,Rss_i, Rss_o] = RA(Xi,Xo,Xc,para);
end

function [Xi,Xo,Xc] = genOriginXgreedy(userNumber, serverNumber,sub_bandNumber,para)
%     根据para中接入H_ASL随机选择一个不为零的置为一
    Xi = zeros(userNumber, serverNumber,sub_bandNumber);
    Xo = zeros(userNumber, serverNumber,sub_bandNumber);
    Xc = zeros(userNumber, serverNumber);
    H_ASL = para.H_ASL;
    userP = para.userP;
    for user_in = 1:userNumber
        user_ASL_vec = H_ASL(user_in,:);
        user_out_ASL_vec = H_ASL(userP(user_in),:);
%         positiveIndices = find(user_ASL_vec > 0);
        sorted_vector = sort(user_ASL_vec, 'descend');
        sorted_vector_out = sort(user_out_ASL_vec, 'descend');
%         nth_largest_element = sorted_vector(n); % 获取第n大的元素
%         nth_largest_element_position = find(vector == nth_largest_element); % 找到第n大元素的位置
        flag = 0;n=1;
        while flag == 0
            nth_largest_element = sorted_vector(n);
            randomServer = find(sorted_vector == nth_largest_element);
            for band = 1:sub_bandNumber
%                 sumResult = sum(Xi(:, randomServer, band));
               if sum(Xi(:, randomServer, band)) == 0
                   Xi(user_in, randomServer, band) = 1;
                   Xc(user_in, randomServer) = 1;
                   flag = 1;
                   break
               end
            end
            n= n+1;
            if n > serverNumber
               break 
            end
        end
        
%         positiveIndices = find(user_out_ASL_vec > 0);
        flag = 0;n=1;
        while flag == 0
            nth_largest_element = sorted_vector_out(n);
            randomServer = find(sorted_vector_out == nth_largest_element);
            for band = 1:sub_bandNumber
               if sum(Xo(:, randomServer, band)) == 0
                   Xo(user_in, randomServer, band) = 1;
                   flag = 1;
                   break
               end
            end
            n = n+1;
            if n > serverNumber
               break 
            end
        end
%         randomCserverNumber = randi(serverNumber);  % 生成1到N区间内的随机数
%         Xc(user_in, randomCserverNumber) = 1;
    end
end