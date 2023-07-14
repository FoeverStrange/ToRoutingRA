function [J1, X1i,X1o,X1c, F1, Rss_i, Rss_o] = optimize_stk_GreedyAcc(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber, ...
    G,shrink,userP ...
    )
    %��װ����
    para.beta_time = beta_time;               %ʱ��/�ܺ�ƫ��
    para.beta_enengy = beta_enengy;
    para.Tu = Tu;                             %��������
    para.W = W;                               %����
    para.H_ASL = H_ASL;                       %ASL�ŵ��������
    para.Ttol = Ttol;                        %ʱ������Լ��
    para.H_ISL = H_ISL;                       %ISL�ŵ��������
    para.Ttol_S = Ttol_S;                       %ISLʱ������Լ��
    para.lamda = lamda;                       %�û����ȼ�����
    para.Sigma_square = Sigma_square;         %��������
    para.Fs = Fs;                             %������������������
    para.RssMax = RssMax;
    para.k = k; 
    para.G = G; 
    para.shrink = shrink;
    para.userP = userP;
    
    [J1, X1i,X1o,X1c, F1,Rss_i, Rss_o] = ta( ...
    userNumber,...              % �û�����
    serverNumber,...            % ����������
    sub_bandNumber,...          % �Ӵ�����
    para ...                    % �������
    );
end

function [J, Xi,Xo,Xc, F,Rss_i, Rss_o] = ta( ...
    userNumber,...              % �û�����
    serverNumber,...            % ����������
    sub_bandNumber,...          % �Ӵ�����
    para...                     % �������
)
%TA Task allocation,��������㷨
[Xi,Xo,Xc] = genOriginXgreedy(userNumber, serverNumber,sub_bandNumber,para);    %�õ���ʼX
[J, F,Rss_i, Rss_o] = RA(Xi,Xo,Xc,para);
end

function [Xi,Xo,Xc] = genOriginXgreedy(userNumber, serverNumber,sub_bandNumber,para)
%     ����para�н���H_ASL���ѡ��һ����Ϊ�����Ϊһ
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
%         nth_largest_element = sorted_vector(n); % ��ȡ��n���Ԫ��
%         nth_largest_element_position = find(vector == nth_largest_element); % �ҵ���n��Ԫ�ص�λ��
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
        end
%         randomCserverNumber = randi(serverNumber);  % ����1��N�����ڵ������
%         Xc(user_in, randomCserverNumber) = 1;
    end
end