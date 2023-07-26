function [J1, X1i,X1o,X1c, F1, Rss_i, Rss_o] = optimize_stk_fixAcc(Fs, Tu, W, RssMax,...
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
[Xi,Xo,Xc] = genOriginXFix(userNumber, serverNumber,sub_bandNumber,para);    %�õ���ʼX
[J, F,Rss_i, Rss_o] = RA(Xi,Xo,Xc,para);
end
function [Xi,Xo,Xc] = genOriginXFix(userNumber, serverNumber,sub_bandNumber,para)
%     ����para�н���H_ASL���ѡ��һ����Ϊ�����Ϊһ
    Xi = zeros(userNumber, serverNumber,sub_bandNumber);
    Xo = zeros(userNumber, serverNumber,sub_bandNumber);
    Xc = zeros(userNumber, serverNumber);
    H_ASL = para.H_ASL;
    userP = para.userP;
    for user_in = 1:userNumber
        user_ASL_vec = H_ASL(user_in,:);
        user_out_ASL_vec = H_ASL(userP(user_in),:);
        positiveIndices = find(user_ASL_vec > 0);
        flag = 0;
        cnt = 0;
        while flag == 0
            cnt = cnt +1;
            if cnt > 5
               break 
            end
            randomServer = positiveIndices(randi(numel(positiveIndices)));
            for band = 1:sub_bandNumber
%                 sumResult = sum(Xi(:, randomServer, band));
               if sum(Xi(:, randomServer, band)) == 0
                   Xi(user_in, randomServer, band) = 1;
                   Xc(user_in, randomServer) = 1;
                   flag = 1;
                   break
               end
            end
        end
        positiveIndices = find(user_out_ASL_vec > 0);
        flag = 0;cnt = 0;
        while flag == 0
            cnt = cnt +1;
            if cnt > 5
               break 
            end
            randomServer = positiveIndices(randi(numel(positiveIndices)));
            for band = 1:sub_bandNumber
               if sum(Xo(:, randomServer, band)) == 0
                   Xo(user_in, randomServer, band) = 1;
                   flag = 1;
                   break
               end
            end
        end
%         randomCserverNumber = randi(serverNumber);  % ����1��N�����ڵ������
%         Xc(user_in, randomCserverNumber) = 1;
    end
end
