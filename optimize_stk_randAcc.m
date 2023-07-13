function [J1, X1i,X1o,X1c, F1, Rss1] = optimize_stk_randAcc(Fs, Tu, W, RssMax,...
    H_ASL, Ttol, H_ISL, Ttol_S,...
    lamda, Sigma_square, beta_time, beta_enengy,...
    k,...
    userNumber, serverNumber, sub_bandNumber, ...
    G ...
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
    
    [J1, X1i,X1o,X1c, F1,Rss1] = ta( ...
    userNumber,...              % �û�����
    serverNumber,...            % ����������
    sub_bandNumber,...          % �Ӵ�����
    para ...                    % �������
    );
end
function [J, Xi,Xo,Xc, F,Rss] = ta( ...
    userNumber,...              % �û�����
    serverNumber,...            % ����������
    sub_bandNumber,...          % �Ӵ�����
    para...                     % �������
)
%TA Task allocation,��������㷨
[Xi,Xo,Xc] = genOriginX(userNumber, serverNumber,sub_bandNumber,para);    %�õ���ʼX
[J, F,Rss] = RA(Xi,Xo,Xc,para);
end

function [Xi,Xo,Xc] = genOriginX(userNumber, serverNumber,sub_bandNumber,para)
%     ����para�н���H_ASL���ѡ��һ����Ϊ�����Ϊһ
    Xi = zeros(userNumber, serverNumber,sub_bandNumber);
    Xo = zeros(userNumber, serverNumber,sub_bandNumber);
    Xc = zeros(userNumber, serverNumber);
    H_ASL = para.H_ASL;
    for user_in = 1:userNumber
        user_ASL_vec = H_ASL(user_in,:);
        user_out_ASL_vec = H_ASL(mod(user_in*2,userNumber)+1,:);
        positiveIndices = find(user_ASL_vec > 0);
        flag = 0;
        while flag == 0
            randomServer = positiveIndices(randi(numel(positiveIndices)));
            for band = 1:sub_bandNumber
%                 sumResult = sum(Xi(:, randomServer, band));
               if sum(Xi(:, randomServer, band)) == 0
                   Xi(user_in, randomServer, band) = 1;
                   flag = 1;
                   break
               end
            end
        end
        positiveIndices = find(user_out_ASL_vec > 0);
        flag = 0;
        while flag == 0
            randomServer = positiveIndices(randi(numel(positiveIndices)));
            for band = 1:sub_bandNumber
               if sum(Xo(:, randomServer, band)) == 0
                   Xo(user_in, randomServer, band) = 1;
                   flag = 1;
                   break
               end
            end
        end
        randomCserverNumber = randi(serverNumber);  % ����1��N�����ڵ������
        Xc(user_in, randomCserverNumber) = 1;
    end
    
    

end