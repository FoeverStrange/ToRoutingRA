function [J1, X1i,X1o,X1c, F1, Rss_i, Rss_o,res_cra,res_comu] = optimize_stk_HAcc2(Fs, Tu, W, RssMax,...
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
    
    [J1, X1i,X1o,X1c, F1,Rss_i, Rss_o,res_cra,res_comu] = ta( ...
    userNumber,...              % �û�����
    serverNumber,...            % ����������
    sub_bandNumber,...          % �Ӵ�����
    para ...                    % �������
    );
end

function [J, Xi,Xo,Xc, F,Rss_i, Rss_o,res_cra,res_comu] = ta( ...
    userNumber,...              % �û�����
    serverNumber,...            % ����������
    sub_bandNumber,...          % �Ӵ�����
    para...                     % �������
)
%TA Task allocation,��������㷨
[Xi,Xo,Xc] = genOriginXFix(userNumber, serverNumber,sub_bandNumber,para);    %�õ���ʼX
[J, F,Rss_i, Rss_o,res_cra,res_comu] = RA(Xi,Xo,Xc,para);
flaj = 0;
iterations = 1;
H_ASL = para.H_ASL;
while(flaj<3)
    genUs = sum(Xc,1);
    genUs(genUs<2) = 0;
    probabilities = genUs / sum(genUs);
%     [~,maxServerp] = max(genUs);
%     indices = find(genUs(:,maxServerp) > 0);
maxServerp = find(rand <= cumsum(probabilities), 1, 'first');
indices = find(Xc(:,maxServerp) > 0);
% �ҵ����Ҫ��user
    user_temp = indices(end);
    Xi_temp = Xi;
    Xc_temp = Xc;
%     �����userռ�õ����ز�����
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
%         �ҵ�ASL_temp ��Ӧ��server�Լ�����subband
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
end
    figure()
    plot(1:iterations-1,picture(:,1));
end

function [Xi,Xo,Xc] = genOriginXH(userNumber, serverNumber,sub_bandNumber,para)
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
            randomServer = find(user_ASL_vec == nth_largest_element);
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
            randomServer = find(user_out_ASL_vec == nth_largest_element);
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
%         randomCserverNumber = randi(serverNumber);  % ����1��N�����ڵ������
%         Xc(user_in, randomCserverNumber) = 1;
    end
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
        while flag == 0
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
%         randomCserverNumber = randi(serverNumber);  % ����1��N�����ڵ������
%         Xc(user_in, randomCserverNumber) = 1;
    end
end