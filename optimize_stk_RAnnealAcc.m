function [J1, X1i,X1o,X1c, F1, Rss_i, Rss_o] = optimize_stk_HRAcc(Fs, Tu, W, RssMax,...
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
% ģ���˻������ʼ��
T = userNumber/10;  %��ʼ�¶�
T_min = 10^-1;
alpha = 0.9;
%TA Task allocation,��������㷨
[Xi,Xo,Xc] = genOriginXH(userNumber, serverNumber,sub_bandNumber,para);    %�õ���ʼX
[J, F,Rss_i, Rss_o] = RA(Xi,Xo,Xc,para);
% ��������Xc���ο�H�㷨
flaj = 0;
iterations = 1;
H_ASL = para.H_ASL;
T = userNumber/10;  %��ʼ�¶�
T_min = 10^-1;
alpha = 0.84;
while(T>T_min)
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
    T = T * alpha;
end
%     figure()
%     hold on
%     plot(picture(:,1));
%     hold off
end
function p = getProbability(delta,t)
    p = exp(delta/t);
end
function [Xi,Xo,Xc] = genOriginXH(userNumber, serverNumber,sub_bandNumber,para)
% ͨ��·�ɾ���Xi,Xo������ΪXcΪ�̶�
G = para.G;
userP = para.userP;
    Xi = zeros(userNumber, serverNumber,sub_bandNumber);
    Xo = zeros(userNumber, serverNumber,sub_bandNumber);
    Xc = zeros(userNumber, serverNumber);
    Xo_bandtemp = zeros(userNumber, serverNumber,sub_bandNumber);
for user = 1:userNumber
    user_out = userP(user);
    flaj1 = 0;
    while flaj1 ==0
        [P,~] = shortestpath(G,serverNumber+user,serverNumber+user_out);
        randomServer_in = P(2);
        randomServer_out = P(end-1);
        user_in = user;
        band_in_full = 0;
        band_out_full = 0;
        for band_in = 1:sub_bandNumber
%                 sumResult = sum(Xi(:, randomServer, band));
           if sum(Xi(:, randomServer_in, band_in)) == 0 && flaj1 == 0
               for band_out = 1:sub_bandNumber
                   if sum(Xo_bandtemp(:, randomServer_out, band_out)) == 0 && flaj1 == 0
                       Xi(user_in, randomServer_in, band_in) = 1;
                       Xc(user_in, randomServer_in) = 1;
                       Xo_bandtemp(user_out, randomServer_out, band_out) = 1;
                       Xo(user_in, randomServer_out, band_out) = 1;
                       flaj1 = 1;
                       break
                   end
                   if band_out == sub_bandNumber && flaj1 == 0
                      band_out_full = 1; 
                   end
               end
               if flaj1 == 1
                  break 
               end
           end
           
           if band_in == sub_bandNumber && flaj1 == 0
              band_in_full = 1; 
           end
        end
        if flaj1 == 0
            if band_in_full == 1
                G = rmedge(G, serverNumber+user, randomServer_in);
            end
            if band_out_full == 1
                G = rmedge(G, serverNumber+user_out, randomServer_out);
            end
        end
    end
    
end

end