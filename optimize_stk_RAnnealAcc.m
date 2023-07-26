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
iterations = 1;
while(T_min<T)
   genUs = sum(Xc,1);
   probabilities = genUs / sum(genUs);
   maxServerp = find(rand <= cumsum(probabilities), 1, 'first');
   indices = find(Xc(:,maxServerp) > 0);
   user_temp = indices(end);
   Xc_temp = Xc;
   Xc_temp(user_temp,maxServerp) = 0;
   Jm = J*(1+0.1);
   Fm = F;
   Rss_im = Rss_i;
   Rss_om = Rss_o;
   for server_temp = 1:serverNumber
       Xc_temp1 = Xc_temp;
       Xc_temp1(user_temp,server_temp) = 1;
       [J_t, F_t,Rss_i_t, Rss_o_t] = RA(Xi,Xo,Xc_temp1,para);
       if J_t < Jm
          Jm = J_t; 
          Fm = F_t;
          Rss_im = Rss_i_t;
          Rss_om = Rss_o_t;
       end
   end
   if Jm < J
       J = Jm;
       F = Fm;
       Rss_i = Rss_im;
       Rss_o = Rss_om;
   else
        delta = J-Jm;
        pro=getProbability(delta,T);
        if(pro>rand)
            J = Jm;
            F = Fm;
            Rss_i = Rss_im;
            Rss_o = Rss_om;
        end
%         flaj = flaj + 1;
   end
    T = T * alpha;
    picture(iterations,1) = J;
    picture(iterations,2) = T;
    iterations = iterations +1;
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