clear all; 
close all;
M = 61;
SNR = -30:5:0;
snr  = 10.^(SNR/10);
N = 256;                                                                   %ofdm调制把总带宽分配成了N个子信道
fc = 60 * 10^9;
W = 0.02*fc;
T = N/W;
B = W/(N + 1);
N0 = 1;
D = N/4;
Ts = 256/W/D;
fs = W;
totolP = N0*snr*N;
C = 100;     
a = 1/M;
b_geshu = 4;                                                               %波束个数
c = 3.0*10^8;
wavelength = fc/c;
d = wavelength/2;
%% ******************计算信干噪比时用到*********************
S_BLL = zeros(1,N);
Lamda = zeros(1,N);
BL = zeros(N,N);
%% *****************************************
fl = zeros(1,N);                                                           %载波
i_power = zeros(1,M);                                                      %每个波束的能量
LL = 8; %路径数                                                                   %路径数
W1 = zeros(1,b_geshu);
K = -30;
K1= zeros(1,M);                                                            %K1的取值
D1 = N;   %将整个带宽划分成D1分；
for l = 1:N
     fl(1,l) = -W/2 + B + (l-1)*B;%子载波频率
end
%% ****************************选择波束************************
D_f = zeros(D1+1,N);
XX = zeros(D1+1,N);
f_k = -W/2:W/D1:W/2;%数值积分
S2 = zeros(b_geshu,D1+1);%在数值计算中记录每一个数值在不同K1（天线）时信道的求和的值
spectral_efficiency2 = zeros(1,length(snr));
spectral_efficiency5 = zeros(1,length(snr));
for cl = 1: C
    alphal = (randn(1,LL)+1j*randn(1,LL))/sqrt(2);                                                  %高斯分布
    delay_L1= rand(1,LL)* D * Ts;                                         %路径的延时服从均匀分布
    siTaL = rand(1,LL) * 2 * pi;
    arriva_anagle1 = d/wavelength*sin(siTaL);
    for i = 1:M
        SH = 0;
        for L=1:LL
            H = 0;
            for K1 = 0:M-1
                H = alphal(1,L)*((1/sqrt(M)).*exp(-1j*2*pi*K1*(arriva_anagle1(1,L)-(i-31)*a)))+H;
            end
            SH = SH + H;
        end
        i_power(1,i) =(abs(SH).^2)*W/D1;
    end
    [R,Y] = sort(i_power);
    W1 = Y(M-b_geshu+1:M);
    %%  **********************************************
    for k = 1:b_geshu
         SH = 0;
         for L=1:LL
             H = 0;
             for K1 = 0:M-1
               H = alphal(1,L)*((1/sqrt(M)).*exp(-1j*2*pi*K1*(arriva_anagle1(1,L).*(f_k./fc+1)-(W1(1,k)-31)*a)).*exp(-1j*2*pi.*f_k*delay_L1(1,L)))+H;
             end
             SH = SH + H;
         end
          S2(k,:) = abs(SH).^2;
    end
    % *********************************计算lamda******************
    for l = 1:N
         D_f(:,l)= 2*pi*fl(1,l)*1j+2*pi*f_k*1j;
         www = exp(-T.*D_f(:,l))-1;
         XX(:,l) = -1/sqrt(T)./D_f(:,l).*www; 
         HL = find(abs(f_k - fl(1,l)) <= B*10^(-6));
         for L4 = 1:length(HL)
             CC = HL(1,L4);
             XX(CC,l) = sqrt(T);
         end
    end 
    BLL2 = zeros(N,N);
    BL = zeros(D1+1,N);
    S_BL = zeros(N,N);
    for k = 1:b_geshu
        for C1 = 1:N
            for CC1 = 1:N
                BL(:,CC1) = conj(XX(:,C1)).*XX(:,CC1);
                S_BL(C1,CC1) = sum(BL(:,CC1).*S2(k,:)')*W/D1;
            end
        end
        BLL2 = BLL2 + S_BL;
    end
    BL2 = diag(BLL2);
    BL = BLL2 - diag(diag(BLL2));
    BL= abs(BL).^2;
    for l = 1:N
        S_BLL(1,l) = norm(sum(BL(l,:)));
    end
    Lamda = (BL2').^2;
    sg = BL2'*N0/2;
    %% *******************************功率分配**************************

    PI = zeros(length(snr),N);
    for i =1:length(snr)
           PI(i,:) = totolP(1,i)/N/N0;
    end
 %% ************************************************* 
 
    for i =1:length(snr)
        S = 0;
        S5 = 0;
        for a1 = 1:N
             S = S + log2(1+ PI(i,a1)* Lamda(1,a1)/((PI(i,a1))*S_BLL(1,a1) + sg(1,a1)));
             S5 = S5 + log2(1+ PI(i,a1)* Lamda(1,a1)/sg(1,a1));
        end
        spectral_efficiency2(1,i) = S/N + spectral_efficiency2(1,i);
        spectral_efficiency5(1,i) = S5/N + spectral_efficiency5(1,i);
    end
end
spectral_efficiency2 = spectral_efficiency2/C;
optimal = spectral_efficiency5/C;
R =zeros(1:length(snr));
R = log2(1 + 2*snr.*61*4*1);
plot(SNR,spectral_efficiency2,'k-s','Markersize',7,'Linewidth',2);
hold on;
plot(SNR,optimal,'k-o','Markersize',7,'Linewidth',2);
plot(SNR,R,'k->','Markersize',7,'Linewidth',2);
xlabel('SNR[dB]'),ylabel('Spectral efficiency[bps/Hz]');
legend('with interference','without interference','approximation(34)');
hold on
grid on
