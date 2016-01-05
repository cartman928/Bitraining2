%% Initialize Parameters

clc
clear

alpha = 0;    %coefficient for block fading model
beta = 0.8^2;  % Attenuation loss from non-direct antennas
w1 = 1;
w2 = 1;
n0 = 10^(-2);    %noise variance

iternums = 1:2; % number of iterations
N_realization = 500; % Number of times to run simulation

averagerateu = zeros(N_realization, length(iternums));
averageratem = zeros(N_realization, length(iternums));
averagerateu_MaxSINR = zeros(N_realization, length(iternums));
averageratem_MaxSINR = zeros(N_realization, length(iternums));
Eu = zeros(N_realization, length(iternums));
Em = zeros(N_realization, length(iternums));
Eu_MaxSINR = zeros(N_realization, length(iternums));
Em_MaxSINR = zeros(N_realization, length(iternums));

%% Training Length
for traininglength = [10 20] % traininglength 2M
        traininglength
%% Start Loop
for realization_idx = 1 : N_realization
        realization_idx
    H11 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
    H22 = (randn(2,2)+1i*randn(2,2))/sqrt(2); 
    H12 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H21 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
 
    M = traininglength/2;
    
    %% one iteration per block
    g1 = rand(2, 1) + 1i*rand(2, 1);    
    g2 = rand(2, 1) + 1i*rand(2, 1);
    g1/norm(g1);
    g2/norm(g2);
 
    v11 = zeros(2, 1); 
    v12 = zeros(2, 1);
    v21 = zeros(2, 1); 
    v22 = zeros(2, 1);
    
    for numiters = 1:length(iternums)
        x1_f = sign(randn(1,M));    
        x2_f = sign(randn(1,M));
        x1_b = sign(randn(1,M));    
        x2_b = sign(randn(1,M));  
        
        %% bi-directional training
            %%Backward Training: sudo-LS Algorithm
            [v11, v12] = S_LS_User1(H11, H12, H21, H22, g1, g2, v21, v22, M, n0, x1_b, x2_b, w1, w2); 
            %[v21, v22] = S_LS_User2(H11, H12, H21, H22, g1, g2, v11, v12, M, n0, x1_b, x2_b, w1, w2); 

            %%Forward Training: LS Algorithm
            %[g1, g2] = LS(H11, H12, H21, H22, v11, v12, v21, v22, M, n0, x1_f, x2_f, w1, w2);
        
        %{    
        averagerateu(realization_idx, numiters, traininglength) = calculate_rateu(Z, n0, Gu, Vu, Gm, upower, mpower);
        averageratem(realization_idx, numiters, traininglength) = calculate_ratem(Z, n0, Gm, Vm, Gu, upower, mpower);
        averagerateu_MaxSINR(realization_idx, numiters) = calculate_rateu(Z, n0, Gu_w, Vu_w, Gm_w, upower, mpower);
        averageratem_MaxSINR(realization_idx, numiters) = calculate_ratem(Z, n0, Gm_w, Vm_w, Gu_w, upower, mpower);
        Eu(realization_idx, numiters,traininglength) = MSEu(Z, Vu, Vm, Gu, Gm, n0, upower, mpower);
        Em(realization_idx, numiters,traininglength) = MSEm(Z, Vu, Vm, Gu, Gm, n0, upower, mpower);
        Eu_MaxSINR(realization_idx, numiters) = MSEu(Z, Vu_w, Vm_w, Gu_w, Gm_w, n0, upower, mpower);
        Em_MaxSINR(realization_idx, numiters) = MSEm(Z, Vu_w, Vm_w, Gu_w, Gm_w, n0, upower, mpower);
        %}
        
            
    end
            
    
end

end


%{
%% Plot C(bits/channel)
figure
subplot(2,1,1);
hold on

p1=plot(iternums, mean(averagerateu(:,:,10))+mean(averageratem(:,:,10)),'Color',[0,0.4470,0.7410]);
p2=plot(iternums,mean(averageratem(:,:,10)),'Color',[0,0.4470,0.7410],'Marker','o');
p3=plot(iternums, mean(averagerateu(:,:,10)),'Color',[0,0.4470,0.7410],'Marker','*');

%{
p4=plot(iternums, mean(averagerateu(:,:,16))+mean(averageratem(:,:,16)),'Color',[0.6350,0.0780,0.1840]);
p5=plot(iternums,mean(averageratem(:,:,16)),'Color',[0.6350,0.0780,0.1840],'Marker','o');
p6=plot(iternums, mean(averagerateu(:,:,16)),'Color',[0.6350,0.0780,0.1840],'Marker','*');
%}

p7=plot(iternums, mean(averagerateu(:,:,20))+mean(averageratem(:,:,20)),'Color',[0.8500,0.3250,0.0980]);
p8=plot(iternums,mean(averageratem(:,:,20)),'Color',[0.8500,0.3250,0.0980],'Marker','o');
p9=plot(iternums, mean(averagerateu(:,:,20)),'Color',[0.8500,0.3250,0.0980],'Marker','*');

p10=plot(iternums, mean(averagerateu_MaxSINR)+mean(averageratem_MaxSINR),'k');
p11=plot(iternums, mean(averageratem_MaxSINR),'k','Marker','o');
p12=plot(iternums, mean(averagerateu_MaxSINR),'k','Marker','*');

legend([p10,p1,p7],'C(Max-SINR)',...
                   'C(Bi-Directional Training);2M=10',...    %'C(Bi-Directional Training);2M=16',...
                   'C(Bi-Directional Training);2M=20')
                  
xlabel('Number of iterations')
ylabel('C(bits/channel)')
%title('2 Users;2X2 MIMO Channel;\sigma^2=10^{-2};1000 Realizations;No-Coop')
%title('2 Users;2X2 MIMO Channel;\sigma^2=10^{-2};1000 Realizations;Coop')
title('2 Users;2X2 MIMO Channel;\sigma^2=10^{-2};1000 Realizations')
axis([1 numiters 0 20])

%% Plot MSE
subplot(2,1,2);
hold on

p13=plot(iternums,mean(Em(:,:,10)),'Color',[0,0.4470,0.7410],'Marker','o');
p14=plot(iternums, mean(Eu(:,:,10)),'Color',[0,0.4470,0.7410],'Marker','*');

p15=plot(iternums,mean(Em(:,:,20)),'Color',[0.8500,0.3250,0.0980],'Marker','o');
p16=plot(iternums, mean(Eu(:,:,20)),'Color',[0.8500,0.3250,0.0980],'Marker','*');

p17=plot(iternums, mean(Em_MaxSINR),'k','Marker','o');
p18=plot(iternums, mean(Eu_MaxSINR),'k','Marker','*');

xlabel('Number of iterations')
ylabel('MSE')

%}