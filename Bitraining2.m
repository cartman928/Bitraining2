%% Initialize Parameters

clc
clear

alpha = 0;    %coefficient for block fading model
beta = 0.8^2;  % Attenuation loss from non-direct antennas
w1 = 1;
w2 = 1;
n0 = 10^(-2);    %noise variance

iternums = 1:10; % number of iterations
N_realization = 500; % Number of times to run simulation

C1 = zeros(N_realization, length(iternums));
C2 = zeros(N_realization, length(iternums));

%% Training Length
for traininglength = [10 20 30] % traininglength 2M
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
            v11_old = v11;
            v12_old = v12;
            v21_old = v21;
            v22_old = v22;

            N1 = sqrt(n0)*(randn(2,M)+1i*randn(2,M))/sqrt(2);
            N2 = sqrt(n0)*(randn(2,M)+1i*randn(2,M))/sqrt(2);
            
            %%Backward Training: sudo-LS Algorithm
            [v11, v12] = S_LS_User1_Brutal(H11, H12, H21, H22, g1, g2, v21_old, v22_old, M, n0, N1, N2, x1_b, x2_b, w1, w2); 
            [v21, v22] = S_LS_User2_Brutal(H11, H12, H21, H22, g1, g2, v11_old, v12_old, M, n0, N1, N2, x1_b, x2_b, w1, w2); 
            %norm(v11)^2+norm(v12)^2


            %%Forward Training: LS Algorithm
            [g1, g2] = LS(H11, H12, H21, H22, v11, v12, v21, v22, M, n0, x1_f, x2_f);
            %norm(g1)^2
            %norm(g2)^2

            SINR1 = norm(g1'*(H11*v11+H12*v21))^2/(norm(g1'*(H11*v12+H12*v22))^2+n0*g1'*g1);
            SINR2 = norm(g2'*(H21*v12+H22*v22))^2/(norm(g2'*(H21*v11+H22*v21))^2+n0*g2'*g2);
            C1(realization_idx, numiters, traininglength) = abs(log2(1+SINR1));
            C2(realization_idx, numiters, traininglength) = abs(log2(1+SINR2));
        
            
    end
            
    
end

end



%% Plot C(bits/channel)
figure
hold on

p1=plot(iternums, mean(C1(:,:,10))+mean(C2(:,:,10)),'--');
p2=plot(iternums, mean(C2(:,:,20))+mean(C2(:,:,20)),'o');
p3=plot(iternums, mean(C1(:,:,30))+mean(C2(:,:,30)));


axis([1 numiters 0 12])

