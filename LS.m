function [g1, g2] = LS(H11, H12, H21, H22, v11, v12, v21, v22, M, n0, x1_f, x2_f)
%update filters by least square algorithm 

%Received Signals
y1 = H11*(v11*x1_f+v12*x2_f)+H12*(v21*x1_f+v22*x2_f) + sqrt(n0)*(randn(2,M)+1i*randn(2,M))/sqrt(2);
y2 = H21*(v11*x1_f+v12*x2_f)+H22*(v21*x1_f+v22*x2_f) + sqrt(n0)*(randn(2,M)+1i*randn(2,M))/sqrt(2);

%LS Algorithm
g1 = (y1*y1')\y1*x1_f';
g2 = (y2*y2')\y2*x2_f';

%Normalize
g1 = g1/norm(g1);
g2 = g2/norm(g2);
end

    
    