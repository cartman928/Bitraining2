function  [v21, v22] = S_LS_User2_Brutal(H11, H12, H21, H22, g1, g2, v11_old, v12_old, M, n0, N1, N2,x1_b, x2_b, w1, w2)
%update filters by sudo-LS algorithm 

    %Power Constraint
    P = 1;

    %received signal vectors
    y1 = H11'*g1*x1_b+H21'*g2*x2_b + N1;
    y2 = H12'*g1*x1_b+H22'*g2*x2_b + N2;

    %estimations
    Ex1y1 = mean([x1_b;x1_b].*y1,2);
    Ex2y1 = mean([x2_b;x2_b].*y1,2);
    Ex1y2 = mean([x1_b;x1_b].*y2,2);
    Ex2y2 = mean([x2_b;x2_b].*y2,2);


    %brutal-force search for lambda1 
    Range = 1000;
    Precision = 0.01;
    for n = -Range:Range

        v21 = ( 2*Ex1y2*Ex1y2'*w1 + 2*Ex2y2*Ex2y2'*w2 + 2*(n'*Precision)*eye(2))\...
          ( 2*Ex1y2*w1 - Ex1y1'*v11_old*Ex1y2*w1...
            -Ex1y2*v11_old'*Ex1y1*w1...
            -Ex2y1'*v11_old*Ex2y2*w2...
            -Ex2y2*v11_old'*Ex2y2*w2);

        v22 = ( 2*Ex1y2*Ex1y2'*w1 + 2*Ex2y2*Ex2y2'*w2 + 2*(n'*Precision)*eye(2))\...
          ( 2*Ex2y2*w2 - Ex1y1'*v12_old*Ex1y2*w1...
            -Ex1y2*v12_old'*Ex1y1*w1...
            -Ex2y1'*v12_old*Ex2y2*w2...
            -Ex2y2*v12_old'*Ex2y1*w2);

        W(n+1+Range) = abs(norm(v11_old)^2+norm(v12_old)^2-P);
    end
    
    [M,I] = min(W);
    lambda2 = (I-1-Range)*Precision;
    
    v21 = ( 2*Ex1y1*Ex1y1'*w1 + 2*Ex2y1*Ex2y1'*w2 + 2*lambda2'*eye(2))\...
          ( 2*Ex1y1*w1 - Ex1y2'*v21*Ex1y1*w1...
            -Ex1y1*v21'*Ex1y2*w1...
            -Ex2y2'*v21*Ex2y1*w2...
            -Ex2y1*v21'*Ex2y2*w2);

    v22 = ( 2*Ex1y1*Ex1y1'*w1 + 2*Ex2y1*Ex2y1'*w2 + 2*lambda2'*eye(2))\...
          ( 2*Ex2y1*w2 - Ex1y2'*v22*Ex1y1*w1...
            -Ex1y1*v22'*Ex1y2*w1...
            -Ex2y2'*v22*Ex2y1*w2...
            -Ex2y1*v22'*Ex2y2*w2);
    
end

    
    
