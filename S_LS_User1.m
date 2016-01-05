function  [v11, v12] = S_LS_User1(H11, H12, H21, H22, g1, g2, v21, v22, M, n0, x1_b, x2_b, w1, w2)
%update filters by sudo-LS algorithm 

P = 1;
upperbound = 1000;
lowerbound = -1000;
lambda1 = (upperbound+lowerbound)/2;

%received signal vectors
y1 = H11'*g1*x1_b+H21'*g2*x2_b;
y2 = H12'*g1*x1_b+H22'*g2*x2_b;

%estimations
Ex1y1 = mean([x1_b;x1_b].*y1,2);
Ex2y1 = mean([x2_b;x2_b].*y1,2);
Ex1y2 = mean([x1_b;x1_b].*y2,2);
Ex2y2 = mean([x2_b;x2_b].*y2,2);

lambda1 = -10;
%filters
v11 = ( 2*Ex1y1*Ex1y1'*w1 + 2*Ex2y1*Ex2y1'*w2 + 2*lambda1'*eye(2))\...
      ( 2*Ex1y1*w1 - Ex1y2'*v21*Ex1y1*w1...
        -Ex1y1*v21'*Ex1y2*w1...
        -Ex2y2'*v21*Ex2y1*w2...
        -Ex2y1*v21'*Ex2y2*w2);
    
v12 = ( 2*Ex1y1*Ex1y1'*w1 + 2*Ex2y1*Ex2y1'*w2 + 2*lambda1'*eye(2))\...
      ( 2*Ex2y1*w2 - Ex1y2'*v22*Ex1y1*w1...
        -Ex1y1*v22'*Ex1y2*w1...
        -Ex2y2'*v22*Ex2y1*w2...
        -Ex2y1*v22'*Ex2y2*w2);
    
P0 = norm(v11)^2+norm(v12)^2


%binary search for lambda1 

while abs(  norm(v11)^2+norm(v12)^2-P   ) > 10^(-5)
    eeee = norm(v11)^2+norm(v12)^2-P 
    %Power Too Large (lambda too small)
    if norm(v11)^2+norm(v12)^2  > P
        %set lower bound
        lowerbound = lambda1;
        %solve for new lambda
        lambda1 = (upperbound+lowerbound)/2
        
        v11 = ( 2*Ex1y1*Ex1y1'*w1 + 2*Ex2y1*Ex2y1'*w2 + 2*lambda1'*eye(2))\...
        ( 2*Ex1y1*w1 - Ex1y2'*v21*Ex1y1*w1...
        -Ex1y1*v21'*Ex1y2*w1...
        -Ex2y2'*v21*Ex2y1*w2...
        -Ex2y1*v21'*Ex2y2*w2);
    
        v12 = ( 2*Ex1y1*Ex1y1'*w1 + 2*Ex2y1*Ex2y1'*w2 + 2*lambda1'*eye(2))\...
        ( 2*Ex2y1*w2 - Ex1y2'*v22*Ex1y1*w1...
        -Ex1y1*v22'*Ex1y2*w1...
        -Ex2y2'*v22*Ex2y1*w2...
        -Ex2y1*v22'*Ex2y2*w2);
    
        P1 = norm(v11)^2+norm(v12)^2
        
        if norm(v11)^2+norm(v12)^2  < P
            upperbound = lambda1;
        else
        end
 
    %Power Too Small (lambda too large)
    else
        %set lower bound
        upperbound = lambda1;
        %solve for new lambda
        lambda1 = (upperbound+lowerbound)/2
        
        v11 = ( 2*Ex1y1*Ex1y1'*w1 + 2*Ex2y1*Ex2y1'*w2 + 2*lambda1'*eye(2))\...
        ( 2*Ex1y1*w1 - Ex1y2'*v21*Ex1y1*w1...
        -Ex1y1*v21'*Ex1y2*w1...
        -Ex2y2'*v21*Ex2y1*w2...
        -Ex2y1*v21'*Ex2y2*w2);
    
        v12 = ( 2*Ex1y1*Ex1y1'*w1 + 2*Ex2y1*Ex2y1'*w2 + 2*lambda1'*eye(2))\...
        ( 2*Ex2y1*w2 - Ex1y2'*v22*Ex1y1*w1...
        -Ex1y1*v22'*Ex1y2*w1...
        -Ex2y2'*v22*Ex2y1*w2...
        -Ex2y1*v22'*Ex2y2*w2);
    
        P2 = norm(v11)^2+norm(v12)^2
        
        if norm(v11)^2+norm(v12)^2  > P
            lowerbound = lambda1;
        else
        end
        
    end
    
end

    
    
