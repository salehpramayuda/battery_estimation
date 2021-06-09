function x_k  = ekf(u_k,y_k,x0, P0, W_k, Z_k, Delta)
    % Extended Kalman Filter for estimating Capacity and U_L
    
    % Initialisierung Prädiktion bei k = 0
    persistent P_k_p x_k_p
    if isempty(P_k_p)
        P_k_p = P0;
    end
    if isempty(x_k_p)
        x_k_p = x0;
    end
    
    % y_hat_k, C, K, x_k, P_k
    y_hat_k = x_k_p(3)*(x_k_p(4)-x_k_p(5)-u_k);
    C_k = [0, 0, x_k_p(4)-x_k_p(5)-u_k, x_k_p(3), -x_k_p(3)];
    K_k = P_k_p*C_k'/(Z_k + C_k*P_k_p*C_k');
    x_k = x_k_p + K_k*(y_k - y_hat_k);
    P_k = P_k_p - K_k*C_k*P_k_p;
    
    % Calculate A, prediction of x and P
    A_k = [[eye(4), zeros(4,1)]; [Delta*(x_k(3)*(x_k(4)-u_k)-x_k(5)*...
        (x_k(2)+x_k(3))), -Delta*x_k(1)*x_k(5), Delta*x_k(1)*...
        (x_k(4)-u_k-x_k(5)), Delta*x_k(1)*x_k(3),1-Delta*x_k(1)*(x_k(2)+x_k(3))]]; 
    x_kp1_p = [x_k(1);x_k(2);x_k(3);x_k(4);x_k(5)+Delta*x_k(1)*...
        (x_k(3)*(x_k(4)-u_k)-x_k(5)*(x_k(2)+x_k(3)))];

    P_kp1_p = A_k*P_k*A_k' + W_k;
    
    % Update persistent variable
    x_k_p = x_kp1_p;
    P_k_p = P_kp1_p;
end