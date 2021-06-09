function x_k = ekf_2(u_k, y_k, x0, P0, W_k, Z_k, Delta, a1, a2, Q_e, a3)
    
    persistent P_pk x_pk
    if isempty(x_pk)
        x_pk = x0; %x[k|k-1]
    end
    if isempty(P_pk)
        P_pk = P0;
    end
    
    if(u_k<0)
        eta = 1;
    else
        eta = 1;
    end
    
    % Kalman-Gain and real Covariance Matrix
    % States (1/C, 1/R, 1/R0, SOC, Uc)
    y_hat_k = a1*x_pk(4)^2+a2*x_pk(4)+a3+x_pk(5)+u_k/x_pk(3);
    C_k = [0, 0, -u_k/x_pk(3)^2, 2*a1*x_pk(4)+a2, 1];
    K_k = P_pk*C_k'/(Z_k + C_k*P_pk*C_k');
    x_k = x_pk + K_k*(y_k - y_hat_k);
    P_k = P_pk - K_k*C_k*P_pk;

    % Calculate A, prediction of x and P
    A_k = [[eye(4); zeros(1,4)], [Delta*(u_k-x_k(5)*x_k(3));0;...
        -Delta*x_k(1)*x_k(5);0;1-Delta*x_k(1)*x_k(3)]];
    x_pkp1 = [x_k(1);x_k(2);x_k(3);x_k(4)+Delta*u_k/Q_e;...
    x_k(5)-Delta*x_k(1)*(x_k(5)*x_k(2)-u_k)];
    
    P_pkp1 = A_k*P_k*A_k' + W_k;
    
    % Update persistent variable
    x_pk = x_pkp1;
    P_pk = P_pkp1;
end


