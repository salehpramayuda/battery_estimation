function x_k = ekf_3(u_k, y_k,x0, P0, W_k, Z_k, Delta, a1, a2, a3)
    
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
        eta = 0.95;
    end
    
     % Version with Q_0 as parameter/state
     % Kalman-Gain and real Covariance Matrix
    y_hat_k = a1*x_pk(5)^2+a2*x_pk(5)+a3+x_pk(6)+u_k/x_pk(3);
    C_k = [0, 0, -u_k*x_pk(3)^-2, 0, 2*a1*x_pk(5)+a2, 1];
    K_k = P_pk*C_k'/(Z_k + C_k*P_pk*C_k');
    x_k = x_pk + K_k*(y_k - y_hat_k);
    P_k = P_pk - K_k*C_k*P_pk;

    % Calculate A, prediction of x and P
    A_k = [[eye(4); zeros(2,4)], [0;0;0;Delta*eta*u_k;1;0],...
        Delta*[x_k(2)*x_k(6)-u_k;x_k(1)*x_k(6);0;0;0;1/Delta-x_k(1)*x_k(2)]]; 
    x_pkp1 = [x_k(1);x_k(2);x_k(3);x_k(4);x_k(5)-Delta*x_k(4)*u_k;...
        x_k(6)-Delta*x_k(1)*(x_k(6)*x_k(2)-u_k)];
    
    P_pkp1 = A_k*P_k*A_k' + W_k;

    % Update persistent variable
    x_pk = x_pkp1;
    P_pk = P_pkp1;
end
