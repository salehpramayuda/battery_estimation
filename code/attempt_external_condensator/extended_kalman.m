function [x, y_hat] = extended_kalman(u_k, y_k, Delta, x0, P0, a, Z_k, W_k, Qe)

    persistent x_p_k P_p_k
    
    if isempty(x_p_k)
        x_p_k = x0;
    end
    if isempty(P_p_k)
        P_p_k = P0;
    end
    
    % calculate states
    C = [0, 0, -u_k/x_p_k(3)^2, 2*a(1)*x_p_k(4)^2+a(2), 1];
    

    
end
