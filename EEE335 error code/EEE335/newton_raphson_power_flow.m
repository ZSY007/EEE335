function [V, theta, iter] = newton_raphson_power_flow()
    % System Parameters
    Ybus = [11.333-1i*22.667, -4+1i*8, -3.333+1i*6.667;
            -4+1i*8, 8-1i*16, -4+1i*8;
            -3.333+1i*6.667, -4+1i*8, 7.333-1i*14.667];
    
    % Bus Data
    Psch = [0; -1.0; -0.8];        % Scheduled P
    Qsch = [0; -0.5; -0.4];        % Scheduled Q
    V = [1.05; 1.0; 1.0];          % Initial voltage magnitudes
    theta = [0; 0.0; 0.0];         % Initial angles in radians
    
    % Newton-Raphson Parameters
    tol = 1e-4;
    max_iter = 50;
    iter = 0;
    
    % Main iteration loop
    while iter < max_iter
        iter = iter + 1;
        
        % Initialize P and Q as column vectors
        P = zeros(3,1);
        Q = zeros(3,1);
        
        % Calculate P and Q
        for i = 1:3
            for j = 1:3
                % Calculate angle difference
                angle_diff = theta(i) - theta(j);
                % Calculate terms using admittance matrix
                P(i) = P(i) + V(i)*V(j)*(real(Ybus(i,j))*cos(angle_diff) + ...
                       imag(Ybus(i,j))*sin(angle_diff));
                Q(i) = Q(i) + V(i)*V(j)*(real(Ybus(i,j))*sin(angle_diff) - ...
                       imag(Ybus(i,j))*cos(angle_diff));
            end
        end
        
        % Calculate power mismatches (exclude slack bus)
        dP = Psch(2:3) - P(2:3);
        dQ = Qsch(2:3) - Q(2:3);
        
        % Check convergence
        if max(abs([dP; dQ])) < tol
            break;
        end
        
        % Initialize Jacobian submatrices
        J11 = zeros(2,2);
        J12 = zeros(2,2);
        J21 = zeros(2,2);
        J22 = zeros(2,2);
        
        % Form Jacobian
        for i = 2:3
            m = i-1;  % Index for Jacobian matrices
            for j = 2:3
                n = j-1;  % Index for Jacobian matrices
                if i == j  % Diagonal elements
                    % Calculate self terms
                    sumP = 0;
                    sumQ = 0;
                    for k = 1:3
                        if k ~= i
                            angle_diff = theta(i) - theta(k);
                            sumP = sumP + V(k)*(real(Ybus(i,k))*sin(angle_diff) - ...
                                   imag(Ybus(i,k))*cos(angle_diff));
                            sumQ = sumQ + V(k)*(real(Ybus(i,k))*cos(angle_diff) + ...
                                   imag(Ybus(i,k))*sin(angle_diff));
                        end
                    end
                    
                    % Diagonal elements of Jacobian submatrices
                    J11(m,n) = V(i)*sumP;
                    J12(m,n) = 2*V(i)*real(Ybus(i,i)) + sumQ;
                    J21(m,n) = -V(i)*sumQ;
                    J22(m,n) = 2*V(i)*imag(Ybus(i,i)) - sumP;
                else  % Off-diagonal elements
                    angle_diff = theta(i) - theta(j);
                    J11(m,n) = V(i)*V(j)*(real(Ybus(i,j))*sin(angle_diff) - ...
                              imag(Ybus(i,j))*cos(angle_diff));
                    J12(m,n) = V(i)*(real(Ybus(i,j))*cos(angle_diff) + ...
                              imag(Ybus(i,j))*sin(angle_diff));
                    J21(m,n) = -V(i)*V(j)*(real(Ybus(i,j))*cos(angle_diff) + ...
                              imag(Ybus(i,j))*sin(angle_diff));
                    J22(m,n) = V(i)*(real(Ybus(i,j))*sin(angle_diff) - ...
                              imag(Ybus(i,j))*cos(angle_diff));
                end
            end
        end
        
        % Form complete Jacobian
        J = [J11 J12; J21 J22];
        
        % Add small value to diagonal to improve conditioning
        J = J + eye(size(J))*1e-10;
        
        % Solve using more stable method
        dx = pinv(J)*[dP; dQ];
        
        % Limit the size of corrections
        dx = sign(dx).*min(abs(dx), 0.5);
        
        % Update state variables
        dtheta = dx(1:2);
        dV_V = dx(3:4);
        
        % Update voltage angles (except slack bus)
        theta(2:3) = theta(2:3) + dtheta;
        % Update voltage magnitudes (except slack bus)
        V(2:3) = V(2:3).*(1 + dV_V);
        
        % Prevent voltage magnitudes from going too low
        V = max(V, 0.5);
    end
    
    % Convert angles to degrees for output
    theta = theta * 180/pi;
end
