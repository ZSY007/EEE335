% 输入阻抗值
Z12 = 0.04  + 0.08j;
Z23 = 0.05 + 0.1j;
Z13 = 0.06 + 0.12j;

% 计算导纳值 Y = 1/Z
Y12 = -1 / Z12;
Y23 = -1 / Z23;
Y13 = -1 / Z13;

% 计算对角线元素 Yii
Y11 = -Y12 - Y13;
Y22 = -Y12 - Y23;
Y33 = -Y13 - Y23;

% 构造导纳矩阵
Y_matrix = [
    Y11, Y12, Y13;
    Y12, Y22, Y23;
    Y13, Y23, Y33
];

% 打印导纳矩阵
disp('导纳矩阵 (Y):');
disp(Y_matrix);

% 假设初始电压 (V)
V1 = 1.05 + 0j; % Slack bus
V2 = 1.0 + 0j;  % PQ bus
V3 = 1.0 + 0j;  % PQ bus
V = [V1; V2; V3];

% 初始化功率
P = zeros(3, 1);
Q = zeros(3, 1);

% 计算功率平衡方程
for i = 1:3
    power = V(i) * conj(sum(Y_matrix(i, :) .* V.'));
    P(i) = real(power);
    Q(i) = -imag(power);
end

disp('有功功率注入 (P):');
disp(P);
disp('无功功率注入 (Q):');
disp(Q);

% 定义计算雅可比矩阵的函数
function [J11, J12, J21, J22] = calculate_jacobian(Y_matrix, V)
    n = length(V);
    G = real(Y_matrix);
    B = imag(Y_matrix);
    V_mag = abs(V);
    V_angle = angle(V);
    
    J11 = zeros(n - 1, n - 1);
    J12 = zeros(n - 1, n - 1);
    J21 = zeros(n - 1, n - 1);
    J22 = zeros(n - 1, n - 1);

    for i = 2:n
        for j = 2:n
            theta_ij = V_angle(i) - V_angle(j);
            if i ~= j
                J11(i - 1, j - 1) = V_mag(i) * V_mag(j) * (G(i, j) * sin(theta_ij) - B(i, j) * cos(theta_ij));
                J12(i - 1, j - 1) = V_mag(i) * (G(i, j) * cos(theta_ij) + B(i, j) * sin(theta_ij));
                J21(i - 1, j - 1) = -V_mag(i) * V_mag(j) * (G(i, j) * cos(theta_ij) + B(i, j) * sin(theta_ij));
                J22(i - 1, j - 1) = V_mag(i) * (G(i, j) * sin(theta_ij) - B(i, j) * cos(theta_ij));
            else
                sum_J11 = 0; sum_J12 = 0; sum_J21 = 0; sum_J22 = 0;
                for k = 1:n
                    if k ~= i
                        theta_ik = V_angle(i) - V_angle(k);
                        sum_J11 = sum_J11 + V_mag(k) * (G(i, k) * sin(theta_ik) - B(i, k) * cos(theta_ik));
                        sum_J12 = sum_J12 + V_mag(k) * (G(i, k) * cos(theta_ik) + B(i, k) * sin(theta_ik));
                        sum_J21 = sum_J21 + V_mag(k) * (G(i, k) * cos(theta_ik) + B(i, k) * sin(theta_ik));
                        sum_J22 = sum_J22 + V_mag(k) * (G(i, k) * sin(theta_ik) - B(i, k) * cos(theta_ik));
                    end
                end
                J11(i - 1, j - 1) = -V_mag(i)^2 * B(i, i) - V_mag(i) * sum_J11;
                J12(i - 1, j - 1) = V_mag(i) * G(i, i) + sum_J12;
                J21(i - 1, j - 1) = V_mag(i)^2 * G(i, i) - V_mag(i) * sum_J21;
                J22(i - 1, j - 1) = V_mag(i) * B(i, i) + sum_J22;
            end
        end
    end
end

% 计算雅可比矩阵
[J11, J12, J21, J22] = calculate_jacobian(Y_matrix, V);

disp('雅可比矩阵 J11 (∂P/∂θ):');
disp(J11);
disp('雅可比矩阵 J12 (∂P/∂|V|):');
disp(J12);
disp('雅可比矩阵 J21 (∂Q/∂θ):');
disp(J21);
disp('雅可比矩阵 J22 (∂Q/∂|V|):');
disp(J22);
