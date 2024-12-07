% 输入阻抗值
Z12 = 0.04 + 0.08j;
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

% 初始电压值 (V)
V = [1.05 + 0j; 1.0 + 0j; 1.0 + 0j];

% 初始功率需求
S_specified = [0; -0.5 - 0.2j; -0.6 - 0.3j]; % 负荷功率需求
P_specified = real(S_specified);
Q_specified = imag(S_specified);

% 收敛参数
tolerance = 1e-4; % 收敛标准
max_iter = 50;    % 最大迭代次数
iter = 0;

while true
    % 计算注入功率
    P_injected = zeros(3, 1);
    Q_injected = zeros(3, 1);

    for i = 1:3
        power = V(i) * conj(sum(Y_matrix(i, :) .* V.'));
        P_injected(i) = real(power);
        Q_injected(i) = -imag(power);
    end

    % 计算功率误差
    delta_P = P_specified(2:end) - P_injected(2:end);
    delta_Q = Q_specified(2:end) - Q_injected(2:end);
    delta = [delta_P; delta_Q];

    % 判断收敛
    if max(abs(delta)) < tolerance
        break;
    end

    % 计算雅可比矩阵
    [J11, J12, J21, J22] = calculate_jacobian(Y_matrix, V);
    J = [J11, J12; J21, J22];

    % 修正值
    correction = J \ delta;

    % 更新电压
    delta_theta = correction(1:2);
    delta_V = correction(3:4);

    % 更新幅值和相角
    V_angle = angle(V);
    V_mag = abs(V);

    V_angle(2:end) = V_angle(2:end) + delta_theta;
    V_mag(2:end) = V_mag(2:end) + delta_V;

    V = V_mag .* exp(1j * V_angle);

    % 更新迭代次数
    iter = iter + 1;
    if iter >= max_iter
        error('Newton-Raphson 方法未收敛');
    end
end

% 输出结果
disp('收敛结果：');
disp(['迭代次数: ', num2str(iter)]);
disp('节点电压 (幅值和相角):');
for i = 1:3
    fprintf('V%d: %.4f ∠ %.4f°\n', i, abs(V(i)), rad2deg(angle(V(i))));
end

disp('注入功率:');
disp('有功功率 P (pu):');
disp(P_injected);
disp('无功功率 Q (pu):');
disp(Q_injected);

% 定义雅可比矩阵计算函数
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
