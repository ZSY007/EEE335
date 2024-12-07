%% 初始化网络参数
% 阻抗矩阵 (Z)
Z12 = 0.04 + 0.08j;
Z23 = 0.05 + 0.1j;
Z13 = 0.06 + 0.12j;

% 导纳矩阵 Y
Y12 = -1 / Z12;
Y23 = -1 / Z23;
Y13 = -1 / Z13;

Y11 = -Y12 - Y13;
Y22 = -Y12 - Y23;
Y33 = -Y13 - Y23;

Y_matrix = [
    Y11, Y12, Y13;
    Y12, Y22, Y23;
    Y13, Y23, Y33
];

% 负荷 (P, Q) 和基准电压
P_load = [0; -1; -0.5]; % 负荷有功功率 (单位: pu)
Q_load = [0; -0.5; -0.3]; % 负荷无功功率 (单位: pu)
V_slack = 1.05 + 0j; % 平衡节点电压 (Slack Bus)

% 初始电压
V = [V_slack; 1.0 + 0j; 1.0 + 0j]; % 初始猜测电压

% 收敛标准
tolerance = 1e-4;
max_iter = 100;

%% Newton-Raphson 方法求解
n = length(V); % 节点数
PQ_indices = 2:n; % PQ 节点索引
num_PQ = length(PQ_indices); % PQ 节点数量

error = inf; % 初始误差
iter = 0; % 迭代计数

while error > tolerance && iter < max_iter
    iter = iter + 1;
    % 保存当前电压
    V_mag = abs(V);
    V_angle = angle(V);
    
    % 计算功率误差
    P_calc = zeros(n, 1);
    Q_calc = zeros(n, 1);
    for i = 1:n
        S = V(i) * conj(sum(Y_matrix(i, :) .* V.'));
        P_calc(i) = real(S);
        Q_calc(i) = -imag(S);
    end
    % 功率不平衡
    dP = P_load - P_calc;
    dQ = Q_load - Q_calc;
    mismatch = [dP(2:end); dQ(2:end)]; % 移除平衡节点
    
    % 计算雅可比矩阵
    [J11, J12, J21, J22] = calculate_jacobian(Y_matrix, V);
    J = [J11, J12; J21, J22]; % 整体雅可比矩阵

    % 更新步长
    delta = J \ mismatch; % 解线性方程组 J * delta = mismatch
    
    % 更新电压幅值和相角
    delta_theta = delta(1:num_PQ);
    delta_Vmag = delta(num_PQ+1:end);
    
    V_angle(PQ_indices) = V_angle(PQ_indices) + delta_theta;
    V_mag(PQ_indices) = V_mag(PQ_indices) + delta_Vmag;
    
    % 更新电压向量
    V(PQ_indices) = V_mag(PQ_indices) .* exp(1j * V_angle(PQ_indices));
    
    % 计算最大误差
    error = max(abs(mismatch));
end

%% 输出结果
% 计算最终功率注入
P_inj = zeros(n, 1);
Q_inj = zeros(n, 1);
for i = 1:n
    S = V(i) * conj(sum(Y_matrix(i, :) .* V.'));
    P_inj(i) = real(S);
    Q_inj(i) = -imag(S);
end

% 打印结果
disp('潮流计算结果:');
fprintf('迭代次数: %d\n', iter);
disp('节点电压 (V):');
disp(V);
disp('有功功率注入 (P):');
disp(P_inj);
disp('无功功率注入 (Q):');
disp(Q_inj);

%% 雅可比矩阵计算函数
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
