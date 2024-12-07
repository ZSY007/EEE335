% 输入阻抗值
Z12 = 0.04 + 0.08j;
Z23 = 0.05 + 0.1j;
Z13 = 0.06 + 0.12j;

% 导纳计算
Y12 = -1 / Z12;
Y23 = -1 / Z23;
Y13 = -1 / Z13;
Y11 = -Y12 - Y13;
Y22 = -Y12 - Y23;
Y33 = -Y13 - Y23;

% 导纳矩阵
Y_matrix = [Y11, Y12, Y13;
            Y12, Y22, Y23;
            Y13, Y23, Y33];

disp('导纳矩阵 (Y):');
disp(Y_matrix);

% 节点功率
P_load = [-1.0; -0.8];  % 有功负荷
Q_load = [-0.5; -0.4];  % 无功负荷

% 初始电压幅值和相角
V0 = [1.05; 1.02; 1.01];  % 初始电压幅值
theta0 = [0; 0; 0];       % 初始相角
x0 = [theta0(2:end); V0(2:end)];  % 初始变量 (PQ节点)

% 定义潮流方程组，支持雅可比矩阵计算
function [F, J] = power_flow_eq_with_jacobian(x, Y_matrix, P_load, Q_load)
    theta = [0; x(1:2)];  % 保持平衡节点相角为 0
    V = [1.05; x(3:4)];   % 保持平衡节点电压为 1.05

    % 初始化误差向量
    P_calc = zeros(3, 1);
    Q_calc = zeros(3, 1);
    J = zeros(4, 4); % 雅可比矩阵（2x有功 + 2x无功）

    % 计算有功功率和无功功率
    for i = 1:3
        for j = 1:3
            P_calc(i) = P_calc(i) + abs(V(i)) * abs(V(j)) * abs(Y_matrix(i, j)) * cos(theta(j) - theta(i) + angle(Y_matrix(i, j)));
            Q_calc(i) = Q_calc(i) + abs(V(i)) * abs(V(j)) * abs(Y_matrix(i, j)) * sin(theta(j) - theta(i) + angle(Y_matrix(i, j)));
        end
    end

    % 生成误差向量 (仅计算 PQ 节点的误差)
    delta_P = P_calc(2:3) - P_load;  % PQ节点的有功功率误差
    delta_Q = Q_calc(2:3) - Q_load;  % PQ节点的无功功率误差
    F = [delta_P; delta_Q]; % 总误差向量

    % 计算雅可比矩阵
    for i = 2:3
        for j = 2:3
            % dP/dtheta, dP/dV
            J(i-1, j-1) = -abs(V(i)) * abs(V(j)) * abs(Y_matrix(i, j)) * sin(theta(j) - theta(i) + angle(Y_matrix(i, j)));
            J(i-1, j+1) = abs(V(j)) * abs(Y_matrix(i, j)) * cos(theta(j) - theta(i) + angle(Y_matrix(i, j)));
            % dQ/dtheta, dQ/dV
            J(i+1, j-1) = -abs(V(i)) * abs(V(j)) * abs(Y_matrix(i, j)) * cos(theta(j) - theta(i) + angle(Y_matrix(i, j)));
            J(i+1, j+1) = -abs(V(j)) * abs(Y_matrix(i, j)) * sin(theta(j) - theta(i) + angle(Y_matrix(i, j)));
        end
    end
end

% 自定义迭代监控函数
global iter_count
iter_count = 0;

function stop = iter_monitor(~, optimValues, state)
    global iter_count
    if strcmp(state, 'iter')
        iter_count = optimValues.iteration; % 更新迭代次数
        fprintf('当前迭代次数: %d, 函数调用次数: %d\n', iter_count, optimValues.funccount);
    end
    stop = false; % 不停止迭代
end

% 设置 fsolve 选项，启用雅可比矩阵和迭代监控
options = optimoptions('fsolve', 'Display', 'iter', 'Jacobian', 'on', ...
                       'OutputFcn', @iter_monitor, 'TolFun', 1e-6, 'MaxIterations', 100);

% 使用 fsolve 求解潮流方程
[x, fval, exitflag] = fsolve(@(x) power_flow_eq_with_jacobian(x, Y_matrix, P_load, Q_load), x0, options);

% 提取结果
theta = [0; x(1:2)];     % 平衡节点相角为 0
V_mag = [1.05; x(3:4)];  % 平衡节点电压为 1.05
V = V_mag .* exp(1j * theta);  % 计算复数形式的节点电压

% 计算节点注入功率
S = zeros(3, 1);  % 节点注入功率
for i = 1:3
    for j = 1:3
        S(i) = S(i) + V(i) * conj(Y_matrix(i, j) * V(j));
    end
end

% 输出结果
disp(['总迭代次数: ', num2str(iter_count)]);
disp('最终节点电压 (V):');
disp(V);

disp('节点注入功率 (S):');
disp(S);
