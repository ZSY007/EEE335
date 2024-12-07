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

%% Gauss-Seidel 方法求解
n = length(V); % 节点数
error = inf; % 初始误差
iter = 0; % 迭代计数

while error > tolerance && iter < max_iter
    iter = iter + 1;
    V_old = V; % 保存上一轮电压值，用于计算误差
    
    % 遍历每个 PQ 节点
    for i = 2:n
        % 计算注入电流
        I_inj = conj((P_load(i) + 1j * Q_load(i)) / V(i));
        
        % 更新电压值
        V(i) = (I_inj - sum(Y_matrix(i, :) .* V.') + Y_matrix(i, i) * V(i)) / Y_matrix(i, i);
    end
    
    % 计算最大误差
    error = max(abs(V - V_old));
end


% 打印结果
disp('潮流计算结果:');
fprintf('迭代次数: %d\n', iter);
disp('节点电压 (V):');
disp(V);
