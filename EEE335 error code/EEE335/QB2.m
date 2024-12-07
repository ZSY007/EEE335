% MATLAB代码：使用Gauss-Seidel方法求解潮流问题

% 阻抗值输入
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
Y_matrix = [Y11, Y12, Y13;
            Y12, Y22, Y23;
            Y13, Y23, Y33];

% 初始节点电压
V1 = 1.05;  % Slack bus, 电压固定
V2 = 1.0;   % PQ bus, 初始猜测电压
V3 = 1.0;   % PQ bus, 初始猜测电压

% 设置功率负荷
P2 = -1.0;  % Bus 2 的有功功率负荷 (负值表示负荷)
Q2 = -0.5;  % Bus 2 的无功功率负荷
P3 = -0.8;  % Bus 3 的有功功率负荷
Q3 = -0.4;  % Bus 3 的无功功率负荷

% 收敛参数
tolerance = 1e-6; % 收敛精度
max_iter = 100;   % 最大迭代次数
iter = 0;         % 初始化迭代次数

% 初始电压向量
V = [V1; V2; V3];

% Gauss-Seidel 迭代
while iter < max_iter
    V_old = V;  % 保存旧电压值
    
    % Slack bus (V1) 不变，跳过
    % PQ bus: 更新 V2
    V(2) = (1 / Y_matrix(2, 2)) * ...
           ((P2 - 1j * Q2) / conj(V(2)) - ...
           Y_matrix(2, 1) * V(1) - ...
           Y_matrix(2, 3) * V(3));
    
    % PQ bus: 更新 V3
    V(3) = (1 / Y_matrix(3, 3)) * ...
           ((P3 - 1j * Q3) / conj(V(3)) - ...
           Y_matrix(3, 1) * V(1) - ...
           Y_matrix(3, 2) * V(2));
    
    % 检查收敛条件
    if max(abs(V - V_old)) < tolerance
        disp('潮流计算收敛！');
        break;
    end
    
    iter = iter + 1;
end

% 如果未收敛
if iter >= max_iter
    disp('潮流计算未收敛！');
end

% 输出结果
disp('迭代次数:');
disp(iter);
disp('最终节点电压 (V):');
disp(V);

% 计算节点功率
S = zeros(3, 1); % 节点注入功率
for i = 1:3
    S(i) = V(i) * conj(sum(Y_matrix(i, :) .* V.'));
end
disp('节点注入功率 (S):');
disp(S);
