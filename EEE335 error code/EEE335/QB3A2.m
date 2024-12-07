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

% 定义符号变量
syms V2 V3 theta2 theta3 real

% 电压幅值和相角的符号变量
V = [1.05; V2; V3];  % 电压幅值
theta = [0; theta2; theta3];  
% 定义潮流方程的符号变量
X = [theta1; theta2; V2; V3];  % PQ节点电压幅值和相角

% 计算功率误差
delta_P = P_calc(2:3) - P_load;
delta_Q = Q_calc(2:3) - Q_load;

F = [delta_P; delta_Q];

% 计算潮流方程的雅可比矩阵 (对X求导)
J = jacobian(F, X);
disp('雅可比矩阵：');
disp(J);

% 使用fsolve求解潮流方程
function F = power_flow_eq(x, Y_matrix, P_load, Q_load)
    % 替换符号变量为数值
    theta = [0; x(1:2)];  % 相角
    V = [1.05; x(3:4)];   % 电压幅值

    % 计算有功功率和无功功率
    P_calc = zeros(3, 1);
    Q_calc = zeros(3, 1);
    for i = 1:3
        for j = 1:3
            % 计算P和Q
            P_calc(i) = P_calc(i) + abs(V(i)) * abs(V(j)) * abs(Y_matrix(i, j)) * cos(angle(Y_matrix(i, j)) - angle(V(i)) + angle(Y_matrix(i, j)));
            Q_calc(i) = Q_calc(i) + abs(V(i)) * abs(V(j)) * abs(Y_matrix(i, j)) * sin(angle(Y_matrix(i, j)) - angle(V(i)) + angle(Y_matrix(i, j)));
        end
    end

    % 计算功率误差
    delta_P = P_calc(2:3) - P_load;  % PQ节点的有功功率误差
    delta_Q = Q_calc(2:3) - Q_load;  % PQ节点的无功功率误差
    F = [delta_P; delta_Q];
end

% 使用fsolve求解潮流方程
options = optimoptions('fsolve', 'Display', 'none', 'TolFun', 1e-6, 'MaxIterations', 100);
[x, fval, exitflag] = fsolve(@(x) power_flow_eq(x, Y_matrix, P_load, Q_load), x0, options);

% 提取结果
theta = [0; x(1:2)];     % 平衡节点相角为 0
V_mag = [1.05; x(3:4)];  % 平衡节点电压为 1.05
V = V_mag .* exp(1j * theta);  % 计算复数形式的节点电压

% 输出最终结果
disp('最终节点电压 (V):');
disp(V);

% 计算节点注入功率
S = zeros(3, 1);  % 节点注入功率
for i = 1:3
    for j = 1:3
        % 计算节点注入功率
        S(i) = S(i) + V(i) * conj(Y_matrix(i, j) * V(j));
    end
end

disp('节点注入功率 (S):');
disp(S);