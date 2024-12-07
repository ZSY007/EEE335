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
P_load = [-1.0; -0.8];  % 节点2和3的有功负荷
Q_load = [-0.5; -0.4];  % 节点2和3的无功负荷

% 初始电压幅值和相角
V0 = [1.05; 1.02; 1.01];  % 初始电压幅值
theta0 = [0; 0; 0];       % 初始相角
x0 = [theta0(2:end); V0(2:end)];  % 初始变量 (PQ节点)

% Newton-Raphson 方法参数
tolerance = 1e-4;  % 误差停止标准
max_iterations = 100;  % 最大迭代次数

% Newton-Raphson 迭代
x = x0;  % 初始值
n_PQ = length(P_load);
for iter = 1:max_iterations
    % 当前电压和相角
    theta = [0; x(1:n_PQ)];  % 节点1为参考点
    V = [V0(1); x(n_PQ+1:end)];
    
    % 计算功率
P_calc = zeros(3, 1);  % 初始化为列向量
Q_calc = zeros(3, 1);

for i = 1:3
    % 电压和导纳矩阵正确匹配计算
    S = V(i) * conj(sum(Y_matrix(i, :) .* (V.')));
    P_calc(i) = real(S);  % 提取有功功率
    Q_calc(i) = -imag(S); % 提取无功功率
end

    
    % 功率误差
    dP = P_load - P_calc(2:3);
    dQ = Q_load - Q_calc(2:3);
    mismatch = [dP; dQ];
    
    % 判断收敛性
    if max(abs(mismatch)) < tolerance
        fprintf('潮流计算在第 %d 次迭代后收敛。\n', iter);
        break;
    end
    
    % 构造雅可比矩阵
    J11 = zeros(n_PQ, n_PQ);  % ∂P/∂θ
    J12 = zeros(n_PQ, n_PQ);  % ∂P/∂|V|
    J21 = zeros(n_PQ, n_PQ);  % ∂Q/∂θ
    J22 = zeros(n_PQ, n_PQ);  % ∂Q/∂|V|
    
    for i = 2:3
        for j = 2:3
            if i == j
                J11(i-1, j-1) = -imag(V(i) * conj(sum(Y_matrix(i, :) .* V .* exp(1j * theta))));
                J12(i-1, j-1) = abs(V(i)) * real(sum(Y_matrix(i, :)));
                J21(i-1, j-1) = real(V(i) * conj(sum(Y_matrix(i, :) .* V .* exp(1j * theta))));
                J22(i-1, j-1) = -abs(V(i)) * imag(sum(Y_matrix(i, :)));
            else
                theta_diff = theta(i) - theta(j);
                J11(i-1, j-1) = -abs(V(i)) * abs(V(j)) * ...
                    (real(Y_matrix(i, j)) * sin(theta_diff) - imag(Y_matrix(i, j)) * cos(theta_diff));
                J12(i-1, j-1) = -abs(V(j)) * ...
                    (real(Y_matrix(i, j)) * cos(theta_diff) + imag(Y_matrix(i, j)) * sin(theta_diff));
                J21(i-1, j-1) = -abs(V(i)) * abs(V(j)) * ...
                    (-real(Y_matrix(i, j)) * cos(theta_diff) - imag(Y_matrix(i, j)) * sin(theta_diff));
                J22(i-1, j-1) = -abs(V(j)) * ...
                    (real(Y_matrix(i, j)) * sin(theta_diff) - imag(Y_matrix(i, j)) * cos(theta_diff));
            end
        end
    end
    
    % 雅可比矩阵
    J = [J11, J12; J21, J22];
    
    % 修正变量
    delta = J \ mismatch;
    x = x + delta;  % 更新
    
end

if iter == max_iterations
    disp('潮流计算未能收敛到设定误差范围内。');
end

% 结果输出
disp('最终节点电压 (幅值和相角)：');
for i = 1:3
    fprintf('节点 %d: |V| = %.4f, θ = %.4f°\n', i, abs(V(i)), angle(V(i)) * 180/pi);
end
