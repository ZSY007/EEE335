% 牛顿-拉夫森潮流计算

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


% 初始电压值
V = [1.05; 1.02; 1.01];  % 初始电压幅值
theta = [0; 0; 0];       % 初始相角

% 节点功率
P_load = [-1.0; -0.8];  % 有功负荷
Q_load = [-0.5; -0.4];  % 无功负荷

% 牛顿-拉夫森参数
tolerance = 1e-6;
max_iter = 100;

% 开始迭代
for iter = 1:max_iter


    % 判断收敛
    if norm(mismatch) < tolerance
        disp('潮流计算收敛');
        break;
    end

    % 构造雅可比矩阵
    J = zeros(4);
    for i = 2:3
        for j = 2:3
            if i == j
                % 对角线元素
                J(i-1, j-1) = -imag(V(i) * conj(Y_matrix(i, :) * V));
                J(i+1, j+1) = real(V(i) * conj(Y_matrix(i, :) * V));
            else
                % 非对角线元素
                J(i-1, j-1) = -abs(V(i)) * abs(V(j)) * imag(Y_matrix(i, j));
                J(i+1, j+1) = -abs(V(i)) * abs(V(j)) * real(Y_matrix(i, j));
            end
        end
    end

    % 修正电压
    delta = J \ mismatch;
    theta(2:end) = theta(2:end) + delta(1:2);
    V(2:end) = V(2:end) + delta(3:4);
end

if iter == max_iter
    disp('潮流计算未收敛');
else
    disp(['潮流计算在 ', num2str(iter), ' 次迭代后收敛']);
end
% 输出迭代次数
disp(['总迭代次数: ', num2str(iter)]);
% 输出最终节点电压
disp('最终节点电压 (V):');
for i = 1:length(V)
    fprintf('   %6.4f %+6.4fi\n', real(V(i)), imag(V(i)));
end

% 计算节点注入功率 (S)
S_injected = V .* conj(Y_matrix * V);

% 输出节点注入功率
disp('节点注入功率 (S):');
for i = 1:length(S_injected)
    fprintf('   %6.4f %+6.4fi\n', real(S_injected(i)), imag(S_injected(i)));
end