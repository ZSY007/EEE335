% MATLAB代码计算三节点系统的导纳矩阵、功率平衡方程和雅可比矩阵

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

disp('导纳矩阵 (Y):');
disp(Y_matrix);

% 初始电压假设 (V)
V1 = 1.05 + 0j;  % Slack bus
V2 = 1.0 + 0j;   % PQ bus
V3 = 1.0 + 0j;   % PQ bus
V = [V1; V2; V3];  % 电压向量

% 初始化功率
P = zeros(3, 1);  % 有功功率
Q = zeros(3, 1);  % 无功功率

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

% 计算雅可比矩阵的各部分
n = length(V);
J11 = zeros(n - 1, n - 1);  % ∂P/∂θ
J12 = zeros(n - 1, n - 1);  % ∂P/∂|V|
J21 = zeros(n - 1, n - 1);  % ∂Q/∂θ
J22 = zeros(n - 1, n - 1);  % ∂Q/∂|V|

for i = 2:n  % 从第2个节点开始（跳过slack bus）
    for j = 2:n
        if i == j
            % 对角线元素
            J11(i-1, j-1) = -imag(V(i) * conj(sum(Y_matrix(i, :) .* V.')));
            J12(i-1, j-1) = abs(V(i)) * real(sum(Y_matrix(i, :)));
            J21(i-1, j-1) = real(V(i) * conj(sum(Y_matrix(i, :) .* V.')));
            J22(i-1, j-1) = -abs(V(i)) * imag(sum(Y_matrix(i, :)));
        else
            % 非对角线元素
            theta_diff = angle(V(i)) - angle(V(j));
            J11(i-1, j-1) = -abs(V(i)) * abs(V(j)) * (real(Y_matrix(i, j)) * sin(theta_diff) - imag(Y_matrix(i, j)) * cos(theta_diff));
            J12(i-1, j-1) = -abs(V(j)) * (real(Y_matrix(i, j)) * cos(theta_diff) + imag(Y_matrix(i, j)) * sin(theta_diff));
            J21(i-1, j-1) = -abs(V(i)) * abs(V(j)) * (-real(Y_matrix(i, j)) * cos(theta_diff) - imag(Y_matrix(i, j)) * sin(theta_diff));
            J22(i-1, j-1) = -abs(V(j)) * (real(Y_matrix(i, j)) * sin(theta_diff) - imag(Y_matrix(i, j)) * cos(theta_diff));
        end
    end
end

disp('雅可比矩阵 J11 (∂P/∂θ):');
disp(J11);
disp('雅可比矩阵 J12 (∂P/∂|V|):');
disp(J12);
disp('雅可比矩阵 J21 (∂Q/∂θ):');
disp(J21);
disp('雅可比矩阵 J22 (∂Q/∂|V|):');
disp(J22);
