% 定义系统参数
% 线路阻抗
Z12 = 0.04 + 0.08j;
Z23 = 0.05 + 0.1j;
Z13 = 0.06 + 0.12j;

% 计算导纳
Y12 = 1/Z12;
Y23 = 1/Z23;
Y13 = 1/Z13;

% 构建导纳矩阵 Y
Y = zeros(3,3);
Y(1,1) = Y12 + Y13;
Y(1,2) = -Y12;
Y(1,3) = -Y13;
Y(2,1) = -Y12;
Y(2,2) = Y12 + Y23;
Y(2,3) = -Y23;
Y(3,1) = -Y13;
Y(3,2) = -Y23;
Y(3,3) = Y13 + Y23;

% 初始条件
V1 = 1.05;
theta1 = 0;
V = [V1; 1; 1];  % 初始电压幅值
theta = [theta1; 0; 0];  % 初始相角(弧度)

% 负荷功率
P2 = -1.0; Q2 = -0.5;  % 负号表示负荷
P3 = -0.8; Q3 = -0.4;


% Gauss-Seidel方法
fprintf('\nGauss-Seidel方法求解:\n')
max_iter = 100;
epsilon = 1e-4;
iter = 0;

while iter < max_iter
    iter = iter + 1;
    V_old = V;
    
    % 更新母线2电压
    sum2 = 0;
    for k = [1,3]
        sum2 = sum2 + Y(2,k)*V(k);
    end
    V(2) = (1/Y(2,2))*((P2-1j*Q2)/conj(V(2)) - sum2);
    
    % 更新母线3电压
    sum3 = 0;
    for k = [1,2]
        sum3 = sum3 + Y(3,k)*V(k);
    end
    V(3) = (1/Y(3,3))*((P3-1j*Q3)/conj(V(3)) - sum3);
    
    % 检查收敛性
    if max(abs(abs(V) - abs(V_old))) < epsilon
        break;
    end
end

fprintf('迭代次数: %d\n', iter)
fprintf('节点电压 (V):\n')
for i = 1:3
    fprintf('   %.4f + %.4fi\n', real(V(i)), imag(V(i)))
end
