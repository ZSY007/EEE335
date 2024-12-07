clc;
clear;

% 导纳矩阵 (Y)
Y = [
    8.33 - 16.67j, -5 + 10j, -3.33 + 6.67j;
    -5 + 10j, 9 - 18j, -4 + 8j;
    -3.33 + 6.67j, -4 + 8j, 7.33 - 14.67j
];

% PQ节点功率
P_spec = [-1.0; -0.8]; % PQ节点有功功率
Q_spec = [-0.5; -0.4]; % PQ节点无功功率

% 初始电压和角度
V = [1.05; 1.0; 1.0]; % 节点电压幅值
theta = [0; 0; 0]; % 节点电压相角
tolerance = 1e-4; % 收敛准则
max_iter = 100; % 最大迭代次数
iter = 0; % 迭代计数器

% 牛顿-拉夫森迭代
while iter < max_iter
    iter = iter + 1;
    
    % 初始化功率误差
    P_calc = zeros(2, 1); % PQ节点的有功功率
    Q_calc = zeros(2, 1); % PQ节点的无功功率
    
    % 计算PQ节点的有功和无功功率
    for k = 2:3 % 仅计算PQ节点
        for m = 1:3
            P_calc(k-1) = P_calc(k-1) + V(k) * V(m) * ...
                (real(Y(k,m)) * cos(theta(k) - theta(m)) + imag(Y(k,m)) * sin(theta(k) - theta(m)));
            Q_calc(k-1) = Q_calc(k-1) + V(k) * V(m) * ...
                (real(Y(k,m)) * sin(theta(k) - theta(m)) - imag(Y(k,m)) * cos(theta(k) - theta(m)));
        end
    end
    
    % 计算功率误差
    delta_P = P_spec - P_calc;
    delta_Q = Q_spec - Q_calc;
    mismatch = [delta_P; delta_Q];
    
    % 检查收敛条件
    if max(abs(mismatch)) < tolerance
        break;
    end
    
    % 构造雅可比矩阵
    J = zeros(4, 4); % 初始化雅可比矩阵
    for k = 2:3 % PQ节点
        for m = 2:3
            % J1: ∂P/∂θ
            if k == m
                J(k-2+1, m-2+1) = -V(k)^2 * imag(Y(k,k)) - ...
                    sum(V(k) * V .* imag(Y(k,:)) .* cos(theta(k) - theta(:)'));
            else
                J(k-2+1, m-2+1) = V(k) * V(m) * ...
                    (real(Y(k,m)) * sin(theta(k) - theta(m)) - imag(Y(k,m)) * cos(theta(k) - theta(m)));
            end
            % J2: ∂P/∂V
            if k == m
                J(k-2+1, m-1+2) = V(k) * sum(real(Y(k,:)) .* cos(theta(k) - theta(:)')) + real(Y(k,k));
            else
                J(k-2+1, m-1+2) = V(k) * ...
                    (real(Y(k,m)) * cos(theta(k) - theta(m)) + imag(Y(k,m)) * sin(theta(k) - theta(m)));
            end
            % J3: ∂Q/∂θ
            if k == m
                J(k+1, m-2+1) = V(k)^2 * real(Y(k,k)) - ...
                    sum(V(k) * V .* real(Y(k,:)) .* sin(theta(k) - theta(:)'));
            else
                J(k+1, m-2+1) = V(k) * V(m) * ...
                    (-real(Y(k,m)) * cos(theta(k) - theta(m)) - imag(Y(k,m)) * sin(theta(k) - theta(m)));
            end
            % J4: ∂Q/∂V
            if k == m
                J(k+1, m-1+2) = -V(k) * sum(imag(Y(k,:)) .* cos(theta(k) - theta(:)')) - imag(Y(k,k));
            else
                J(k+1, m-1+2) = V(k) * ...
                    (imag(Y(k,m)) * cos(theta(k) - theta(m)) - real(Y(k,m)) * sin(theta(k) - theta(m)));
            end
        end
    end
    
    % 求解修正值
    delta = -J \ mismatch;
    
    % 更新电压幅值和相角
    theta(2:3) = theta(2:3) + delta(1:2);
    V(2:3) = V(2:3) + delta(3:4);
end

% 输出结果
disp('牛顿-拉夫森法收敛电压:');
disp(V);
disp('收敛相角 (弧度):');
disp(theta);
disp(['收敛迭代次数: ', num2str(iter)]);
