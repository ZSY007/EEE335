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
% 打印导纳矩阵
disp('导纳矩阵 (Y):');
disp(Y);

syms P [1 3]  % 定义有功功率P为符号向量
syms Q [1 3]  % 定义无功功率Q为符号向量
% 计算功率平衡方程
for i = 1:3
    power = V(i) * conj(sum(Y_matrix(i, :) .* V.'));
    P(i) = real(power);
    Q(i) = -imag(power);
end

% 将有功功率和无功功率转换为小数形式
P_decimal = double(P);  % 转为小数
Q_decimal = double(Q);  % 转为小数

disp('有功功率注入 (P):');
disp(P_decimal);
disp('无功功率注入 (Q):');
disp(Q_decimal);

% 计算雅可比矩阵
    J = zeros(4,4);
   
    % 计算雅可比矩阵的元素
    for i = 2:3
        for j = 2:3
            if i == j
                % 对角元素
                J(i-1,j-1) = -Q(i) - V(i)^2*imag(Y(i,i));  % dP/dθ
                J(i+1,j+1) = (P(i) + V(i)^2*real(Y(i,i)))/V(i);  % dQ/dV
            else
                % 非对角元素
                J(i-1,j-1) = V(i)*V(j)*abs(Y(i,j))*sin(theta(i)-theta(j)-angle(Y(i,j)));  % dP/dθ
                J(i+1,j+1) = V(i)*abs(Y(i,j))*cos(theta(i)-theta(j)-angle(Y(i,j)));  % dQ/dV
            end
            % 交叉项
            J(i-1,j+1) = V(i)*abs(Y(i,j))*cos(theta(i)-theta(j)-angle(Y(i,j)));  % dP/dV
            J(i+1,j-1) = -V(i)*V(j)*abs(Y(i,j))*cos(theta(i)-theta(j)-angle(Y(i,j)));  % dQ/dθ
        end
    end
    % 打印雅可比矩阵
disp('雅可比矩阵 (J):');
disp(J);
