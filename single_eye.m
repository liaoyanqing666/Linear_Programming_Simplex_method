% writen by 21122890黎思源
format rational

c = [1, 3, 1, 0, 0];
b = [15; 20];
A = [2, -3, 2, 1, 0;
    1/3, 1, 5, 0, 1;];

[x_opt, fx_opt, iter] = Simplex_eye(A,b,c);
if fx_opt == Inf
    disp("无界解")
else
    disp("解为")
    disp(x_opt)
    fprintf("最优函数值为%d\n\n", fx_opt)
end


function [x_opt, fx_opt, iter] = Simplex_eye(A,b,c) 
% 其中x_opt为最优解，fx_opt为最优函数值，iter为迭代次数
    m = size(A, 1); % 行数（秩）
    n = size(A, 2); % 列数

    % 寻找单位矩阵并计算初始值
    row = Find_eye(n, m, A); % 基的下标
    cb = c(row);
    check = c - cb * A;
    iter = 0;
    
    % 开始迭代
    while ~all(check <= 0) % 检验数不全小于0
        iter = iter + 1;
        [max_num, max_index] = max(check); % 入基变量为max_index
        theta_index = 0;
        theta = Inf;
        for i = 1:m
            if A(i, row(i)) > 0 && b(i) / A(i, row(i)) <= theta
                theta_index = i;
                theta = b(i) / A(i, row(i));
            end
        end
        if theta_index == 0 % 此时存在无界解
            x_opt = [];
            fx_opt = Inf;
            return
        end
        row(theta_index) =  max_index; % 换基迭代
        
        % 线性变换
        A(theta_index, :) = A(theta_index, :) / A(theta_index, max_index);
        b(theta_index) = b(theta_index) / A(theta_index, max_index);
        for i = 1:m
            if i == theta_index
                continue
            end
            temp = A(i, max_index);
            A(i, :) = A(i, :) - temp * A(theta_index, :);
            b(i) = b(i) - temp * b(theta_index);
        end

        % 更新数据
        cb = c(row);
        check = c - cb * A;
    end
    x_opt = zeros(1, n);
    x_opt(row) = b;
    x_opt = x_opt';
    fx_opt = c * x_opt;
end


function row = Find_eye(n, m, A) % 寻找单位矩阵
    row = zeros(1, m);
    for i = 1:m
        F = zeros(m,1);
        F(i) = 1;
        for j = 1:n
            if A(:,j) == F
                row(i) = j;
                break;
            end
        end
    end
end