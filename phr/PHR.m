function [x, v, arr] = PHR(f, h, g, x0, ep, detail)
    %PHR - Description
    %
    % Syntax: x = PHR(f, h, g, x0, ep, detail)
    %
    % Long description
    % f 目标函数
    % h 等式约束
    % g 不等式约束
    % x0 初值
    % ep 精度
    % detail 是否输出中间迭代过程

    x = symvar(f); % 提取目标函数的符号
    x = x0; % 迭代初值

    lambda = zeros(size(g)); % 拉格朗日乘子lambda

    r = ones(size(h)); % 等式约束罚函数系数r
    mu = 2; % 罚函数因子mu
    k = 0; % 迭代次数k
    max_K = 100; % 最大迭代次数

    alpha = 2; % 变尺度因子
    beta = 0.5; % 检验系数
    phi_aux = 100; % phi值
    phi_pre = 100; % 前一次迭代的phi值

    % 精度要求
    % 迭代次数限制
    while phi_aux > ep & k < max_K

        % 保留上次迭代的值
        x_pre = x;

        % 转化为无约束子问题
        t = mu * g + lambda; % 临时保存
        % max转化为绝对值
        uncon_f = f + sum(r .* h) + (mu / 2) * sum(h.^2) + ...
            (1/2 / mu) * sum(((t + abs(-t) / 2)).^2 - lambda.^2);

        % 调用无约束优化拟牛顿法dfp算法
        x = double(dfp(uncon_f, x, ep, max_K));

        % 转化为cell数组
        x_aux = num2cell(x);

        % 得到约束方程的数值结果
        h_aux = double(h(x_aux{:}));
        g_aux = double(g(x_aux{:}));

        % 保存前一次迭代的phi值
        phi_pre = phi_aux;

        % 更新phi值
        phi = norm(h);
        phi_aux = double(phi(x_aux{:}));

        % 如果两次phi之比大于beta，则更新mu值

        if phi_aux / phi_pre > beta
            mu = alpha * mu;
        end

        % 更新lambda
        for i = 1:size(g)
            lambda = mu * g_aux(i) + lambda(i);
        end

        % 更新r
        for i = 1:size(h)
            r = r + mu * h_aux(i);
        end

        % 更新迭代次数k
        k = k + 1;

        % 输出中间迭代过程
        if detail == true
            result = ['Iteration ', num2str(k), ': '];
            % disp(result);
            % disp(x);
            fprintf(result);
            disp(x');
        end 

        % 输出目标函数值
        n = length(x);
        str = 'f(';

        for i = 1:n
            xx(i) = x(i);
        end

        for i = 1:n - 1
            str = [str, 'xx(', num2str(i), '),'];
        end

        str = [str, 'xx(', num2str(n), '))'];
        v = double(eval(str));
        % 保存每一次的函数值
        arr(k) = v;
    end

end
