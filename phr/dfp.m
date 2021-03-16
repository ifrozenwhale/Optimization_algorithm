function x = dfp(f, x0, ep, max_k)
    %dfp - 无约束优化算法（拟牛顿法）
    %
    % Syntax: x = dpf(f, x0, ep, max_k)
    % f 目标函数
    % x0 初始解
    % ep 精度
    % max_k 最大迭代次数
    %
    % Initializing variables

    aux = sym('aux'); % Auxiliar symbolic variable for lambda
    x = symvar(f); % 目标函数变量x
    grad(x) = gradient(f); % 目标函数梯度
    [~, m] = size(x); % 向量长度
    x = x0; % 初始迭代点
    k = 0; % 迭代次数
    e = 100; % 范数
    x_aux = num2cell(x); % 转化为cell数组
    g = double(grad(x_aux{:})); % 梯度初值
    H = eye(m, m); % 一般初始化H为单位阵

    % 精度要求
    % 迭代次数限制
    while e > ep && k < max_k
        % 保存前值
        x_pre = x;
        g_pre = g;

        % 计算方向
        d = -H * g;

        % 一维搜索
        x = x_pre + aux * d;
        x_aux = num2cell(x);
        f_aux(aux) = f(x_aux{:});
        % 得到最优的lambda k
        % 使用0.618法一维搜索优化
        lambda = golden_search(f_aux, 0, 1, 1e-6);

        % 计算x(k+1)
        x = double(x_pre + lambda * d);
        x_aux = num2cell(x);
        g = double(grad(x_aux{:}));

        % 计算范数 ||df(x)||
        e = abs(norm(g));
        % If not satisfied, then update H
        if e > ep
            dx = double(x - x_pre);
            dg = double(g - g_pre);

            if dx' * dg <= 0
                H = eye(m, m);
            else
                H = H + ((dx * dx') / (dx' * dg)) - ((H * dg * dg' * H) / (dg' * H * dg));
            end

        end

        % 更新迭代
        k = k + 1;

    end

end