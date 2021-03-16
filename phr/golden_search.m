function lambda = golden_search(func, low, up, ep)
    % golden_search - 黄金分割一维搜索
    %
    % Syntax: lambda = golden_search(func, low, up)
    %
    % func 目标函数
    % low 下界
    % up 上届
    %
    func = matlabFunction(func); % 转化为函数句柄
    tow = (3 - sqrt(5)) / 2;
    N = ceil(log(ep) / log(1 - tow) + 3);

    x1 = (1 - tow) * low + tow * up;
    x2 = tow * low + (1 - tow) * up;
    v1 = func(x1) * 1.0;
    v2 = func(x2) * 1.0;

    for k = 4:N

        if v1 > v2
            low = x1;
            x1 = x2;
            v1 = v2;
            x2 = tow * low + (1 - tow) * up;
            v2 = func(x2);
        else
            up = x2;

            x2 = x1;
            v2 = v1;
            x1 = (1 - tow) * low + tow * up;
            v1 = func(x1);
        end

    end

    lambda = (x1 + x2) / 2;
end