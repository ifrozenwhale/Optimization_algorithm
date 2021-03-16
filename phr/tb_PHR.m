ep = 1e-4;
x = sym('x', [1 3]);
f(x) = 0.2 * (x(1))^2 + 0.2 * (x(2))^2 + 0.2 * (x(3))^2 + 310 * x(1) + 305 * x(2) + 300 * x(3) - 1000;
h = symfun([240 - x(1) - x(2) - x(3)], x);
g = symfun([60 - x(1); 140 - x(1) - x(2)], x);
x0 = [0; 0; 0];

fprintf("\nResult by PHR method:\n")
[minx, min_value, arr] = PHR(f, h, g, x0, ep, true)    

% plot
xp = 1:length(arr);
figure()
plot(xp, arr, '-p');
legend('目标函数值')
title('目标函数值收敛曲线')
print(gcf, '-r600', '-dpng', 'opt.png');