# PHR算法

## 背景

约束优化 (Constrained Optimization)，即约束优化问题，是优化问题的分支。它是在一系列约
束条件下，寻找一组参数值，使某个或某一组函数的目标值达到最优。其中约束条件既可以是等
式约束也可以是不等式约束。寻找这一组参数值的关键可是：满足约束条件和目标值要达到最
优。求解约束问题的方法可分为传统方法和进化算法。

本次实验采用传统算法 PHR 算法。混合约束问题的乘子法，该方法是 Rockafeller 在 PH 算法
的基础上提出的。方法结合了拉格朗日乘子和罚函数，构造无约束的辅助函数。

使用matlab实现。

## 结构

### 文件说明

| 代码编号 | 文件名          | 功能描述                              |
| -------- | --------------- | ------------------------------------- |
| 1        | PHR.m           | 求解约束混合优化问题的 PHR 算法函数   |
| 2        | dfp.m           | 求解无约束混合优化问题的 dpf 算法函数 |
| 3        | golden_search.m | 一维搜索优化问题的黄金分割法函数      |
| 4        | tb_PHR.m        | 问题输入测试                          |

### 输入说明

| 变量   | 意义             | 说明       | 示例                               |
| ------ | ---------------- | ---------- | ---------------------------------- |
| f      | 目标函数         | 符号函数   | f(x) = (x(1) - 2)^2 + (x(2) - 3)^2 |
| x0     | 初始解           | 1×n维数组  | x0 = [0; 0]                        |
| ep     | 精度             | double类型 | ep = 1e-4                          |
| max_k  | 最大迭代次数     | int类型    | max_k = 200                        |
| detail | 输出中间迭代过程 | bool类型   | detail = true                      |

### 输出说明

| 变量 | 意义   | 说明       |
| ---- | ------ | ---------- |
| x    | 最优解 | 1×n维数组  |
| v    | 最优值 | double类型 |



## 实现

> 截自自己latex完成的文档

![image-20210317000203899](https://frozenwhale.oss-cn-beijing.aliyuncs.com/img/image-20210317000203899.png)

![image-20210317000215173](https://frozenwhale.oss-cn-beijing.aliyuncs.com/img/image-20210317000215173.png)

![image-20210317000251885](https://frozenwhale.oss-cn-beijing.aliyuncs.com/img/image-20210317000251885.png)
