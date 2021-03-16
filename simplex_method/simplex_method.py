
import math
import numpy as np


class Simplex:
    def __init__(self, m, P):
        """
        构造函数
        第一行为目标函数系数矩阵，目标值位于最后一列。
        下面是约束方程系数矩阵，全部为等式约束，最后一列常数项。

        :param m: 系数矩阵
        :param P: 初始基可行解
        """
        self.cn = m.shape[1]
        self.bn = m.shape[0]
        self.mat = m * 1.0
        self.P = P
        self.X = [0] * self.cn
        print("init")

    def display(self):
        print("Simple table")
        for i in range(self.mat.shape[0]):
            for j in range(self.mat.shape[1]):
                print("%10s " % str(np.round(self.mat[i, j], 4)), end='')
            print()
        print("result z = " + str(-np.round(self.mat[0, self.cn - 1], 4)))

    def update_x(self):
        """
        更新可行解列表 X
        :return: null
        """
        for j in range(self.cn):
            s = 0
            px = 0
            for i in range(self.bn):
                s += self.mat[i, j]
                if self.mat[i, j] == 1:  # 寻找元素1
                    px = i
            if s == 1:  # 整列只有一个元素为1
                self.X[j] = np.round(self.mat[px, self.cn - 1], 5)

    def pivot(self):
        """
        寻找转换变量
        :return: 新的基变量
        """
        x = y = 0
        c_max = -math.inf
        C = self.mat[0]
        B = []
        for i in range(self.bn):
            B.append(self.mat[i, self.cn - 1])

        # 寻找最大检验数，找到入基变量
        for i in range(C.shape[1] - 1):
            # if (c_max < C[0, i]) and (i not in self.P):
            if c_max < C[0, i]:
                c_max = C[0, i]
                y = i

        # 更新可行解列表
        self.update_x()
        if c_max <= 0:
            return 0

        # 判断是否无最优解
        for i in range(self.cn):
            flag = 0
            for j in range(1, self.bn):
                if C[0, j] > 0 and self.mat[j, i] <= 0:
                    flag += 1
            if flag == self.bn:
                return -1

        # 寻找离基变量
        b_min = math.inf
        for i in range(1, self.bn):
            # if B[i] <= 0:
            #     continue
            if self.mat[i, y] > 0:
                tmp = B[i] / self.mat[i, y]
                if b_min > tmp:
                    b_min = tmp
                    x = i

        # 新的基变量
        p = [x, y]
        # 移除离基变量
        for i in self.P:
            if self.mat[x, i] != 0:
                self.P.remove(i)
                break

        self.P.add(y)
        return p

    def Gaussian(self, p):
        """
        高斯消元法
        :param p: 新的基变量
        :return:
        """
        x = p[0]
        y = p[1]
        norm = self.mat[x, y]
        for i in range(self.cn):
            self.mat[x, i] /= norm  # 主行归一
        for i in range(self.bn):
            if i == x:
                continue
            if self.mat[i, y] != 0:
                tmp = self.mat[i, y]
                for j in range(self.cn):
                    self.mat[i, j] -= tmp * self.mat[x, j]

    def dual_pivot(self):
        """
        对偶单纯形法
        :return: 新的入基变量
        """
        x = y = 0
        B = []
        for i in range(0, self.bn):
            B.append(self.mat[i, self.cn - 1])
        b_min = min(B[1:])
        x = np.argmin(B)

        # update X
        self.update_x()

        if b_min >= 0:
            return 0

        for i in range(1, self.bn):
            flag = 0
            for j in range(self.cn):
                if B[i] < 0 and self.mat[i, j] >= 0:
                    flag += 1
            if flag == self.cn:
                return -1

        c_min = math.inf
        C = self.mat[0]
        for i in range(0, self.cn - 1):
            if self.mat[x, i] < 0:
                tmp = C[0, i] / self.mat[x, i]
                if c_min > tmp:
                    c_min = tmp
                    y = i
        p = [x, y]

        for i in self.P:
            if self.mat[x, i] != 0:
                self.P.remove(i)
                break
        self.P.add(y)
        return p

    def init_P(self):
        for p in P:
            for i in range(self.bn):
                if self.mat[i, p] == 1:
                    self.Gaussian([i, p])
                    break

    def judge_type(self):
        """
        判断是否使用对偶单纯形法
        如果所有的检验数都小于0，则将系数矩阵置相反
        :return: True为单纯形法， False为对偶单纯形法
        """
        if np.max(mat[0]) <= 0:
            for i in range(1, self.bn):
                if self.mat[i, self.cn - 1] > 0:
                    self.mat[i, :] *= -1
            return False
        return True

    def simple_simplex(self):
        """
        单纯形法的求解
        :return:
        """
        t = [0] * 2
        mat_type = self.judge_type()
        # 将基变量高斯消元
        self.init_P()

        while True:
            self.display()
            self.X = [0] * self.cn
            # 判断类型
            if mat_type:
                t = self.pivot()
            else:
                t = self.dual_pivot()

            # 无最优解的情况
            if t == -1:
                print("no solution")
                return

            if t == 0:
                print("result x = " + str(self.X))
                return
            print()
            print("pivot p = ", end='')
            print(str(t[0]) + " " + str(t[1]))
            print("base P = ", end='')
            for i in self.P:
                print(i, end=' ')
            print()
            self.Gaussian(t)

    def solve(self):
        self.simple_simplex()


if __name__ == '__main__':
    P = set()
    M = 10000
    # 当从键盘输入时

    # cn, bn = map(int, input().split())
    # mat = np.mat(np.zeros((bn, cn)))
    # for i in range(bn):
    #     tmp = list(map(int, input().split()))
    #     for j in range(cn):
    #         mat[i, j] = tmp[j]

    # 当矩阵从这里构造时

    # 测试对偶单纯性法，求min的时候
    # mat = np.mat([[-1.0, -2, 0, 0, 0],
    #               [1, 1, -1, 0, 1],
    #               [1, 1, 0, 1, 2]])

    # 一般测试
    # mat = np.mat([[-5, 0, -21, 0, 0, 0],
    #               [-1, 1, -6, 1, 0.0, -2],
    #               [-1, -1, -2, 0, 1, -1]])

    # 测试钢管下料
    # mat = np.mat([[-0, -0.1, -0.2, -0.3, -0.8, 0, 0, 0, 0],
    #               [1, 2, 0, 1, 0, -1, 0, 0,  100],
    #               [0, 0, 2, 2, 1, 0, -1, 0,  100],
    #               [3, 1, 2, 0, 3, 0, 0, -1,  100]])

    # 测试大M法求解的时候
    mat = np.mat([[2., 3, 0, 0, -M, 0],
                  [1, 3, 1, 0, 0, 5],
                  [1, 1, 0, -1, 1, 2]])

    # 添加初始基可行解
    P.add(2)
    P.add(4)

    # 构造并执行
    simplex = Simplex(mat, P)
    simplex.solve()
