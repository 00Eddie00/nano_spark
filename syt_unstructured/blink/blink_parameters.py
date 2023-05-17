import numpy as np

# Jdiffusion扩散项
DCAJSR = 3.5 * 10 ** 8  # 终池jSR 自由Ca离子扩散系数
DF = 2 * 10 ** 7  # 荧光指示剂染料 扩散系数

# Jbuff
BCSQ = 14.0  # 方程中的B，所有肌钙集蛋白（肌浆网SR上调控钙储存的主要结合蛋白，调节SR钙释放）的量
KDCSQ = 0.63  # 方程中的K，肌钙集蛋白的Ca离子解离常数

# Jrefill,Jryr两项
# 扩散系数
DCARYR = 0.25 * 6.5 * 10 ** 7  # RyR通道 Ca离子扩散系数 KRyR 6.5 * 10 ** 7
DCAFSR = 0.7854 * 10 ** 6  # 连接的fSR通道 Ca离子扩散系数 KfSR
# 浓度
CCAMYO = 0.0001  # [Ca2+]cyto  胞浆的Ca浓度
CCAFSR = 1.0  # [Ca2+]fSR  fSR通道中Ca离子浓度

# Jdye染料
F = 0.1  # [F]T   fluo-3总浓度
K1 = 48800  # KF+  Ca离子和fluo-3的结合速率
K2 = 19520  # KF-  Ca离子和fluo-3的解离速率
# K1 = 0
# K2 = 0

# 关闭RyR通道用到
DT = 2 * 10 ** -6  # dt
RELEASE_TIMES = 2 * 10 ** -2  # 0.02s,这个就是ryr通道开放的时间，后80毫秒是恢复的

# 点的属性
B_INNER = 0  # 内部
B_INFLOW = 2  # 入流
B_WALL = 1  # 壁面
B_OUTFLOW = 4  # 出流
B_SYMMETRY = 8
B_SOURCE = 16

grid_file_name = "../config/blink/gridt.dat"
nod_file_name = "../config/blink/nod.dat"
npoch_file_name = "../config/blink/npoch.dat"
grid = np.loadtxt(grid_file_name)
nod = np.loadtxt(nod_file_name) - 1
npoch = np.loadtxt(npoch_file_name)
NP = len(grid)
NE = len(nod)
