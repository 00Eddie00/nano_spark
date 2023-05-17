import numpy as np

# 时间相关，单位：秒
DT = 2 * 10 ** -6  # DT 时间间隔
RELEASE_TIME = 0.02  # 释放时间，即到此刻RYR通道关闭
STOP_TIME = 0.1  # 结束时间

INITIAL_C_CA = 0.0001  # 肌质中钙离子初始浓度
D_CA = 3.5 * 10 ** 8
K_RYR = 1.22522 * 10 ** 10
V = 0.5 * np.pi * 5 ** 2
CA_JSR = 1
K_F_PLUS = 27000
K_F_MINUS = 17
F_T = 0.01



K_Calmodulin_PLUS = 100000
K_Calmodulin_MINUS = 31
Calmodulin_T = 0.036

K_TroponinC_PLUS = 125000
K_TroponinC_MINUS = 250
TroponinC_T = 0.07

K_SR_PLUS = 115000
K_SR_MINUS = 100
SR_T = 0.047

K_SL_PLUS = 115000
K_SL_MINUS = 1000
SL_T = 1.124
