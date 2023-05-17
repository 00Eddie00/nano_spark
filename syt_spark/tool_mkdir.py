import os


def mkdir(path):

    path = path.strip()  # 去除首位空格
    path = path.rstrip("\\")  # 去除尾部 \ 符号

    is_exists = os.path.exists(path)     # 判断路径是否存在，存在：True，不存在：False

    # 判断结果
    if not is_exists:
        # 如果不存在则创建目录
        os.makedirs(path)  # 创建目录操作函数
        print(path + ' 创建成功')
        return True
    else:
        # 如果目录存在则不创建，并提示目录已存在
        print(path + ' 目录已存在')
        return False
