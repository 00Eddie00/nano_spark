import numpy as np
import os


def shrink_files(position_con, multiple):

    new_position_con = [position_con[i] for i in range(0, len(position_con.shape[1]), multiple)]
    return np.array(new_position_con)
