import numpy as np


def integral(disp: np.array, force: np.array) -> np.array:
    ret = np.zeros(disp.size)

    total_volume: float = 0
    for i in range(disp.size - 1):
        volume = (disp[i+1] - disp[i]) * 0.5 * (force[i+1] + force[i])
        total_volume += volume
        ret[i+1] = total_volume

    return ret

def derivative(x_data:np.array, y_data:np.array) -> np.array:
    return np.diff(y_data) / np.diff(x_data)

