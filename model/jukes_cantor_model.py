import numpy as np

def estimate_jukes_cantor_distance(p_hat: float) -> float:
    """estimate Jukes-Cantor distance from match probability"""
    if p_hat >= 1:
        p_hat = 0.999  # prevent log(0) when p_hat = 1
    try:
        d = -3/4 * np.log(1 - 4/3 * (1 - p_hat))
        return max(0, d)  # prevent negative distances
    except ValueError:
        return 0  # return 0 if the formula fails