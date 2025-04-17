import numpy as np
from utils.logger import setup_logger

# 初始化日志记录器
logger = setup_logger()

def estimate_jukes_cantor_distance(p_hat: float) -> float:
    """estimate Jukes-Cantor distance from match probability"""
    logger.info(f"Estimating Jukes-Cantor distance for p_hat: {p_hat}")
    if p_hat >= 1:
        p_hat = 0.999  # prevent log(0) when p_hat = 1
        logger.warning("p_hat was >= 1, adjusted to 0.999 to prevent log(0)")
        print("p_hat was >= 1, adjusted to 0.999 to prevent log(0)")
    try:
        d = -3/4 * np.log(1 - 4/3 * (1 - p_hat))
        logger.info(f"Calculated distance: {d}")
        return max(0, d)  # prevent negative distances
    except ValueError as e:
        logger.error(f"Error calculating distance: {e}", exc_info=True)
        return 0  # return 0 if the formula fails