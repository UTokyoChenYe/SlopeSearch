from model.F_k_function import F_k_funtion
from model.evolution_models import estimate_jukes_cantor_distance

def compute_distance(seq1, seq2, args):
    """
    Compute the distance between two sequences using the upper model.
    """
    # 1. calculate p_hat
    F_k_function_obj = F_k_funtion(seq1, seq2, args)
    p_hat = F_k_function_obj.calculate_p_hat()

    # 2. calculate distance
    distance = estimate_jukes_cantor_distance(p_hat)

    return distance
    