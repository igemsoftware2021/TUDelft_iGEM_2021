
def micromolar_conc_to_math_exp(conc: float, decimals: int):
    """Function sets a float value that with micromolar units into a correct scientific math expression

    Parameters
    ----------
    conc: float
        concentration in micromolar
    decimals: int
        how many decimals you want

    Returns
    -------
    math_exp: str
        resulting mathematical expression

    Notes
    -----

    Example:
    concentration = 5.5*10**-3
    -> 5.5 x 10^-3 \mu M


    """
    if not isinstance(decimals, int):
        raise TypeError("decimals variable should be an integer")

    count = 0
    temp_conc = conc  # temporary variable to store a possible
    while temp_conc < 1 and temp_conc != 0:
        temp_conc = temp_conc * 10
        count += 1
    if count == 0:
        math_exp = f"${temp_conc:.1f}$ $\mu M$"
    elif count == 1:
        math_exp = f"${conc:.2f}$ $\mu M$"
    else:
        math_exp = f"${temp_conc:.{decimals}f} \\times 10^{{-{count}}}$ $\mu M$"
    return math_exp
