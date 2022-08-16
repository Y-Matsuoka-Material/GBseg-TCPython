import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from tc_python import *

# Configurations-------------------------------------------------------------------------------------------

# Definition of the target system
ELEMENTS = ["Co", "Cr", "Fe", "Ni", "Mn"]
COMPOSITIONS = [-1, 0.2, 0.2, 0.2, 0.2] # Composition in atomic fraction. The composition of the first element is dependent variable.

# Phases used in the calculation of chemical potentials in grain interior.
GRAIN_PHASES = ["FCC_L12"] # Set grain interior as FCC-single phase.
# GRAIN_PHASES = []  # Use default phases defined in the database.

# The database used in the calculation.
DATABASE = "TCHEA5"

# Temperature range of the calculation.
TEMPERATURES = np.linspace(600, 1500, 101)[::-1]

# Visualization of intermediate results during calculation.
VISUALIZATION = True
# --------------------------------------------------------------------------------------------------------

COMPOSITIONS[0] = 1 - sum(COMPOSITIONS[1:])
N_ELEMENTS = len(ELEMENTS)


# Penalty function
def J(comp_gb, mu_grain, calc, T):
    try:
        calc_result = calc.set_condition("T", T)
        for e, cc in zip(ELEMENTS[1:], comp_gb):
            calc_result = calc_result.set_condition(f"X({e})", cc)

        calc_result = calc_result.calculate()

        mu_gb = [calc_result.get_value_of(f"MU({e})") for e in ELEMENTS]

        result = 0
        for i in range(N_ELEMENTS - 1):
            result += ((mu_grain[0] - mu_grain[i + 1]) - (mu_gb[0] - mu_gb[i + 1])) ** 2
        return result
    except:
        return np.inf


with TCPython() as start:
    # Calculation of chemical potentials of grain interior----------------------------------------------------
    seg_calc = (
        start.
        select_database_and_elements(DATABASE, ELEMENTS)
    )

    if GRAIN_PHASES != []:
        seg_calc = seg_calc.without_default_phases()
        for p in GRAIN_PHASES:
            seg_calc = seg_calc.select_phase(p)

    seg_calc = (seg_calc.
                get_system().
                with_single_equilibrium_calculation().
                set_condition("P", 100000.0).
                set_condition("N", 1.0)
                )

    for e, c in zip(ELEMENTS[1:], COMPOSITIONS[1:]):
        seg_calc = seg_calc.set_condition(f"X({e})", c)

    mu_grain = []
    for T in TEMPERATURES:
        calc_result = seg_calc.set_condition("T", T).calculate()

        mu_grain.append([calc_result.get_value_of(f"MU({e})") for e in ELEMENTS])
    # --------------------------------------------------------------------------------------------------------

    print("First step completed!")

    # Calculation of grain boundary composition---------------------------------------------------------------
    seg_calc = (
        start.
        select_database_and_elements(DATABASE, ELEMENTS).
        without_default_phases().
        select_phase("LIQUID").
        get_system().
        with_single_equilibrium_calculation().
        set_condition("P", 100000.0).
        set_condition("N", 1.0).
        disable_global_minimization()
    )

    comps_gb = []
    prev_success = False
    for T, mu in zip(TEMPERATURES, mu_grain):
        # Search for grain boundary composition by Nelder-Mead method.
        if prev_success:
            opt_res = fmin(J, comp_result, (mu, seg_calc, T), full_output=True, disp=False)
        else:
            opt_res = fmin(J, COMPOSITIONS[1:], (mu, seg_calc, T), full_output=True, disp=False)
        comp_result, j_result = opt_res[0], opt_res[1]

        # Validation of the calculation result.
        if j_result > 1.0:
            print(f"At T = {T}, Residual error is too high (J = {j_result})")
            comps_gb.append(np.full_like(comp_result, np.nan))
            prev_success = False
        else:
            print(f"At T = {T}, J = {j_result}")
            comps_gb.append(comp_result)
            prev_success = True

        # Visualization of the intermediate results of a calculation.
        if VISUALIZATION:
            comps_gb_visualize = np.array(comps_gb)

            plt.clf()
            plt.plot(TEMPERATURES[:comps_gb_visualize.shape[0]], (1 - np.sum(comps_gb_visualize, axis=1)) * 100, label=ELEMENTS[0])
            for i in range(N_ELEMENTS - 1):
                plt.plot(TEMPERATURES[:comps_gb_visualize.shape[0]], comps_gb_visualize[:, i] * 100, label=ELEMENTS[i + 1])

            plt.legend()
            plt.xlim(TEMPERATURES.min(), TEMPERATURES.max())
            plt.ylim(0, 100)

            plt.xlabel("Temperature (K)")
            plt.ylabel("Concentration in GBP (at.%)")

            plt.tight_layout()
            plt.pause(0.001)
    # --------------------------------------------------------------------------------------------------------

    print("Second step completed!")
    comps_gb = np.array(comps_gb)

# Visualization of the calculation results----------------------------------------------------------------
plt.clf()
plt.plot(TEMPERATURES, (1 - np.sum(comps_gb, axis=1)) * 100, label=ELEMENTS[0])
for i in range(N_ELEMENTS - 1):
    plt.plot(TEMPERATURES, comps_gb[:, i] * 100, label=ELEMENTS[i + 1])

plt.legend()
plt.xlim(TEMPERATURES.min(), TEMPERATURES.max())
plt.ylim(0, 100)

plt.xlabel("Temperature (K)")
plt.ylabel("Concentration in GBP (at.%)")

plt.tight_layout()
plt.show()
# --------------------------------------------------------------------------------------------------------
