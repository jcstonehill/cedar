import cedar
import numpy as np


class UC_ZrC_C(cedar.base.Material):
    """
    Composite NERVA Fuel

    Reference:

    F. P. Durham, “Nuclear Engine Definition Study PreliminaryReport, Volume 1 -
    Engine Description,” Los Alamos National Laboratory, Report LA-5044-MS Vol 1,
    Los Alamos, NM, Sept. 1972.
    """

    _k_T_data = np.array(
        [0, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400,
         1600, 1800, 2000, 2200, 2400, 2600, 5000],
        dtype=float
    )
    _k_data = np.array(
        [90, 90, 71, 60, 53, 48, 44, 38, 36, 30,
         30, 30, 30, 30, 30, 30, 30],
        dtype=float
    )

    _cp_T_data = np.array(
        [0, 300, 400, 500, 600, 700, 800, 900, 1000, 1200,
         1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 5000],
        dtype=float
    )

    _cp_data = np.array(
        [448.194, 448.194, 576.464, 661.4, 714.034, 757.114, 790.406,
         809.458, 825.716, 857.043, 884.296, 905.349, 928.836,
         953.694, 980.986, 1007.792, 1033.661, 1060.81, 1060.81],
        dtype=float
    )

    def rho_rt(self):
        return 6700
    
    def k(self, T):
        return np.interp(T, self._k_T_data, self._k_data)

    def cp(self, T):
        return np.interp(T, self._cp_T_data, self._cp_data)