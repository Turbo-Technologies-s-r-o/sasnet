"""Gas properties for secondary air system calculations.

v0.1: calorically perfect gas (air by default) with Sutherland viscosity.
Real-gas properties via CoolProp can be plugged in later through the same
interface (R, gamma, cp, mu(T)).
"""
from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class PerfectGas:
    """Calorically perfect gas.

    Defaults are dry air: R = 287.05 J/(kg K), gamma = 1.4,
    Sutherland constants mu_ref = 1.716e-5 Pa s at T_ref = 273.15 K, S = 110.4 K.
    """

    R: float = 287.05
    gamma: float = 1.4
    mu_ref: float = 1.716e-5
    T_ref: float = 273.15
    S: float = 110.4

    @property
    def cp(self) -> float:
        return self.gamma * self.R / (self.gamma - 1.0)

    def mu(self, T: float) -> float:
        """Dynamic viscosity by Sutherland's law, Pa s."""
        return self.mu_ref * (T / self.T_ref) ** 1.5 * (self.T_ref + self.S) / (T + self.S)

    def rho(self, p: float, T: float) -> float:
        return p / (self.R * T)

    def critical_pressure_ratio(self) -> float:
        g = self.gamma
        return (2.0 / (g + 1.0)) ** (g / (g - 1.0))

    def flow_function(self, pr: float) -> float:
        """Isentropic nozzle flow function.

        Returns m * sqrt(R * T0) / (A * p0) for static-to-total pressure
        ratio pr, limited at the critical (choked) ratio.
        """
        g = self.gamma
        pr = max(pr, self.critical_pressure_ratio())
        pr = min(pr, 1.0)
        val = 2.0 * g / (g - 1.0) * (pr ** (2.0 / g) - pr ** ((g + 1.0) / g))
        return math.sqrt(max(val, 0.0))


AIR = PerfectGas()
