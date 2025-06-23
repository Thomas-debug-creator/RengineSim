
def thrust_equation(mass_change_rate, velocity_exhaust, pressure_exhaust, pressure_outside, area_exit):
    # Rocket thrust equation
    return mass_change_rate * velocity_exhaust + (pressure_exhaust - pressure_outside) * area_exit