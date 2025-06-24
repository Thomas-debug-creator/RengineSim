import math

def compute_thrust(mass_change_rate, velocity_exhaust, pressure_exhaust, pressure_outside, area_exit):
    # Thrust = mdot * C
    return mass_change_rate * compute_effective_velocity(mass_change_rate, velocity_exhaust, pressure_exhaust, pressure_outside, area_exit)

def compute_effective_velocity(mass_change_rate, velocity_exhaust, pressure_exhaust, pressure_outside, area_exit):
    # C = ve + (Pe - Po) * Ae / mdot
    return velocity_exhaust + (pressure_exhaust - pressure_outside) * area_exit / mass_change_rate

def compute_specific_impulse(effective_velocity, gravity):
    # Isp = C / g
    return effective_velocity / gravity

def compute_weight_flow_rate(thrust, specific_impulse):
    # Wdot = Thrust / Isp
    return thrust / specific_impulse

def compute_delta_v(mass_full, mass_empty, gravity, specific_impulse, burn_time):
    # Rocket equation
    return math.log(mass_full / mass_empty) * gravity * specific_impulse - gravity * burn_time
    
