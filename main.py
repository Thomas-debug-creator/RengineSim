from utils.grain_sim_functions import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import math

tolerance = 0.01

def main():
    ## Prepare image for the animated gif
    ims = []
    fig, ax = plt.subplots()

    ## Initialize simulation space
    print("Initializing simulation space")

    N = 31 # WARNING: N MUST BE ODD
    Nx = N
    Ny = N
    centerx, centery, x, y = initialize_sim_space(Nx, Ny)

    ## Propellant casing (non burning surface)
    print("Initializing grain casing")
    grain_radius = 0.9
    prop_indices = create_grain_casing(x, y, centerx, centery, grain_radius)
    

    ## Combustion parameters
    max_combustion = 10
    initial_combustion = 1
    combustion_rate = 1
    combustion_propagation_rate = 0.2
    ignition_threshold = 2

    # Propellant initial grain configuration (star)
    print("Creating initial grain configuration")
    grain_config = "internal_tube"
    prop = create_grain_config(Nx, Ny, centerx, centery, x, y, initial_combustion, grain_config)

    ## Plot initial condition
    im = plot_propellant_surface(ax, prop, max_combustion, "Initial condition")
    plot_propellant_surface(ax, prop, max_combustion, "Initial condition", False)
    ims.append([im])

    ## Simulation
    num_steps = 15
    delta_t = 1
    current_time = 0

    # Compute initial number of burning elements for final graph
    nb_burning_elements = []
    nb_burning_elements.append(compute_nb_burning_elements(prop, Nx, Ny, max_combustion))

    # Compute ignition radius
    print("Computing ignition radius")
    ref_ignition_indices = compute_ignition_radius(combustion_propagation_rate, delta_t, Nx, Ny, centerx, centery, x, y)
    
    # Simulate steps
    print("Simulating")

    for n in range(num_steps):
        print(f'Simulating step {n}')
        current_time += delta_t
        prop = simulate_next_step(prop, max_combustion, initial_combustion, combustion_rate, ignition_threshold, ref_ignition_indices, prop_indices)

        im = plot_propellant_surface(ax, prop, max_combustion, f"Propellant surface evolution over {num_steps} steps")
        ims.append([im])


        nb_burning_elements.append(compute_nb_burning_elements(prop, Nx, Ny, max_combustion))

    
    ## Plots
    process_gif(fig, ims, grain_config)
    plot_combustion_surface(num_steps, nb_burning_elements, grain_config)







if __name__ == '__main__':
    main()


