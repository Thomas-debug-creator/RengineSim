from utils.equations import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import math

tolerance = 0.01
def main():
    ims = []
    fig, ax = plt.subplots()

    ## Create background grid
    print("Creating background grid")

    # N*N grid, N odd
    N = 31
    Nx = N
    Ny = N

    centerx = math.floor(Nx/2)
    centery = math.floor(Ny/2)

    x = np.linspace(-1,1,Nx)
    y = np.linspace(-1,1,Ny)

    prop = np.zeros((Nx,Ny))

    ## Propellant casing (non burning surface)
    non_combustion_radius = 0.9

    prop_indices = []

    for i in range(Nx):
        for j in range(Ny):
            distance_to_center = math.sqrt((x[centerx] - x[i])**2 + (y[centery] - y[j])**2)

            if distance_to_center <= non_combustion_radius:
                prop_indices.append([i,j])
    

    ## Combustion parameters
    max_combustion = 10
    initial_combustion = 1
    combustion_rate = 1
    combustion_propagation_rate = 0.2
    combustion_threshold = 2

    # Propellant shape (star)
    print("Creating propellant shape")
    n_star = 7
    for m in range(n_star):
        prop[centerx + m, centery] = initial_combustion
        prop[centerx, centery - m] = initial_combustion
        prop[centerx - m, centery] = initial_combustion
        prop[centerx, centery + m] = initial_combustion

    for m in range(n_star - 2):
        prop[centerx + m, centery + m] = initial_combustion
        prop[centerx + m, centery - m] = initial_combustion
        prop[centerx - m, centery - m] = initial_combustion
        prop[centerx - m, centery + m] = initial_combustion

    ## Plot initial condition
    im = plot_propellant_surface(ax, prop, max_combustion, "Initial condition")
    plot_propellant_surface(ax, prop, max_combustion, "Initial condition", False)
    ims.append([im])

    ## Simulation
    num_steps = 15
    delta_t = 1
    current_time = 0

    nb_burning_elements = [sum(prop[i,j] > tolerance for i in range(Nx) for j in range(Ny))]

    # Compute combustion radius
    print("Computing ignition radius")
    combustion_radius = combustion_propagation_rate * delta_t

    ref_ignition_indices = [[u - centerx,v - centery] for u in range(Nx) for v in range(Ny) if (math.sqrt((x[centerx] - x[u])**2 + (y[centery] - y[v])**2) <= combustion_radius)]

    # Simulate steps
    print("Simulating")

    for n in range(num_steps):
        print(f'Simulating step {n}')
        current_time += delta_t
        prop = simulate_next_step(prop, max_combustion, initial_combustion, combustion_rate, combustion_threshold, ref_ignition_indices, prop_indices)

        im = plot_propellant_surface(ax, prop, max_combustion, f"Propellant surface at time t = {current_time}")
        plot_propellant_surface(ax, prop, max_combustion, f"Propellant surface at time t = {current_time}")
        ims.append([im])


        nb_burning_elements.append(sum(prop[i,j] > tolerance for i in range(Nx) for j in range(Ny) if (prop[i,j] < max_combustion)))

    
    ## Animated gif
    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                    repeat_delay=1000)
    ani.save("movie.mp4")

    plt.show()

    plt.plot(nb_burning_elements)
    plt.show()



def plot_propellant_surface(ax, prop, max_combustion, title = "", animated = True):
    
    im = ax.imshow(prop, interpolation="bilinear", cmap = "hot", vmin = 0, vmax = max_combustion, animated = animated)
    ax.set_title(title)
    
    ax.set_yticklabels([])
    ax.set_xticklabels([])

    if animated:
        return im

def simulate_next_step(prop, max_combustion, initial_combustion, combustion_rate, combustion_threshold, ref_ignition_indices, prop_indices):
    new_prop = np.zeros_like(prop)
    burning_elements_indices = np.argwhere(prop > tolerance)

    for [i,j] in burning_elements_indices:
        p = prop[i,j]
    
        ## Burn more
        new_prop[i,j] = min(p + combustion_rate, max_combustion)
        p = new_prop[i,j]
        
        ## Ignite other elements if burning enough
        if p >= combustion_threshold:                    
            # Ignite neighbours in the area of effect
            neighbours = compute_non_burning_neighbours(prop, i, j, ref_ignition_indices, prop_indices)
            for [u,v] in neighbours:
                new_prop[u,v] =  initial_combustion


    return new_prop


def compute_non_burning_neighbours(prop, i, j, ref_ignition_indices, prop_indices):
    neighbours = [] 

    ignition_indices = []

    for [u,v] in ref_ignition_indices:
        new_u = np.clip(i + u, 0, prop.shape[0] - 1)
        new_v = np.clip(j + v, 0, prop.shape[1] - 1)

        if ([new_u, new_v] in prop_indices):
            ignition_indices.append([new_u, new_v])

    for [u,v] in ignition_indices:
        if prop[u,v] < tolerance:
            neighbours.append([u,v])

    return neighbours




if __name__ == '__main__':
    main()


