from utils.equations import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import math

tolerance = 0.01
def main():
    ims = []
    fig, ax = plt.subplots()

    ## Create grid + solid propellant shape

    # N*N grid, N odd
    N = 51
    Nx = N
    Ny = N

    centerx = math.floor(Nx/2)
    centery = math.floor(Ny/2)

    prop = np.zeros((Nx,Ny))

    # Create initial propellant shape
    max_combustion = 10
    initial_combustion = 1

    # Star shape
    n_star = 15
    for m in range(n_star):
        prop[centerx + m, centery] = initial_combustion
        prop[centerx, centery - m] = initial_combustion
        prop[centerx - m, centery] = initial_combustion
        prop[centerx, centery + m] = initial_combustion

    ## Plot initial condition
    im = plot_propellant_surface(ax, prop, max_combustion, "Initial condition")
    plot_propellant_surface(ax, prop, max_combustion, "Initial condition", False)
    ims.append([im])

    ## Define combustion and simulation parameters
    num_steps = 100
    combustion_rate = 1
    burning_threshold = 8
    delta_t = 1

    ## Simulate
    current_time = 0

    for _ in range(num_steps):
        current_time += delta_t
        prop = simulate_next_step(prop, delta_t, combustion_rate, max_combustion, burning_threshold, initial_combustion)
        im = plot_propellant_surface(ax, prop, max_combustion, f"Propellant surface at time t = {current_time}")
        ims.append([im])

    

    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                    repeat_delay=1000)
    ani.save("movie.mp4")

    plt.show()



def plot_propellant_surface(ax, prop, max_combustion, title = "", animated = True):
    
    im = ax.imshow(prop, interpolation="bilinear", cmap = "hot", vmin = 0, vmax = max_combustion, animated = animated)
    ax.set_title(title)
    
    ax.set_yticklabels([])
    ax.set_xticklabels([])

    if animated:
        return im

def simulate_next_step(prop, delta_t, combustion_rate, max_combustion, burning_threshold, initial_combustion):
    new_prop = np.zeros_like(prop)
    for i in range(prop.shape[0]):
        for j in range(prop.shape[1]):
            p = prop[i,j]
    
            if p > tolerance: # Burning element
                # # Burn more
                new_prop[i,j] = min(p + combustion_rate, max_combustion)

                if p >= burning_threshold:
                    # Compute area of effect
                    combustion_radius = math.floor(combustion_rate * delta_t) 
                    
                    # Ignite neighbours in the area of effect
                    neighbours = compute_non_burning_neighbours(prop, i, j, combustion_radius)
                    for [u,v] in neighbours:
                        new_prop[u,v] =  initial_combustion


    return new_prop


def compute_non_burning_neighbours(prop, i, j, combustion_radius):
    neighbours = []
    min_i = max(i - combustion_radius, 0)
    max_i = min(i + combustion_radius, prop.shape[0] - 1)
    range_i = range(min_i, max_i+1)

    min_j = max(j - combustion_radius, 0)
    max_j = min(j + combustion_radius, prop.shape[1] - 1)
    range_j = range(min_j, max_j+1)


    for u in range_i:
        for v in range_j:
            if ((abs(i - u) + abs(j - v)) <= combustion_radius):
                if prop[u,v] < tolerance:
                    neighbours.append([u,v])

    return neighbours




if __name__ == '__main__':
    main()


