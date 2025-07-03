
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import math

tolerance = 0.01

output_folder = 'output/'
def initialize_sim_space(Nx, Ny):
    "Nx and Ny must be odd!"
    centerx = math.floor(Nx/2)
    centery = math.floor(Ny/2)

    x = np.linspace(-1,1,Nx)
    y = np.linspace(-1,1,Ny)

    return centerx, centery, x, y

def create_grain_casing(x, y, centerx, centery, grain_radius = 0.9):
    "Return indices of the actual grain inside the casing assuming a circular shape"
    "Elements inside these indices can burn, while those outside are considered out of the casing and are therefore not considered in the simulation"

    prop_indices = []
    non_prop_indices = []
    default_config = np.zeros((x.shape[0], y.shape[0]))

    for i in range(x.shape[0]):
        for j in range(y.shape[0]):
            distance_to_center = math.sqrt((x[centerx] - x[i])**2 + (y[centery] - y[j])**2)

            if distance_to_center <= grain_radius:
                prop_indices.append([i,j])
            else:
                non_prop_indices.append([i,j])
                default_config[i,j] = -1
                
    return prop_indices, non_prop_indices, default_config

def create_grain_config(Nx, Ny, centerx, centery, prop_indices, non_prop_indices, x, y, max_combustion, grain_config = "star_8"):
    prop = np.zeros((Nx,Ny))

    if grain_config == "internal_tube_slots":
        tube_radius = 0.25
        star_radius = 0.55
        max_diag_diff = 2
        create_grain_config_internal_tube_slots(prop, centerx, centery, prop_indices, x, y, max_combustion, star_radius, tube_radius, max_diag_diff)
    elif grain_config == "internal_tube":
        tube_radius = 0.3
        create_grain_config_internal_tube(prop, centerx, centery, prop_indices, x, y, max_combustion, tube_radius)
    elif grain_config == "external_burning_rod":
        inner_tube_radius = 0.7
        create_grain_config_external_rod(prop, centerx, centery, prop_indices, x, y, max_combustion, inner_tube_radius)
    elif grain_config == "rod_and_tube":
        inner_tube_radius = 0.47
        outer_tube_radius = 0.53
        create_grain_config_rod_and_tube(prop, centerx, centery, prop_indices, x, y, max_combustion, inner_tube_radius, outer_tube_radius)
    elif grain_config == "double_anchor":
        anchor_inner_radius = 0.5
        anchor_outer_radius = 0.6
        anchor_pieces_max_angle = math.pi / 3
        max_center_diff = 1
        central_tube_radius = 0.1
        create_grain_config_double_anchor(prop, centerx, centery, prop_indices, x, y, max_combustion, anchor_inner_radius, anchor_outer_radius, anchor_pieces_max_angle, max_center_diff, central_tube_radius)


    for [i,j] in non_prop_indices:
        prop[i,j] = -1
    return prop

def create_grain_config_internal_tube_slots(prop, centerx, centery, prop_indices, x, y, max_combustion, star_radius, tube_radius, max_diag_diff):
    # Internal center tube
    for [i,j] in prop_indices:
        distance_to_center = math.sqrt((x[centerx] - x[i])**2 + (y[centery] - y[j])**2)

        if distance_to_center <= tube_radius:
            prop[i,j] = max_combustion

        # (i - centerx == j - centery or i - centerx == -(j - centery)):
            

        
        if distance_to_center <= star_radius:
            diag_diff_1 = abs(i - centerx - (j - centery))
            diag_diff_2 = abs(i - centerx + (j - centery))
            if (abs(i - centerx) <= max_diag_diff or abs(j - centery) <= max_diag_diff or diag_diff_1 <= max_diag_diff  or diag_diff_2 <= max_diag_diff):
                prop[i,j] = max_combustion

    

def create_grain_config_internal_tube(prop, centerx, centery, prop_indices, x, y, max_combustion, tube_radius):
    
    for [i,j] in prop_indices:
        distance_to_center = math.sqrt((x[centerx] - x[i])**2 + (y[centery] - y[j])**2)

        if distance_to_center <= tube_radius:
            prop[i,j] = max_combustion

def create_grain_config_external_rod(prop, centerx, centery, prop_indices, x, y, max_combustion, inner_tube_radius):
    for [i,j] in prop_indices:
        distance_to_center = math.sqrt((x[centerx] - x[i])**2 + (y[centery] - y[j])**2)

        if distance_to_center > inner_tube_radius:
            prop[i,j] = max_combustion

def create_grain_config_rod_and_tube(prop, centerx, centery, prop_indices, x, y, max_combustion, inner_tube_radius, outer_tube_radius):
    for [i,j] in prop_indices:
        distance_to_center = math.sqrt((x[centerx] - x[i])**2 + (y[centery] - y[j])**2)

        if distance_to_center > inner_tube_radius and distance_to_center <= outer_tube_radius:
            prop[i,j] = max_combustion

def create_grain_config_double_anchor(prop, centerx, centery, prop_indices, x, y, max_combustion, anchor_inner_radius, anchor_outer_radius, anchor_pieces_max_angle, max_center_diff, central_tube_radius):
    
    for [i,j] in prop_indices:
        distance_to_center = math.sqrt((x[centerx] - x[i])**2 + (y[centery] - y[j])**2)
        angle_to_center = math.atan((y[centery] - y[j])/(x[centerx] - x[i]))

        # Central tube
        if distance_to_center <= central_tube_radius:
            prop[i,j] = max_combustion

        # Straight line
        if distance_to_center <= anchor_inner_radius and abs(j - centery) <= max_center_diff:
            prop[i,j] = max_combustion
        
        # Angle pieces
        if distance_to_center > anchor_inner_radius and distance_to_center <= anchor_outer_radius:
            if abs(angle_to_center) <= anchor_pieces_max_angle:
                prop[i,j] = max_combustion





def compute_nb_burning_elements(prop, Nx, Ny, max_combustion):
    return sum(prop[i,j] > tolerance for i in range(Nx) for j in range(Ny) if prop[i,j] < max_combustion)

def compute_ignition_radius(combustion_propagation_rate, delta_t, Nx, Ny, centerx, centery, x, y):
    "Compute the local indices which can be ignited by a burning element"
    combustion_radius = combustion_propagation_rate * delta_t

    ref_ignition_indices = [[u - centerx,v - centery] for u in range(Nx) for v in range(Ny) if (math.sqrt((x[centerx] - x[u])**2 + (y[centery] - y[v])**2) <= combustion_radius)]

    return ref_ignition_indices

def simulate_next_step(prop, max_combustion, initial_combustion, combustion_rate, ignition_threshold, ref_ignition_indices, prop_indices, default_config):
    new_prop = default_config.copy()

    burning_elements_indices = np.argwhere(prop > tolerance)

    for [i,j] in burning_elements_indices:
        p = prop[i,j]
    
        ## Burn more
        next_p = min(p + combustion_rate, max_combustion)
        new_prop[i,j] = next_p

        ## Ignite other elements if burning enough
        if next_p >= ignition_threshold:                    
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

""" Plots """

def plot_propellant_surface(ax, prop, max_combustion, grain_config, animated = True):
    
    im = ax.imshow(prop, interpolation="bilinear", cmap = "hot", vmin = -1, vmax = max_combustion, animated = animated)
    ax.set_title(f'Combustion surface for the {grain_config} configuration')
    
    ax.set_yticklabels([])
    ax.set_xticklabels([])

    if animated:
        return im

def process_gif(fig, ims, grain_config):
    ani = animation.ArtistAnimation(fig, ims, interval=200, blit=True,
                                    repeat_delay=1000)
    ani.save(output_folder + f'{grain_config}.mp4')

    plt.show()

def plot_combustion_surface(num_steps, nb_burning_elements, grain_config):
    plt.plot(range(num_steps + 1), nb_burning_elements)
    plt.xlabel("Simulation step")
    plt.ylabel("Number of burning elements")
    plt.title(f'Burning elements for the {grain_config} configuration')

    plt.savefig(output_folder + f'Combustion_surface_{grain_config}')

    plt.show()