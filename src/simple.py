import openmc
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import h5py
from mpl_toolkits.mplot3d import Axes3D

"""

Voxel plots not generating properly, k value needs to be way higher, would like to get better stats

"""


def create_model(rod_position=0.0, plot_geometry=False, output_dir='.'):
    """
    Create a simple model of a block of natural uranium in water with a single cadmium control rod.
    
    Parameters
    ----------
    rod_position : float
        Position of the control rod (0.0 = fully inserted, 1.0 = fully withdrawn)
    plot_geometry : bool, optional
        Whether to create and save geometry plots
    output_dir : str, optional
        Directory to save geometry plots
    
    Returns
    -------
    openmc.model.Model
        OpenMC model
    """
    # Materials
    # Natural uranium 
    uranium = openmc.Material(name='Natural uranium')
    uranium.set_density('g/cm3', 19.0)
    uranium.add_element('U', 1.0)  # Natural uranium
    
    # Water moderator
    water = openmc.Material(name='Water moderator')
    water.set_density('g/cm3', 1.0)
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')  # Thermal scattering data for hydrogen in water
    
    # Cadmium control rod
    cadmium = openmc.Material(name='Cadmium control rod')
    cadmium.set_density('g/cm3', 8.65)
    cadmium.add_element('Cd', 1.0)
    
    # Create a materials collection
    materials = openmc.Materials([uranium, water, cadmium])
    materials.export_to_xml()
    
    # Geometry
    # Dimensions
    pool_radius = 100.0  # cm
    pool_height = 150.0  # cm

    uranium_radius = 40.0  # cm
    uranium_height = 80.0  # cm
    
    control_rod_radius = 5.0  # cm 
    control_rod_height = 90.0  # cm
    
    # Outer boundary of the pool
    min_x = openmc.XPlane(-pool_radius, boundary_type='vacuum')
    max_x = openmc.XPlane(pool_radius, boundary_type='vacuum')
    min_y = openmc.YPlane(-pool_radius, boundary_type='vacuum')
    max_y = openmc.YPlane(pool_radius, boundary_type='vacuum')
    min_z = openmc.ZPlane(-pool_height/2, boundary_type='vacuum')
    max_z = openmc.ZPlane(pool_height/2, boundary_type='vacuum')
    
    pool_region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z
    
    # Simple uranium cylinder
    uranium_cyl = openmc.ZCylinder(r=uranium_radius)
    uranium_min_z = openmc.ZPlane(z0=-uranium_height/2)
    uranium_max_z = openmc.ZPlane(z0=uranium_height/2)
    uranium_region = -uranium_cyl & +uranium_min_z & -uranium_max_z
    
    # Control rod position
    # rod_position=0.0 means fully inserted, 1.0 means fully withdrawn
    rod_z_pos = (control_rod_height/2 + pool_height/2) * rod_position
    
    # Control rod
    control_cyl = openmc.ZCylinder(r=control_rod_radius)
    control_min_z = openmc.ZPlane(z0=-control_rod_height/2 + rod_z_pos)
    control_max_z = openmc.ZPlane(z0=control_rod_height/2 + rod_z_pos)
    control_region = -control_cyl & +control_min_z & -control_max_z
    
    # Create cells
    fuel_cell = openmc.Cell(name='Uranium Block', fill=uranium, region=uranium_region & ~control_region)
    control_cell = openmc.Cell(name='Control Rod', fill=cadmium, region=control_region & uranium_region)
    water_cell = openmc.Cell(name='Water', fill=water, region=pool_region & ~uranium_region)
    
    # Create universe
    root_universe = openmc.Universe(cells=[fuel_cell, control_cell, water_cell])
    
    # Create geometry
    geometry = openmc.Geometry(root_universe)
    geometry.export_to_xml()
    
    # Settings
    settings = openmc.Settings()
    settings.batches = 100
    settings.inactive = 30
    settings.particles = 10000
    settings.max_lost_particles = 10
        
    # Increase fission bank size to avoid warnings
    settings.fission_sites = {'size': 200000}
    
    # Create a source that's distributed in the uranium block
    source = openmc.IndependentSource(space=openmc.stats.Box(
        [-uranium_radius, -uranium_radius, -uranium_height/2],
        [uranium_radius, uranium_radius, uranium_height/2]
    ), energy=openmc.stats.Watt(a=0.988e6, b=2.249e-6))
    settings.source = source
    
    # Add entropy mesh for convergence
    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.lower_left = [-uranium_radius, -uranium_radius, -uranium_height/2]
    entropy_mesh.upper_right = [uranium_radius, uranium_radius, uranium_height/2]
    entropy_mesh.dimension = [15, 15, 15]
    settings.entropy_mesh = entropy_mesh
    
    settings.export_to_xml()
    
    # Tallies
    tallies = openmc.Tallies()
    
    # Mesh tally for flux visualization
    mesh = openmc.RegularMesh()
    mesh.lower_left = [-pool_radius, -pool_radius, -pool_height/2]
    mesh.upper_right = [pool_radius, pool_radius, pool_height/2]
    mesh.dimension = [50, 50, 30]
    
    mesh_filter = openmc.MeshFilter(mesh)
    tally = openmc.Tally(name='flux')
    tally.filters = [mesh_filter]
    tally.scores = ['flux', 'fission', 'absorption']
    tallies.append(tally)
    
    # Create a tally to compute k
    k_tally = openmc.Tally(name='k')
    k_tally.scores = ['nu-fission', 'absorption']
    tallies.append(k_tally)
    
    tallies.export_to_xml()
    
    # Create the model
    model = openmc.model.Model(geometry, materials, settings, tallies)
    
    # Generate geometry plots if requested
    if plot_geometry:
        # Create directory for plots if it doesn't exist
        plots_dir = 'geometry_plots'
        os.makedirs(plots_dir, exist_ok=True)
        
        # Colors for the plots
        colors = {
            water_cell: 'lightblue',
            fuel_cell: 'yellow',
            control_cell: 'red'
        }
        
        # XY plot (top view)
        xy_plot = openmc.Plot()
        xy_plot.filename = os.path.join(plots_dir, 'geometry_xy')
        xy_plot.width = (2*pool_radius, 2*pool_radius)
        xy_plot.pixels = (1000, 1000)
        xy_plot.color_by = 'cell'
        xy_plot.colors = colors
        xy_plot.basis = 'xy'
        xy_plot.origin = (0, 0, 0)
        
        # XZ plot (side view)
        xz_plot = openmc.Plot()
        xz_plot.filename = os.path.join(plots_dir, 'geometry_xz')
        xz_plot.width = (2*pool_radius, pool_height)
        xz_plot.pixels = (1000, 500)
        xz_plot.color_by = 'cell'
        xz_plot.colors = colors
        xz_plot.basis = 'xz'
        xz_plot.origin = (0, 0, 0)
        
        # YZ plot (side view)
        yz_plot = openmc.Plot()
        yz_plot.filename = os.path.join(plots_dir, 'geometry_yz')
        yz_plot.width = (2*pool_radius, pool_height)
        yz_plot.pixels = (1000, 500)
        yz_plot.color_by = 'cell'
        yz_plot.colors = colors
        yz_plot.basis = 'yz'
        yz_plot.origin = (0, 0, 0)
        
        # Material plot (shows materials instead of cells)
        mat_plot = openmc.Plot()
        mat_plot.filename = os.path.join(plots_dir, 'materials_xy')
        mat_plot.width = (2*pool_radius, 2*pool_radius)
        mat_plot.pixels = (1000, 1000)
        mat_plot.color_by = 'material'
        mat_plot.basis = 'xy'
        mat_plot.origin = (0, 0, 0)
        
        # Add voxel plot for 3D visualization
        voxel_plot = openmc.Plot()
        voxel_plot.filename = os.path.join(plots_dir, 'voxel')
        # Focus on the uranium and control rod region
        voxel_plot.width = (2*uranium_radius, 2*uranium_radius, uranium_height*1.2)
        voxel_plot.pixels = (100, 100, 100)
        voxel_plot.color_by = 'cell'
        voxel_plot.colors = colors
        voxel_plot.type = 'voxel'
        voxel_plot.origin = (0, 0, 0)
        
        # Add a cutaway voxel plot to see inside
        voxel_cut_plot = openmc.Plot()
        voxel_cut_plot.filename = os.path.join(plots_dir, 'voxel_cutaway')
        voxel_cut_plot.width = (2*uranium_radius, 2*uranium_radius, uranium_height*1.2)
        voxel_cut_plot.pixels = (100, 100, 100)
        voxel_cut_plot.color_by = 'cell'
        voxel_cut_plot.colors = colors
        voxel_cut_plot.type = 'voxel'
        voxel_cut_plot.origin = (0, 0, 0)
        # Add a quarter cutout to see inside
        voxel_cut_plot.view = [0., 0., 0., 
                              uranium_radius, uranium_radius, uranium_height/2]
        
        # Create and save the plots
        plots = openmc.Plots([xy_plot, xz_plot, yz_plot, mat_plot, voxel_plot, voxel_cut_plot])
        plots.export_to_xml()
        
        # Run the plotting
        openmc.plot_geometry()
        
        print(f"Geometry plots saved to {plots_dir}")
    
    return model

def plot_flux(statepoint_path):
    """Plot the neutron flux distribution"""
    sp = openmc.StatePoint(statepoint_path)
    
    # Get the eigenvalue
    k_eff = sp.keff
    
    # Get the flux tally
    flux_tally = sp.get_tally(name='flux')
    mesh_filter = flux_tally.find_filter(openmc.MeshFilter)
    mesh = mesh_filter.mesh
    
    flux = flux_tally.get_values(scores=['flux'])
    flux.shape = mesh.dimension
    
    # Create plots of the flux
    # XY plane (middle Z)
    z_mid = mesh.dimension[2] // 2
    plt.figure(figsize=(10, 8))
    plt.imshow(np.log(flux[:, :, z_mid]), cmap='plasma', 
              extent=[-100, 100, -100, 100], origin='lower')
    plt.colorbar(label='Log Neutron Flux')
    plt.title(f'Neutron Flux Distribution - XY Plane (k-effective = {k_eff.nominal_value:.5f})')
    plt.xlabel('X [cm]')
    plt.ylabel('Y [cm]')
    plt.savefig(os.path.join(os.path.dirname(statepoint_path), 'flux_map_xy.png'), dpi=300)
    
    # XZ plane (middle Y)
    y_mid = mesh.dimension[1] // 2
    plt.figure(figsize=(12, 6))
    plt.imshow(np.log(flux[:, y_mid, :]).T, cmap='plasma', 
              extent=[-100, 100, -75, 75], origin='lower')
    plt.colorbar(label='Log Neutron Flux')
    plt.title(f'Neutron Flux Distribution - XZ Plane (k-effective = {k_eff.nominal_value:.5f})')
    plt.xlabel('X [cm]')
    plt.ylabel('Z [cm]')
    plt.savefig(os.path.join(os.path.dirname(statepoint_path), 'flux_map_xz.png'), dpi=300)
    
    return k_eff.nominal_value

def visualize_voxel_h5(h5_file, output_file, title="3D Voxel Plot", downsample=2, data_key=None):
    """
    Visualize a voxel plot HDF5 file and save it as a PNG image.
    
    Parameters
    ----------
    h5_file : str
        Path to the voxel HDF5 file
    output_file : str
        Path to save the visualization
    title : str, optional
        Title for the plot
    downsample : int, optional
        Factor to downsample the voxel data for faster rendering
    data_key : str, optional
        The HDF5 dataset key to use for voxel data
    """
    try:
        # Read the voxel data from the HDF5 file
        with h5py.File(h5_file, 'r') as f:
            # Use the provided data_key if given
            if data_key and data_key in f:
                voxel_data = f[data_key][()]
            # Otherwise try common dataset names
            elif 'id' in f:
                voxel_data = f['id'][()]
            elif 'cell' in f:
                voxel_data = f['cell'][()]
            elif 'material' in f:
                voxel_data = f['material'][()]
            else:
                # Print all available keys for debugging
                print(f"Could not find appropriate data in {h5_file}")
                print(f"Available keys: {list(f.keys())}")
                return
                
            # Get dimensions if available
            if 'dimension' in f:
                dims = f['dimension'][()]
            elif 'dimensions' in f:
                dims = f['dimensions'][()]
            else:
                print("No dimension information found, assuming (100, 100, 100)")
                dims = (100, 100, 100)
            
        # Reshape if needed
        if len(voxel_data.shape) == 1:
            voxel_data = voxel_data.reshape(dims)
            
        # Downsample for performance
        voxel_downsampled = voxel_data[::downsample, ::downsample, ::downsample]
        
        # Create the 3D plot
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Get unique cell/material IDs
        cell_ids = np.unique(voxel_downsampled)
        
        # Create indices for the downsampled array
        x, y, z = np.indices(voxel_downsampled.shape)
        
        # Define colors for each cell/material
        colors = {
            1: 'yellow',    # Uranium
            2: 'red',       # Control rod
            3: 'skyblue'    # Water
        }
        
        # Plot each cell with a different color
        for cell_id in cell_ids:
            # Skip background (usually 0)
            if cell_id == 0:
                continue
                
            # Create a mask for this cell/material
            mask = (voxel_downsampled == cell_id)
            
            # Get the color for this cell/material
            color = colors.get(cell_id, 'gray')  # Default to gray if not in colors dict
            
            # Plot the voxels
            ax.voxels(mask, facecolor=color, edgecolor=None, 
                     alpha=0.7 if cell_id == 3 else 0.9)  # Water more transparent
        
        # Set labels and title
        ax.set_title(title)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # Improve viewing angle
        ax.view_init(elev=30, azim=45)
        
        # Remove axis ticks for cleaner look
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        
        # Save the figure
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Saved 3D visualization to {output_file}")
        
    except Exception as e:
        print(f"Error visualizing voxel data: {e}")

def plot_geometry_only(rod_positions=[0.0, 1.0]):
    """
    Generate geometry plots to visualize the model without running simulations.
    
    Parameters
    ----------
    rod_positions : list, optional
        List of rod positions to generate geometry plots for
    """
    output_dir = "simple-output"
    os.makedirs(output_dir, exist_ok=True)
    
    for pos in rod_positions:
        print(f"\nGenerating geometry plots for rod position {pos:.2f}...")
        position_dir = os.path.join(output_dir, f"geometry_pos_{pos:.2f}")
        os.makedirs(position_dir, exist_ok=True)
        
        # Change to the position directory
        os.chdir(position_dir)
        
        try:
            # Create model and plots
            create_model(rod_position=pos, plot_geometry=True, output_dir=position_dir)
        finally:
            # Return to the parent directory
            os.chdir('../..')
    
    print(f"\nGeometry plots complete. Check the {output_dir}/geometry_pos_* directories.")

def run_rod_positions():
    """Run simulations for different rod positions and plot results"""
    output_dir = "simple-output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Array of rod positions to simulate (0.0 = fully inserted, 1.0 = fully withdrawn)
    rod_positions = np.linspace(0.0, 1.0, 5)
    k_values = []
    
    for i, pos in enumerate(rod_positions):
        print(f"\nRunning simulation with rod position {pos:.2f}...")
        position_dir = os.path.join(output_dir, f"pos_{pos:.2f}")
        os.makedirs(position_dir, exist_ok=True)
        
        # Change to the position directory
        os.chdir(position_dir)
        
        try:
            # Create and run the model (with geometry plots for the first position)
            plot_geom = (i == 0)  # Only plot geometry for the first position
            model = create_model(rod_position=pos, plot_geometry=plot_geom, output_dir=position_dir)
            statepoint_path = f"statepoint.100.h5"
            
            # Only run if the statepoint doesn't already exist
            if not os.path.exists(statepoint_path):
                model.run()
            
            # Plot flux and get k-effective
            k_eff = plot_flux(statepoint_path)
            k_values.append(k_eff)
            print(f"Rod position {pos:.2f}: k-effective = {k_eff:.5f}")
            
        finally:
            # Return to the parent directory
            os.chdir('../..')
    
    # Plot k-effective vs rod position
    plt.figure(figsize=(10, 6))
    plt.plot(rod_positions, k_values, 'o-', linewidth=2, markersize=8)
    plt.xlabel('Control Rod Position (0.0 = fully inserted, 1.0 = fully withdrawn)')
    plt.ylabel('k-effective')
    plt.title('Effect of Control Rod Position on Criticality')
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, 'k_vs_position.png'), dpi=300)
    plt.close()
    
    # Save the data to a CSV file
    with open(os.path.join(output_dir, 'k_results.csv'), 'w') as f:
        f.write('Rod Position,k-effective\n')
        for pos, k in zip(rod_positions, k_values):
            f.write(f'{pos:.2f},{k:.5f}\n')
    
    print("\nSimulations complete. Results saved to the 'simple-output' directory.")

if __name__ == '__main__':
    import sys
    
    # Default is to run simulations
    if len(sys.argv) <= 1:
        run_rod_positions()
    else:
        # Check command line arguments
        if sys.argv[1] == 'plot-geometry':
            # Just plot the geometry (faster for visualization)
            print("Generating geometry plots only...")
            # Show fully inserted and fully withdrawn positions
            plot_geometry_only([0.0, 0.5, 1.0])
        elif sys.argv[1] == 'run':
            # Run the full simulation
            run_rod_positions()
        else:
            print(f"Unknown command: {sys.argv[1]}")
            print("Available commands:")
            print("  plot-geometry - Only create geometry plots (no simulation)")
            print("  run - Run full simulations")
            print("  (no args) - Same as 'run'")