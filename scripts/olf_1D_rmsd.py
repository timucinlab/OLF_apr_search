import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import matplotlib.pyplot as plt

# File paths for system 1
psf_file_1 = "path_to_psf1"
dcd_file_1 = "path_to_dcd1"

# File paths for system 2
psf_file_2 = "path_to_psf2"
dcd_file_2 = "path_to_dcd2"

# Function to calculate and return RMSD
def get_rmsd(psf_file, dcd_file, color):
    # Load the trajectory
    u = mda.Universe(psf_file, dcd_file)

    # Select the protein atoms for RMSD calculation
    protein = u.select_atoms("protein")

    # RMSD calculation (set the first frame as the reference)
    rmsd = RMSD(protein, u, select="protein", ref_frame=0)
    rmsd.run()

    # Extract RMSD values
    frame_numbers = rmsd.rmsd[:, 0]  # Frame numbers
    rmsd_values = rmsd.rmsd[:, 2]  # RMSD values in Å

    # Convert frame numbers to time in ns (assuming 0.1 ns per frame)
    time_ns = frame_numbers * 0.1  # Convert to ns

    return time_ns, rmsd_values, color

# Get RMSD for system 1 and system 2
time_ns_1, rmsd_values_1, color_1 = get_rmsd(psf_file_1, dcd_file_1, '#822EFF')  
time_ns_2, rmsd_values_2, color_2 = get_rmsd(psf_file_2, dcd_file_2, '#F433FF')  

# Plot RMSD for both systems
plt.figure(figsize=(18, 4))

# Add vertical lines to given frames
plt.axvline(x=4400 * 0.1, color='#BCC6CC', linestyle='--', linewidth=2, label='Frame 9431')  
plt.axvline(x=12500 * 0.1, color='#BCC6CC', linestyle='--', linewidth=2, label='Frame 17903') 
plt.axvline(x=9900 * 0.1, color='black', linestyle='--', linewidth=2, label='Frame 257')  
plt.axvline(x=0 * 0.1, color='#52595D', linestyle='--', linewidth=2, label='Frame 10132')  

plt.plot(time_ns_1, rmsd_values_1, label='System 1', color=color_1, linewidth=2)
plt.plot(time_ns_2, rmsd_values_2, label='System 2', color=color_2, linewidth=2)

# Add labels and title
plt.xlabel('time (ns)', fontsize=30)
plt.ylabel('RMSD (Å)', fontsize=30)
plt.title('', fontsize=32)
plt.ylim(0,5)
plt.tick_params(axis='both', which='major', labelsize=28)

# Save the plot
output_path = "path_to_output.png"
plt.tight_layout()
plt.savefig(output_path, format='png', dpi=600)

# Show the plot
plt.show()
