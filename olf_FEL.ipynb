import numpy as np
from scipy.ndimage import gaussian_filter
from sklearn.decomposition import PCA
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj
import plotly.graph_objs as go

# Load PDB and traj
mypdbfile = "path_to_pdb_file"
mydcdfile = "path_to_dcd_file"
u = mda.Universe(mypdbfile, mydcdfile)

# Select atoms and align traj
ca_atoms = u.select_atoms("name CA")
ref = mda.Universe(mypdbfile)
aligner = AlignTraj(u, ref, select="name CA", in_memory=True)
aligner.run()

aligned_coords = np.array([ca_atoms.positions for ts in u.trajectory])
n_frames = aligned_coords.shape[0]

# Perform PCA
xyz_flattened = aligned_coords.reshape(n_frames, -1)
pca = PCA(n_components=2)
pca_scores = pca.fit_transform(xyz_flattened)

pc1, pc2 = pca_scores[:, 0], pca_scores[:, 1]
padding_factor = 0.2  # Increased padding factor (30% of the range)
pc1_min, pc1_max = pc1.min(), pc1.max()
pc2_min, pc2_max = pc2.min(), pc2.max()

pc1_padding = (pc1_max - pc1_min) * padding_factor
pc2_padding = (pc2_max - pc2_min) * padding_factor
pc1_min_expanded, pc1_max_expanded = pc1_min - pc1_padding, pc1_max + pc1_padding
pc2_min_expanded, pc2_max_expanded = pc2_min - pc2_padding, pc2_max + pc2_padding

# Expand 2D grid
nbins = 100
X_expanded = np.linspace(pc1_min_expanded, pc1_max_expanded, nbins)
Y_expanded = np.linspace(pc2_min_expanded, pc2_max_expanded, nbins)
X_mesh, Y_mesh = np.meshgrid(X_expanded, Y_expanded)

# Histogram for expanded grid
hist, x_edges, y_edges = np.histogram2d(
    pc1, pc2, bins=nbins, 
    range=[[pc1_min_expanded, pc1_max_expanded], [pc2_min_expanded, pc2_max_expanded]], 
    density=True
)

hist = hist + 1e-10
kb = 0.0083145  # Boltzmann constant in kJ/(mol⋅K)
temp = 300  # Temperature in Kelvin
F_expanded = -kb * temp * np.log(hist / np.max(hist))
free_energy_smoothed = gaussian_filter(F_expanded.T, sigma=0.5)

# Create 3D surface plot
trace = go.Surface(
    x=X_mesh,
    y=Y_mesh,
    z=free_energy_smoothed,
    colorscale='rainbow',
    opacity=1,
    contours=dict(
        z=dict(
            show=True,
            usecolormap=True,
            highlightcolor='black',
            project=dict(x=True, y=True, z=True)
        )
    )
)

layout = go.Layout(
    title='Enlarged Free Energy Landscape',
    scene=dict(
        xaxis=dict(
            showgrid=False,  # Hide grid lines
            title='PC1',
            zeroline=False,  
            showline=True,  # Show axis lines
            linecolor='black',  
            linewidth=3, 
            range=[pc1_min_expanded, pc1_max_expanded],
            showticklabels=False, 
            titlefont=dict(size=24, color='black'), 
            backgroundcolor="white" 
        ),
        yaxis=dict(
            showgrid=False, 
            title='PC2',
            zeroline=False, 
            showline=True, 
            linecolor='black', 
            linewidth=3,  
            range=[pc2_min_expanded, pc2_max_expanded],
            showticklabels=False,  
            titlefont=dict(size=24, color='black'), 
            backgroundcolor="white" 
        ),
        zaxis=dict(
            showgrid=False, 
            title='Free Energy',
            zeroline=False,  
            showline=True,  
            linecolor='black',  
            linewidth=3,  # 
            showticklabels=False,  
            titlefont=dict(size=24, color='black'), 
            backgroundcolor="white"  
        )
    ),
    width=1100,
    height=900,
    margin=dict(l=50, r=50, b=50, t=50)
)

fig = go.Figure(data=[trace], layout=layout)
fig.show()
