import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj
import plotly.graph_objs as go

# ----------------------------
# INPUTS
# ----------------------------
mypdbfile = "path_to_pdb"
mydcdfile = "path_to_dcd"

# ----------------------------
# LOAD & ALIGN TRAJECTORY
# ----------------------------
u = mda.Universe(mypdbfile, mydcdfile)
ca_atoms = u.select_atoms("name CA")
ref = mda.Universe(mypdbfile)
AlignTraj(u, ref, select="name CA", in_memory=True).run()

aligned_coords = np.array([ca_atoms.positions for _ in u.trajectory])
n_frames = aligned_coords.shape[0]

# ----------------------------
# PCA ON ALIGNED COORDINATES
# ----------------------------
xyz_flattened = aligned_coords.reshape(n_frames, -1)
pca = PCA(n_components=2)
pca_scores = pca.fit_transform(xyz_flattened)
pc1, pc2 = pca_scores[:, 0], pca_scores[:, 1]

# ----------------------------
# K-MEANS (k=5) ON PCA SCORES
# ----------------------------
k = 3
kmeans = KMeans(n_clusters=k, n_init=20, random_state=42)
labels = kmeans.fit_predict(pca_scores)

# ----------------------------
# DEFINE EXPANDED GRID
# ----------------------------
padding_factor = 0.2
pc1_min, pc1_max = pc1.min(), pc1.max()
pc2_min, pc2_max = pc2.min(), pc2.max()
pc1_pad = (pc1_max - pc1_min) * padding_factor
pc2_pad = (pc2_max - pc2_min) * padding_factor
pc1_min_e, pc1_max_e = pc1_min - pc1_pad, pc1_max + pc1_pad
pc2_min_e, pc2_max_e = pc2_min - pc2_pad, pc2_max + pc2_pad

nbins = 100
X_e = np.linspace(pc1_min_e, pc1_max_e, nbins)
Y_e = np.linspace(pc2_min_e, pc2_max_e, nbins)

# ----------------------------
# 2D HISTOGRAM → FREE ENERGY
# ----------------------------
hist, _, _ = np.histogram2d(
    pc1, pc2, bins=nbins,
    range=[[pc1_min_e, pc1_max_e], [pc2_min_e, pc2_max_e]],
    density=True
)
hist = hist + 1e-10
kb = 0.0083145  # kJ/mol·K
temp = 300.0    # K
F = -kb * temp * np.log(hist / np.max(hist))

# Smooth and clamp to 0–25 kJ/mol
F_smooth = gaussian_filter(F.T, sigma=0.5)
zmin, zmax = 0.0, 25.0
F_smooth = np.clip(F_smooth, zmin, zmax)

# ----------------------------
# CUSTOM COLORSCALE: top (1.0) = white
# ----------------------------
blackbody_to_white = [
    [0.00, 'rgb(0,0,0)'],
    [0.35, 'rgb(180,0,0)'],
    [0.60, 'rgb(255,130,0)'],
    [0.80, 'rgb(255,220,90)'],
    [1.00, 'rgb(255,255,255)']
]

# ----------------------------
# BASE LAYER: FREE ENERGY
# ----------------------------
fe_layer = go.Contour(
    x=X_e, y=Y_e, z=F_smooth,
    colorscale=blackbody_to_white,
    reversescale=False,
    zmin=zmin, zmax=zmax,
    contours=dict(coloring='heatmap', showlines=False, start=zmin, end=zmax, size=1.0),
    showscale=False,
    connectgaps=True,
    name='Free Energy',
    showlegend=False
)

# ----------------------------
# OVERLAY: BLACK DOTTED CONVEX HULLS ONLY (thicker)
# ----------------------------
hull_traces = []
for k_idx in range(k):
    mask = (labels == k_idx)
    if np.sum(mask) < 3:
        continue
    pts = np.column_stack((pc1[mask], pc2[mask]))
    try:
        hull = ConvexHull(pts)
        hull_pts = pts[hull.vertices]
        hx = np.append(hull_pts[:, 0], hull_pts[0, 0])  # close loop
        hy = np.append(hull_pts[:, 1], hull_pts[0, 1])

        hull_traces.append(
            go.Scatter(
                x=hx, y=hy,
                mode='lines',
                line=dict(color='black', width=4, dash='dot'),
                showlegend=False,
                hoverinfo='skip'
            )
        )
    except Exception:
        pass  # skip degenerate hulls

# ----------------------------
# LAYOUT (legend disabled)
# ----------------------------
layout = go.Layout(
    title=f'Free Energy Landscape (Top View) with K-Means Hulls (k={k})',
    xaxis=dict(
        title='PC1',
        showgrid=False, zeroline=False,
        showline=True, linecolor='black', linewidth=3,
        range=[pc1_min_e, pc1_max_e],
        showticklabels=False
    ),
    yaxis=dict(
        title='PC2',
        showgrid=False, zeroline=False,
        showline=True, linecolor='black', linewidth=3,
        scaleanchor='x', scaleratio=1,
        range=[pc2_min_e, pc2_max_e],
        showticklabels=False
    ),
    width=1100, height=900,
    margin=dict(l=50, r=50, b=50, t=50),
    plot_bgcolor='red', paper_bgcolor='red',
    showlegend=False
)

fig = go.Figure(data=[fe_layer] + hull_traces, layout=layout)
fig.show()
