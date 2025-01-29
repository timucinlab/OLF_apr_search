import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align
import matplotlib.pyplot as plt
import numpy as np

u = mda.Universe("path_to_psf",
                 "path_to_dcd")

# Align the trajectory on protein
aligner = align.AlignTraj(u, u, select='protein', in_memory=True).run()

# Calculate the distance matrix for protein 
matrix = diffusionmap.DistanceMatrix(u, select='protein').run()
matrix = matrix.dist_matrix

fig = plt.imshow(matrix, cmap='OrRd', vmin=0, vmax=10)
plt.xlabel('', weight='bold', fontsize=12)
plt.ylabel('system_name', weight='bold', fontsize=12)

# Add color bar
cbar = plt.colorbar(label='RMSD (Å)')
cbar.ax.tick_params(labelsize=16) 
cbar.set_label('RMSD (Å)', fontsize=16)  

ax = plt.gca()
ax.set_xticks([])
ax.set_yticks([])

plt.title("Protein (ref:protein)", fontsize=18)

# Save the plot as a PDF 
output_path = 'path_to_output_prot-prot.pdf'
plt.savefig(output_path, format='pdf',bbox_inches="tight")
plt.close()

# Align the trajectory on beta-strands
aligner = align.AlignTraj(u, u, select='protein and resid 248 to 251 255 to 259 268 to 271 285 to 289 297 to 301 314 to 317 328 to 330 333 to 337 343 to 348 353 to 359 380 to 384 387 to 392 400 to 406 413 to 422 428 to 432 435 to 439 447 to 454 460 to 467 474 to 480 485 to 490 493 to 501', in_memory=True).run()
print(u)

# Calculate the distance matrix for beta-strands
matrix = diffusionmap.DistanceMatrix(u, select='protein and resid 248 to 251 255 to 259 268 to 271 285 to 289 297 to 301 314 to 317 328 to 330 333 to 337 343 to 348 353 to 359 380 to 384 387 to 392 400 to 406 413 to 422 428 to 432 435 to 439 447 to 454 460 to 467 474 to 480 485 to 490 493 to 501').run()
matrix = matrix.dist_matrix

fig = plt.imshow(matrix, cmap='OrRd', vmin=0, vmax=10)
plt.xlabel('', weight='bold', fontsize=12)
plt.ylabel('system_name', weight='bold', fontsize=12)

# Add color bar
cbar = plt.colorbar(label='RMSD (Å)')
cbar.ax.tick_params(labelsize=16)  
cbar.set_label('RMSD (Å)', fontsize=16) 

ax = plt.gca()
ax.set_xticks([])
ax.set_yticks([])

plt.title("β-strands (ref:β-strands)", fontsize=18)

# Save the plot as a PDF 
output_path = 'path_to_output_beta-beta.pdf''
plt.savefig(output_path, format='pdf',bbox_inches="tight") ##kenarlara oturması için bbox inches ekledim
plt.close()

# Align the trajectory on beta-strands
aligner = align.AlignTraj(u, u, select='protein and resid 248 to 251 255 to 259 268 to 271 285 to 289 297 to 301 314 to 317 328 to 330 333 to 337 343 to 348 353 to 359 380 to 384 387 to 392 400 to 406 413 to 422 428 to 432 435 to 439 447 to 454 460 to 467 474 to 480 485 to 490 493 to 501', in_memory=True).run()
print(u)

# Calculate the distance matrix for protein
matrix = diffusionmap.DistanceMatrix(u, select='protein').run()
matrix = matrix.dist_matrix

fig = plt.imshow(matrix, cmap='OrRd', vmin=0, vmax=10)
plt.xlabel('', weight='bold', fontsize=12)
plt.ylabel('system_name', weight='bold', fontsize=12)

# Add color bar 
cbar = plt.colorbar(label='RMSD (Å)')
cbar.ax.tick_params(labelsize=16) 
cbar.set_label('RMSD (Å)', fontsize=16)  

ax = plt.gca()
ax.set_xticks([])
ax.set_yticks([])

plt.title("Protein (ref:β-strands)", fontsize=18)

# Save the plot as a PDF 
output_path = 'path_to_output_beta-prot.pdf'
plt.savefig(output_path, format='pdf', bbox_inches="tight")  # Correct alignment
plt.close()
