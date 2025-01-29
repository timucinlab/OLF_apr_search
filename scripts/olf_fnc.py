import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis.tests.datafiles import PSF, DCD

# Load the traj
u = mda.Universe("path_to_psf",
                 "path_to_dcd")

sel_0 = "protein and resid 268 to 317"
sel_1 = "protein and resid 328 to 359"
sel_2 = "protein and resid 380 to 422"
sel_3 = "protein and resid 428 to 467"
sel_4 = "protein and resid 474 to 501 255 to 259"

# Select atoms for each blade
sel0 = u.select_atoms(sel_0)
sel1 = u.select_atoms(sel_1)
sel2 = u.select_atoms(sel_2)
sel3 = u.select_atoms(sel_3)
sel4 = u.select_atoms(sel_4)

# Set up the analysis with soft-cut
def soft_cut_contacts(r, r0=6.0, beta=5.0, lambda_constant=1.8):
    return contacts.soft_cut_q(r, r0=r0, beta=beta, lambda_constant=lambda_constant)

ca1 = contacts.Contacts(
    u, select=(sel_0, sel_1),
    refgroup=(sel0, sel1),
    method=soft_cut_contacts  # Use soft-cut function
)
ca2 = contacts.Contacts(
    u, select=(sel_1, sel_2),
    refgroup=(sel1, sel2),
    method=soft_cut_contacts  # Use soft-cut function
)
ca3 = contacts.Contacts(
    u, select=(sel_2, sel_3),
    refgroup=(sel2, sel3),
    method=soft_cut_contacts  # Use soft-cut function
)
ca4 = contacts.Contacts(
    u, select=(sel_3, sel_4),
    refgroup=(sel3, sel4),
    method=soft_cut_contacts  # Use soft-cut function
)
ca5 = contacts.Contacts(
    u, select=(sel_4, sel_0),
    refgroup=(sel4, sel0),
    method=soft_cut_contacts  # Use soft-cut function
)

# Run the contact analysis
ca1.run()
ca2.run()
ca3.run()
ca4.run()
ca5.run()

# Calculate and print the average number of contacts for each pair
average_contacts = np.mean(ca1.timeseries[:, 1])
print('Average contacts between sel0 and sel1: {:.2f}'.format(average_contacts))
average_contacts = np.mean(ca2.timeseries[:, 1])
print('Average contacts between sel1 and sel2: {:.2f}'.format(average_contacts))
average_contacts = np.mean(ca3.timeseries[:, 1])
print('Average contacts between sel2 and sel3: {:.2f}'.format(average_contacts))
average_contacts = np.mean(ca4.timeseries[:, 1])
print('Average contacts between sel3 and sel4: {:.2f}'.format(average_contacts))
average_contacts = np.mean(ca5.timeseries[:, 1])
print('Average contacts between sel4 and sel0: {:.2f}'.format(average_contacts))

# Save the contact number matrices for each frame to the specified path
output_path = "path_to_output"

contact_matrices = []
for i, ca in enumerate([ca1, ca2, ca3, ca4, ca5], start=1):
    contact_matrix = np.zeros((len(ca.timeseries), 1))
    contact_matrix[:, 0] = ca.timeseries[:, 1]  # Store the contact values in the matrix
    contact_matrices.append(contact_matrix)
    # Save each matrix to a file in the specified path
    matrix_filename = f"{output_path}contact_matrix_{i}.txt"
    np.savetxt(matrix_filename, contact_matrix)
    print(f"Contact matrix for pair {i} saved as '{matrix_filename}'")
