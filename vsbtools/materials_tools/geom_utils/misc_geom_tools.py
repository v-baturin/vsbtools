import numpy as np
from ase.io import read, write
from pathlib import Path

def center_and_align_xyz(xyz_file):
    # Read the xyz file
    xyz_file = Path(xyz_file)
    atoms = read(xyz_file)

    # Center the structure
    atoms.center()

    # Compute the moments of inertia tensor
    _, eigenvectors = atoms.get_moments_of_inertia(vectors=True)
    principal_axes = eigenvectors.T

    # Rotate the structure so that the principal axes align with x, y, and z axes
    rotation_matrix = np.linalg.inv(principal_axes)
    atoms.set_positions(atoms.positions.dot(eigenvectors))

    # Write the modified structure to a new xyz file
    new_xyz_file = xyz_file.parent / f"aligned_{xyz_file.name}"
    write(new_xyz_file, atoms, format='xyz')

    print(f"Aligned structure saved to {new_xyz_file}")

# Example usage
if __name__ == "__main__":
    center_and_align_xyz("/home/vsbat/SYNC/00__WORK/20240308_Charpentier/AlN6H18.xyz")