import numpy as np
from ase.io import read, write
from ase.build import make_supercell, stack
from ase import Atoms




def stack_poscars(in_atoms, multipliers=None, inversions=None, centering_atoms=None, write2file=None):
    """

    @param in_atoms:
    @param multipliers:
    @param inversions:
    @param centering_atoms: three numbers of atoms for choosing the origin, e.g. [0,1,3] means that x-coordinate of 0th
                            atom, y-coordinate of 1st atom and z-coordinate of 3rd atom are all zero
    @param write2file:
    @return:
    """
    # infile = '/home/vsbat/SYNC/00__WORK/20200209_CoPc_project/type5_POSCAR'
    # write2file = '/home/vsbat/mnt/arkuda_VBATURIN/PROJECT_CoPc/Polymers/POLYMER_t3_1bond/step0_Relax_lowprec/POSCAR'
    if multipliers is None:
        multipliers = [2, 1, 1]
    if inversions is None:
        inversions = [False, False, False]

    if isinstance(in_atoms, str):
        in_atoms = read(in_atoms, format='vasp')

    if centering_atoms is not None:
        assert len(centering_atoms) == 3
        center = np.zeros((1, 3))
        for dim in range(3):
            center[0, dim] = in_atoms.get_scaled_positions()[centering_atoms[dim]][dim]
        in_atoms.set_scaled_positions(in_atoms.get_scaled_positions() - center + np.array([[0.5, 0.5, 0.5]]))

    assert all([m_i >= 1 for m_i in multipliers])

    in_atoms.wrap()
    out_atoms = in_atoms
    for i_mult, mult_i in enumerate(multipliers):
        if mult_i > 1:
            second_cell = out_atoms.copy()
            second_cell.set_scaled_positions(
                second_cell.get_scaled_positions() * np.array([[-1 if x else 1 for x in inversions]]))
            out_atoms = stack(out_atoms, second_cell, axis=i_mult)
    if write2file:
        write(write2file, out_atoms, format='vasp', vasp5=True, direct=True)

    return out_atoms


def vacuum_dress(ats, vacuum_sizes):
    """
    @param ats: ase.Atoms object
    @param vacuum_sizes: ordered container of vacuum distances along cartesian coordinates
    @return:  new cell parameters
    """
    # circum_box = np.max(ats.positions, 0) - np.min(ats.positions, 0)
    ats.set_cell(ats.cell + vacuum_sizes * np.eye(3))
    ats.pbc = [True, True, True]
    ats.center()
    return ats.cell


def exclude_atoms(nparray1d, del_ind):
    for i in sorted(del_ind, reverse=True):
        del nparray1d[i]


def mol2poscars(init_mol, rm_atoms=None, basis_atoms=None, vacuums=None, res_filename=None):
    """
    Turns a molecule into a periodic ASE.Atoms object.
    @param init_mol: initial ase.Atoms object of a molecule or a name of a file containing it
    @param rm_atoms: list of atoms to delete in order to make a polymer
    @param basis_atoms: two pairs of zero-based numbers of atoms, first pair to create new X axis,
                        second - for a direction of Y axis (orthogonal transform, so shape is conserved)
    @param vacuums: vacuum sizes to dress a molecule. One of components should be close to interatomic distance if a
                    Polymer is to be considered
    @param res_filename: a filename to write a result
    @return: ase.Atoms object
    """
    if isinstance(init_mol, str):
        out_mol = read(init_mol)
    else:
        out_mol = init_mol.copy()
    coords = out_mol.positions

    if basis_atoms is not None:
        assert np.array(basis_atoms).shape == (2, 2)
        e_hor = coords[basis_atoms[0][1]] - coords[basis_atoms[0][0]]
        a = e_hor / np.linalg.norm(e_hor)
        e_vert = coords[basis_atoms[1][1]] - coords[basis_atoms[1][0]]
        b = e_vert - a * (np.dot(e_vert, a))
        b /= np.linalg.norm(b)
        c = np.cross(a, b)
        rot_mat = np.vstack((a, b, c))
        new_coords = np.dot(coords, rot_mat.T)
        out_mol.set_positions(new_coords)

    exclude_atoms(out_mol, rm_atoms)

    if vacuums is not None:
        vacuum_dress(out_mol, vacuums)
    out_mol.center()

    if res_filename is not None:
        write(res_filename, out_mol, format='vasp', vasp5=True)

    return out_mol

def sphere_cut(in_poscar, center_frac, radius, out_poscar=None):
    primitive_atoms = read(in_poscar)
    multipliers = np.array([np.ceil(2 * radius * np.linalg.norm(v)) + 1 for v in primitive_atoms.get_reciprocal_cell()])
    circum_cell = make_supercell(primitive_atoms, np.diag(multipliers))

    sphere_origin = (np.floor(multipliers / 2) + center_frac) @ primitive_atoms.cell
    sphere_filter = np.linalg.norm(circum_cell.positions - sphere_origin, axis=1) <= radius
    sphere_atoms = circum_cell[sphere_filter]

    if out_poscar:
        write(out_poscar, sphere_atoms, vasp5=True, direct=True)

    return sphere_atoms
