from itertools import chain
from nanome.api import structure
from scipy.spatial import KDTree


__all__ = ['extract_residues', 'merge_complexes']


def get_neighboring_atoms(target_reference: structure.Complex, selected_atoms: list, site_size=6):
    """Use KDTree to find target atoms within site_size radius of selected atoms."""
    mol = next(
        mol for i, mol in enumerate(target_reference.molecules)
        if i == target_reference.current_frame)
    ligand_positions = [atom.position.unpack() for atom in selected_atoms]
    target_atoms = chain(*[ch.atoms for ch in mol.chains if not ch.name.startswith("H")])
    target_tree = KDTree([atom.position.unpack() for atom in target_atoms])
    target_point_indices = target_tree.query_ball_point(ligand_positions, site_size)
    near_point_set = set()
    for point_indices in target_point_indices:
        for point_index in point_indices:
            near_point_set.add(tuple(target_tree.data[point_index]))
    neighbor_atoms = []
    for targ_atom in mol.atoms:
        if targ_atom.position.unpack() in near_point_set:
            neighbor_atoms.append(targ_atom)
    return neighbor_atoms
