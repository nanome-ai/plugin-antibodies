from enum import Enum
from itertools import chain
from nanome.api import structure
from scipy.spatial import KDTree

from nanome.util import Color


__all__ = ['IMGTCDRColorScheme', 'extract_residues', 'merge_complexes']


class IMGTCDRColorScheme(Enum):
    """
    Source: https://www.imgt.org/IMGTScientificChart/RepresentationRules/colormenu.php#h1_26
    """
    HEAVY_CDR1 = Color(200, 0, 0)
    HEAVY_CDR2 = Color(255, 169, 0)
    HEAVY_CDR3 = Color(156, 65, 215)
    LIGHT_CDR1 = Color(96, 96, 228)
    LIGHT_CDR2 = Color(70, 213, 0)
    LIGHT_CDR3 = Color(63, 157, 63)
    # Added by us
    FR = Color.White()


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
