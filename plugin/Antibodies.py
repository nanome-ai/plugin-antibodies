import functools
import itertools
import tempfile
import time
import nanome
from abnumber import Chain as AbChain
from abnumber.exceptions import ChainParseError
from Bio import SeqIO, SeqUtils
from concurrent.futures import ThreadPoolExecutor
from nanome.api import ui, structure
from nanome.util import async_callback, Color, Logs, enums

from enum import Enum
from .utils import get_neighboring_atoms

protein_letters_3to1 = SeqUtils.IUPACData.protein_letters_3to1
run_btn = enums.PluginListButtonType.run


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


class Antibodies(nanome.AsyncPluginInstance):

    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.menu = ui.Menu()
        self.menu.root.layout_orientation = enums.LayoutTypes.horizontal
        self.menu.title = "Antibody Regions"
        self.integration.structure_prep = self.integration_request

    @async_callback
    async def integration_request(self, request):
        complexes = request.get_args()
        comp = complexes[0]
        self.prep_antibody_complex(comp)
        request.send_response([comp])

    @async_callback
    async def on_run(self):
        # Get selected antibody complex
        Logs.debug("Loading Complex")
        self.set_plugin_list_button(run_btn, 'Loading Complex', False)
        comp_list = await self.request_complex_list()
        shallow_comp = next((cmp for cmp in comp_list if cmp.get_selected()), None)
        if not shallow_comp:
            self.send_notification(enums.NotificationTypes.error, "Please select an antibody")
            self._reset_run_btn()
            return
        comp = (await self.request_complexes([shallow_comp.index]))[0]
        if not self.validate_antibody(comp):
            self.send_notification(enums.NotificationTypes.error, "Selected complex is not an antibody")
            return
        self.set_plugin_list_button(run_btn, 'Finding CDR Loops...', False)
        self.prep_antibody_complex(comp)
        Logs.debug("Updating Structures.")
        self.update_structures_deep([comp])
        self.set_plugin_list_button(run_btn, 'Building menu...', True)
        self.build_menu(comp)
        self.menu.enabled = True
        self.update_menu(self.menu)
        self._reset_run_btn()

    def build_menu(self, comp: structure.Complex):
        self.menu.root.children = []
        for chain in comp.chains:
            seq_str = self.get_sequence_from_struct(chain)
            try:
                abchain = AbChain(seq_str, scheme='imgt')
            except ChainParseError as e:
                Logs.debug(f"Could not parse Chain {chain.name}")
                continue
            self.add_menu_chain_column(self.menu, chain, abchain)

    @classmethod
    def prep_antibody_complex(cls, comp):
        start_time = time.time()
        comp.set_all_selected(False)
        # Loop through chain and color cdr loops
        Logs.debug("Processing Chains.")
        chains_to_color = []
        for chain in comp.chains:
            Logs.debug(f"Chain {chain.name}")
            seq_str = cls.get_sequence_from_struct(chain)
            if not seq_str:
                Logs.debug(f"Unable to sequence chain {chain.name}")
                continue
            try:
                abchain = AbChain(seq_str, scheme='imgt')
            except ChainParseError as e:
                Logs.debug(f"Could not parse Chain {chain.name}")
                continue
            else:
                chains_to_color.append((chain, abchain))

        with ThreadPoolExecutor(max_workers=len(chains_to_color)) as executor:
            for ch, abch in chains_to_color:
                executor.submit(cls.format_chain, ch, abch)
        # Log data about run
        end_time = time.time()
        elapsed_time = round(end_time - start_time, 2)
        log_extra = {'elapsed_time': elapsed_time, 'residue_count': len(list(comp.residues))}
        Logs.message(f"Complex prepped in {elapsed_time} seconds", extra=log_extra)

    @classmethod
    def format_chain(cls, chain, abchain):
        """Color CDR loops and add labels."""
        # Make entire complex Grey.
        try:
            cdr1_residues = cls.get_cdr1_residues(chain)
            cdr2_residues = cls.get_cdr2_residues(chain)
            cdr3_residues = cls.get_cdr3_residues(chain)
            fr1_residues = cls.get_fr1_residues(chain)
            fr2_residues = cls.get_fr2_residues(chain)
            fr3_residues = cls.get_fr3_residues(chain)
            fr4_residues = cls.get_fr4_residues(chain)
        except ChainParseError as e:
            Logs.warning(f"Could find cdr loops for Chain {chain.name}")

        fr1_color = IMGTCDRColorScheme.FR.value
        fr2_color = IMGTCDRColorScheme.FR.value
        fr3_color = IMGTCDRColorScheme.FR.value
        fr4_color = IMGTCDRColorScheme.FR.value
        chain_type = abchain.chain_type
        if chain_type == 'H':
            cdr1_color = IMGTCDRColorScheme.HEAVY_CDR1.value
            cdr2_color = IMGTCDRColorScheme.HEAVY_CDR2.value
            cdr3_color = IMGTCDRColorScheme.HEAVY_CDR3.value
        else:
            cdr1_color = IMGTCDRColorScheme.LIGHT_CDR1.value
            cdr2_color = IMGTCDRColorScheme.LIGHT_CDR2.value
            cdr3_color = IMGTCDRColorScheme.LIGHT_CDR3.value

        cdr1_res_indices = [res.index for res in cdr1_residues]
        cdr2_res_indices = [res.index for res in cdr2_residues]
        cdr3_res_indices = [res.index for res in cdr3_residues]
        fr1_res_indices = [res.index for res in fr1_residues]
        fr2_res_indices = [res.index for res in fr2_residues]
        fr3_res_indices = [res.index for res in fr3_residues]
        fr4_res_indices = [res.index for res in fr4_residues]

        Logs.debug("Coloring cdr loops and framework")
        for residue in chain.residues:
            current_color = Color.Grey()
            use_wire_rendering = False
            res_index = residue.index
            if res_index in cdr1_res_indices:
                current_color = cdr1_color
                use_wire_rendering = True
            elif res_index in cdr2_res_indices:
                current_color = cdr2_color
                use_wire_rendering = True
            elif res_index in cdr3_res_indices:
                current_color = cdr3_color
                use_wire_rendering = True
            elif res_index in fr1_res_indices:
                current_color = fr1_color
            elif res_index in fr2_res_indices:
                current_color = fr2_color
            elif res_index in fr3_res_indices:
                current_color = fr3_color
            elif res_index in fr4_res_indices:
                current_color = fr4_color
            residue.ribbon_color = current_color
            for atom in residue.atoms:
                atom.atom_color = current_color
                if use_wire_rendering:
                    atom.set_visible(True)
                    atom.atom_mode = atom.AtomRenderingMode.Wire

        cls.label_residue_set(cdr1_residues, f'CDR{chain_type}1')
        cls.label_residue_set(cdr2_residues, f'CDR{chain_type}2')
        cls.label_residue_set(cdr3_residues, f'CDR{chain_type}3')
        cls.label_residue_set(fr1_residues, 'FR1')
        cls.label_residue_set(fr2_residues, 'FR2')
        cls.label_residue_set(fr3_residues, 'FR3')
        cls.label_residue_set(fr4_residues, 'FR4')

        comp = chain.complex
        cls.display_neighboring_atoms(comp, cdr1_residues)
        cls.display_neighboring_atoms(comp, cdr2_residues)
        cls.display_neighboring_atoms(comp, cdr3_residues)

    @classmethod
    def display_neighboring_atoms(cls, comp, residue_list, distance=3.0):
        # Make sure all atoms near cdr loop are in wire mode
        # This makes viewing interactions easier.
        Logs.debug("Making neighboring atoms wires")
        cdr_atoms = itertools.chain(*[res.atoms for res in residue_list])
        neighbor_atoms = get_neighboring_atoms(comp, cdr_atoms)
        for atom in neighbor_atoms:
            atom.set_visible(True)
            atom.atom_mode = atom.AtomRenderingMode.Wire
        return comp

    def on_chain_btn_pressed(self, residue_list, btn):
        self.zoom_on_structures(residue_list)

    def on_cdr_btn_pressed(self, residue_list, btn):
        # Select atoms
        for atom in itertools.chain(*[res.atoms for res in residue_list]):
            atom.selected = btn.selected
        self.update_structures_deep(residue_list)

    def add_menu_chain_column(self, menu: ui.Menu, chain: structure.Chain, abchain: AbChain):
        ln_chain = menu.root.create_child_node()
        ln_chain_btn = ln_chain.create_child_node()
        chain_btn = ln_chain_btn.add_new_button(f'{abchain.chain_type}')

        cdr1_residues = list(self.get_cdr1_residues(chain))
        cdr2_residues = list(self.get_cdr2_residues(chain))
        cdr3_residues = list(self.get_cdr3_residues(chain))
        cdr_residues = cdr1_residues + cdr2_residues + cdr3_residues
        chain_btn.register_pressed_callback(
            functools.partial(self.on_chain_btn_pressed, cdr_residues))

        ln_cdr1 = ln_chain.create_child_node()
        cdr1_btn = ln_cdr1.add_new_button(f"CDR{abchain.chain_type}1")
        cdr1_btn.toggle_on_press = True
        cdr1_btn.register_pressed_callback(
            functools.partial(self.on_cdr_btn_pressed, cdr1_residues))

        ln_cdr2 = ln_chain.create_child_node()
        cdr2_btn = ln_cdr2.add_new_button(f"CDR{abchain.chain_type}2")
        cdr2_btn.toggle_on_press = True
        cdr2_btn.register_pressed_callback(
            functools.partial(self.on_cdr_btn_pressed, cdr2_residues))

        ln_cdr3 = ln_chain.create_child_node()
        cdr3_btn = ln_cdr3.add_new_button(f"CDR{abchain.chain_type}3")
        cdr3_btn.toggle_on_press = True
        cdr3_btn.register_pressed_callback(
            functools.partial(self.on_cdr_btn_pressed, cdr3_residues))

    @staticmethod
    def label_residue_set(residue_list, label_text):
        """Add label to middle residue in residue list."""
        middle_residue = residue_list[(len(residue_list) // 2)]
        middle_residue.labeled = True
        middle_residue.label_text = label_text

    def validate_antibody(self, comp):
        # Make sure at least one chain can be parsed with ABChain
        for chain in comp.chains:
            seq_str = self.get_sequence_from_struct(chain)
            try:
                abchain = AbChain(seq_str, scheme='imgt')
                if abchain:
                    return True
            except ChainParseError as e:
                continue
            return False

    @classmethod
    def get_cdr1_residues(cls, struc, abchain=None):
        cdr_name = 'cdr1'
        residues = cls._get_region_residues(struc, cdr_name, abchain=abchain)
        return residues

    @classmethod
    def get_cdr2_residues(cls, struc, abchain=None):
        cdr_name = 'cdr2'
        residues = cls._get_region_residues(struc, cdr_name, abchain=abchain)
        return residues

    @classmethod
    def get_cdr3_residues(cls, struc, abchain=None):
        cdr_name = 'cdr3'
        residues = cls._get_region_residues(struc, cdr_name, abchain=abchain)
        return residues

    @classmethod
    def get_fr1_residues(cls, struc, abchain=None):
        fr_name = 'fr1'
        residues = cls._get_region_residues(struc, fr_name, abchain=abchain)
        return residues

    @classmethod
    def get_fr2_residues(cls, struc, abchain=None):
        fr_name = 'fr2'
        residues = cls._get_region_residues(struc, fr_name, abchain=abchain)
        return residues

    @classmethod
    def get_fr3_residues(cls, struc, abchain=None):
        fr_name = 'fr3'
        residues = cls._get_region_residues(struc, fr_name, abchain=abchain)
        return residues

    @classmethod
    def get_fr4_residues(cls, struc, abchain=None):
        fr_name = 'fr4'
        residues = cls._get_region_residues(struc, fr_name, abchain=abchain)
        return residues

    @classmethod
    def _get_region_residues(cls, chain, region_name: str, abchain=None):
        """Get nanome residues corresponding to provided cdr name.

        valid region names are 'cdr1', 'cdr2', 'cdr3', 'fr1', 'fr2', 'fr3', 'fr4'
        """
        valid_region_choices = ['cdr1', 'cdr2', 'cdr3', 'fr1', 'fr2', 'fr3', 'fr4']
        if region_name not in valid_region_choices:
            raise ValueError(f"Invalid region name: {region_name}. Valid choices are {', '.join(valid_region_choices)}")
        if not abchain:
            seq_str = cls.get_sequence_from_struct(chain)
            abchain = AbChain(seq_str, scheme='imgt')

        seq_attr_name = f'{region_name}_seq'
        cdr_seq = getattr(abchain, seq_attr_name)
        cdr_residues = []
        try:
            chain_seq = ''.join([protein_letters_3to1[res.name.title()] for res in chain.residues])
        except KeyError:
            return []
        if cdr_seq in chain_seq:
            residues = list(chain.residues)
            # Find the subset of the chain that comprises CDR Loop
            for idx in range(len(residues) - len(cdr_seq) + 1):
                residue_sublist = residues[idx:idx + len(cdr_seq)]
                residue_sublist_seq = ''.join([
                    protein_letters_3to1[res.name.title()] for res in residue_sublist])
                if residue_sublist_seq == cdr_seq:
                    cdr_residues = residue_sublist
                    break
        return filter(lambda res: res in cdr_residues, chain.residues)

    @staticmethod
    def get_sequence_from_pdb(pdb_filepath):
        with open(pdb_filepath) as handle:
            sequence = next(SeqIO.parse(handle, "pdb-atom"))
        seq = str(sequence.seq)
        return seq

    @staticmethod
    def get_sequence_from_struct(struct):
        try:
            chain_seq = ''.join([
                protein_letters_3to1[res.name.title()]
                for res in struct.residues
            ])
        except KeyError:
            chain_seq = ''
        return chain_seq

    def _reset_run_btn(self):
        self.set_plugin_list_button(run_btn, 'Run', True)


def main():
    name = 'Antibodies'
    description = "Select antibody in entry list, then run plugin to add IMGT color scheme and highlight CDR loops."
    plugin = nanome.Plugin(
        name, description, 'other', False,
        integrations=[enums.Integrations.structure_prep])
    plugin.set_plugin_class(Antibodies)
    plugin.run()


if __name__ == '__main__':
    main()
