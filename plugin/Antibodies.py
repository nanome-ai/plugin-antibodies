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
        self.menu = self.build_menu(comp)
        self.menu.enabled = True
        self.update_menu(self.menu)
        self._reset_run_btn()

    @classmethod
    def prep_antibody_complex(cls, comp):
        start_time = time.time()
        # comp.set_all_selected(False)
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

        fr_color = IMGTCDRColorScheme.FR.value
        chain_type = abchain.chain_type
        if chain_type == 'H':
            cdr1_color = IMGTCDRColorScheme.HEAVY_CDR1.value
            cdr2_color = IMGTCDRColorScheme.HEAVY_CDR2.value
            cdr3_color = IMGTCDRColorScheme.HEAVY_CDR3.value
        else:
            cdr1_color = IMGTCDRColorScheme.LIGHT_CDR1.value
            cdr2_color = IMGTCDRColorScheme.LIGHT_CDR2.value
            cdr3_color = IMGTCDRColorScheme.LIGHT_CDR3.value

        cdr1_res_serials = [res.serial for res in cdr1_residues]
        cdr2_res_serials = [res.serial for res in cdr2_residues]
        cdr3_res_serials = [res.serial for res in cdr3_residues]
        fr_res_serials = [res.serial for res in itertools.chain(fr1_residues, fr2_residues, fr3_residues, fr4_residues)]

        Logs.debug("Coloring cdr loops and framework")
        for residue in chain.residues:
            current_color = Color.Grey()
            use_wire_rendering = False
            res_serial = residue.serial
            if res_serial in cdr1_res_serials:
                current_color = cdr1_color
                use_wire_rendering = True
            elif res_serial in cdr2_res_serials:
                current_color = cdr2_color
                use_wire_rendering = True
            elif res_serial in cdr3_res_serials:
                current_color = cdr3_color
                use_wire_rendering = True
            elif res_serial in fr_res_serials:
                current_color = fr_color

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
        Logs.debug("Making neighboring atoms wires")
        cls.display_neighboring_atoms(comp, cdr1_residues)
        cls.display_neighboring_atoms(comp, cdr2_residues)
        cls.display_neighboring_atoms(comp, cdr3_residues)

    @classmethod
    def display_neighboring_atoms(cls, comp, residue_list):
        # Make sure all atoms near cdr loop are in wire mode
        # This makes viewing interactions easier.
        cdr_atoms = itertools.chain(*[res.atoms for res in residue_list])
        neighbor_atoms = get_neighboring_atoms(comp, cdr_atoms)
        for atom in neighbor_atoms:
            is_framework = atom.atom_color.rgb == IMGTCDRColorScheme.FR.value.rgb
            if not is_framework:
                atom.set_visible(True)
                atom.atom_mode = atom.AtomRenderingMode.Wire
        return comp

    def build_menu(self, comp: structure.Complex):
        menu = ui.Menu()
        menu.root.layout_orientation = enums.LayoutTypes.horizontal
        menu.title = f"Antibody Regions {comp.full_name}"
        comp.register_selection_changed_callback(self.on_selection_changed)
        for chain in comp.chains:
            seq_str = self.get_sequence_from_struct(chain)
            try:
                abchain = AbChain(seq_str, scheme='imgt')
            except ChainParseError as e:
                Logs.debug(f"Could not parse Chain {chain.name}")
                continue
            self.add_menu_chain_column(menu, chain, abchain)
        self.update_cdr_btns(menu, comp)
        return menu

    def add_menu_chain_column(self, menu: ui.Menu, chain: structure.Chain, abchain: AbChain):
        ln_chain = menu.root.create_child_node()
        ln_chain.chain_index = chain.index
        ln_chain_btn = ln_chain.create_child_node()
        chain_btn = ln_chain_btn.add_new_button(f'{abchain.chain_type}')

        cdr1_residues = self.get_cdr1_residues(chain)
        cdr2_residues = self.get_cdr2_residues(chain)
        cdr3_residues = self.get_cdr3_residues(chain)
        cdr_residues = cdr1_residues + cdr2_residues + cdr3_residues
        chain_btn.register_pressed_callback(
            functools.partial(self.on_chain_btn_pressed, cdr_residues))

        ln_cdr1 = ln_chain.create_child_node()
        cdr1_btn = ln_cdr1.add_new_button(f"CDR{abchain.chain_type}1")
        cdr1_btn.toggle_on_press = True
        cdr1_btn.cdr_residues = cdr1_residues
        cdr1_btn.register_pressed_callback(
            functools.partial(self.on_cdr_btn_pressed, cdr1_residues))

        ln_cdr2 = ln_chain.create_child_node()
        cdr2_btn = ln_cdr2.add_new_button(f"CDR{abchain.chain_type}2")
        cdr2_btn.toggle_on_press = True
        cdr2_btn.cdr_residues = cdr2_residues
        cdr2_btn.register_pressed_callback(
            functools.partial(self.on_cdr_btn_pressed, cdr2_residues))

        ln_cdr3 = ln_chain.create_child_node()
        cdr3_btn = ln_cdr3.add_new_button(f"CDR{abchain.chain_type}3")
        cdr3_btn.toggle_on_press = True
        cdr3_btn.cdr_residues = cdr3_residues
        cdr3_btn.register_pressed_callback(
            functools.partial(self.on_cdr_btn_pressed, cdr3_residues))

    def on_selection_changed(self, comp):
        """Update the selection in the plugin."""
        Logs.debug(f"Selection changes for {comp.full_name}")
        self.update_cdr_btns(self.menu, comp)
        self.update_menu(self.menu, comp)
        Logs.debug("Finished updating selections")

    def update_cdr_btns(self, menu, comp):
        btns_to_update = []
        for ln in menu.root.get_children():
            # Get most up to date chain selections
            chain_index = ln.chain_index
            comp_chain = next(ch for ch in comp.chains if ch.index == chain_index)
            cdr1_btn = ln.get_children()[1].get_content()
            cdr2_btn = ln.get_children()[2].get_content()
            cdr3_btn = ln.get_children()[3].get_content()
            
            abchain = AbChain(self.get_sequence_from_struct(comp_chain), scheme='imgt')
            cdr1_residues = self.get_cdr1_residues(comp_chain, abchain=abchain)
            cdr2_residues = self.get_cdr2_residues(comp_chain, abchain=abchain)
            cdr3_residues = self.get_cdr3_residues(comp_chain, abchain=abchain)
            cdr1_btn.selected = any([
                atom.selected for atom in itertools.chain(*[res.atoms for res in cdr1_residues])])
            cdr2_btn.selected = any([
                atom.selected for atom in itertools.chain(*[res.atoms for res in cdr2_residues])])
            cdr3_btn.selected = any([
                atom.selected for atom in itertools.chain(*[res.atoms for res in cdr3_residues])])
            btns_to_update.extend([cdr1_btn, cdr2_btn, cdr3_btn])
        self.update_content(*btns_to_update)

    def on_cdr_btn_pressed(self, residue_list, btn):
        # Select atoms
        for atom in itertools.chain(*[res.atoms for res in residue_list]):
            atom.selected = btn.selected
        self.update_structures_deep(residue_list)

    def on_chain_btn_pressed(self, residue_list, btn):
        self.zoom_on_structures(residue_list)

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
        return list(filter(lambda res: res in cdr_residues, chain.residues))

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
        # integrations=[enums.Integrations.structure_prep]
    )
    plugin.set_plugin_class(Antibodies)
    plugin.run()


if __name__ == '__main__':
    main()
