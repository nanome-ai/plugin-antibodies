import itertools
import tempfile
import time
import nanome
from abnumber import Chain as AbChain
from abnumber.exceptions import ChainParseError
from Bio import SeqUtils
from concurrent.futures import ThreadPoolExecutor
from nanome.util import async_callback, Color, Logs, enums

from .utils import get_neighboring_atoms, IMGTCDRColorScheme
from .menu import RegionMenu, SettingsMenu

protein_letters_3to1 = SeqUtils.IUPACData.protein_letters_3to1
run_btn = enums.PluginListButtonType.run


class Antibodies(nanome.AsyncPluginInstance):
    current_menu_index = 1  # incremented to support multiple menus

    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.integration.structure_prep = self.integration_request
        self.menus = {}
        self.settings_menu = SettingsMenu(self)

    def on_stop(self):
        self.temp_dir.cleanup()

    def on_advanced_settings(self):
        self.settings_menu.render()

    @async_callback
    async def on_run(self):
        start_time = time.time()
        # Get selected antibody complex
        Logs.debug("Loading Complexes")
        self.set_plugin_list_button(run_btn, 'Loading...', False)
        comp_list = await self.request_complex_list()
        shallow_comps = (cmp for cmp in comp_list if cmp.get_selected())
        if not shallow_comps:
            self.send_notification(enums.NotificationTypes.error, "Please select an antibody")
            self._reset_run_btn()
            return
        comps = (await self.request_complexes([cmp.index for cmp in shallow_comps]))
        for i, comp in enumerate(comps):
            counter_str = f'({i + 1}/{len(comps)}) ' if len(comps) > 1 else ''
            self.set_plugin_list_button(run_btn, f'{counter_str}Coloring...', False)
            if not self.validate_antibody(comp):
                self.send_notification(enums.NotificationTypes.warning, f"{comp.full_name} is not an antibody")
                continue

            self.prep_antibody_complex(comp)
            # self.set_plugin_list_button(run_btn, f'{counter_str}Building menu...', False)
            Logs.debug("Building Menu...")
            new_menu = RegionMenu(self)
            new_menu.build_menu(comp)
            self.menus[new_menu.index] = new_menu
            self.current_menu_index += 1

        Logs.debug("Updating Structures.")
        self.set_plugin_list_button(run_btn, 'Updating..', False)
        self.update_structures_deep(comps)
        for menu in self.menus.values():
            menu.enable()
        self._reset_run_btn()
        end_time = time.time()
        elapsed_time = round(end_time - start_time, 2)
        log_extra = {
            'elapsed_time': elapsed_time,
            'complex_count': len(comps),
            'residue_count': sum([len(list(comp.residues)) for comp in comps])}
        Logs.message(f"Complexes updated in {elapsed_time} seconds", extra=log_extra)
        return comps

    @async_callback
    async def integration_request(self, request):
        # TODO Reactivate when we better handle swappable structure prep plugins.
        complexes = request.get_args()
        for comp in complexes:
            self.prep_antibody_complex(comp)
        request.send_response(complexes)
        return complexes

    @classmethod
    def prep_antibody_complex(cls, comp):
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

    @classmethod
    def format_chain(cls, chain, abchain):
        """Color CDR loops and add labels."""
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

    @staticmethod
    def label_residue_set(residue_list, label_text):
        """Add label to middle residue in residue list."""
        middle_residue = residue_list[(len(residue_list) // 2)]
        middle_residue.labeled = True
        middle_residue.label_text = label_text

    @classmethod
    def validate_antibody(cls, comp):
        # Make sure at least one chain can be parsed with ABChain
        for chain in comp.chains:
            seq_str = cls.get_sequence_from_struct(chain)
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
    name = 'Antibody Representation'
    description = "Select antibody in entry list, then run plugin to add IMGT color scheme and highlight CDR loops."
    category = 'other'
    has_advanced_options = True
    plugin = nanome.Plugin(
        name, description, category, has_advanced_options,
        # integrations=[enums.Integrations.structure_prep]
    )
    plugin.set_plugin_class(Antibodies)
    plugin.run()


if __name__ == '__main__':
    main()
