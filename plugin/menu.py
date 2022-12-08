
import functools
import itertools
import os
from abnumber import Chain as AbChain
from abnumber.exceptions import ChainParseError
from nanome.api import ui, structure
from nanome.util import enums, Logs

from .utils import IMGTCDRColorScheme

ASSETS_FOLDER = os.path.join(os.getcwd(), 'plugin', 'assets')
CHAIN_BTN_JSON = os.path.join(ASSETS_FOLDER, 'chain_btn.json')
ZOOM_ICON_PNG = os.path.join(ASSETS_FOLDER, 'ZoomIcon.png')

class RegionMenu:

    def __init__(self, plugin):
        self._menu = ui.Menu()
        self._plugin = plugin

    @property
    def root(self):
        return self._menu.root

    @property
    def chain_btn_sets(self):
        row_lns = self.root.get_children()
        for row_ln in row_lns:
            for chain_btn_set in row_ln.get_children():
                yield chain_btn_set

    def build_menu(self, comp: structure.Complex):
        self._menu = ui.Menu()
        self._menu.root.layout_orientation = enums.LayoutTypes.vertical
        self._menu.title = f"{comp.full_name} Regions "
        comp.register_selection_changed_callback(self.on_selection_changed)
        antibody_chains = []
        for chain in comp.chains:
            seq_str = self._plugin.get_sequence_from_struct(chain)
            try:
                abchain = AbChain(seq_str, scheme='imgt')
            except ChainParseError as e:
                continue
            else:
                antibody_chains.append((chain, abchain))

        row_count = max(len(antibody_chains) // 6, 1)
        cols_per_row = len(antibody_chains) // row_count

        start_index = 0
        for _ in range(0, row_count):
            row_ln = self._menu.root.create_child_node()
            row_ln.layout_orientation = enums.LayoutTypes.horizontal
            for chain, abchain in antibody_chains[start_index:start_index + cols_per_row]:
                self.add_menu_chain_column(row_ln, chain, abchain)
            start_index += cols_per_row
        self.update_cdr_btns(self._menu, comp)

    def format_region_btn(self, region_name, mesh_color, cdr_residues):
        if not hasattr(self, '__prefab_btn'):
            json_path = os.path.join(os.getcwd(), 'plugin', 'assets', 'region_btn.json')
            self.__prefab_btn = ui.LayoutNode.io.from_json(json_path)
        ln_cdr = self.__prefab_btn.clone()
        cdr_btn = ln_cdr.get_children()[0].get_content()
        cdr_btn.text.value.set_all(region_name)
        cdr_mesh = ln_cdr.get_children()[0].get_children()[0].get_content()
        cdr_mesh.mesh_color = mesh_color
        cdr_btn.toggle_on_press = True
        cdr_btn.cdr_residues = cdr_residues
        cdr_btn.register_pressed_callback(
            functools.partial(self.on_cdr_btn_pressed, cdr_residues))
        return ln_cdr

    def format_chain_zoom_btn(self, chain_type: str, cdr_residues: list):
        if not hasattr(self, '__prefab_chain_btn'):
            self.__prefab_chain_btn = ui.LayoutNode.io.from_json(CHAIN_BTN_JSON)

        ln_chain_btn = self.__prefab_chain_btn.clone()
        chain_btn = ln_chain_btn.get_children()[0].get_content()
        chain_btn_icon = ln_chain_btn.get_children()[0].get_children()[0].get_content()
        chain_btn_label = ln_chain_btn.get_children()[0].get_children()[1].get_content()
        chain_btn_label.text_value = chain_type
        chain_btn_icon.file_path = ZOOM_ICON_PNG
        chain_btn.register_pressed_callback(
            functools.partial(self.on_chain_btn_pressed, cdr_residues))

        return ln_chain_btn

    def add_menu_chain_column(self, row_ln: ui.LayoutNode, chain: structure.Chain, abchain: AbChain):
        ln_chain = row_ln.create_child_node()
        ln_chain.chain_index = chain.index
        cdr1_residues = self._plugin.get_cdr1_residues(chain)
        cdr2_residues = self._plugin.get_cdr2_residues(chain)
        cdr3_residues = self._plugin.get_cdr3_residues(chain)
        cdr_residues = cdr1_residues + cdr2_residues + cdr3_residues
        ln_chain_btn = self.format_chain_zoom_btn(abchain.chain_type, cdr_residues)
        ln_chain.add_child(ln_chain_btn)

        chain_type = abchain.chain_type
        if chain_type == 'H':
            cdr1_color = IMGTCDRColorScheme.HEAVY_CDR1.value
            cdr2_color = IMGTCDRColorScheme.HEAVY_CDR2.value
            cdr3_color = IMGTCDRColorScheme.HEAVY_CDR3.value
        else:
            cdr1_color = IMGTCDRColorScheme.LIGHT_CDR1.value
            cdr2_color = IMGTCDRColorScheme.LIGHT_CDR2.value
            cdr3_color = IMGTCDRColorScheme.LIGHT_CDR3.value

        cdr1_region_name = f"CDR{abchain.chain_type}1"
        ln_cdr1 = self.format_region_btn(cdr1_region_name, cdr1_color, cdr1_residues)
        ln_chain.add_child(ln_cdr1)

        cdr2_region_name = f"CDR{abchain.chain_type}2"
        ln_cdr2 = self.format_region_btn(cdr2_region_name, cdr2_color, cdr2_residues)
        ln_chain.add_child(ln_cdr2)

        cdr3_region_name = f"CDR{abchain.chain_type}3"
        ln_cdr3 = self.format_region_btn(cdr3_region_name, cdr3_color, cdr3_residues)
        ln_chain.add_child(ln_cdr3)

    def on_cdr_btn_pressed(self, residue_list, btn):
        """When cdr button pressed, select all atoms in the residue_list."""
        for atom in itertools.chain(*[res.atoms for res in residue_list]):
            atom.selected = btn.selected
        self._plugin.update_structures_deep(residue_list)

    def on_chain_btn_pressed(self, residue_list, btn):
        self._plugin.zoom_on_structures(residue_list)

    def on_selection_changed(self, comp):
        """Update the selection in the plugin."""
        Logs.debug(f"Selection changes for {comp.full_name}")
        self.update_cdr_btns(self._menu, comp)
        Logs.debug("Finished updating selections")
        self._plugin.update_menu(self._menu)

    def update_cdr_btns(self, menu, comp):
        """Update the CDR buttons to reflect the current selections."""
        for ln in self.chain_btn_sets:
            # Get most up to date chain selections
            chain_index = ln.chain_index
            comp_chain = next(ch for ch in comp.chains if ch.index == chain_index)
            cdr1_btn = ln.get_children()[1].get_children()[0].get_content()
            cdr2_btn = ln.get_children()[2].get_children()[0].get_content()
            cdr3_btn = ln.get_children()[3].get_children()[0].get_content()
            seq_str = self._plugin.get_sequence_from_struct(comp_chain)
            try:
                abchain = AbChain(seq_str, scheme='imgt')
            except ChainParseError as e:
                Logs.debug(f"Could not parse Chain {comp_chain.name}")
                continue
            cdr1_residues = self._plugin.get_cdr1_residues(comp_chain, abchain=abchain)
            cdr2_residues = self._plugin.get_cdr2_residues(comp_chain, abchain=abchain)
            cdr3_residues = self._plugin.get_cdr3_residues(comp_chain, abchain=abchain)
            cdr1_btn.selected = any([
                atom.selected for atom in itertools.chain(*[res.atoms for res in cdr1_residues])])
            cdr2_btn.selected = any([
                atom.selected for atom in itertools.chain(*[res.atoms for res in cdr2_residues])])
            cdr3_btn.selected = any([
                atom.selected for atom in itertools.chain(*[res.atoms for res in cdr3_residues])])
