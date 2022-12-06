
import functools
import itertools
import os
from abnumber import Chain as AbChain
from abnumber.exceptions import ChainParseError
from nanome.api import ui, structure
from nanome.util import enums, Logs

from .utils import IMGTCDRColorScheme


class RegionMenu:

    def __init__(self, plugin):
        self._menu = ui.Menu()
        self._plugin = plugin

    def build_menu(self, comp: structure.Complex):
        self._menu = ui.Menu()
        self._menu.root.layout_orientation = enums.LayoutTypes.horizontal
        self._menu.title = f"{comp.full_name} Regions "
        comp.register_selection_changed_callback(self.on_selection_changed)
        for chain in comp.chains:
            seq_str = self._plugin.get_sequence_from_struct(chain)
            try:
                abchain = AbChain(seq_str, scheme='imgt')
            except ChainParseError as e:
                Logs.debug(f"Could not parse Chain {chain.name}")
                continue
            self.add_menu_chain_column(self._menu, chain, abchain)
        # self.update_cdr_btns(self._menu, comp)

    def add_menu_chain_column(self, menu: ui.Menu, chain: structure.Chain, abchain: AbChain):
        ln_chain = menu.root.create_child_node()
        ln_chain.chain_index = chain.index
        ln_chain_btn = ln_chain.create_child_node()
        chain_btn = ln_chain_btn.add_new_button(f'{abchain.chain_type}')

        cdr1_residues = self._plugin.get_cdr1_residues(chain)
        cdr2_residues = self._plugin.get_cdr2_residues(chain)
        cdr3_residues = self._plugin.get_cdr3_residues(chain)
        cdr_residues = cdr1_residues + cdr2_residues + cdr3_residues
        chain_btn.register_pressed_callback(
            functools.partial(self.on_chain_btn_pressed, cdr_residues))

        chain_type = abchain.chain_type
        if chain_type == 'H':
            cdr1_color = IMGTCDRColorScheme.HEAVY_CDR1.value
            cdr2_color = IMGTCDRColorScheme.HEAVY_CDR2.value
            cdr3_color = IMGTCDRColorScheme.HEAVY_CDR3.value
        else:
            cdr1_color = IMGTCDRColorScheme.LIGHT_CDR1.value
            cdr2_color = IMGTCDRColorScheme.LIGHT_CDR2.value
            cdr3_color = IMGTCDRColorScheme.LIGHT_CDR3.value

        ln_cdr1 = ui.LayoutNode.io.from_json(f"{os.getcwd()}/plugin/region_btn.json")
        cdr1_btn = ln_cdr1.get_children()[0].get_content()
        cdr1_btn.text.value.set_all(f"CDR{abchain.chain_type}1")
        cdr1_mesh = ln_cdr1.get_children()[0].get_children()[0].get_content()
        cdr1_mesh.mesh_color = cdr1_color
        ln_chain.add_child(ln_cdr1)
        cdr1_btn.toggle_on_press = True
        cdr1_btn.cdr_residues = cdr1_residues
        cdr1_btn.register_pressed_callback(
            functools.partial(self.on_cdr_btn_pressed, cdr1_residues))

        ln_cdr2 = ui.LayoutNode.io.from_json(f"{os.getcwd()}/plugin/region_btn.json")
        cdr2_btn = ln_cdr2.get_children()[0].get_content()
        cdr2_mesh = ln_cdr2.get_children()[0].get_children()[0].get_content()
        cdr2_mesh.mesh_color = cdr2_color
        cdr2_btn.text.value.set_all(f"CDR{abchain.chain_type}2")
        ln_chain.add_child(ln_cdr2)
        cdr2_btn.toggle_on_press = True
        cdr2_btn.cdr_residues = cdr2_residues
        cdr2_btn.register_pressed_callback(
            functools.partial(self.on_cdr_btn_pressed, cdr2_residues))

        ln_cdr3 = ui.LayoutNode.io.from_json(f"{os.getcwd()}/plugin/region_btn.json")
        cdr3_btn = ln_cdr3.get_children()[0].get_content()
        cdr3_mesh = ln_cdr3.get_children()[0].get_children()[0].get_content()
        cdr3_mesh.mesh_color = cdr3_color
        cdr3_btn.text.value.set_all(f"CDR{abchain.chain_type}3")
        ln_chain.add_child(ln_cdr3)
        cdr3_btn.toggle_on_press = True
        cdr3_btn.cdr_residues = cdr3_residues
        cdr3_btn.register_pressed_callback(
            functools.partial(self.on_cdr_btn_pressed, cdr3_residues))

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
        # self.update_cdr_btns(self._menu, comp)
        # Logs.debug("Finished updating selections")
        # self._plugin.update_menu(self._menu)
    
    def update_cdr_btns(self, menu, comp):
        """Update the CDR buttons to reflect the current selections."""
        for ln in menu.root.get_children():
            # Get most up to date chain selections
            chain_index = ln.chain_index
            comp_chain = next(ch for ch in comp.chains if ch.index == chain_index)
            cdr1_btn = ln.get_children()[1].get_content()
            cdr2_btn = ln.get_children()[2].get_content()
            cdr3_btn = ln.get_children()[3].get_content()
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
