
import functools
import itertools
import os
from abnumber import Chain as AbChain
from abnumber.exceptions import ChainParseError
from nanome.api import ui, structure
from nanome.util import enums, Logs, Color, async_callback

from .utils import IMGTCDRColorScheme

ASSETS_FOLDER = os.path.join(os.getcwd(), 'plugin', 'assets')
CHAIN_BTN_JSON = os.path.join(ASSETS_FOLDER, 'chain_btn.json')
CHECKMARK_PNG = os.path.join(ASSETS_FOLDER, 'checkmark.png')


class RegionMenu:

    def __init__(self, plugin):
        self._menu = ui.Menu()
        self._plugin = plugin
        self.index = plugin.current_menu_index
        self._menu.register_closed_callback(self.close_menu)
        self._menu.root.set_padding(0.01, 0.01, 0.01, 0.01)

    @property
    def root(self):
        return self._menu.root

    @property
    def index(self):
        return self._menu.index

    @index.setter
    def index(self, value):
        self._menu.index = value

    def enable(self):
        self._menu.enabled = True
        self._plugin.update_menu(self._menu)

    @property
    def chain_btn_sets(self):
        """Parse menu to get button sets for each chain."""
        row_lns = self.root.get_children()
        for row_ln in row_lns:
            for chain_btn_set in row_ln.get_children():
                yield chain_btn_set

    @property
    def region_btns(self):
        """Parse menu to get all region buttons (chain/cdrs)."""
        row_lns = self.root.get_children()
        for row_ln in row_lns:
            for chain_btn_set in row_ln.get_children():
                for ln_btn in chain_btn_set.get_children():
                    btn = ln_btn.get_children()[0].get_content()
                    yield btn

    @property
    def chain_btns(self):
        """Parse menu to get all chain buttons."""
        row_lns = self.root.get_children()
        for row_ln in row_lns:
            for region_btn_set in row_ln.get_children():
                # Chain button is always first child of region button set
                chain_ln_btn = region_btn_set.get_children()[0]
                btn = chain_ln_btn.get_children()[0].get_content()
                assert isinstance(btn, ui.Button)
                assert hasattr(btn, 'chain_type')
                yield btn

    def build_menu(self, comp: structure.Complex):
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
        self.update_cdr_btns(comp)

    def format_cdr_btn(self, region_name, mesh_color, cdr_residues):
        if not hasattr(self, '__prefab_btn'):
            json_path = os.path.join(os.getcwd(), 'plugin', 'assets', 'region_btn.json')
            self.__prefab_btn = ui.LayoutNode.io.from_json(json_path)
        ln_cdr = self.__prefab_btn.clone()
        cdr_btn = ln_cdr.get_children()[0].get_content()
        cdr_btn.icon.value.set_all(CHECKMARK_PNG)
        cdr_btn.text.value.set_all(region_name)
        cdr_btn.text.value.unusable = f"{region_name}..."
        cdr_btn.disable_on_press = True
        cdr_btn.icon.color.set_all(Color.White())
        cdr_mesh = ln_cdr.get_children()[0].get_children()[0].get_content()
        cdr_mesh.mesh_color = mesh_color
        cdr_btn.toggle_on_press = True
        cdr_btn.cdr_residues = cdr_residues
        cdr_btn.register_pressed_callback(
            functools.partial(self.on_cdr_btn_pressed, cdr_residues))
        return ln_cdr

    def format_chain_btn(self, chain: structure.Chain, chain_type: str):
        if not hasattr(self, '__prefab_chain_btn'):
            self.__prefab_chain_btn = ui.LayoutNode.io.from_json(CHAIN_BTN_JSON)

        ln_chain_btn = self.__prefab_chain_btn.clone()
        chain_btn = ln_chain_btn.get_children()[0].get_content()
        chain_index = chain.index
        chain_btn.chain_index = chain_index
        chain_btn.chain_type = chain_type
        chain_btn.toggle_on_press = True
        chain_btn.icon.value.set_all(CHECKMARK_PNG)
        btn_text = chain_type
        chain_btn.text.value.set_all(btn_text)
        chain_btn.text.value.unusable = f"{btn_text}..."
        chain_btn.disable_on_press = True
        chain_btn.register_pressed_callback(
            functools.partial(self.on_chain_btn_pressed, chain))
        chain_btn.selected = all([atom.selected for atom in chain.atoms])
        chain_btn.icon.active = chain_btn.selected
        return ln_chain_btn

    def add_menu_chain_column(self, row_ln: ui.LayoutNode, chain: structure.Chain, abchain: AbChain):
        ln_chain = row_ln.create_child_node()
        ln_chain.chain_index = chain.index
        cdr1_residues = self._plugin.get_cdr1_residues(chain)
        cdr2_residues = self._plugin.get_cdr2_residues(chain)
        cdr3_residues = self._plugin.get_cdr3_residues(chain)
        ln_chain_btn = self.format_chain_btn(chain, abchain.chain_type)
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
        ln_cdr1 = self.format_cdr_btn(cdr1_region_name, cdr1_color, cdr1_residues)
        ln_chain.add_child(ln_cdr1)

        cdr2_region_name = f"CDR{abchain.chain_type}2"
        ln_cdr2 = self.format_cdr_btn(cdr2_region_name, cdr2_color, cdr2_residues)
        ln_chain.add_child(ln_cdr2)

        cdr3_region_name = f"CDR{abchain.chain_type}3"
        ln_cdr3 = self.format_cdr_btn(cdr3_region_name, cdr3_color, cdr3_residues)
        ln_chain.add_child(ln_cdr3)

    @async_callback
    async def on_cdr_btn_pressed(self, residue_list, cdr_btn):
        """When cdr button pressed, select all atoms in the residue_list."""
        cdr_name = cdr_btn.text.value.selected
        Logs.message(f"CDR button {cdr_name} {'Selected' if cdr_btn.selected else 'Deselected'}")
        chain = residue_list[0].chain
        latest_chain = await self.get_latest_chain(chain)
        latest_residues = [res for res in latest_chain.residues if res.index in [res.index for res in residue_list]]
        for atom in itertools.chain.from_iterable(residue.atoms for residue in latest_residues):
            atom.selected = cdr_btn.selected
        cdr_btn.icon.active = cdr_btn.selected
        self._plugin.update_content(cdr_btn)
        self._plugin.update_structures_deep(latest_residues)

    async def get_latest_chain(self, chain):
        # Get latest version of chain from plugin, which should contain most recent colors/representation.
        comp = chain.complex
        [latest_comp] = await self._plugin.request_complexes([comp.index])
        latest_chain = next(ch for ch in latest_comp.chains if ch.index == chain.index)
        return latest_chain

    @async_callback
    async def on_chain_btn_pressed(self, chain: structure.Chain, chain_btn):
        chain_type = chain_btn.chain_type
        chain_btn.icon.active = chain_btn.selected
        Logs.message(f"Chain button {chain_type} {'Selected' if chain_btn.selected else 'Deselected'}")
        latest_chain = await self.get_latest_chain(chain)
        for atom in latest_chain.atoms:
            atom.selected = chain_btn.selected
        self._plugin.update_structures_deep([latest_chain])
        self._plugin.update_content(chain_btn)

    def on_selection_changed(self, comp):
        """Update the region buttons in the plugin when selection changed."""
        btns_selected = [btn.selected for btn in self.region_btns]
        Logs.debug(f"Selection changes on complex")
        self.update_cdr_btns(comp)
        # Update chain buttons
        for chain_btn in self.chain_btns:
            chain_index = chain_btn.chain_index
            chain = next((ch for ch in comp.chains if ch.index == chain_index), None)
            if not chain:
                continue
            chain_btn.selected = all([atom.selected for atom in chain.atoms])
            chain_btn.icon.active = chain_btn.selected
        updated_btns = [btn for btn in self.region_btns]
        changed_btns = [b for a, b in zip(btns_selected, updated_btns) if a != b.selected]

        if changed_btns:
            changed_btn_names = [btn.text.value.selected for btn in changed_btns]
            Logs.debug(f"Updating {', '.join(changed_btn_names)} buttons")
            self._plugin.update_content(changed_btns)

    def update_cdr_btns(self, comp):
        """Update the CDR buttons to reflect the current selections."""
        for ln in self.chain_btn_sets:
            # Get most up to date chain selections
            chain_index = ln.chain_index
            comp_chain = next((ch for ch in comp.chains if ch.index == chain_index), None)
            if not comp_chain:
                # Chain was probably deleted. Skip it.
                continue
            cdr1_btn = ln.get_children()[1].get_children()[0].get_content()
            cdr2_btn = ln.get_children()[2].get_children()[0].get_content()
            cdr3_btn = ln.get_children()[3].get_children()[0].get_content()
            seq_str = self._plugin.get_sequence_from_struct(comp_chain)
            try:
                abchain = AbChain(seq_str, scheme='imgt')
            except ChainParseError:
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
            cdr1_btn.icon.active = cdr1_btn.selected
            cdr2_btn.icon.active = cdr2_btn.selected
            cdr3_btn.icon.active = cdr3_btn.selected

    def close_menu(self, menu):
        """Delete menu when closed."""
        Logs.message(f"Closing menu {menu.index}")
        menu.enabled = False
        if menu.index in self._plugin.menus:
            del self._plugin.menus[menu.index]
