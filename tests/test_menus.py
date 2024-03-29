import asyncio
import nanome
import os
import unittest
from nanome.api import structure, PluginInstance
from plugin.Antibodies import Antibodies
from random import randint
from unittest.mock import MagicMock
from plugin.menu import RegionMenu, SettingsMenu


fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')


class RegionMenuTestCase(unittest.IsolatedAsyncioTestCase):

    def setUp(self):
        self.plugin = Antibodies()
        PluginInstance._instance = self.plugin
        nanome._internal.network.plugin_network.PluginNetwork._instance = MagicMock()
        self.menu = RegionMenu(self.plugin)
        self.plugin.start()

        self.pdb_file = os.path.join(fixtures_dir, '2q8b.pdb')
        self.complex = structure.Complex.io.from_pdb(path=self.pdb_file)
        self._set_indices(self.complex)

    def test_build_menu(self):
        """Validate that the menu is built properly."""
        self.assertEqual(len(list(self.menu.chain_btn_sets)), 0)
        self.menu.build_menu(self.complex)
        self.assertEqual(len(list(self.menu.chain_btn_sets)), 2)

    def test_update_cdr_btns(self):
        # Validate button states change when a CDR is selected or deselected.
        # Assert None of the chains are selected, and all buttons are deselected.
        atoms_selected = any(atom.selected for atom in self.complex.atoms)
        self.assertFalse(atoms_selected)
        for ln_chain_col in self.menu._menu.root.get_children():
            for ln_btn in ln_chain_col.get_children()[1:]:
                btn = ln_btn.get_children()[0].get_content()
                self.assertFalse(btn.selected)
        heavy_chain = next(ch for ch in self.complex.chains if ch.name == 'H')
        for atom in heavy_chain.atoms:
            atom.selected = True

        self.menu.update_cdr_btns(self.complex)
        for ln_chain_col in self.menu.root.get_children():
            top_btn = ln_chain_col.get_children()[0].get_content()
            expected_selected = top_btn.text.active == 'H'
            for ln_btn in ln_chain_col.get_children()[1:]:
                btn = ln_btn.get_content()
                self.assertEqual(btn.selected, expected_selected)

    def test_on_selection_changed(self):
        # Validate that the menu is updated when a CDR is selected or deselected.
        self.plugin.update_content = MagicMock()
        self.menu.build_menu(self.complex)
        # Assert that content is not updated, because buttons are in sync with current selection.
        self.menu.on_selection_changed(self.complex)
        self.plugin.update_content.assert_not_called()
        # Select all atoms, and assert that content is updated.
        for atom in self.complex.atoms:
            atom.selected = True
        self.menu.on_selection_changed(self.complex)
        self.plugin.update_content.assert_called_once()

    async def test_on_cdr_btn_pressed(self):
        """Validate that the menu is updated when a CDR is selected or deselected."""
        self.plugin._network = MagicMock()
        # request_complexes is called when a CDR is selected.
        fut = asyncio.Future()
        fut.set_result([self.complex])
        self.plugin.request_complexes = MagicMock(return_value=fut)

        self.menu.build_menu(self.complex)
        btn = next(self.menu.chain_btn_sets).get_children()[1].get_children()[0].get_content()
        for atom in self.complex.atoms:
            atom.selected = False
        self.assertEqual(len(list(atom for atom in self.complex.atoms if atom.selected)), 0)
        i = 0
        res_set = set()
        for atom in self.complex.atoms:
            if i == 10:
                break
            res_set.add(atom.residue)
            i += 1
        btn.selected = True
        await self.menu.on_cdr_btn_pressed(list(res_set), btn)
        for res in res_set:
            for atom in res.atoms:
                self.assertTrue(atom.selected)

    def _set_indices(self, comp):
        """Set random indices for residues and atoms.

        Prevents issues with duplicate indices.
        """
        min_index = 100000
        max_index = 999999
        comp.index = randint(min_index, max_index)
        for chain in comp.chains:
            chain.index = randint(min_index, max_index)
            for residue in chain.residues:
                residue.index = randint(min_index, max_index)
                for atom in residue.atoms:
                    atom.index = randint(min_index, max_index)


class SettingsMenuTestCase(unittest.TestCase):

    def setUp(self):
        self.plugin = MagicMock()
        self.settings_menu = SettingsMenu(self.plugin)
        self.settings_menu.render()

    def test_get_numbering_scheme(self):
        """Validate that the menu is built properly."""
        expected_default_value = 'imgt'
        default_value = self.settings_menu.numbering_scheme
        self.assertEqual(default_value, expected_default_value)

        dd_numbering_scheme = self.settings_menu.dd_numbering_scheme
        valid_numbering_schemes = ['imgt', 'chothia', 'kabat']
        for scheme in valid_numbering_schemes:
            for item in dd_numbering_scheme.items:
                item.selected = item.name.lower() == scheme
            dd_numbering_scheme = self.settings_menu.dd_numbering_scheme
            self.assertEqual(self.settings_menu.numbering_scheme, scheme)
