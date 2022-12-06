import asyncio
import itertools
import nanome
import os
import unittest
from nanome.api import structure, PluginInstance
from plugin.Antibodies import Antibodies
from random import randint
from unittest.mock import MagicMock
from plugin.menu import RegionMenu

fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')


def run_awaitable(awaitable, *args, **kwargs):
    loop = asyncio.get_event_loop()
    if loop.is_running:
        loop = asyncio.new_event_loop()
    result = loop.run_until_complete(awaitable(*args, **kwargs))
    loop.close()
    return result


class RegionMenuTestCase(unittest.TestCase):

    def setUp(self):
        self.plugin = Antibodies()
        PluginInstance._instance = self.plugin
        nanome._internal._network._plugin_network.PluginNetwork._instance = MagicMock()
        self.menu = RegionMenu(self.plugin)

        self.pdb_file = os.path.join(fixtures_dir, '2q8b.pdb')
        self.complex = structure.Complex.io.from_pdb(path=self.pdb_file)

    def test_build_menu(self):
        """Validate that the menu is built properly."""
        self.menu.build_menu(self.complex)
        self.assertEqual(len(self.menu.root.get_children()), 2)

    def test_update_cdr_btns(self):
        # Validate button states change when a CDR is selected or deselected.
        # Assert None of the chains are selected, and all buttons are deselected.
        atoms_selected = any(atom.selected for atom in self.complex.atoms)
        self.assertFalse(atoms_selected)
        for ln_chain_col in self.menu._menu.root.get_children():
            for ln_btn in ln_chain_col.get_children()[1:]:
                breakpoint()
                btn = ln_btn.get_children()[0].get_content()
                self.assertFalse(btn.selected)
        heavy_chain = next(ch for ch in self.complex.chains if ch.name == 'H')
        for atom in heavy_chain.atoms:
            atom.selected = True

        self.menu.update_cdr_btns(self.menu._menu, self.complex)
        for ln_chain_col in self.menu.root.get_children():
            top_btn = ln_chain_col.get_children()[0].get_content()
            expected_selected = top_btn.text.active == 'H'
            for ln_btn in ln_chain_col.get_children()[1:]:
                btn = ln_btn.get_content()
                self.assertEqual(btn.selected, expected_selected)

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
