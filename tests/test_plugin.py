import asyncio
import os
import unittest
from abnumber import Chain as AbChain
from abnumber.exceptions import ChainParseError
from nanome.api import structure, PluginInstance
from plugin.Antibodies import Antibodies, IMGTCDRColorScheme
from unittest.mock import MagicMock, patch
fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')


def run_awaitable(awaitable, *args, **kwargs):
    loop = asyncio.get_event_loop()
    if loop.is_running:
        loop = asyncio.new_event_loop()
    result = loop.run_until_complete(awaitable(*args, **kwargs))
    loop.close()
    return result


class AntibodiesPluginTestCase(unittest.TestCase):

    def setUp(self):
        self.plugin = MagicMock()
        PluginInstance._instance = self.plugin
        self.pdb_file = os.path.join(fixtures_dir, '3chn.pdb')
        self.complex = structure.Complex.io.from_pdb(path=self.pdb_file)

    def test_color_complex_chains(self):
        """Validate that the complex is colored by component and chain."""
        comp = self.complex
        Antibodies.color_complex_chains(comp)
        # Get abchain to validate colors for each region.
        for chain in comp.chains:
            seq_str = Antibodies.get_sequence_from_struct(chain)
            if not seq_str:
                continue
            try:
                abchain = AbChain(seq_str, scheme='imgt')
            except ChainParseError as e:
                continue
            cdr1_residues = Antibodies.get_cdr1_residues(chain, abchain)
            cdr2_residues = Antibodies.get_cdr2_residues(chain, abchain)
            cdr3_residues = Antibodies.get_cdr3_residues(chain, abchain)
            fr1_residues = Antibodies.get_fr1_residues(chain, abchain)
            fr2_residues = Antibodies.get_fr2_residues(chain, abchain)
            fr3_residues = Antibodies.get_fr3_residues(chain, abchain)
            fr4_residues = Antibodies.get_fr4_residues(chain, abchain)
            
            # Assert each cdr region is colored properly.
            if abchain.chain_type == 'H':
                self.assertEqual(
                    set(res.ribbon_color for res in cdr1_residues),
                    set([IMGTCDRColorScheme.HEAVY_CDR1.value]))
                self.assertEqual(
                    set(res.ribbon_color for res in cdr2_residues),
                    set([IMGTCDRColorScheme.HEAVY_CDR2.value]))
                self.assertEqual(
                    set(res.ribbon_color for res in cdr3_residues),
                    set([IMGTCDRColorScheme.HEAVY_CDR3.value]))
            else:
                self.assertEqual(
                    set(res.ribbon_color for res in cdr1_residues),
                    set([IMGTCDRColorScheme.LIGHT_CDR1.value]))
                self.assertEqual(
                    set(res.ribbon_color for res in cdr2_residues),
                    set([IMGTCDRColorScheme.LIGHT_CDR2.value]))
                self.assertEqual(
                    set(res.ribbon_color for res in cdr3_residues),
                    set([IMGTCDRColorScheme.LIGHT_CDR3.value]))
            # Validate Framework color matches the chain type.
            self.assertEqual(
                set(res.ribbon_color for res in fr1_residues),
                set([IMGTCDRColorScheme.FR.value]))
            self.assertEqual(
                set(res.ribbon_color for res in fr2_residues),
                set([IMGTCDRColorScheme.FR.value]))
            self.assertEqual(
                set(res.ribbon_color for res in fr3_residues),
                set([IMGTCDRColorScheme.FR.value]))
            self.assertEqual(
                set(res.ribbon_color for res in fr4_residues),
                set([IMGTCDRColorScheme.FR.value]))
