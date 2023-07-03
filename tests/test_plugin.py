import asyncio
import nanome
import os
import unittest
from abnumber import Chain as AbChain
from abnumber.exceptions import ChainParseError
from nanome.api import structure, PluginInstance
from plugin.Antibodies import Antibodies, IMGTCDRColorScheme
from random import randint
from unittest.mock import MagicMock


fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')


class AntibodiesPluginTestCase(unittest.IsolatedAsyncioTestCase):

    def setUp(self):
        self.plugin = Antibodies()
        PluginInstance._instance = self.plugin
        nanome._internal.network.plugin_network.PluginNetwork._instance = MagicMock()

        # Mock args that are passed to setup plugin instance networking
        session_id = plugin_network = pm_queue_in = pm_queue_out = custom_data = \
            log_pipe_conn = original_version_table = permissions = MagicMock()
        self.plugin._setup(
            session_id, plugin_network, pm_queue_in, pm_queue_out, log_pipe_conn,
            original_version_table, custom_data, permissions
        )
        self.plugin.start()

        self.pdb_file = os.path.join(fixtures_dir, '2q8b.pdb')
        self.complex = structure.Complex.io.from_pdb(path=self.pdb_file)

    def tearDown(self) -> None:
        self.plugin.on_stop()
        return super().tearDown()

    async def test_on_run(self):
        """Validate that the plugin starts properly."""
        
        comp = self.complex
        comp._selected = True
        fut = asyncio.Future()
        fut.set_result([comp])

        presenter_fut = asyncio.Future()
        presenter_fut.set_result(MagicMock())

        self.plugin.request_complex_list = MagicMock(return_value=fut)
        self.plugin.request_complexes = MagicMock(return_value=fut)
        self.plugin.request_presenter_info = MagicMock(return_value=presenter_fut)
        self.plugin.update_structures_deep = MagicMock()
        comps = await self.plugin.on_run()
        self._validate_complex_coloring(comps)

    async def test_integration(self):
        """Validate that the plugin starts properly."""        
        comp = self.complex
        request = MagicMock()
        request.get_args.return_value = [comp]
        # self.plugin.update_structures_deep = MagicMock()
        modified_comps = await self.plugin.integration_request(request)
        self._validate_complex_coloring(modified_comps)

    def test_prep_antibody_complex(self):
        """Validate that the complex is colored by component and chain."""
        comp = self.complex
        scheme = 'imgt'
        Antibodies.prep_antibody_complex(comp, scheme)
        self._validate_complex_coloring([comp])

    def _validate_complex_coloring(self, comps):
        # Get abchain to validate colors for each region.
        for comp in comps:
            for chain in comp.chains:
                seq_str = Antibodies.get_sequence_from_struct(chain)
                if not seq_str:
                    continue
                try:
                    abchain = AbChain(seq_str, scheme='imgt')
                except ChainParseError:
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
                    # Validate heavy chain CDRs
                    self.assertEqual(
                        set(res.ribbon_color.rgb for res in cdr1_residues),
                        set([IMGTCDRColorScheme.HEAVY_CDR1.value.rgb]))
                    self.assertEqual(
                        set(res.ribbon_color.rgb for res in cdr2_residues),
                        set([IMGTCDRColorScheme.HEAVY_CDR2.value.rgb]))
                    self.assertEqual(
                        set(res.ribbon_color.rgb for res in cdr3_residues),
                        set([IMGTCDRColorScheme.HEAVY_CDR3.value.rgb]))
                else:
                    # Validate light chain CDRs
                    self.assertEqual(
                        set(res.ribbon_color.rgb for res in cdr1_residues),
                        set([IMGTCDRColorScheme.LIGHT_CDR1.value.rgb]))
                    self.assertEqual(
                        set(res.ribbon_color.rgb for res in cdr2_residues),
                        set([IMGTCDRColorScheme.LIGHT_CDR2.value.rgb]))
                    self.assertEqual(
                        set(res.ribbon_color.rgb for res in cdr3_residues),
                        set([IMGTCDRColorScheme.LIGHT_CDR3.value.rgb]))
                # Validate Framework color matches the chain type.
                self.assertEqual(
                    set(res.ribbon_color.rgb for res in fr1_residues),
                    set([IMGTCDRColorScheme.FR.value.rgb]))
                self.assertEqual(
                    set(res.ribbon_color.rgb for res in fr2_residues),
                    set([IMGTCDRColorScheme.FR.value.rgb]))
                self.assertEqual(
                    set(res.ribbon_color.rgb for res in fr3_residues),
                    set([IMGTCDRColorScheme.FR.value.rgb]))
                self.assertEqual(
                    set(res.ribbon_color.rgb for res in fr4_residues),
                    set([IMGTCDRColorScheme.FR.value.rgb]))

    def test_validate_antibody(self):
        """Validate that antibodys and non-antibody structures are recognized."""
        is_antibody = Antibodies.validate_antibody(self.complex)
        self.assertTrue(is_antibody)
        non_antibody = structure.Complex.io.from_pdb(path=os.path.join(fixtures_dir, '1tyl.pdb'))
        is_antibody = Antibodies.validate_antibody(non_antibody)
        self.assertFalse(is_antibody)

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
