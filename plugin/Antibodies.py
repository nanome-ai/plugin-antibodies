import tempfile
import nanome
from nanome.util import async_callback, Logs, enums

from abnumber import Chain as AbChain
from abnumber.exceptions import ChainParseError
from Bio import SeqIO, SeqUtils


class Antibodies(nanome.AsyncPluginInstance):

    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    @async_callback
    async def on_run(self):
        # Get selected antibody complex
        run_btn = enums.PluginListButtonType.run
        Logs.debug("Loading Complex")
        self.set_plugin_list_button(run_btn, 'Loading Complex', False)
        comp_list = await self.request_complex_list()
        shallow_comp = next((cmp for cmp in comp_list if cmp.get_selected()), None)
        if not shallow_comp:
            self.send_notification(enums.NotificationTypes.error, "Please select an antibody")
            self.set_plugin_list_button(run_btn, 'Run', True)
            return
        comp = (await self.request_complexes([shallow_comp.index]))[0]
        comp.set_all_selected(False)
        # Highlight and zoom in on cdr3 loop
        Logs.debug("Finding loops...")
        self.set_plugin_list_button(run_btn, 'Finding loops...', False)
        try:
            cdr1_residues = self.get_cdr1_residues(comp)
            cdr2_residues = self.get_cdr2_residues(comp)
            cdr3_residues = self.get_cdr3_residues(comp)
        except ChainParseError as e:
            msg = "Please make sure selected structure is an antibody."
            self.send_notification(enums.NotificationTypes.error, msg)
            self.set_plugin_list_button(run_btn, 'Run', True)
            return

        # Set color scheme to IMGT after step above confirms complex is an antibody
        Logs.debug("Coloring complex")
        self.set_plugin_list_button(run_btn, 'Coloring...', False)
        self.update_structures_shallow(list(comp.atoms))
        self.apply_color_scheme(
            enums.ColorScheme.IMGT,
            enums.ColorSchemeTarget.All,
            False)

        for i, cdr_residues in enumerate([cdr1_residues, cdr2_residues, cdr3_residues]):
            # Add label to middle residue
            middle_residue = cdr_residues[len(cdr_residues)//2]
            middle_residue.labeled = True
            middle_residue.label_text = f"CDR{i + 1}"
            # Select CDR residues
            Logs.debug("Formatting...")
            self.set_plugin_list_button(run_btn, 'Formatting...', False)
            for res in cdr_residues:
                for atom in res.atoms:
                    atom.selected = True
        await self.update_structures_deep([comp])
        # Zoom to CDR3 loop
        await self.zoom_on_structures(cdr3_residues)
        self.set_plugin_list_button(run_btn, 'Run', True)

    def get_cdr1_residues(self, comp):
        cdr_name = 'cdr1'
        residues = self._get_cdr_residues(comp, cdr_name)
        return residues

    def get_cdr2_residues(self, comp):
        cdr_name = 'cdr2'
        residues = self._get_cdr_residues(comp, cdr_name)
        return residues

    def get_cdr3_residues(self, comp):
        cdr_name = 'cdr3'
        residues = self._get_cdr_residues(comp, cdr_name)
        return residues

    def _get_cdr_residues(self, comp, cdr: str):
        """Get nanome residues corresponding to provided cdr name.
        
        valid cdr names are 'cdr1', 'cdr2', and 'cdr3'
        """
        if cdr not in ['cdr1', 'cdr2', 'cdr3']:
            raise ValueError(f"Invalid cdr name: {cdr}. Valid choices are 'cdr1', 'cdr2', and 'cdr3'")
        pdb_file = tempfile.NamedTemporaryFile(suffix=".pdb", dir=self.temp_dir.name)
        pdb_file = pdb_file.name
        comp.io.to_pdb(path=pdb_file)

        seq_str = self.get_sequence_from_pdb(pdb_file)
        abchain = AbChain(seq_str, scheme='imgt')

        protein_letters_3to1 = SeqUtils.IUPACData.protein_letters_3to1

        seq_attr_name = f'{cdr}_seq'
        cdr_seq = getattr(abchain, seq_attr_name)
        cdr_residues = []
        for chain in comp.chains:
            try:
                chain_seq = ''.join([protein_letters_3to1[res.name.title()] for res in chain.residues])
            except KeyError:
                continue
            if cdr_seq in chain_seq:
                residues = list(chain.residues)
                # Find the subset of the chain that comprises CDR Loop
                for idx in range(len(residues) - len(cdr_seq) + 1):
                    residue_sublist = residues[idx:idx+len(cdr_seq)]
                    residue_sublist_seq = ''.join([
                        protein_letters_3to1[res.name.title()
                    ] for res in residue_sublist])
                    if residue_sublist_seq == cdr_seq:
                        cdr_residues = residue_sublist
                        break
        return cdr_residues

        

    def get_sequence_from_pdb(self, pdb_filepath):
        with open(pdb_filepath) as handle:
            sequence = next(SeqIO.parse(handle, "pdb-atom"))
        seq = str(sequence.seq)
        return seq

def main():
    plugin = nanome.Plugin('Antibodies', 'Visualize and build antibody proteins', 'other', False)
    plugin.set_plugin_class(Antibodies)
    plugin.run()


if __name__ == '__main__':
    main()
