import tempfile
import nanome
from nanome.util import async_callback, Logs, enums

from abnumber import Chain as AbChain
from Bio import SeqIO, SeqUtils


class Antibodies(nanome.AsyncPluginInstance):

    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    @async_callback
    async def on_run(self):
        # Get selected antibody complex
        run_btn = enums.PluginListButtonType.run
        self.set_plugin_list_button(run_btn, 'Loading Complex', False)
        comp_list = await self.request_complex_list()
        shallow_comp = next((cmp for cmp in comp_list if cmp.get_selected()), None)
        if not shallow_comp:
            self.send_notification(enums.NotificationTypes.error, "Please select an antibody")
            self.set_plugin_list_button(run_btn, 'Run', True)
            return
        comp = (await self.request_complexes([shallow_comp.index]))[0]

        # Set color scheme to IMGT
        self.set_plugin_list_button(run_btn, 'Applying Color Scheme...', False)
        self.update_structures_shallow(list(comp.atoms))
        self.apply_color_scheme(
            enums.ColorScheme.IMGT,
            enums.ColorSchemeTarget.All,
            False)

        comp.set_all_selected(False)
        
        # Highlight and zoom in on cdr3 loop
        self.set_plugin_list_button(run_btn, 'Finding CDR3...', False)
        cdr3_residues = self.get_cdr3_residues(comp)
        if cdr3_residues:
            # Add label to middle residue
            middle_residue = cdr3_residues[len(cdr3_residues)//2]
            middle_residue.labeled = True
            middle_residue.label_text = "CDR3"
            # Select CDR3 residues 
            self.set_plugin_list_button(run_btn, 'Zooming...', False)
            for res in cdr3_residues:
                for atom in res.atoms:
                    atom.selected = True
            await self.update_structures_deep([comp])
            # Zoom to CDR3 loop
            await self.zoom_on_structures(cdr3_residues)
            self.set_plugin_list_button(run_btn, 'Run', True)

    def get_cdr3_residues(self, comp):
        pdb_file = tempfile.NamedTemporaryFile(suffix=".pdb", dir=self.temp_dir.name)
        pdb_file = pdb_file.name
        comp.io.to_pdb(path=pdb_file)

        seq_str = self.get_sequence_from_pdb(pdb_file)
        abchain = AbChain(seq_str, scheme='imgt')
        
        protein_letters_3to1 = SeqUtils.IUPACData.protein_letters_3to1
        cdr3_seq = abchain.cdr3_seq
        cdr3_residues = []
        for chain in comp.chains:
            chain_seq = ''.join([protein_letters_3to1[res.name.title()] for res in chain.residues])
            if cdr3_seq in chain_seq:
                residues = list(chain.residues)
                # Find the subset of the chain that comprises CDR3 Loop
                for idx in range(len(residues) - len(cdr3_seq) + 1):
                    residue_sublist = residues[idx:idx+len(cdr3_seq)]
                    residue_sublist_seq = ''.join([
                        protein_letters_3to1[res.name.title()
                    ] for res in residue_sublist])
                    if residue_sublist_seq == cdr3_seq:
                        cdr3_residues = residue_sublist
                        break
        return cdr3_residues

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
