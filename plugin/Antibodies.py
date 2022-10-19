import nanome
from nanome.api.structure import Complex
from nanome.util import async_callback, Logs, enums

from abnumber import Chain as AbChain
from Bio import SeqIO, SeqUtils


class Antibodies(nanome.AsyncPluginInstance):

    @async_callback
    async def on_run(self):
        # input_path = "3chn.pdb"
        # self.send_files_to_load([input_path])
        # time.sleep(2)
        comp_list = await self.request_complex_list()
        shallow_comp = next((cmp for cmp in comp_list if cmp.get_selected()), None)
        if not shallow_comp:
            self.send_notification(enums.NotificationTypes.error, "Please select an antibody")
            return
        comp = (await self.request_complexes([shallow_comp.index]))[0]

        comp.set_all_selected(False)
        # Set color scheme to IMGT
        self.update_structures_shallow(list(comp.atoms))
        self.apply_color_scheme(
            enums.ColorScheme.IMGT,
            enums.ColorSchemeTarget.All,
            False)

        pdb_file = 'test_pdb.pdb'
        comp.io.to_pdb(path=pdb_file)
        # Zoom in on CDR3 loop
        seq_str = self.get_sequence_from_pdb(pdb_file)
        abchain = AbChain(seq_str, scheme='imgt')
        
        protein_letters_3to1 = SeqUtils.IUPACData.protein_letters_3to1
        cdr3_seq = abchain.cdr3_seq
        for chain in comp.chains:
            chain_seq = ''.join([protein_letters_3to1[res.name.title()] for res in chain.residues])
            if cdr3_seq in chain_seq:
                residues = list(chain.residues)
                cdr3_residues = None
                # Find the subset of the chain that comprises CDR3 Loop
                for idx in range(len(residues) - len(cdr3_seq) + 1):
                    residue_sublist = residues[idx:idx+len(cdr3_seq)]
                    residue_sublist_seq = ''.join([
                        protein_letters_3to1[res.name.title()
                    ] for res in residue_sublist])
                    if residue_sublist_seq == cdr3_seq:
                        cdr3_residues = residue_sublist
                        break
                if cdr3_residues:
                    for res in cdr3_residues:
                        for atom in res.atoms:
                            atom.selected = True
                        # map(lambda atom: setattr(atom, 'selected', True), res.atoms)
                    await self.update_structures_deep(cdr3_residues)
                    await self.zoom_on_structures(cdr3_residues)
                    break
                Logs.message("found it!")

        # chain_seqs = {}
        # for chain in comp.chains:
        #     aa_names = [protein_letters_3to1[res.name.title()] for res in chain.residues]
        #     seq = "".join(aa_names)
        #     Logs.message(f"Chain {chain.name}: {seq}")
        #     abchain = AbChain(seq, scheme='imgt')
        #     if abchain.cdr3_dict:
        #         Logs.message('cdr3 found')
        #     chain_seqs[chain.name] = seq
        #     pass
        #     # Logs.debug(str([len(val) for val in abchain.regions.values()]))
        Logs.message('Done')

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
