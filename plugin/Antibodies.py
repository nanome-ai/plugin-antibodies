import itertools
import tempfile
import nanome
from nanome.util import async_callback, Color, Logs, enums

from abnumber import Chain as AbChain
from abnumber.exceptions import ChainParseError
from Bio import SeqIO, SeqUtils

from enum import Enum

protein_letters_3to1 = SeqUtils.IUPACData.protein_letters_3to1


class IMGTColorScheme1(Enum):
    CDR1 = Color(200, 0, 0)
    CDR2 = Color(255, 169, 0)
    CDR3 = Color(156, 65, 215)


class IMGTColorScheme2(Enum):
    CDR1 = Color(96, 96, 228)
    CDR2 = Color(70, 213, 0)
    CDR3 = Color(63, 157, 63)


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
        await self.check_chains(comp)

        # await self.highlight_cdr_loops(comp)
        self.set_plugin_list_button(run_btn, 'Run', True)
        Logs.debug("Done")

    async def check_chains(self, comp):
        for residue in comp.residues:
            residue.ribbon_color = Color.Grey()
            for atom in residue.atoms:
                atom.color = Color.Grey()
        comp.set_all_selected(False)
        await self.update_structures_deep([comp])
        for chain in comp.chains:
            Logs.message(f"Chain {chain.name}")
            seq_str = self.get_sequence_from_struct(chain)
            if not seq_str:
                Logs.warning(f"Unable to sequence chain {chain.name}")
                continue
            abchain = AbChain(seq_str, scheme='imgt')
            cdr3_seq = abchain.cdr3_seq
            if not cdr3_seq:
                Logs.message(f"No CDR3 in chain {chain.name}")
                continue
            
            cdr3_residues = self.get_cdr3_residues(chain)
            cdr3_atoms = itertools.chain(*[res.atoms for res in cdr3_residues])
            
            for res in cdr3_residues:
                if res.ribboned:
                    res.ribbon_color = IMGTColorScheme1.CDR3.value

                for atom in res.atoms:
                    atom.selected = True
                    atom.color =IMGTColorScheme1.CDR3.value

            self.update_structures_deep(cdr3_residues)

    async def highlight_cdr_loops(self, comp):
        Logs.debug("Finding loops...")
        run_btn = enums.PluginListButtonType.run
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

        Logs.debug("Formatting...")
        comp.set_all_selected(False)
        self.set_plugin_list_button(run_btn, 'Formatting...', False)
        color_scheme = IMGTColorScheme1
        for i, cdr_residues in enumerate([cdr1_residues, cdr2_residues, cdr3_residues]):
            # Add label to middle residue
            cdr_val = f"CDR{i + 1}"
            middle_residue = cdr_residues[(len(cdr_residues) // 2) + 1]
            middle_residue.labeled = True
            middle_residue.label_text = cdr_val
            # color CDR residues
            color_hex = getattr(color_scheme, cdr_val)
            for res in cdr_residues:
                for atom in res.atoms:
                    atom.color = Color(whole_num=color_hex)
        cdr_residues = cdr1_residues + cdr2_residues + cdr3_residues
        atom_iterchain = itertools.chain(*[res.atoms for res in cdr_residues])
        self.update_structures_shallow(list(atom_iterchain))

        # Zoom to CDR3 loop
        await self.zoom_on_structures(cdr1_residues)

    def get_cdr1_residues(self, comp):
        cdr_name = 'cdr1'
        residues = self._get_cdr_residues(comp, cdr_name)
        return residues

    def get_cdr2_residues(self, comp):
        cdr_name = 'cdr2'
        residues = self._get_cdr_residues(comp, cdr_name)
        return residues

    def get_cdr3_residues(self, struct):
        cdr_name = 'cdr3'
        residues = self._get_cdr_residues(struct, cdr_name)
        return residues

    def cdr_in_chain(self, chain_seq, cdr) -> bool:
        if cdr not in ['cdr1', 'cdr2', 'cdr3']:
            raise ValueError(f"Invalid cdr value: {cdr}. Valid choices are 'cdr1', 'cdr2', and 'cdr3'")
        abchain = AbChain(chain_seq)
        cdr_name = f'{cdr}_seq'
        cdr_seq = getattr(abchain, cdr_name)
        return cdr_seq in chain_seq

    def _get_cdr_residues(self, chain, cdr: str):
        """Get nanome residues corresponding to provided cdr name.

        valid cdr names are 'cdr1', 'cdr2', and 'cdr3'
        """
        if cdr not in ['cdr1', 'cdr2', 'cdr3']:
            raise ValueError(f"Invalid cdr name: {cdr}. Valid choices are 'cdr1', 'cdr2', and 'cdr3'")
        # pdb_file = tempfile.NamedTemporaryFile(suffix=".pdb", dir=self.temp_dir.name)
        # pdb_file = pdb_file.name
        # comp.io.to_pdb(path=pdb_file)

        seq_str = self.get_sequence_from_struct(chain)
        abchain = AbChain(seq_str, scheme='imgt')

        seq_attr_name = f'{cdr}_seq'
        cdr_seq = getattr(abchain, seq_attr_name)
        cdr_residues = []
        try:
            chain_seq = ''.join([protein_letters_3to1[res.name.title()] for res in chain.residues])
        except KeyError:
            return []
        if cdr_seq in chain_seq:
            residues = list(chain.residues)
            # Find the subset of the chain that comprises CDR Loop
            for idx in range(len(residues) - len(cdr_seq) + 1):
                residue_sublist = residues[idx:idx + len(cdr_seq)]
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
    
    def get_sequence_from_struct(self, struct):
        try:
            chain_seq = ''.join([
                protein_letters_3to1[res.name.title()]
                for res in struct.residues
            ])
        except KeyError:
            chain_seq = ''
        return chain_seq


def main():
    description = "Select antibody in entry list, then run plugin to add IMGT color scheme and highlight CDR loops."
    plugin = nanome.Plugin('Antibodies', description, 'other', False)
    plugin.set_plugin_class(Antibodies)
    plugin.run()


if __name__ == '__main__':
    main()
