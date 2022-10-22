import itertools
import tempfile
import nanome
from nanome.util import async_callback, Color, Logs, enums

from abnumber import Chain as AbChain
from abnumber.exceptions import ChainParseError
from Bio import SeqIO, SeqUtils

from enum import Enum

protein_letters_3to1 = SeqUtils.IUPACData.protein_letters_3to1
run_btn = enums.PluginListButtonType.run


class IMGTCDRColorScheme(Enum):
    """
    Source: https://www.imgt.org/IMGTScientificChart/RepresentationRules/colormenu.php#h1_26
    """
    HEAVY_CDR1 = Color(200, 0, 0)
    HEAVY_CDR2 = Color(255, 169, 0)
    HEAVY_CDR3 = Color(156, 65, 215)
    LIGHT_CDR1 = Color(96, 96, 228)
    LIGHT_CDR2 = Color(70, 213, 0)
    LIGHT_CDR3 = Color(63, 157, 63)
    # Added by us
    FR1 = Color(211, 211, 211)
    FR2 = Color(189, 189, 189)
    FR3 = Color(158, 158, 158)
    FR4 = Color(125, 125, 125)


class Antibodies(nanome.AsyncPluginInstance):

    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    @async_callback
    async def on_run(self):
        # Get selected antibody complex
        Logs.debug("Loading Complex")
        self.set_plugin_list_button(run_btn, 'Loading Complex', False)
        comp_list = await self.request_complex_list()
        shallow_comp = next((cmp for cmp in comp_list if cmp.get_selected()), None)
        if not shallow_comp:
            self.send_notification(enums.NotificationTypes.error, "Please select an antibody")
            self.set_plugin_list_button(run_btn, 'Run', True)
            return
        comp = (await self.request_complexes([shallow_comp.index]))[0]
        await self.highlight_cdr_loops(comp)
        self.set_plugin_list_button(run_btn, 'Run', True)
        Logs.debug("Done")

    async def highlight_cdr_loops(self, comp):
        if not self.validate_antibody(comp):
            self.send_notification(enums.NotificationTypes.error, "Selected complex is not an antibody")
            return
        
        # Make entire complex Grey.
        Logs.debug("Making Complex Grey")
        self.set_plugin_list_button(run_btn, 'Coloring...', False)
        for residue in comp.residues:
            residue.ribbon_color = Color.Grey()
            for atom in residue.atoms:
                atom.color = Color.Grey()

        self.set_plugin_list_button(run_btn, 'Coloring...', False)
        comp.set_all_selected(False)
        self.update_structures_deep([comp])
        # Loop through chain and color cdr loops
        Logs.debug("Processing Chains.")
        self.set_plugin_list_button(run_btn, 'Finding CDR Loops...', False)
        for chain in comp.chains:
            Logs.debug(f"Chain {chain.name}")
            seq_str = self.get_sequence_from_struct(chain)
            if not seq_str:
                Logs.warning(f"Unable to sequence chain {chain.name}")
                continue
            try:
                abchain = AbChain(seq_str, scheme='imgt')
            except ChainParseError as e:
                Logs.warning(f"Could not parse Chain {chain.name}")
                continue
            
            cdr3_seq = abchain.cdr3_seq
            if not cdr3_seq:
                Logs.debug(f"No CDR3 in chain {chain.name}")
                continue

            try:    
                cdr1_residues = self.get_cdr1_residues(chain)
                cdr2_residues = self.get_cdr2_residues(chain)
                cdr3_residues = self.get_cdr3_residues(chain)
                fr1_residues = self.get_fr1_residues(chain)
                fr2_residues = self.get_fr2_residues(chain)
                fr3_residues = self.get_fr3_residues(chain)
                fr4_residues = self.get_fr4_residues(chain)
            except ChainParseError as e:
                Logs.warning(f"Could find cdr loops for Chain {chain.name}")
                continue
            fr1_color = IMGTCDRColorScheme.FR1.value
            fr2_color = IMGTCDRColorScheme.FR2.value
            fr3_color = IMGTCDRColorScheme.FR3.value
            fr4_color = IMGTCDRColorScheme.FR4.value

            chain_type = abchain.chain_type
            if chain_type == 'H':
                cdr1_color = IMGTCDRColorScheme.HEAVY_CDR1.value
                cdr2_color = IMGTCDRColorScheme.HEAVY_CDR2.value
                cdr3_color = IMGTCDRColorScheme.HEAVY_CDR3.value
            else:
                cdr1_color = IMGTCDRColorScheme.LIGHT_CDR1.value
                cdr2_color = IMGTCDRColorScheme.LIGHT_CDR2.value
                cdr3_color = IMGTCDRColorScheme.LIGHT_CDR3.value
            
            cdr_residue_lists = [cdr1_residues, cdr2_residues, cdr3_residues]
            fr_residue_lists = [fr1_residues, fr2_residues, fr3_residues, fr4_residues]
            cdr_colors = [cdr1_color, cdr2_color, cdr3_color]
            fr_colors = [fr1_color, fr2_color, fr3_color, fr4_color]

            residue_lists = cdr_residue_lists + fr_residue_lists
            chain_colors = cdr_colors +  fr_colors
            
            i = 0
            for residue_list, res_color in zip(residue_lists, chain_colors):
                # Add label to middle residue
                # The last 4 values are Fr1, Fr2, Fr3, Fr4
                if i < len(list(residue_lists)) - 4:
                    label_val = f"CDR{chain_type}{i + 1}"
                else:
                    label_val = f"FR{chain_type}{i + 1}"
                middle_residue = residue_list[(len(residue_list) // 2)]
                middle_residue.labeled = True
                middle_residue.label_text = label_val
                # Set residue and atom colors
                for res in residue_list:
                    res.ribbon_color = res_color
                    for atom in res.atoms:
                        atom.color = res_color
                i += 1
                self.update_structures_deep(residue_list)
            # Get residues from list of lists
        self.set_plugin_list_button(run_btn, 'Done', False)

    def validate_antibody(self, comp):
        # Make sure at least one chain can be parsed with ABChain
        for chain in comp.chains:
            seq_str = self.get_sequence_from_struct(chain)
            try:
                abchain = AbChain(seq_str, scheme='imgt')
                if abchain:
                    return True
            except ChainParseError as e:
                continue
            return False

    def get_cdr1_residues(self, struc):
        cdr_name = 'cdr1'
        residues = self._get_cdr_residues(struc, cdr_name)
        return residues

    def get_cdr2_residues(self, struc):
        cdr_name = 'cdr2'
        residues = self._get_cdr_residues(struc, cdr_name)
        return residues

    def get_cdr3_residues(self, struc):
        cdr_name = 'cdr3'
        residues = self._get_cdr_residues(struc, cdr_name)
        return residues

    def get_fr1_residues(self, struc):
        fr_name = 'fr1'
        residues = self._get_cdr_residues(struc, fr_name)
        return residues

    def get_fr2_residues(self, struc):
        fr_name = 'fr2'
        residues = self._get_cdr_residues(struc, fr_name)
        return residues
    
    def get_fr3_residues(self, struc):
        fr_name = 'fr3'
        residues = self._get_cdr_residues(struc, fr_name)
        return residues
    
    def get_fr4_residues(self, struc):
        fr_name = 'fr4'
        residues = self._get_cdr_residues(struc, fr_name)
        return residues

    def _get_cdr_residues(self, chain, cdr: str):
        """Get nanome residues corresponding to provided cdr name.

        valid cdr names are 'cdr1', 'cdr2', and 'cdr3'
        """
        if cdr not in ['cdr1', 'cdr2', 'cdr3', 'fr1', 'fr2', 'fr3', 'fr4']:
            raise ValueError(f"Invalid cdr name: {cdr}. Valid choices are 'cdr1', 'cdr2', and 'cdr3'")
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
