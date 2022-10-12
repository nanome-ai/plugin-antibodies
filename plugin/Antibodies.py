from turtle import color
import nanome
from nanome.util import async_callback, Logs, enums

from abnumber import Chain as AbChain
from Bio import SeqIO


class Antibodies(nanome.AsyncPluginInstance):

    @async_callback
    async def on_run(self):
        input_path = "3chn.pdb"
        self.send_files_to_load([input_path])
        seq = self.get_sequence_from_pdb(input_path)
        # abchain = AbChain(seq, scheme='imgt')
        color_scheme = enums.ColorScheme.IMGT
        self.apply_color_scheme(color_scheme, enums.ColorSchemeTarget.All, False)
        Logs.message('Done')

    def get_sequence_from_pdb(self, pdb_filepath):
        with open(pdb_filepath) as handle:
            sequence = next(SeqIO.parse(handle, "pdb-atom"))
        seq = sequence.seq.__str__()
        return seq


def main():
    plugin = nanome.Plugin('Antibodies', 'Visualize and build antibody proteins', 'other', False)
    plugin.set_plugin_class(Antibodies)
    plugin.run()


if __name__ == '__main__':
    main()
