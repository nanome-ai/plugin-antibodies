import nanome
from nanome.api.ui import Menu
from nanome.util import async_callback, Logs

from ABlooper import CDR_Predictor


class Antibodies(nanome.AsyncPluginInstance):

    def start(self):
        self.menu = Menu()
        self.menu.title = 'Antibodies'
        self.menu.width = 1
        self.menu.height = 1

        msg = 'Hello Nanome!'
        node = self.menu.root.create_child_node()
        self.label = node.add_new_label(msg)
        Logs.message(msg)

    @async_callback
    async def on_run(self):
        # Print the number of complexes in the workspace
        # to the Menu.
        # comps = await self.request_complex_list()
        # msg = f'{len(comps)} Complex(es) in Workspace'
        # Logs.message(msg)
        input_path = "3chn.pdb"
        output_path = "ABlooper_model.pdb"

        pred = CDR_Predictor(input_path, model='imgt')  # chains = ("H", "L")
        pred.write_predictions_in_pdb_format(output_path)
        # self.label.text_value = msg
        # self.menu.enabled = True
        # self.update_menu(self.menu)


def main():
    plugin = nanome.Plugin('Antibodies', 'Visualize and build antibody proteins', 'other', False)
    plugin.set_plugin_class(Antibodies)
    plugin.run()


if __name__ == '__main__':
    main()
