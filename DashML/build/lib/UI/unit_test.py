import unittest
from unittest.mock import patch
import pandas as pd
import configparser

from DashML.UI.DT_CLI import DashingTurtleCLI

class TestDashingTurtleCLI(unittest.TestCase):

    def setUp(self):
        self.cli = DashingTurtleCLI()
        # Provide a fake config
        self.cli.config = configparser.ConfigParser()
        self.cli.config.add_section('defaults')
        self.cli.config.set('defaults', 'temp', '37.0')
        self.cli.config.set('defaults', 'run', '1')
        self.cli.config.set('defaults', 'list_type', 'library')

    @patch("rich.console.Console.print")
    def test_do_config_prints_config(self, mock_print):
        self.cli.do_config("")
        mock_print.assert_called()

    @patch("DashML.Database_fx.Select_DB.select_library_full")
    @patch("rich.console.Console.print")
    def test_handle_seq_list(self, mock_rich_print, mock_db_select):
        # Mock dataframe
        df = pd.DataFrame({
            "ID": [1],
            "is_modified": [0],
            "contig": ["contig1"],
            "sequence_len": [10],
            "complex": [0],
            "sequence": ["AUGCUAUGC"],
            "timestamp": ["2022-01-01"],
            "secondary": [""],
            "name": ["TestSeq"]
        })
        mock_db_select.return_value = df

        self.cli.handle_seq_list([])
        mock_rich_print.assert_called()

    @patch("DashML.Database_fx.Insert_DB.insert_library")
    @patch("rich.console.Console.print")
    def test_handle_seq_add(self, mock_rich_print, mock_insert):
        mock_insert.return_value = 1
        args = ["-s", "AGCUAGCUA", "-n", "TestName", "-t1", "TypeA", "-t2", "TypeB"]
        self.cli.handle_seq_add(args)
        mock_rich_print.assert_called()

    @patch("rich.console.Console.print")
    def test_do_predict(self, mock_print):
        args = "-u 1,2 -l 3,4 -v T"
        self.cli.do_predict(args)
        mock_print.assert_called_with("Predicting with: unmod_lids=1,2, lids=3,4, vienna=T")

    @patch("rich.console.Console.print")
    def test_do_create_landscape(self, mock_print):
        args = "-l 3 -u 2 -o T"
        self.cli.do_create_landscape(args)
        mock_print.assert_called_with("Creating landscape with: lid=3, unmod_lid=2, optimize_clusters=T")

    @patch("rich.console.Console.print")
    def test_print_manual_known(self, mock_print):
        self.cli.print_manual("seq")
        mock_print.assert_called()

    @patch("rich.console.Console.print")
    def test_print_manual_unknown(self, mock_print):
        self.cli.print_manual("unknown")
        mock_print.assert_called_with("No manual entry for 'unknown'. Available: seq, load, predict, create_landscape, config, exit")

    def test_resolve_subcommand_unique(self):
        cmd = self.cli.resolve_subcommand("-l", ["-list", "-add"])
        self.assertEqual(cmd, "-list")

    @patch("rich.console.Console.print")
    def test_resolve_subcommand_ambiguous(self, mock_print):
        self.cli.resolve_subcommand("-a", ["-add", "-align"])
        mock_print.assert_called()

    @patch("rich.console.Console.print")
    def test_resolve_subcommand_unknown(self, mock_print):
        self.cli.resolve_subcommand("-x", ["-list", "-add"])
        mock_print.assert_called()

if __name__ == "__main__":
    unittest.main()
