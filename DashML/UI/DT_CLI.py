import sys, re
import traceback
import readline
import cmd
import configparser
import os
import pandas as pd
from colorama import init, Fore, Style
from rich.console import Console
from rich.table import Table
import argparse
import DashML.Basecall.run_basecall as run_basecall
import DashML.Landscape.Cluster.run_landscape as landscape
import DashML.Predict.run_predict as predict
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_DB as dbsel


pd.options.mode.chained_assignment = None

try:
    import readline
except ImportError:
    print("Module readline not available.")
else:
    import rlcompleter
    readline.parse_and_bind("tab: complete")
    readline.parse_and_bind("set editing-mode emacs")


# Configure libedit on macOS (if needed)
if 'libedit' in readline.__doc__:
    readline.parse_and_bind("bind ^I rl_complete")
    readline.parse_and_bind("bind ^A beginning-of-line")
    readline.parse_and_bind("bind ^E end-of-line")
    readline.parse_and_bind("bind ^K kill-line")
    readline.parse_and_bind("bind ^Y yank")
    readline.parse_and_bind("bind ^P previous-history")
    readline.parse_and_bind("bind ^N next-history")
else:
    # GNU readline (Linux or correctly installed on macOS)
    readline.parse_and_bind("tab: complete")

init(autoreset=True)

CONFIG_LOCATIONS = [
    '/etc/dashingrc',
    os.path.expanduser('~/.dashingrc'),
    './dashing_config.ini'
]

def handle_cli_errors(func):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            print("❌ Command failed!")
            print(f"Details: {e}")
            # Uncomment if you want full traceback
            # traceback.print_exc()
    return wrapper

class DashingTurtleCLI(cmd.Cmd):
    intro = Fore.GREEN + "Welcome to DashingTurtle! Type 'help' or '?' to list commands.\n"
    prompt = Fore.YELLOW + "🐢> " + Style.RESET_ALL
    VALID_LIST_TYPES = {'library', 'unmodified', 'modified'}

    @handle_cli_errors
    def db_check(self, lid, unmod_lid=None):
        """
        Check if a library ID exists in the database.
        Returns:
            check (bool), df (DataFrame or None)
        """
        errors = []
        check = False
        df = None

        try:
            lid_int = int(lid)
        except ValueError:
            errors.append(f"LID '{lid}' is not a valid integer.")

        if unmod_lid is not None:
            try:
                unmod_lid_int = int(unmod_lid)
            except ValueError:
                errors.append(f"Unmodified LID '{unmod_lid}' is not a valid integer.")

        # If integer conversion errors occurred, print errors and return early
        if errors:
            print(Fore.RED + "\n".join(errors))
            return check, None

        if unmod_lid is None:
            df = dbsel.select_librarybyid(lid_int)
            if len(df) == 0:
                errors.append(f"LID {lid} not found in database. Please run seq -list to verify.")
            else:
                check = True
        else:
            df = dbsel.select_library_full()
            df_mod = df[df['ID'] == lid_int]
            df_unmod = df[df['ID'] == unmod_lid_int]

            # Check if each LID exists individually
            if df_mod.empty:
                errors.append(f"LID {lid} not found in database. Please run seq -list to verify.")

            if df_unmod.empty:
                errors.append(f"LID {unmod_lid} not found in database. Please run seq -list to verify.")

            # Only continue if both exist
            if not df_mod.empty and not df_unmod.empty:
                # Check modification status of mod
                if df_mod['type1'].unique()[0] == df_mod['type2'].unique()[0]:
                    errors.append(f"LID {lid} is listed as NOT modified. Please run seq -list to verify.")

                # Check modification status of unmod
                if df_unmod['type1'].unique()[0] != df_unmod['type2'].unique()[0]:
                    errors.append(f"LID {unmod_lid} is listed as modified. Please run seq -list to verify.")

                # Check contig match
                contig_mod = df_mod['contig'].unique()[0]
                contig_unmod = df_unmod['contig'].unique()[0]
                if contig_mod != contig_unmod:
                    errors.append(
                        f"Contig {contig_mod} and {contig_unmod} are not the same sequence. Please run seq -list to verify."
                    )
            else:
                # Already added specific "not found" errors above, so no generic combined error needed
                pass
        if not errors:
            check = True

        if errors:
            print(Fore.RED + "\n".join(errors))

        return check, df

    def preloop(self):
        self.commands = ['seq', 'load', 'predict', 'create_landscape', 'exit', 'help', 'man', 'config']
        self.seq_subcommands = {
            '-list': self.handle_seq_list,
            '-add': self.handle_seq_add,
        }
        self.config = self.load_config()
        self.validate_config()

    def load_config(self):
        config = configparser.ConfigParser()
        for path in CONFIG_LOCATIONS:
            if os.path.exists(path):
                config.read(path)
        return config

    def get_config_default(self, section, key, fallback=None):
        return self.config.get(section, key, fallback=fallback)

    def validate_config(self):
        if not self.config.sections():
            print(Fore.YELLOW + "No config file loaded; skipping validation.\n")
            return

        errors = []
        if not self.config.has_section('defaults'):
            errors.append("Missing [defaults] section in config.")

        temp = self.config.get('defaults', 'temp', fallback=None)
        if temp is None:
            errors.append("Missing 'temp' in [defaults].")
        else:
            try:
                float(temp)
            except ValueError:
                errors.append(f"Invalid 'temp': must be a number, got '{temp}'.")

        run = self.config.get('defaults', 'run', fallback=None)
        if run is None:
            errors.append("Missing 'run' in [defaults].")
        else:
            try:
                int(run)
            except ValueError:
                errors.append(f"Invalid 'run': must be an integer, got '{run}'.")

        list_type = self.config.get('defaults', 'list_type', fallback=None)
        if list_type and list_type not in self.VALID_LIST_TYPES:
            errors.append(f"Invalid 'list_type': must be one of {self.VALID_LIST_TYPES}, got '{list_type}'.")

        if errors:
            print(Fore.RED + "Configuration validation failed:")
            for err in errors:
                print(Fore.RED + f" - {err}")
            sys.exit(1)
        else:
            print(Fore.GREEN + "Configuration loaded and validated successfully.\n")

    @handle_cli_errors
    def do_config(self, arg):
        "Show loaded configuration values."
        if not self.config.sections():
            print(Fore.YELLOW + "No config loaded.")
            return
        for section in self.config.sections():
            print(Fore.CYAN + f"[{section}]")
            for key, value in self.config[section].items():
                print(f"  {key} = {value}")
        print()

    def print_dataframe_rich(self, df):
        console = Console()
        table = Table(show_header=True, header_style="bold magenta")

        for column in df.columns:
            table.add_column(column, min_width=10, overflow="ellipsis")

        for _, row in df.iterrows():
            table.add_row(*[str(x) for x in row])

        console.print(table)
    def handle_seq_list(self, args):
        list_type = args[0] if args and args[0] in self.VALID_LIST_TYPES else self.get_config_default("defaults", "list_type", "library")
        print(Fore.CYAN + f"Listing sequences of type '{list_type}'")
        df = dbsel.select_library_full()
        df.rename(columns={'ID': 'LID'}, inplace=True)
        df.drop(columns=['contig', 'sequence_len', 'complex', 'sequence', 'timestamp', 'secondary'], inplace=True)
        if list_type == 'unmodified':
            df = df.loc[df['is_modified'] == 0]
            df.drop(columns=['is_modified'], inplace=True)
        elif list_type =='modified':
            df = df.loc[df['is_modified'] == 1]
            df.drop(columns=['is_modified'], inplace=True)
        df.drop(columns=['is_putative', 'experiment'], inplace=True)
        df.rename(columns={'sequence_name': 'Seq Name', 'temp': 'Temp', 'type1':'Condition1',
                           'type2':'Condition2', 'run':'Run'}, inplace=True)
        self.print_dataframe_rich(df)
    def handle_seq_add(self, args):
        parser = argparse.ArgumentParser(prog="seq -add", description="Add a new sequence")
        parser.add_argument("-s", "--sequence", required=True, help="RNA sequence")
        parser.add_argument("-sec", "--secondary", default='', help="Secondary structure")
        parser.add_argument("-e", "--experiment", default='', help="Control Experiment type")
        parser.add_argument("-n", "--name", required=True, help="Sequence name")
        parser.add_argument("-t", "--temp", default=self.get_config_default("defaults", "temp", "37.0"),
                            help="Temperature (°C)")
        parser.add_argument("-t1", "--type1", required=True, help="First type label")
        parser.add_argument("-t2", "--type2", required=True, help="Second type label")
        parser.add_argument("-r", "--run", default=self.get_config_default("defaults", "run", "1"), help="Run ID. Default 1.")

        try:
            # Create a console locally
            console = Console()
            opts = parser.parse_args(args)
            print(
                Fore.CYAN + f"Adding sequence with: {opts.sequence}, {opts.secondary}, {opts.experiment}, {opts.name}, {opts.temp}, {opts.type1}, {opts.type2}, {opts.run}")

            sequence = opts.sequence
            # Validate that the sequence contains only valid RNA characters
            if not re.fullmatch(r"^[ACGU]+$", sequence.upper()):
                console = Console()
                console.print("[red]Error: Sequence must contain only A, C, G, or U characters.[/red]")
                return

            secondary = opts.secondary
            # Check that secondary structure is valid dot-bracket notation
            if not re.fullmatch(r"^[().]+$", secondary):
                console.print("[red]Error: Secondary structure must contain only '.', '(', and ')' characters.[/red]")
                return

            # Validate secondary structure length
            if secondary:

                # Check secondary structure length
                if len(secondary) != len(sequence):
                    console.print(
                        f"[red]Error: Secondary structure length ({len(secondary)}) does not match sequence length ({len(sequence)}).[/red]"
                    )
                    return

                # Check that experiment is not blank
                if not opts.experiment.strip():
                    console.print(
                        "[red]Error: Experiment ID must be provided when a secondary structure is specified.[/red]"
                    )
                    return

            # Create DataFrame
            df = pd.DataFrame([{
                "contig": opts.name,
                "sequence": sequence,
                "secondary": secondary,
                "sequence_len": len(sequence),
                "experiment": opts.experiment,
                "sequence_name": opts.name,
                "temp": opts.temp,
                "type1": opts.type1,
                "type2": opts.type2,
                "is_modified": (0 if opts.type1 == opts.type2 else 1),
                "complex":0,
                "run": opts.run
            }])
            lid = dbins.insert_library(df)
            #formmatting
            df = pd.DataFrame([{
                "LID": lid,
                "Sequence Name": opts.name,
                "Temp": opts.temp,
                "Type1": opts.type1,
                "Type2": opts.type2,
                "Run": opts.run
            }])
            # Print nicely
            self.print_dataframe_rich(df)
        except Exception as err:
            print(f"Error adding sequence str({err})")
        except SystemExit:
            # argparse calls sys.exit on error, catch to prevent exit
            print()
            pass

    import argparse
    import os
    from colorama import Fore

    def handle_load(self, args):
        parser = argparse.ArgumentParser(prog="load", description="Load data (signal or basecall)")
        subparsers = parser.add_subparsers(dest="subcommand", help="Subcommands: signal or basecall")

        # --- Signal subcommand ---
        signal_parser = subparsers.add_parser("signal", help="Load nanopore signal-level data")
        signal_parser.add_argument("-l", "--lid", required=True, help="Library ID")
        signal_parser.add_argument("-p", "--path", required=True,
                                   help="Path to the signal file (Nanopolish tab separated txt file)")

        # --- Basecall subcommand ---
        basecall_parser = subparsers.add_parser("basecall", help="Load nanopore basecall data")
        basecall_parser.add_argument("-l", "--lid", required=True, help="Library ID")
        basecall_parser.add_argument("-p", "--path", required=True, help="Path to directory containing alignment data")
        basecall_parser.add_argument("--plot", dest="plot", action="store_true", help="Generate plots (default: True)")
        basecall_parser.add_argument("--no-plot", dest="plot", action="store_false", help="Disable plot generation")
        basecall_parser.set_defaults(plot=True)

        try:
            opts = parser.parse_args(args)
            check, df_lid = self.db_check(lid=opts.lid)

            # --- Check paths ---
            if opts.subcommand == "signal":
                if not os.path.isfile(opts.path):
                    print(Fore.RED + f"Signal file path does not exist: {opts.path}")
                    return

                if not check:
                    return
                else:
                    print(
                        Fore.CYAN + f"Loading signal data with: lid={opts.lid}, path={opts.path}")

                    # Select only current contig for upload from tx
                    cols = ['contig', 'position', 'reference_kmer', 'read_index',
                            'event_level_mean', 'event_length', 'event_stdv']
                    df = pd.read_csv(opts.path, sep='\t', usecols=cols)
                    contig = df_lid['contig'].unique()[0]
                    df = df[df['contig'] == contig.strip()]
                    df['type1'] = df_lid['type1'].unique()[0]
                    df['type2'] = df_lid['type2'].unique()[0]

                    if len(df) <= 0:
                        print(
                            Fore.RED + f"Sequence name '{contig}' not found in file path={opts.path}. "
                                        f"Please be sure sequence name entered matches sequence "
                                        f"name in FASTA file and in Nanopolish file.")
                        return
                    df['LID'] = int(opts.lid)
                    dbins.insert_signal(df)
                    print(Fore.CYAN + "Load Complete.")

            elif opts.subcommand == "basecall":
                if not os.path.isdir(opts.path):
                    print(Fore.RED + f"Basecall directory path does not exist: {opts.path}")
                    return

                if not check:
                    return
                else:
                    contig = df_lid['contig'].unique()[0]
                    print(
                        Fore.CYAN + f"Loading basecall data with: lid={opts.lid},path={opts.path}, plot={opts.plot}")

                    try:
                        run_basecall.get_modification(opts.lid, contig, opts.path,
                                                  df_lid['type2'].unique()[0].capitalize(), plot=opts.plot)
                    except Exception as err:
                        print(Fore.RED + 'Basecall Error: ' + str(err))
                        return

            else:
                print(Fore.RED + "Please specify either 'signal' or 'basecall' as subcommand.")

        except SystemExit:
            # argparse throws SystemExit on parse errors; we ignore to avoid exiting CLI
            print()
            pass

    @handle_cli_errors
    def do_seq(self, arg):
        args = arg.strip().split()
        if not args:
            print("Usage: seq [-list [library|unmodified|modified]] | [-add ...]")
            return

        subcmd_input = args[0]
        resolved_subcmd = self.resolve_subcommand(subcmd_input, self.seq_subcommands.keys())
        if not resolved_subcmd:
            return

        handler = self.seq_subcommands.get(resolved_subcmd)
        if handler:
            handler(args[1:])

    @handle_cli_errors
    def do_load(self, arg):
        args = arg.strip().split()
        if not args:
            print("Usage: load <basecall|signal> [parameters]")
            return

        # Remove optional dash prefix if present
        if args[0].startswith("-"):
            args[0] = args[0][1:]

        self.handle_load(args)

    @handle_cli_errors
    def do_predict(self, arg):
        parser = argparse.ArgumentParser(prog="predict", description="Predict RNA modification reactivity scores")
        parser.add_argument("-u", "--unmod_lid", required=True, help="Unmodified library ID. check seq -list umodified")
        parser.add_argument("-l", "--mod_lid", required=True, help="Modified library ID. check seq -list modified")
        parser.add_argument("-v", "--vienna", action="store_true", help="Enable Vienna base pairing prediction")

        try:
            opts = parser.parse_args(arg.strip().split())
            print(Fore.CYAN + f"Predicting with: unmod_lid={opts.unmod_lid}, lid={opts.mod_lid}, vienna={opts.vienna}")
            check, df = self.db_check(lid=opts.mod_lid, unmod_lid=opts.unmod_lid)
            if not check:
                return
            else:
                predict.run_predict(unmod_lids=opts.unmod_lid, lids=opts.mod_lid, continue_reads=False,
                                    vienna=opts.vienna)
        except SystemExit:
            print()
            pass

    @handle_cli_errors
    def do_create_landscape(self, arg):
        parser = argparse.ArgumentParser(prog="create_landscape", description="Create reactivity or mutational landscape")
        parser.add_argument("-l", "--mod_lid", required=True, help="Modified library ID")
        parser.add_argument("-u", "--unmod_lid", required=True, help="Unmodified library ID")
        parser.add_argument(
            "-o",
            "--optimize",
            action="store_true",
            help="Optimize clusters numbers and plot (default: False)"
        )
        try:
            opts = parser.parse_args(arg.strip().split())
            print(Fore.GREEN + f"Creating landscape with: lid={opts.mod_lid}, unmod_lid={opts.unmod_lid}, optimize_clusters={opts.optimize}")
            check, df = self.db_check(lid=opts.mod_lid, unmod_lid=opts.unmod_lid)
            if not check:
                return
            else:
                images = landscape.run_landscape(unmod_lid=opts.unmod_lid, lid=opts.mod_lid,
                                                 optimize_clusters=opts.optimize)
        except SystemExit:
            print()
            pass

    @handle_cli_errors
    def do_exit(self, arg):
        print(Fore.RED + "Exiting DashingTurtle.")
        return True

    @handle_cli_errors
    def do_man(self, arg):
        if not arg:
            print("Usage: man <command>")
        else:
            self.print_manual(arg.strip())

    def print_manual(self, command):
        manuals = {
            "seq": (
                "Usage: seq <subcommand> [options]\n\n"
                "Subcommands:\n"
                "  -list [library|unmodified|modified]\n"
                "      List sequences stored. Default is 'library'.\n\n"
                "Parameters for -list:\n"
                "  library           List all library sequences (default)\n"
                "  unmodified        List unmodified sequences\n"
                "  modified          List modified sequences\n\n"
                "  -add <sequence> <secondary> <experiment> <name> <temp> <type1> <type2> <run>\n"
                "      Add a new sequence.\n\n"
                "Parameters for -add:\n"
                "  -s, --sequence       RNA sequence (required)\n"
                "  -sec, --secondary    Secondary structure if experimentally confirmed. (optional)\n"
                "  -e, --experiment     Control experiment if using a verified secondary structure (optional)\n"
                "  -n, --name           Sequence name. Must be identical to library entry in FASTA. (required)\n"
                "  -t, --temp           Temperature in °C Default is 37.(optional)\n"
                "  -t1, --type1         Unmodified Condition (required)\n"
                "  -t2, --type2         Modified Condition (same as t1 if unmodified) (required)\n"
                "  -r, --run            Run ID. Default 1. Multiple identical experiments should be assigned different runs.\n"
            ),
            "load": (
                "Usage: load <subcommand> [parameters]\n\n"
                "Subcommands:\n"
                " basecall -l <lid> -p <path> -plot <FlAG>\n"
                "      Load basecalling alignment data.\n\n"
                "   Parameters for basecall (can be used in any order):\n"
                "   -l, --lid            Library ID\n"
                "   -p, --path           Path to the basecall file\n"
                "   -plot                Plot alignments flag \n\n"
                " signal -l <lid> -p <path>\n"
                "      Load nanopore signal-level data.\n\n"
                "   Parameters for signal (can be used in any order):\n"
                "   -l, --lid          Library ID of sequence (seq -l)\n"
                "   -p, --path         Path to the signal file (Nanopolish tab separated txt file)\n"
            ),
            "predict": (
                "Usage: predict -l <lid> -u <unmodified lid> -v <FLAG>\n\n"
                "Predict RNA modification reactivity scores. Validates sequence matches.\n\n"
                "   Parameters (can be used in any order):\n"
                "   -u, --unmod_lid     Umodified library ID. Check seq -list unmodified\n"
                "   -l, --mod_lid       Modified library ID. Check seq -list modified\n"
                "   -v, --vienna        Use Vienna base pairing probability in prediction FLAG\n"
            ),
            "create_landscape": (
                "Usage: create_landscape -l <lid> -u <unmodified lid> -o <FLAG>\n\n"
                "Generate reactivity or mutational landscape visualization.\n\n"
                "Parameters (can be used in any order):\n"
                "  -l, --lid                Modified library ID\n"
                "  -u, --unmod_lid          Unmodified library ID\n"
                "  -o, --optimize           Optimize clusters for visualization using Elbow Method (T or F)\n"
            ),
            "config": (
                "Usage: config\n\n"
                "Show currently loaded configuration values and defaults."
            ),
            "exit": (
                "Usage: exit\n\n"
                "Exit the CLI session."
            )
        }

        manual = manuals.get(command)
        if manual:
            print(f"\nManual for '{command}':\n")
            print(manual)
        else:
            print(f"No manual entry for '{command}'. Available: {', '.join(manuals.keys())}")

    def resolve_subcommand(self, input_cmd, valid_cmds):
        matches = [cmd for cmd in valid_cmds if cmd.startswith(input_cmd)]
        if len(matches) == 1:
            return matches[0]
        elif len(matches) > 1:
            print(Fore.RED + f"*** Ambiguous subcommand '{input_cmd}'. Matches: {matches}")
        else:
            print(Fore.RED + f"*** Unknown subcommand '{input_cmd}'")
        return None

    @handle_cli_errors
    def do_help(self, arg):
        if not arg:
            print("""
DashingTurtle
--------------

Sequence Management:
  seq -list [library|unmodified|modified]  
      List stored sequences. Default type is 'library'.

  seq -add -s <sequence> -sec <secondary> -e <experiment> -n <name> -t <temp> -t1 <type1> -t2 <type2> -r <run>
      Add a new sequence.

Loading Data:
  load basecall -l <lid> -p <path> -plot <FLAG>
      Load basecall alignment data.

  load signal -l <lid>  -p <path>
      Load signal-level data.

Prediction & Landscape:
  predict -u <unmod_lids> -l <lids> -v <FLAG>
      Run modification reactivity prediction.

  create_landscape -l <lid> -u <unmod_lid> -o <FLAG>
      Generate landscape visualization.

General:
  config         Show configuration
  man <command>  Detailed manual entry for a command
  exit           Quit the interface
""")
        else:
            self.do_man(arg)

    def completenames(self, text, *ignored):
        return [c for c in self.commands if c.startswith(text)]

    def complete_seq(self, text, *args):
        return [sc for sc in self.seq_subcommands if sc.startswith(text)]

    def complete_man(self, text, *args):
        return [c for c in self.commands if c.startswith(text)]


def main():
    print("Launching CLI...")
    DashingTurtleCLI().cmdloop()

if __name__ == '__main__':
    main()
