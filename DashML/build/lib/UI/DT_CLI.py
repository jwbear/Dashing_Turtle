import sys
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

class DashingTurtleCLI(cmd.Cmd):
    intro = Fore.GREEN + "Welcome to DashingTurtle! Type 'help' or '?' to list commands.\n"
    prompt = Fore.YELLOW + "🐢> " + Style.RESET_ALL
    VALID_LIST_TYPES = {'library', 'unmodified', 'modified'}

    def db_check(self, lid, unmod_lid=None):
        """
        Check if a library ID exists in the database.
        """
        check = False
        if unmod_lid is None:
            df = dbsel.select_librarybyid(int(lid))
            if (len(df) > 0):
                check = True
            else:
                print(Fore.RED + f"LID {lid} not found in database. Please run seq -list to verify.")
        else:
            df = dbsel.select_library_full()
            df = df[(df['ID'] == int(lid)) | (df['ID']==int(unmod_lid))]
           #print(df)
            if len(df) != 2:
                print(Fore.RED + f"LID {lid} or {unmod_lid} not found in database or are the same sequence. Please run seq -list to verify.")
                return check, df
            #check types
            dt = df.loc[(df['ID'] == int(lid))]
            contig1 = dt['contig'].unique()[0]
            if dt['type1'].unique()[0] == dt['type2'].unique()[0]:
                print(Fore.RED + f"LID {lid} is listed as NOT modified. Please run seq -list to verify.")
                return check, df
            dt = df.loc[(df['ID'] == int(unmod_lid))]
            contig2 = dt['contig'].unique()[0]
            if dt['type1'].unique()[0] != dt['type2'].unique()[0]:
                print(Fore.RED + f"LID {lid} is listed as modified. Please run seq -list to verify.")
                return check, df
            if contig1 != contig2:
                print(Fore.RED + f"Contig {contig1} and {contig2} are not the same sequence. Please run seq -list to verify.")
                return check, df
            check = True

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
        parser.add_argument("-r", "--run", default=self.get_config_default("defaults", "run", "1"), help="Run ID")

        try:
            opts = parser.parse_args(args)
            print(
                Fore.CYAN + f"Adding sequence with: {opts.sequence}, {opts.secondary}, {opts.experiment}, {opts.name}, {opts.temp}, {opts.type1}, {opts.type2}, {opts.run}")

            sequence = opts.sequence
            secondary = opts.secondary

            # Validate secondary structure length
            if secondary:
                # Create a console locally
                console = Console()

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
            # test
            # seq -add -s GGAUCGAUCG -sec .......... -e EXP001 -n TestSeq1 -t 37 -t1 TypeA -t2 TypeB -r 1
            # seq -add -s GGAUCGAUCG -sec ..........-- -e EXP001 -n BadSeq -t 37 -t1 TypeA -t2 TypeB -r 1
            # seq -add -s AGCUAGCUA -n TestName -t1 TypeA -t2 TypeB
            # seq -add -s AGCUAGCUA -sec ((((....)) -n TestName -t1 TypeA -t2 TypeB
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

    def do_load(self, arg):
        args = arg.strip().split()
        if not args:
            print("Usage: load <basecall|signal> [parameters]")
            return

        # Remove optional dash prefix if present
        if args[0].startswith("-"):
            args[0] = args[0][1:]

        self.handle_load(args)

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

    def do_exit(self, arg):
        print(Fore.RED + "Exiting DashingTurtle.")
        return True

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
                "  -sec, --secondary    Secondary structure (optional)\n"
                "  -e, --experiment     Experiment ID (optional)\n"
                "  -n, --name           Sequence name (required)\n"
                "  -t, --temp           Temperature in °C (optional)\n"
                "  -t1, --type1         First type label (required)\n"
                "  -t2, --type2         Second type label (required)\n"
                "  -r, --run            Run ID (optional)\n"
            ),
            "load": (
                "Usage: load <subcommand> [parameters]\n\n"
                "Subcommands:\n"
                "  basecall -l <lid> -c <contig> -p <path> -m <modification> -plot <T|F>\n"
                "      Load basecalling alignment data.\n\n"
                "Parameters for basecall (can be used in any order):\n"
                "  -l, --lid            Library ID\n"
                "  -p, --path           Path to the basecall file\n"
                "  -plot                Plot alignments (T or F). Default is plot \n\n"
                "  signal -l <lid> -p <path>\n"
                "      Load nanopore signal-level data.\n\n"
                "Parameters for signal (can be used in any order):\n"
                "  -l, --lid          Library ID of sequence (seq -l)\n"
                "  -p, --path         Path to the signal file (Nanopolish tab separated txt file)\n"
            ),
            "predict": (
                "Usage: predict [parameters]\n\n"
                "Predict RNA modification reactivity scores. Validates sequence matches.\n\n"
                "Parameters (can be used in any order):\n"
                "  -u, --unmod_lid     Umodified library ID. Check seq -list unmodified\n"
                "  -l, --mod_lid       Modified library ID. Check seq -list modified\n"
                "  -v, --vienna        Use Vienna base pairing in prediction FLAG\n"
            ),
            "create_landscape": (
                "Usage: create_landscape [parameters]\n\n"
                "Generate reactivity or mutational landscape visualization.\n\n"
                "Parameters (can be used in any order):\n"
                "  -l, --lid                Modified library ID\n"
                "  -u, --unmod_lid          Unmodified library ID\n"
                "  -o, --optimize           Optimize clusters for visualization (T or F)\n"
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
