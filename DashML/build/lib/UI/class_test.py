import pexpect
import sys
import os
import re
import pytest

# pytest -v class_test.py

# Utility to strip ANSI color codes (if needed in assertions)
def strip_ansi(text):
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    return ansi_escape.sub('', text)


# Regex for prompt (adjust if your prompt changes)
PROMPT_REGEX = r'\x1b\[33m🐢> \x1b\[0m'


@pytest.fixture(scope="module")
def cli_process():
    """Start the CLI process and yield it, then close at the end."""
    script_path = os.path.join(os.path.dirname(__file__), "DT_CLI.py")
    script_path = os.path.abspath(script_path)

    child = pexpect.spawn(sys.executable, [script_path], encoding='utf-8', timeout=20)

    # Enable for debugging
    # child.logfile = sys.stdout

    child.expect("Welcome to DashingTurtle!")
    child.expect(PROMPT_REGEX)
    yield child

    if not child.closed:
        child.sendline("exit")
        child.expect("Exiting DashingTurtle.")
        child.close()


def send_and_check(child, command, expected):
    """Helper to send a command, expect specific output, and wait for prompt."""
    child.sendline(command)
    child.expect(expected)
    # After expected output, always expect prompt to ensure CLI is ready
    child.expect(PROMPT_REGEX)
    # Get complete output before prompt
    output = strip_ansi(child.before + child.after)
    assert expected in output


def test_help(cli_process):
    send_and_check(cli_process, "help", "exit           Quit the interface")


def test_config(cli_process):
    send_and_check(cli_process, "config", "Defaults from configuration")


def test_man_seq(cli_process):
    send_and_check(cli_process, "man seq", "Manual for 'seq'")


def test_seq_no_args(cli_process):
    send_and_check(cli_process, "seq", "Available subcommands for seq")


def test_seq_unknown(cli_process):
    send_and_check(cli_process, "seq -unknown", "Unknown subcommand")


def test_predict(cli_process):
    send_and_check(cli_process, "predict -u 1,2 -l 3,4 -v T", "Predicting with: unmod_lids=1,2, lids=3,4, vienna=T")


def test_create_landscape(cli_process):
    send_and_check(cli_process, "create_landscape -l 3 -u 2 -o T",
                   "Creating landscape with: lid=3, unmod_lid=2, optimize_clusters=T")


@pytest.mark.parametrize(
    "command, expected",
    [
        (
                "load basecall -l 1 -n my_seq -t modA -p /path/to/alignment --plot",
                "Loading basecall with lid=1, name=my_seq, type=modA, path=/path/to/alignment, plot=True"
        ),
        (
                "load signal -m 2 -u 1 -p /path/to/data.tsv",
                "Loading signal with mod_lid=2, unmod_lid=1, path=/path/to/data.tsv"
        ),
    ]
)
def test_load_commands(cli_process, command, expected):
    send_and_check(cli_process, command, expected)


def test_exit(cli_process):
    cli_process.sendline("exit")
    cli_process.expect("Exiting DashingTurtle.")
    output = strip_ansi(cli_process.before + cli_process.after)
    assert "Exiting DashingTurtle." in output
    cli_process.close()
