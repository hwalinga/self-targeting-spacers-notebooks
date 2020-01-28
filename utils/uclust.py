#!/usr/bin/env python3

import shutil
import subprocess
import tempfile
from functools import partial
from os import path


def create_fasta(sequences):
    start = 1
    for seq in sequences:
        yield ">%d\n" % start
        start += 1
        yield "%s\n" % seq


def execute(cmd: str):
    exe = partial(
        subprocess.run,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    command = cmd.split()[0]
    if not shutil.which(command):
        raise OSError(f"{command} not found")
    return exe(cmd)


def uclust(sequences, match_id=0.9):  # list of sequences
    with tempfile.TemporaryDirectory() as tmpdirname:

        seq = path.join(tmpdirname, "seq")
        with open(seq, "w") as seq_file:
            seq_file.writelines(create_fasta(sequences))

        seq_sorted = path.join(tmpdirname, "seqsorted")
        execute("uclust --sort %s --output %s" % (seq, seq_sorted))

        seq_clustered = path.join(tmpdirname, "seqclustered")
        execute(
            "uclust --input %s --uc %s --id %.2f"
            % (seq_sorted, seq_clustered, match_id)
        )

        with open(seq_clustered) as seq_clustered_file:
            for line in seq_clustered_file:
                if line.startswith("C"):
                    yield int(line.split()[2])
