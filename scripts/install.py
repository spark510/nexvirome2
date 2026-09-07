#!/usr/bin/env python3
"""Install the checkout on Linux using only conda-forge and bioconda."""

import argparse
import os
from pathlib import Path
import platform
import shutil
import subprocess


CHANNEL_SETTINGS = {
    "CONDA_CHANNELS": "conda-forge,bioconda",
    "CONDA_ALLOWLIST_CHANNELS": "conda-forge,bioconda",
    "CONDA_CHANNEL_PRIORITY": "strict",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", action="store_true", help="Install only the core assembly environment.")
    parser.add_argument("--dry-run", action="store_true", help="Solve without creating an environment.")
    args = parser.parse_args()
    if platform.system() != "Linux" or platform.machine().lower() not in {"x86_64", "amd64"}:
        parser.error("The external bioinformatics tools require Linux x86_64 (including WSL2).")
    conda = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if not conda:
        parser.error("Conda is not available. Initialize Conda before running this installer.")
    root = Path(__file__).resolve().parents[1]
    spec = root / ("environment.yml" if args.base else "environments/extended.yaml")
    name = "nexvirome2" if args.base else "nexvirome2-extended"
    env = dict(os.environ, **CHANNEL_SETTINGS)
    command = [conda, "env", "create", "--no-default-packages", "--file", str(spec)]
    if args.dry_run:
        command.append("--dry-run")
    subprocess.run(command, cwd=root, env=env, check=True)
    if args.dry_run:
        return
    # Persist priority only. A persistent allowlist can break unrelated Conda
    # operations when a user's global channels include defaults (lists merge).
    activation_env = dict(env)
    activation_env.pop("CONDA_ALLOWLIST_CHANNELS", None)
    subprocess.run(
        [conda, "env", "config", "vars", "set", "--name", name]
        + ["CONDA_CHANNEL_PRIORITY=strict"],
        cwd=root,
        env=activation_env,
        check=True,
    )
    print(f"Installation complete. Run: conda activate {name}")
    print("Then run: nexvirome2 doctor" + ("" if args.base else " --extended"))
    print("For later package additions: conda install --override-channels -c conda-forge -c bioconda PACKAGE")


if __name__ == "__main__":
    main()
