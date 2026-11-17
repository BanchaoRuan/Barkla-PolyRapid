import argparse
import fcntl
import os
from pathlib import Path
import re
import subprocess


def file_stamp(path):
    stat = path.stat()
    return stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def assign_forcefield(directory, structure, template):
    with (directory / ".dl_field.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if any(character.isspace() for character in structure.name):
            raise ValueError("The XYZ filename must not contain whitespace.")
        with structure.open() as stream:
            atom_count = int(stream.readline())
        if atom_count <= 0:
            raise ValueError("The XYZ atom count must be positive.")

        previous = {
            path: file_stamp(path)
            for path in directory.glob("gromacs*")
            if path.is_file()
        }
        control = template.replace("polymer_optimized.xyz", structure.name, 1)
        (directory / "dl_field.control").write_text(control, encoding="utf-8")
        log = directory / "dl_field_run.log"
        with log.open("w", encoding="utf-8") as stream:
            result = subprocess.run(
                ["dl_field"], cwd=directory, stdin=subprocess.DEVNULL,
                stdout=stream, stderr=subprocess.STDOUT, check=False,
            )
        if result.returncode != 0:
            raise RuntimeError(f"DL_FIELD exited with code {result.returncode}; see {log.name}.")
        if "Program executed successfully." not in log.read_text(errors="replace"):
            raise RuntimeError(f"DL_FIELD did not report success; see {log.name}.")

        topology = directory / "gromacs.top"
        coordinates = directory / "gromacs.gro"
        includes = re.findall(
            r'^\s*#include\s+"([^"]+\.itp)"',
            topology.read_text(), flags=re.MULTILINE,
        )
        if not includes:
            raise RuntimeError("gromacs.top does not include a molecular topology.")
        for path in [topology, coordinates, *(directory / name for name in includes)]:
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError(f"Missing or empty output: {path.name}")
            if file_stamp(path) == previous.get(path):
                raise RuntimeError(f"Output was not refreshed by this run: {path.name}")
        with coordinates.open() as stream:
            stream.readline()
            output_atoms = int(stream.readline())
        if output_atoms != atom_count:
            raise RuntimeError(f"Atom count changed from {atom_count} to {output_atoms}.")


def main():
    parser = argparse.ArgumentParser(description="Assign PolyRapid force fields on Barkla.")
    parser.add_argument("output", nargs="?", type=Path, default=Path("Output"),
                        help="Directory containing one subdirectory per material (default: Output).")
    args = parser.parse_args()
    root = args.output.resolve()
    if not root.is_dir():
        parser.error(f"Material directory does not exist: {root}")
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Submit with assign_forcefield_barkla.sh to run inside Slurm.")
    template = Path(__file__).with_name("polymer.control").read_text(encoding="utf-8")
    failures = 0
    completed = 0
    for directory in sorted(path for path in root.iterdir() if path.is_dir()):
        structures = sorted(directory.glob("*.xyz"))
        if not structures:
            print(f"SKIP {directory.name}: no XYZ file", flush=True)
            continue
        print(f"RUN {directory.name}", flush=True)
        try:
            if len(structures) != 1:
                raise ValueError("Expected exactly one XYZ file in the material directory.")
            assign_forcefield(directory, structures[0], template)
        except (OSError, ValueError, RuntimeError) as error:
            print(f"FAIL {directory.name}: {error}", flush=True)
            failures += 1
        else:
            print(f"OK {directory.name}", flush=True)
            completed += 1
    print(f"Completed: {completed}; failed: {failures}", flush=True)
    return int(failures > 0 or completed == 0)


if __name__ == "__main__":
    raise SystemExit(main())
