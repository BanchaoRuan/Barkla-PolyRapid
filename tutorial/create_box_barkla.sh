#!/bin/bash -l
#SBATCH --job-name=create_box
#SBATCH --partition=troisi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=5-00:00:00
#SBATCH --output=create_box_%j.out

set -eo pipefail

if [[ ${1:-} != --run ]]; then
    if (( $# > 1 )); then
        echo 'Usage: bash create_box_barkla.sh [Output]' >&2
        exit 2
    fi
    script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
    cd -P -- "${1:-Output}"
    shopt -s nullglob
    materials=(W01*/)
    if (( ${#materials[@]} == 0 )); then
        echo "No W01 material directories in $PWD." >&2
        exit 1
    fi
    # Submit once from the real Output directory, using the original tutorial on shared storage.
    exec sbatch --chdir="$PWD" "$script_dir/create_box_barkla.sh" --run "$script_dir/create_box.py"
fi

driver=${2:?Missing path to create_box.py}
export TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" LMOD_IGNORE_CACHE=yes
module load gromacs/2024.3-gcc14.2.0
plain_gmx=$(type -P gmx)
module swap gromacs/2024.3-gcc14.2.0 gromacs/2024.3-openmpi5.0.8-gcc14.2.0
export PATH="${plain_gmx%/*}:$PATH"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:?Submit this job with Slurm}" FI_PROVIDER=opx
gmx --version
gmx_mpi --version

# Hold the same file locks as the force-field step while the tutorial edits the real topologies.
srun --mpi=none --ntasks=1 --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores \
    python3 -Bu - "$driver" <<'PY'
from contextlib import ExitStack
from pathlib import Path
import fcntl
import subprocess
import sys

materials = sorted(path for path in Path.cwd().glob('W01*') if path.is_dir())
if not materials:
    sys.exit('No W01 material directories at job start.')

with ExitStack() as locks:
    for material in materials:
        lock = locks.enter_context((material / '.dl_field.lock').open('a'))
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            sys.exit(f'Another job is using {material}.')
        if (any(material.glob('*.oldgro')) or any(p.name != 'gromacs.gro' for p in material.glob('*.gro'))
                or len(list(material.glob('*.top'))) > 1):
            sys.exit(f'Check {material} for a previous box run or ambiguous input files before submitting again.')

    result = subprocess.run([sys.executable, '-Bu', sys.argv[1]])
    if result.returncode:
        sys.exit(result.returncode)
    for material in materials:
        # The tutorial renames enmin.gro to enmin.oldgro only after the insertion command succeeds.
        if (material / 'enmin.oldgro').is_file():
            box = material / 'enmin_box.gro'
            if not box.is_file() or not box.stat().st_size:
                sys.exit(f'Missing or empty packed coordinates: {box}')
PY
