#!/bin/bash -l
#SBATCH --job-name=assign_forcefield
#SBATCH --partition=troisi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=assign_forcefield_%j.out

set -eo pipefail

# Slurm receives absolute paths because its copy of this script lives in a spool directory.
if [[ "${1:-}" == "--run" ]]; then
    module load dl_field/4.13
    exec python3 -u "$2" "$3"
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
output_dir=$(cd -- "${1:-Output}" && pwd -P)
exec sbatch --chdir="$output_dir" "$script_dir/assign_forcefield_barkla.sh" \
    --run "$script_dir/assign_forcefield.py" "$output_dir"
