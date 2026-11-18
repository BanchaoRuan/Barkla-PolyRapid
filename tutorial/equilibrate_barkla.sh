#!/bin/bash -l
#SBATCH --job-name=equilibrate
#SBATCH --partition=troisi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=42
#SBATCH --time=5-00:00:00
#SBATCH --output=equilibrate_%j.out

# From a batch directory: bash ../tutorial/equilibrate_barkla.sh [--maxcycles N] [Output | Output/PID]
# PolyRapid source: https://github.com/paolacarbone/PolyRapid/tree/879145fbe91997bfb033a61089aaa32c3fed7979/tutorial
set -eo pipefail
shopt -s nullglob dotglob

die() { printf '%s\n' "$*" >&2; exit 1; }

check_box() {
    # The upstream driver chooses the first *.gro and *.top in its working directory.
    local gro=( *.gro ) top=( *.top ) steps=( Step[0-9]* )
    [[ ${#gro[@]} == 1 && -s enmin_box.gro && ${#top[@]} == 1 &&
       -s gromacs.top && -s gromacs1.itp ]] ||
        die "$PWD: expected one packed enmin_box.gro, one gromacs.top and its gromacs1.itp."
    [[ ! -e Enmin && ${#steps[@]} == 0 && ! -e rdf_groups.ndx &&
       ! -e All_properties.csv && ! -e All_convergence.csv ]] ||
        die "$PWD: existing equilibration output. Inspect the previous run before restarting."
}

check_cycles() {
    [[ $1 =~ ^[1-9][0-9]*$ ]] || die '--maxcycles needs a positive integer.'
}

# The job receives the tutorial path explicitly because Slurm spools this script.
# Direct single-material sbatch with the tutorial directory as its argument also works.
if [[ ${1:-} == --run || ( -n ${SLURM_JOB_ID:-} && -n ${1:-} && -f $1/equilibrate.py ) ]]; then
    if [[ $1 == --run ]]; then shift; fi
    (( $# == 1 || $# == 2 )) || die 'Expected a tutorial directory and optional cycle limit.'
    cpus=${SLURM_CPUS_PER_TASK:?Submit this computation through Slurm}
    tutorial=$1
    maxcycles=${2:-6}
    check_cycles "$maxcycles"
    python_bin=${POLYRAPID_PYTHON:-"$HOME/.conda/envs/polyrapid/bin/python"}
    [[ -x $python_bin && -s $tutorial/equilibrate.py ]] || die 'Missing polyrapid Python or equilibrate.py.'

    [[ -s enmin_box.gro ]] || die "$PWD: missing enmin_box.gro."
    exec 9>>.dl_field.lock
    flock -n 9 || die "$PWD: another workflow job is using this material."
    check_box

    export TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" LMOD_IGNORE_CACHE=yes
    module purge
    module load gromacs/2024.3-gcc14.2.0
    type -P gmx >/dev/null || die 'The single-node gmx executable is missing.'
    if type -P gmx_mpi >/dev/null; then
        die 'gmx_mpi is still on PATH. Check the single-node GROMACS environment.'
    fi

    for mdp in enmin.mdp npt.mdp anneal.mdp; do
        [[ -s $tutorial/$mdp ]] || die "Missing template $tutorial/$mdp."
        if [[ -e $mdp ]]; then
            # create_box.py writes the same enmin.mdp without the final newline.
            [[ "$(<"$mdp")" == "$(<"$tutorial/$mdp")" ]] || die "$PWD/$mdp differs from the tutorial."
        else
            cp -- "$tutorial/$mdp" "$mdp"
        fi
    done

    export OMP_NUM_THREADS=7 FI_PROVIDER=opx
    printf 'Material %s | CPUs %s | OpenMP threads per rank %s | max cycles %s\n' \
        "${PWD##*/}" "$cpus" "$OMP_NUM_THREADS" "$maxcycles"
    gmx --version
    srun --mpi=none --ntasks=1 --cpus-per-task="$cpus" --cpu-bind=cores \
        "$python_bin" -I -B -u "$tutorial/equilibrate.py" --deltardf 0.02 --maxcycles "$maxcycles"
    exit
fi

maxcycles=6
if [[ ${1:-} == --maxcycles ]]; then
    maxcycles=${2:?Supply the number of cycles}
    check_cycles "$maxcycles"
    shift 2
fi
(( $# <= 1 )) || die 'Usage: bash equilibrate_barkla.sh [--maxcycles N] [Output | Output/PID]'
tutorial=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
for source in equilibrate.py enmin.mdp npt.mdp anneal.mdp; do
    [[ -s $tutorial/$source ]] || die "Missing $tutorial/$source."
done
target=$(cd -- "${1:-Output}" && pwd -P)
if [[ -f $target/enmin_box.gro || -f $target/gromacs.top ]]; then
    materials=("$target")
    batch_dir=$(dirname -- "$target")
else
    batch_dir=$target
    materials=("$target"/W01*/)
    (( ${#materials[@]} )) || die "No W01 material directories under $target."
fi

# Serialize only the submission pass. The compute jobs run independently.
exec 8>>"$batch_dir/.equilibrate-submit.lock"
flock 8 || die "Could not lock submissions in $batch_dir."
export TMPDIR="$batch_dir" TMP="$batch_dir" TEMP="$batch_dir"
active_dirs=$(squeue --me --name=equilibrate --noheader --format=%Z)

submit_material() (
    cd -P -- "$1" || exit 1
    if grep -Fxq -- "$PWD" <<< "$active_dirs"; then
        printf 'Skipping %s (already in Slurm).\n' "$PWD"
        exit 0
    fi
    if [[ -s All_properties.csv && -s All_convergence.csv ]]; then
        printf 'Skipping %s (result tables exist; see the previous log for convergence).\n' "$PWD"
        exit 0
    fi
    check_box
    receipt=$(sbatch --parsable --chdir="$PWD" "$tutorial/equilibrate_barkla.sh" \
        --run "$tutorial" "$maxcycles") || die "Submission failed for $PWD."
    [[ $receipt =~ ^[0-9]+(\;[[:alnum:]_.-]+)?$ ]] || die "$PWD: unexpected receipt '$receipt'. Check squeue."
    printf '%s -> Slurm %s\n' "${PWD##*/}" "$receipt"
)

status=0
for material in "${materials[@]}"; do
    submit_material "$material" || status=1
done
exit "$status"
