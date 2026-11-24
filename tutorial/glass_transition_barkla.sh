#!/bin/bash -l
#SBATCH --job-name=glass_transition
#SBATCH --partition=troisi
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=7
#SBATCH --time=5-00:00:00
#SBATCH --output=glass_transition_%j.out

# bash glass_transition_barkla.sh [--list] [Output | Output/PID ...]
# --list only reports eligibility.
set -eo pipefail
shopt -s nullglob dotglob

die() { printf '%s\n' "$*" >&2; exit 1; }

equil_job() {
    local file job latest=0
    for file in "$1"/equilibrate_*.out; do
        job=${file##*/equilibrate_}
        job=${job%.out}
        if [[ $job =~ ^[1-9][0-9]*$ ]] && (( job > latest )); then latest=$job; fi
    done
    printf '%s\n' "$latest"
}

# Print the one coordinate file that the upstream driver will actually select.
final_coordinates() (
    cd -P -- "$1" || exit 1
    grep -Fxq 'Box is equilibrated' "equilibrate_$2.out" || die 'no convergence message in the latest equilibration log'
    [[ -s All_properties.csv && -s All_convergence.csv ]] || die 'missing equilibration result tables'
    last=$(awk -F, 'NR>1 && NF>1 {step=$1} END {print step}' All_properties.csv)
    pair=$(awk -F, 'NR>1 && NF>1 {pair=$1} END {print pair}' All_convergence.csv)
    [[ $last =~ ^[1-9][0-9]*$ && $pair == "$((last - 1))->$last" ]] || die 'inconsistent final steps in the equilibration tables'
    for step in Step*/; do
        name=${step%/}
        if [[ $name =~ ^Step([1-9][0-9]*)$ ]]; then
            (( BASH_REMATCH[1] <= last )) || die "$name is newer than the final equilibration tables"
        elif [[ -d $step/anneal ]]; then
            die "ambiguous annealing directory $name"
        fi
    done
    gro=( "Step$last/anneal/"*.gro )
    top=( *.top )
    expected="Step$last/anneal/step${last}_anneal.gro"
    [[ ${#gro[@]} == 1 && -s $expected ]] || die "expected a single final coordinate file at $expected"
    [[ ${#top[@]} == 1 && -s gromacs.top && -s gromacs1.itp ]] || die 'expected gromacs.top and its gromacs1.itp'
    if [[ -e glasstransition.mdp || -L glasstransition.mdp ]]; then
        cmp -s glasstransition.mdp "$tutorial/glasstransition.mdp" || die 'existing glasstransition.mdp differs from the tutorial'
    fi
    printf '%s/%s\n' "$PWD" "$expected"
)

if [[ ${1:-} == --run ]]; then
    (( $# == 3 )) || die 'Internal worker arguments are missing.'
    tutorial=$2
    previous_job=$3
    material=$(pwd -P)
    : "${SLURM_JOB_ID:?Submit through Slurm}"
    ranks=${SLURM_NTASKS:?Missing Slurm task count}
    threads=${SLURM_CPUS_PER_TASK:?Missing Slurm CPU allocation}
    python_bin=${POLYRAPID_PYTHON:-"$HOME/.conda/envs/polyrapid/bin/python"}
    [[ -x $python_bin ]] || die "Missing Python: $python_bin"
    exec 9>>.dl_field.lock
    flock -n 9 || die "$material: another workflow is using this material."
    [[ $(equil_job "$material") == "$previous_job" ]] || die 'The equilibration run changed after submission.'
    coordinates=$(final_coordinates "$material" "$previous_job") || exit 1
    if [[ -e GlassTransition || -L GlassTransition ]]; then
        die 'GlassTransition already exists; inspect the previous run.'
    fi

    export TMPDIR="${SLURM_TMPDIR:-$material}" LMOD_IGNORE_CACHE=yes
    export TMP="$TMPDIR" TEMP="$TMPDIR" MPLBACKEND=Agg MPLCONFIGDIR="$PWD/.mplconfig"
    module purge
    module load gromacs/2024.3-gcc14.2.0
    plain_gmx=$(type -P gmx) || die 'Missing single-node gmx.'
    module swap gromacs/2024.3-gcc14.2.0 gromacs/2024.3-openmpi5.0.8-gcc14.2.0
    export PATH="${plain_gmx%/*}:$PATH"
    type -P gmx_mpi >/dev/null || die 'Missing MPI GROMACS.'
    real_mpirun=$(type -P mpirun) || die 'Missing Open MPI.'
    export FI_PROVIDER=opx OMP_NUM_THREADS="$threads" OMP_PLACES=cores OMP_PROC_BIND=close
    # The upstream driver calls bare mpirun. Supply binding options on the actual launch command.
    launcher_dir=$(mktemp -d "$TMPDIR/.polyrapid-mpi.XXXXXXXX")
    trap 'rm -f -- "$launcher_dir/mpirun"; rmdir -- "$launcher_dir"' EXIT
    cat > "$launcher_dir/mpirun" <<'LAUNCH'
#!/bin/sh
exec "$POLYRAPID_MPIRUN" -n "$SLURM_NTASKS" --map-by "slot:PE=$SLURM_CPUS_PER_TASK" --bind-to core "$@"
LAUNCH
    chmod 700 "$launcher_dir/mpirun"
    export POLYRAPID_MPIRUN="$real_mpirun" PATH="$launcher_dir:$PATH"
    # The upstream energy command opens its XVG twice. A positive backup count keeps stdout separate.
    export GMX_MAXBACKUP=99
    printf 'Material %s | equilibration %s | %s | MPI ranks %s x threads %s\n' "$material" "$previous_job" "$coordinates" "$ranks" "$threads"
    if [[ ! -e glasstransition.mdp ]]; then cp -- "$tutorial/glasstransition.mdp" .; fi
    "$python_bin" -I -B -u "$tutorial/glass_transition.py"
    "$python_bin" -I -B -u "$tutorial/glass_transition_weighted.py" "$material"
    exit
fi

mode=run
case ${1:-} in
    --list) mode=list; shift ;;
    --*) die "Unknown option: $1" ;;
esac
tutorial=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
for source in glass_transition.py glasstransition.mdp glass_transition_weighted.py; do
    [[ -s $tutorial/$source ]] || die "Missing $tutorial/$source."
done
if (( $# == 0 )); then set -- Output; fi
materials=()
for arg in "$@"; do
    path=$(cd -P -- "$arg" && pwd -P) || die "Cannot open $arg."
    if [[ ${path##*/} == Output ]]; then
        children=( "$path"/W01*/ )
        (( ${#children[@]} )) || die "No W01 materials under $path."
        for child in "${children[@]}"; do materials+=("${child%/}"); done
    else
        materials+=("$path")
    fi
done

submit_batch() (
    batch=$1
    if [[ $mode != list ]]; then
        exec 8>>"$batch/.tg-submit.lock" || exit 1
        flock 8 || die "Could not lock submissions in $batch."
    fi
    active=$(squeue --me --name=equilibrate,glass_transition --noheader --format=%Z) || die 'squeue failed.'
    batch_materials=()
    job_ids=()
    for material in "${unique_materials[@]}"; do
        [[ ${material%/*} == "$batch" ]] || continue
        batch_materials+=("$material")
        job_ids+=("$(equil_job "$material")")
    done
    ids=$(printf '%s\n' "${job_ids[@]}" | awk '$1 != 0 && !seen[$1]++' | paste -sd , -)
    completed=''
    if [[ -n $ids ]]; then
        completed=$(sacct -X -n -P -j "$ids" -o JobIDRaw,State,ExitCode |
            awk -F'|' '$2=="COMPLETED" && $3=="0:0" {print $1}') || die 'sacct failed.'
    fi
    status=0
    for i in "${!batch_materials[@]}"; do
        material=${batch_materials[i]}
        job=${job_ids[i]}
        glass="$material/GlassTransition"
        if grep -Fxq -- "$material" <<< "$active"; then
            reason='a material job is pending or running'
        elif [[ -e $glass || -L $glass ]]; then
            if [[ -s $glass/Analysis/Tg_weighted.csv ]]; then reason='weighted Tg output already exists'
            elif [[ -s $glass/Analysis/Tg_results.csv ]]; then reason='author results exist; run glass_transition_weighted.py separately if needed'
            else reason='partial GlassTransition output exists; inspect the previous run'; fi
        elif [[ $job == 0 ]]; then reason='no equilibration job log yet'
        elif ! grep -Fxq -- "$job" <<< "$completed"; then reason="equilibration $job has no COMPLETED/0:0 accounting record"
        elif ! result=$(final_coordinates "$material" "$job" 2>&1); then reason=$result
        else
            if [[ $mode == list ]]; then
                printf 'READY %s (equilibration %s)\n' "$material" "$job"
                continue
            fi
            if receipt=$(sbatch --parsable --chdir="$material" "$tutorial/glass_transition_barkla.sh" --run "$tutorial" "$job") &&
                [[ $receipt =~ ^[0-9]+(\;[[:alnum:]_.-]+)?$ ]]; then
                printf '%s -> Slurm %s\n' "$material" "$receipt"
            else
                printf 'Submission failed or receipt unclear for %s; check squeue.\n' "$material" >&2
                status=1
            fi
            continue
        fi
        printf 'SKIP %s (%s)\n' "$material" "$reason"
    done
    exit "$status"
)

unique_materials=()
while IFS= read -r material; do unique_materials+=("$material"); done < <(printf '%s\n' "${materials[@]}" | LC_ALL=C sort -u)
status=0
while IFS= read -r batch; do
    submit_batch "$batch" || status=1
done < <(printf '%s\n' "${unique_materials[@]%/*}" | LC_ALL=C sort -u)
exit "$status"
