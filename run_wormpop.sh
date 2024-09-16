#!/bin/bash

run_wormpop()
{(
    set -Eeo pipefail
    eval "$(conda shell.bash hook)"
    conda activate simulation_env

    python3 main.py --parameters=constants.json --variants=variants.json --database=database.sqlite
)}

export -f run_wormpop

srun --time=2-00:00:00 --export=ALL --exclusive --job-name=plot_green bash -c 'run_wormpop'
#srun -N1 --exclusive --pty /usr/bin/env python3 simulation.py --parameters=constants.json --variants=variants.json --database=database.sqlite
#srun -N1 --exclusive --pty /usr/bin/env python3 python3 simulation.py --parameters=constants.json --variants=variants.json --database=database.sqlite