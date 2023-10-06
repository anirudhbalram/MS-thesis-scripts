import pathlib

#this function creates a jobscript that will run Despotic on all selected clouds. "max_lines" is the number of cores that will be used. Change mem-per-cpu, mail-user, and the HiPerGator username (I've written my name - "a.ravishankar" here) as necessary!!!

def create_jobscript(param_filename, max_lines, param_file, SBATCH_args={}):
    default_SBATCH_args = {
        "qos": "narayanan-b",
        "nodes": "1",
        "tasks-per-node": "1",
        "cpus-per-task": "1",
        "mem-per-cpu": "4gb",
        "time": "96:00:00",
        "output": "/dev/null",
        "mail-type":"END,FAIL",
        "mail-user":"ranirudh@students.iisertirupati.ac.in"
    }
    # NOTE: this makes no guarantees about the order in which SBATCH args are laid out
    full_SBATCH_args = {**default_SBATCH_args, **SBATCH_args}
    full_SBATCH_args["array"] = f"1-{max_lines}"

    with open("slick_run_jobscript_2.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.writelines(f"#SBATCH --{arg_name}={arg_val}\n" for arg_name, arg_val in full_SBATCH_args.items())
        f.write("module load conda\n")
        f.write("conda activate /blue/narayanan/a.ravishankar/py38\n")
        f.write(f"python {pathlib.Path(__file__).parent.resolve()}/slick_run.py --parameters {param_file} --cloudinfofile {param_filename} --cloudinfoline $SLURM_ARRAY_TASK_ID\n")
