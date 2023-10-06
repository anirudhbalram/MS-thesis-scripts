from sys import argv
from os import mkdir

from create_basic_characteristics_table import create_basic_table
from create_cloudspercore_list import create_cloudspercore_list
from create_jobscript import create_jobscript
from read_config import parse_parameters

def main():
    if len(argv) > 1:
        config_file = argv[1]
    else:
        print("no parameter file given.")
        # TODO: print some sort of usage message here?
        exit(1)

    config = parse_parameters(config_file)    #this parses the parameters provided in the config_file "parameters.ini"
    
    mkdir(config["output_dir"])    #creates a new output directory if it's not already there. Comment out this line of code if you've already created a directory by hand

    # Creates table with basic characteristics for all the clouds
    create_basic_table(config)

    # Creates table with basic characteristics for all the clouds (if mode='total'), or for a sample of them (mode = 'randomize')
    param_filename, max_lines = create_cloudspercore_list(config)    #max_lines will be the number of cores used and will have to be written into slick_run_jobscript.sh

    # Creates slick_run_jobscript.sh
    create_jobscript(param_filename, max_lines, config_file, config["sbatch"])

    if config["skip_run"]:
        # non-0 exit code signals to not do the slick run step
        exit(1)

if __name__ == "__main__":
    main()
