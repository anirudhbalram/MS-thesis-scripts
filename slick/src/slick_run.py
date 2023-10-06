import pandas as pd
import argparse
import glob

from limfunctions import creating_table
from read_config import parse_parameters

parser = argparse.ArgumentParser(prog='lim_code_simple')
parser.add_argument("-p", "--parameters", type=str)
parser.add_argument("--cloudinfofile", type=str)
parser.add_argument("--cloudinfoline", type=int)

args = parser.parse_args()

cloud_info_file = args.cloudinfofile
with open(cloud_info_file, "r") as f:
    [*cloud_list] = [int(x) for x in next(x for i, x in enumerate(f) if i == args.cloudinfoline-1).split()]    #reading the cloud particle indices for running Despotic
    
config = parse_parameters(args.parameters)
df_basic = pd.read_csv("results/Basic_Characteristics_z_2_MS.csv")
# Creates final luminosity table
creating_table(cloud_list, df_basic, config["output_dir"])
print("Slick run completed successfully")
