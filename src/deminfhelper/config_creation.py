import yaml
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--out_dir', type=str, required=True)
parser.add_argument('--vcf_dir', type=str, required=True)
parser.add_argument('--generation_time', type=float, required=True)
parser.add_argument('--mutation_rate', type=float, required=True)
parser.add_argument('--p0', type=str, required=True)
parser.add_argument('--lower_bound', type=str, required=True)
parser.add_argument('--upper_bound', type=str, required=True)
args = parser.parse_args()


d = {}
out_dir = args.out_dir
d['out_dir'] = out_dir
vcf = args.vcf_dir
d['vcf'] = vcf
gen_time = args.generation_time
d['gen_time'] = gen_time
mut_rate = args.mutation_rate # Mutation rate
d['mut_rate'] = mut_rate
p0 = args.p0
d['p0'] = p0
lower_bound = args.lower_bound
d['lower_bound'] = lower_bound
upper_bound = args.upper_bound
d['upper_bound'] = upper_bound

d["name_pop"] = "cerevisiae"
d["out_dir_sfs"] = out_dir 
d["path_to_sfs"] = out_dir + "SFS_" + d["name_pop"] + ".fs"
d["out_dir_stairwayplot2"] = out_dir + "output_stairwayplot2/"
d["summary_file_stw"] = d["out_dir_stairwayplot2"] + "s_" + d["name_pop"] + "/" + d["name_pop"] + ".final.summary"
d["out_dir_dadi"] = out_dir + "output_dadi/"
d["out_dir_msmc2"] = out_dir + "output_msmc2/"
d["out_dir_smcpp"] = out_dir + "output_smcpp/"
d["out_dir_psmc"] = out_dir + "output_psmc/"
d["plot_file_smcpp"] = d["out_dir_smcpp"] + d["name_pop"] + "_inference.csv"
d["out_dir_gq_distrib"] = out_dir + "output_stats/"
d["final_out_dir"] = out_dir + "inferences/"
d["out_dir_stats"] = out_dir + "output_stats/"



def update_yaml_params(yaml_file_path, params, values, output_file_path=None):
    """
    Update parameters in a YAML file with given values.

    Args:
        yaml_file_path (str): Path to the original YAML file.
        params (list of str): List of parameter keys to update.
        values (list): List of new values corresponding to params.
        output_file_path (str, optional): Path to save updated YAML.
                                          If None, overwrite original file.

    Raises:
        ValueError: If params and values lengths do not match.
    """

    if len(params) != len(values):
        raise ValueError("Length of params and values must be the same.")

    # Load YAML
    with open(yaml_file_path, 'r') as f:
        data = yaml.safe_load(f)
   
    # Update parameters
    for param, value in zip(params, values):
        if param in data:
            data[param] = value
        else:
            print(f"Warning: Parameter '{param}' not found in YAML. Adding it.")
            data[param] = value

    # Save updated YAML
    save_path = output_file_path if output_file_path else yaml_file_path
    with open(save_path, 'w') as f:
        yaml.safe_dump(data, f, sort_keys=False)

    print(f"YAML file updated and saved to: {save_path}")

yaml_path = "./test.yml"

update_yaml_params(yaml_path, list(d.keys()), list(d.values()))
