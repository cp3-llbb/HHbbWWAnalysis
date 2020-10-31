import sys
import os
import yaml
import glob
import argparse
from copy import copy
from pprint import pprint

class PseudoData:
    def __init__(self,path,yaml,era):
        config = self.extractYAML(os.path.join(path,yaml))

        self.producePseudoDataYAML(path,config,era)

    @staticmethod
    def extractYAML(path_yaml):
        with open(path_yaml,'r') as handle:
            config = yaml.load(handle,Loader=yaml.FullLoader)
        return config

    @staticmethod
    def producePseudoDataYAML(path,config,era):
        files_dict = {}
        for rootfile in glob.glob(os.path.join(path,'results','*root')):
            sample = os.path.basename(rootfile)
            if '__skeleton__' in sample:
                continue
            if sample in config["files"]:
                files_dict[sample] = copy(config["files"][sample])
                files_dict[sample].update({'stack-index':0})
            else:
                if 'Pseudodata' in sample:
                    base_sample = sample.replace('Pseudodata.root','')
                    for name in config["files"].keys():
                        if base_sample in name:
                            files_dict[sample] = copy(config["files"][name])
                            files_dict[sample].update(
                                    {'group':'pseudodata','type':'mc',
                                     'stack-index':1})
                            break
                    if sample not in files_dict.keys():
                        raise RuntimeError('Could not find an entry associated to %s'%sample)

        config["files"] = files_dict

        config['groups'].update({'pseudodata':{'legend':'pseudo-data',
                                               'drawing-options': 'P',
                                               'marker-color': 1,
                                               'marker-size': 1,
                                               'marker-type': 20}})

        # Save yaml #
        path_yaml = os.path.join(path,'plots_pseudodata.yml')
        with open(path_yaml,'w') as handle:
            yaml.dump(config,handle)
        print ("Saved new YAML in %s"%(path_yaml))

        # Print plotIt command #
        out_path = os.path.join(path,'plots_pseudodata_%s'%era)
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        cmd = "plotIt -i {input} -o {output} -e {era} {yaml}".format(**{'input': path,
                                                                           'output': out_path,
                                                                           'era': era,
                                                                           'yaml': path_yaml})
        print ("PlotIt command (in bamboo env) : ")
        print (cmd)





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Produce plots with MC estimate')
    parser.add_argument('--path', action='store', required=True, type=str,
                        help='Bamboo output path')
    parser.add_argument('--yaml', action='store', required=True, type=str,
                        help='Name of the yaml')
    parser.add_argument('--era', action='store', required=True, type=str,
                        help='Era')

    args = parser.parse_args()

    instance = PseudoData(args.path,args.yaml,args.era)
