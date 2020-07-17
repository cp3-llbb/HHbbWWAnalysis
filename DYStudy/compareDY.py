import sys
import os
import yaml
import argparse
from pprint import pprint

class CompareDY:
    def __init__(self,path,yaml_mc,yaml_dd,era):
        config_mc = self.extractYAML(os.path.join(path,yaml_mc))
        config_dd = self.extractYAML(os.path.join(path,yaml_dd))

        self.produceComparisonYAML(path,config_mc,config_dd,era)

    @staticmethod
    def extractYAML(path_yaml):
        with open(path_yaml,'r') as handle:
            config = yaml.load(handle,Loader=yaml.FullLoader)
        return config

    @staticmethod
    def produceComparisonYAML(path,config_mc,config_dd,era):
        config = {}

        # Mix `files` part #
        config['files'] = {}
        for key,val in config_mc["files"].items():
            if "DY" in key:
                config['files'][key] = val
                config['files'][key]['stack-index'] = 0
        for key,val in config_dd["files"].items():
            if "DYEstimation" in key:
                config['files'][key] = val
                config['files'][key]['stack-index'] = 1

        # Mix `groups` part #
        config['groups'] = {
            'DY':{'line-color' : '#1a83a1',
                  'line-width' : 2,
                  'fill-color' : -1,
                  'legend'     : 'Drell-Yan (MC)',
                  'order'      : 1},
             'DYEstimation':{'fill-color' : '#0e0b6b',
                             'line-width' : 2,
                             'fill-color' : -1,
                             'legend'     : 'Drell-Yan (from data)',
                             'order'      : 2}}
            
        # Copy other elements #
        config['configuration'] = config_dd['configuration']
        config['legend'] = config_dd['legend']
        config['plots'] = config_dd['plots']
        config['systematics'] = config_dd['systematics']

        # Save yaml #
        path_yaml = os.path.join(path,'plots_DYComparison.yml')
        with open(path_yaml,'w') as handle:
            yaml.dump(config,handle)
        print ("Saved new YAML in %s"%(path_yaml))

        # Print plotIt command #
        out_path = os.path.join(path,'plots_DYComparison_%s'%era)
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        cmd = "plotIt -i {input} -o {output} -y -e {era} {yaml}".format(**{'input': path,
                                                                           'output': out_path,
                                                                           'era': era,
                                                                           'yaml': path_yaml})
        print ("PlotIt command (in bamboo env) : ")
        print (cmd)





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge usual bamboo output with the DY estimation from data')
    parser.add_argument('--path', action='store', required=True, type=str,
                        help='Bamboo output path')
    parser.add_argument('--yaml_mc', action='store', required=True, type=str,
                        help='Name of the yaml with classic MC')
    parser.add_argument('--yaml_dd', action='store', required=True, type=str,
                        help='Name of the yaml with data-driven')
    parser.add_argument('--era', action='store', required=True, type=str,
                        help='Era')

    args = parser.parse_args()

    instance = CompareDY(args.path,args.yaml_mc,args.yaml_dd,args.era)
