"""
Generate pdbqt files for docking from the sampled files generated by the model.
"""
import argparse
from os import listdir, makedirs
from os.path import exists, isfile, join
import yaml
import subprocess
from multiprocessing import Pool


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-num_workers",
                        required=False,
                        default=20,
                        help="number of workders to generate the pdbqt files")

    parser.add_argument("-in_dir",
                        required=False,
                        default="../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample/",
                        help="directory of smi files after being selected by ranking frequencies")

    parser.add_argument("-smi_dir",
                        required=False,
                        default="../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_ranked_smi/",
                        help="directory of smi files of sampled SMILES")

    parser.add_argument("-pdbqt_dir",
                        required=False,
                        default="../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_ranked_pdbqt/",
                        help="directory of pdbqt files of sampled SMILES")

    parser.add_argument("-dock_box_dir",
                        required=False,
                        default="../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_ranked_dock_box/",
                        help="directory of pdbqt files of sampled SMILES")

    return parser.parse_args()


class PreparePDBQT:
    def __init__(self, in_dir, smi_dir, pdbqt_dir, dock_box_dir, num_mols=100, num_workers=20):
        if not exists(smi_dir):
            makedirs(smi_dir)
        self.smi_dir = smi_dir

        if not exists(pdbqt_dir):
            makedirs(pdbqt_dir)
        self.pdbqt_dir = pdbqt_dir
        
        if not exists(dock_box_dir):
            makedirs(dock_box_dir)
        self.dock_box_dir = dock_box_dir
        
        # number of molecules selected for each pocket
        self.num_mols = num_mols

        # number of processs for parallelism
        self.num_workers = num_workers
        
        self.eboxsize = './eboxsize.pl'
        
        # input files
        self.in_files = [f for f in listdir(in_dir) if isfile(join(in_dir, f))]

    def prepare_pdbqt(self):
        with Pool(self.num_workers) as p:
            p.map(self.prepare_pdbqt_single_proc, self.in_files)

    def prepare_pdbqt_single_proc(self, in_file):
        """
        Prepare pdbqt files of the molecules for a single pocket, and compute their
        docking box sizes. 
        """
        pocket = in_file.split('_')[0]

        # load the SMILES
        with open(join(in_dir, in_file), 'r') as f:
            smiles_freq = yaml.full_load(f)

        # sort the SMILES according to their frequencies
        smiles = [(freq, smiles) for smiles, freq in smiles_freq.items()]
        smiles.sort(key=lambda x: x[0], reverse=True)
        
        # take the SMILES with highest frequencies
        smi_path = join(self.smi_dir, pocket + '.smi')
        # write them into a .smi file
        self.write_smi_file(smiles, smi_path)

        # convert the .smi file to pdbqt files
        pocket_pdbqt_dir = join(self.pdbqt_dir, pocket)
        self.gen_pdbqt(smi_path, pocket_pdbqt_dir, pocket)

        # compute the docking boxes of the pdbqt files
        self.compute_dock_box(pocket_pdbqt_dir, self.dock_box_dir, pocket)


    def write_smi_file(self, smiles, smi_path):
        with open(smi_path, 'w') as f:
            for s in smiles[0:self.num_mols]:
                # s[0] is the frequency
                f.write(s[1] + '\n')

    @staticmethod
    def gen_pdbqt(smi_path, pdbqt_dir, pocket):
        if not exists(pdbqt_dir):
            makedirs(pdbqt_dir)
        pdbqt_path = join(pdbqt_dir, f'{pocket}-.pdbqt')
        p = subprocess.run(
            'obabel' + 
            f' -ismi {smi_path}' +
            ' -opdbqt' +
            f' -O {pdbqt_path}' +
            ' -h -m --gen3D',
            shell=True,
            capture_output=True
        )

        if p.returncode != 0:
            print('Something went wrong when generating pdbqt file.')
            print('input file: ', smi_path)
            print('output file: ', pdbqt_path)
            return

    def compute_dock_box(self, pocket_pdbqt_dir, dock_box_dir, pocket):
        # list all pdbqt files of this pocket
        pdbqt_files = [f for f in listdir(pocket_pdbqt_dir) if isfile(join(pocket_pdbqt_dir, f))]
        docking_boxes = {}

        # compute size of docking box of each molecule
        for mol in pdbqt_files:
            mol_path = join(pocket_pdbqt_dir, mol)
            p = subprocess.run(
                self.eboxsize + ' {}'.format(mol_path),
                shell=True,
                stdout=subprocess.PIPE,
                text=True
            )

            # when there is no error
            if p.returncode == 0:
                result = p.stdout
                # add to result dict
                docking_boxes[mol.split()[0]] = float(result.strip())
            else:
                print('Something went wrong when computing docking box, ligand: {}'.format(mol_path))

        # save docking box sizes of this pocket
        with open(dock_box_dir + pocket + "_docking_box.yaml", "w") as f:
            yaml.dump(docking_boxes, f)


if __name__ == '__main__':
    args = get_args()
    in_dir = args.in_dir
    smi_dir = args.smi_dir
    pdbqt_dir = args.pdbqt_dir
    dock_box_dir = args.dock_box_dir
    num_workers = args.num_workers

    pp = PreparePDBQT(
        in_dir=in_dir, 
        smi_dir=smi_dir,
        pdbqt_dir=pdbqt_dir,
        dock_box_dir=dock_box_dir,
        num_workers=num_workers
    )

    pp.prepare_pdbqt()

    