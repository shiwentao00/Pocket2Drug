import subprocess


def parse_smina_output(result):
    """parse the output texts of smina and return
    the free energy of the best docking pose"""
    # binary stream to string
    result = result.decode('utf-8')
    result = result.split('\n')

    # line 29 is the result of best docking pose
    return float(result[29].split()[1])


def smina_dock(smina, protein, ligand, pocket_center, docking_box,
               exhaustiveness=8, num_modes=2, energy_range=99):
    x = pocket_center[0]
    y = pocket_center[1]
    z = pocket_center[2]
    box_size = docking_box
    p = subprocess.run(
        smina + ' -r {} -l {}'.format(protein, ligand) +
        ' --center_x {} --center_y {} --center_z {}'.format(x, y, z) +
        ' --size_x {} --size_y {} --size_z {}'.format(box_size, box_size, box_size) +
        # ' -o {}'.format(out_path) +
        ' --exhaustiveness {}'.format(exhaustiveness) +
        ' --num_modes {}'.format(num_modes) +
        ' --energy_range {}'.format(energy_range),
        shell=True,
        capture_output=True
        # stdout=subprocess.PIPE,
        # text=True
    )

    # print(p.stdout)
    # print(p.stderr)
    if p.returncode == 0:
        result = p.stdout
        score = parse_smina_output(result)
        return score, p.returncode
    else:
        return None, p.returncode
