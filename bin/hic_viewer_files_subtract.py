"""

"""

import numpy as np
import os
from optparse import OptionParser



def parse_options():
    parser = OptionParser()
    parser.add_option("-f", "--first_folder", dest="folder1",
                      help="First folder containing viewer files", 
                      metavar="FOLDER1")
    parser.add_option("-s", "--second_folder", dest="folder2",
                      help="Second folder containing viewer files", 
                      metavar="FOLDER2")
    parser.add_option("-o", "--out_folder", dest="outfolder",
                      help="Folder for outputs.", 
                      metavar="OUTFOLDER")

    (options, args) = parser.parse_args()
    return options

def get_files(dir):
    if not os.path.isdir(dir):
        raise ValueError('Not a directory')
    files = os.listdir(dir)
    files_good = []
    for f in files:
        if f[-2:] == 'gz':
            files_good.append(f)
    return files_good

options = parse_options()

folder1 = options.folder1
folder2 = options.folder2
outfolder = options.outfolder

files1 = get_files(folder1)
files2 = get_files(folder2)

for f in files1:
    if f[:2] == '2L':
        continue
    if f in files2:
        mat1 = np.loadtxt(os.path.join(folder1, f))
        mat2 = np.loadtxt(os.path.join(folder2, f))
        mat_div = mat1 / (mat2 + 0.001)
        outfile = os.path.join(outfolder, f)
        np.savetxt(outfile, mat_div)
    else:
        print('File ' + f + ' is not in folder 2.')




