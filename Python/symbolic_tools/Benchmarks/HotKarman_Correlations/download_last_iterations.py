import os
import pandas as pd
from Benchmarks.HotKarman_Correlations.HT_Nu_Correlations import get_Nu_cylinder_by_Churchill_Bernstein
import numpy as np
import pwd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import re
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, strip_folder_name
#######################################################
#                       README
# to see OS environment variables from python script,
# you need to put them in ~/.profile
#
#######################################################

home = pwd.getpwuid(os.getuid()).pw_dir
local_download_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'stuff_pro')
host_folder = os.path.join(f'$FROM_PRO', 'batch_HotKarman3D_CM_HIGHER', 'keep_nu_and_k_sizer*')

if not os.path.exists(local_download_folder):
    os.makedirs(local_download_folder)

cmd = "rsync -zarv  --prune-empty-dirs --include \"*/\"  --include=\"*.csv\" --exclude=\"*\" "\
      + f"\"{host_folder}\""\
      + f" \"{local_download_folder}\""

print(cmd)
# os.system(cmd)

'ssh plgmuaddieb@prometheus.cyfronet.pl ls -l /net/scratch/people/plgmuaddieb/output/'