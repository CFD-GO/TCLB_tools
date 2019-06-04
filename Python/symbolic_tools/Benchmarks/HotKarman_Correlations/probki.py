
import os
import subprocess
import pandas as pd
import numpy as np
import pwd

#
# s = pd.Series([0, 1, 2, 3, 4, 5, 6, 7])
#
# window = 3
# aa = s.rolling(window).mean()
#
#
# y = aa[window-1:]
# x = np.linspace(start=window-1, stop=len(y)+1, num=len(y), endpoint=True)
#
# print(aa)

# cmd = "rsync -zarv  --prune-empty-dirs --include \"*/\"  --include=\"*.csv\" --exclude=\"*\" \"$FROM_PRO/batch_HotKarman3D\"" + f" \"{local_logs_folder}\""

# cmd = f"rsync -zarv  --prune-empty-dirs " \
#     f"--include \"*/\"  " \
#     f"--include=\"*.csv\" " \
#     f"--exclude=\"*\" " \
#     f"\"plgmuaddieb@prometheus.cyfronet.pl:/net/scratch/people/plgmuaddieb/output/batch_HotKarman3D/keep_nuk_old/\" \"{local_logs_folder}\""

cmd ="echo  hoho ${FROM_PRO}"
# print(cmd)

os.system(cmd)
os.system('/bin/bash -c "echo $HOME aa $FROM_PRO"')

subprocess.call(["bash", "-c", cmd])

print('FROM_PRO' in os.environ)
print('HOME' in os.environ)