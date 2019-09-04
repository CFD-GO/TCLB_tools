
from SymbolicCollisions.core.cm_symbols import e_D3Q27, moments_dict
from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_indx_notation
from sympy import Matrix

import re
import os
import pwd

# SETUP
d = 3
q = 27

populations = get_print_symbols_in_indx_notation(q, print_symbol='h')
print_as_vector(populations, outprint_symbol='f')
print()
print_as_vector(populations, outprint_symbol='f', output_order_of_moments=e_D3Q27)
print()
print_as_vector(populations, outprint_symbol='f', output_order_of_moments=moments_dict[f'D{d}Q{q}'])
print("Example: Parsing a file.")

# home = pwd.getpwuid(os.getuid()).pw_dir
# main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'HotBarman3D')
my_dir = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(my_dir, "some_code.txt")
with open(path, "r") as f:
    stuff = f.read()

print(stuff)

# stuff = re.sub(r'f000', 'h[0]', stuff)  # this is rather wrong...
# stuff = re.sub(r'f100', 'h[1]', stuff)
# stuff = re.sub(r'f200', 'h[2]', stuff)
# stuff = re.sub(r'f010', 'h[3]', stuff)
# stuff = re.sub(r'f020', 'h[4]', stuff)
# stuff = re.sub(r'f001', 'h[5]', stuff)
# stuff = re.sub(r'f002', 'h[6]', stuff)
# stuff = re.sub(r'f111', 'h[7]', stuff)
# stuff = re.sub(r'f211', 'h[8]', stuff)
# stuff = re.sub(r'f121', 'h[9]', stuff)
# stuff = re.sub(r'f221', 'h[10]', stuff)
# stuff = re.sub(r'f112', 'h[11]', stuff)
# stuff = re.sub(r'f212', 'h[12]', stuff)
# stuff = re.sub(r'f122', 'h[13]', stuff)
# stuff = re.sub(r'f222', 'h[14]', stuff)
# stuff = re.sub(r'f110', 'h[15]', stuff)
# stuff = re.sub(r'f210', 'h[16]', stuff)
# stuff = re.sub(r'f120', 'h[17]', stuff)
# stuff = re.sub(r'f220', 'h[18]', stuff)
# stuff = re.sub(r'f101', 'h[19]', stuff)
# stuff = re.sub(r'f201', 'h[20]', stuff)
# stuff = re.sub(r'f102', 'h[21]', stuff)
# stuff = re.sub(r'f202', 'h[22]', stuff)
# stuff = re.sub(r'f011', 'h[23]', stuff)
# stuff = re.sub(r'f021', 'h[24]', stuff)
# stuff = re.sub(r'f012', 'h[25]', stuff)
# stuff = re.sub(r'f022', 'h[26]', stuff)

stuff = re.sub(r'f000', 'h[0]', stuff) # this is rather good...
stuff = re.sub(r'f100', 'h[1]', stuff)
stuff = re.sub(r'f010', 'h[2]', stuff)
stuff = re.sub(r'f001', 'h[3]', stuff)
stuff = re.sub(r'f110', 'h[4]', stuff)
stuff = re.sub(r'f101', 'h[5]', stuff)
stuff = re.sub(r'f011', 'h[6]', stuff)
stuff = re.sub(r'f200', 'h[7]', stuff)
stuff = re.sub(r'f020', 'h[8]', stuff)
stuff = re.sub(r'f002', 'h[9]', stuff)
stuff = re.sub(r'f120', 'h[10]', stuff)
stuff = re.sub(r'f102', 'h[11]', stuff)
stuff = re.sub(r'f210', 'h[12]', stuff)
stuff = re.sub(r'f201', 'h[13]', stuff)
stuff = re.sub(r'f012', 'h[14]', stuff)
stuff = re.sub(r'f021', 'h[15]', stuff)
stuff = re.sub(r'f111', 'h[16]', stuff)
stuff = re.sub(r'f220', 'h[17]', stuff)
stuff = re.sub(r'f202', 'h[18]', stuff)
stuff = re.sub(r'f022', 'h[19]', stuff)
stuff = re.sub(r'f211', 'h[20]', stuff)
stuff = re.sub(r'f121', 'h[21]', stuff)
stuff = re.sub(r'f112', 'h[22]', stuff)
stuff = re.sub(r'f122', 'h[23]', stuff)
stuff = re.sub(r'f212', 'h[24]', stuff)
stuff = re.sub(r'f221', 'h[25]', stuff)
stuff = re.sub(r'f222', 'h[26]', stuff)

print(stuff)

path2 = os.path.join(my_dir, "parsed_code.txt")
with open(path2, "w") as f:
    f.write(stuff)

print("BYE")

