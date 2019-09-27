
import os
import re

from SymbolicCollisions.core.cm_symbols import walberla_csys

print("Example: Parsing a file.")

# home = pwd.getpwuid(os.getuid()).pw_dir
# main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'HotBarman3D')
my_dir = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(my_dir, "some_walberla_code.cpp")
with open(path, "r") as f:
    stuff = f.read()

print('=== raw input ===')
print(stuff)

# old parser
# for key, value in walberla_csys.items():
#     stuff = re.sub(key, value, stuff)
#
# print('=== NSEWTB -> mom ===')
# stuff = re.sub('const ', '', stuff)
# print(stuff)
# stuff = re.sub('dst', 'src', stuff)
# print('\n\n')
# print(stuff)

print('=== NSEWTB -> mom ===')


def read_from(src_or_dst):
    for key, value in walberla_csys.items():
        parsed_key = key.rstrip(' ')
        parsed_key = re.sub('v', '', parsed_key)
        print(f"{value} = {src_or_dst}->get(x, y, z, Stencil_T::idx[{parsed_key}]);")


read_from('src')
print('\n\n')
read_from('dst')


print('=== mom -> NSEWTB ===')


def write_to(src_or_dst):
    for key, value in walberla_csys.items():
        parsed_key = key.rstrip(' ')
        parsed_key = re.sub('v', '', parsed_key)
        print(f"{src_or_dst}->get(x, y, z, Stencil_T::idx[{parsed_key}]) = {value};")


write_to('src')
print('\n\n')
write_to('dst')


path2 = os.path.join(my_dir, "parsed_code.txt")
with open(path2, "w") as f:
    f.write(stuff)

print("BYE")

