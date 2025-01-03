import os

flag = '-avchP'
source = '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/'
dest = '/Volumes/ECHO/Research/'
cmd = f'rsync {flag} {source} {dest}'
##### os.system(cmd)
print('\n'*5)
print(cmd)
print('\n'*5)