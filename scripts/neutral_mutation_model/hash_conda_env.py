import hashlib
import os

# Path as viewed from cluster
_env_dir=os.path.realpath("/mnt/grid/siepel/hpc_norepl/home/data/ndukler/.extraINSIGHT_scratch_grch38/.conda")
md5hash = hashlib.md5()
md5hash.update(_env_dir.encode())

f = open("/sonas-hs/siepel/hpc_norepl/home/data/ndukler/.extraINSIGHT_scratch_grch38/.conda", 'rb')
md5hash.update(f.read())
f.close()
h = md5hash.hexdigest()
# Only use 1st eight characters
print(h)
# Resulting has is first 8bb96700f
