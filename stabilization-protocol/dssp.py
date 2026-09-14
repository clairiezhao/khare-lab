import collections
import collections.abc
collections.Callable = collections.abc.Callable


from ssbio.protein.structure.properties.stride import STRIDE

# Instantiating the wrapper handles the parsing mechanics
stride_runner = STRIDE()
results = stride_runner.run(pdb_file='nitaly.pdb')
