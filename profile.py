import cProfile
import os
import sys
import pstats
"""
Import functions
"""

"""
TODO:
- open all data needed here
- run any initial methods that require output to be passed to functions in next section
"""


def test():
    """
    TODO: Call methods you want to profile here
    """


cProfile.run('test()', 'OUTFILE')
p = pstats.Stats('OUTFILE')
p.sort_stats('tottime')
p.print_stats()