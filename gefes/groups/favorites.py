# Built-in modules #

# Internal modules #
import gefes
from gefes.aggretate import Aggretate

# Third party modules #

###############################################################################
# Test #
proj = gefes.projects['test']
test_agg  = Aggretate(proj[0:3],  'test_agg')

# The depth profile #
proj = gefes.projects['alinen']
alinen_epi  = Aggretate(proj[0:4],  'alinen_epi')
alinen_meta = Aggretate(proj[5:8],  'alinen_meta')
alinen_hypo = Aggretate(proj[8:12], 'alinen_hypo')