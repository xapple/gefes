# Built-in modules #

# Internal modules #
import gefes
from gefes.groups.aggregates import Aggregate

# Third party modules #

###############################################################################
# Test #
proj = gefes.projects['test']
test_agg  = Aggregate( 'test_agg', proj[0:3])

# The depth profile #
proj = gefes.projects['alinen']
alinen_epi  = Aggregate('alinen_epi',  proj[0:4])
alinen_meta = Aggregate('alinen_meta', proj[5:8])
alinen_hypo = Aggregate('alinen_hypo', proj[8:12])