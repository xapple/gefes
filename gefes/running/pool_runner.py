# Built-in modules #

# Internal modules #
from gefes.running import Runner
from gefes.common.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
class PoolRunner(Runner):
    """Will run stuff on a pool"""
    default_time = '1-00:00:00'

    default_steps = [
        # Initial steps #
        {'clean_reads':           {}},
        {'make_plots':            {}},
        # Run after the assembly is done #
        {'map_reads':             {}},
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.pool = parent, parent
        # Default variables #
        self.job_start_time = None
        self.job_end_time = None
        self.job_runtime = None

    def run_slurm(self, steps=None, **kwargs):
        # Script #
        command = """steps = %s
                     pool = [p for p in gefes.pools if str(p)=='%s'][0]
                     pool.runner(steps)""" % (steps, self.pool)
        # Params #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        if 'dependency' not in kwargs: kwargs['dependency'] = 'singleton'
        # Test case #
        if self.pool.project.name == 'test':
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'
        # Send it #
        job_name = "gefes_%s" % self.pool
        self.pool.slurm_job = SLURMJob(command, self.pool.p.logs_dir, job_name=job_name, **kwargs)
        return self.pool.slurm_job.launch()