# Built-in modules #

# Internal modules #
from gefes.running import Runner
from gefes.common.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
class AggregateRunner(Runner):
    default_time = '6:00:00'

    default_steps = [
        {'assemble':     {}},
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.project = parent, parent
        self.pools = parent.pools

    def run_slurm(self, steps=None, **kwargs):
        # Make script #
        command = """steps = %s
                     exp = [exp for exp in XXXXXXXXX if exp.name=='%s'][0]
                     exp.run(steps)""" % (steps, self.experiment.name)
        # Test case #
        if 'test' in self.experiment.name:
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'
        # Send it #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        self.slurm_job = SLURMJob(command, self.experiment.p.logs_dir, job_name="humic_" + self.experiment.name, **kwargs)
        self.slurm_job.launch()
