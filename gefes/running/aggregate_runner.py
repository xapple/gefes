# Built-in modules #

# Internal modules #
from plumbing.runner import Runner
from plumbing.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
class AggregateRunner(Runner):
    """Will run stuff on an aggregate"""
    default_time = '2-00:00:00'

    default_steps = [
        {'assemble':     {}},
        {'make_plots':   {}},
    ]

    def __init__(self, parent):
        self.parent, self.project = parent, parent
        self.samples = parent.samples

    def run_slurm(self, steps=None, **kwargs):
        # Test case #
        if self.cluster.name == 'test':
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'
        # Make script #
        command =  ["steps = %s" % steps]
        command += ["name = '%s'" % self.cluster.name]
        command += ["cluster = getattr(illumitag.clustering.favorites, name)"]
        command += ["cluster.run(steps)"]
        # Send it #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        if 'dependency' not in kwargs: kwargs['dependency'] = 'singleton'
        job_name = "agg_%s" % self.cluster.name
        self.slurm_job = SLURMJob(command, self.parent.p.logs_dir, job_name=job_name, **kwargs)
        return self.slurm_job.run()