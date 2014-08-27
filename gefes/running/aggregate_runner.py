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
        {'assemble.run':     {}},
    ]

    def __init__(self, parent):
        self.parent, self.project = parent, parent
        self.samples = parent.samples

    @property
    def job_name(self): return "agg_%s" % self.parent.name

    def command(self, steps):
        command =  ["steps = %s" % steps]
        command += ["name = '%s'" % self.parent.name]
        command += ["cluster = getattr(gefes.groups.favorites, name)"]
        command += ["cluster.run(steps)"]
        return command

    def run_slurm(self, steps=None, **kwargs):
        # Test case #
        if self.parent.name == 'test':
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'
        # Extra params #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        if 'dependency' not in kwargs: kwargs['dependency'] = 'singleton'
        # Send it #
        self.slurm_job = SLURMJob(self.command(steps), self.parent.p.logs_dir, job_name=self.job_name, **kwargs)
        return self.slurm_job.run()