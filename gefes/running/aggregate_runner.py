# Built-in modules #

# Internal modules #
import gefes
from plumbing.runner import Runner
from plumbing.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
class AggregateRunner(Runner):
    """Will run stuff on an aggregate"""
    default_time = '2-00:00:00'

    default_steps = [
        {'assembly.run':     {}},
    ]

    def __init__(self, parent):
        self.parent, self.aggregate = parent, parent
        self.samples = parent.samples

    @property
    def job_name(self): return "agg_%s" % self.parent.name

    def command(self, steps):
        from gefes.groups.projects import Project
        command =  ["steps = %s" % steps]
        command += ["name = '%s'" % self.parent.name]
        if isinstance(self.parent, Project): command += ["aggregate = gefes.projects[name]"]
        else: command += ["aggregate = getattr(gefes.groups.favorites, name)"]
        command += ["aggregate.load()"]
        command += ["aggregate.runner.run(steps)"]
        return command

    def run_slurm(self, steps=None, **kwargs):
        # Test case #
        if self.parent.name == 'test':
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = 'short'
            kwargs['email'] = '/dev/null'
        # Extra params #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        if 'dependency' not in kwargs: kwargs['dependency'] = 'singleton'
        # Send it #
        if 'job_name' not in kwargs: kwargs['job_name'] = self.job_name
        self.slurm_job = SLURMJob(self.command(steps), self.parent.p.logs_dir, gefes, **kwargs)
        return self.slurm_job.run()