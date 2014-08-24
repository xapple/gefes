# Built-in modules #

# Internal modules #
from gefes.running import Runner
from plumbing.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
class ProjectRunner(Runner):
    """Will run stuff on a project"""
    default_time = '7-00:00:00'

    default_steps = [
        {'assemble':     {}},
        {'make_plots':   {}},
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.project = parent, parent
        self.samples = parent.samples

    def run_slurm(self, steps=None, **kwargs):
        # Make script #
        command = """steps = %s
                     proj = [proj for proj in gefes.projects if proj.name=='%s'][0]
                     proj.runner(steps)""" % (steps, self.project.name)
        # Test case #
        if 'test' in self.project.name:
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'
            if kwargs.get('cluster') == 'halvan':
                kwargs['time'] = '1-00:00:00'
                kwargs.pop('qos')
        # Send it #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        job_name = "gefes_%s" % self.project
        self.slurm_job = SLURMJob(command, self.project.p.logs_dir, job_name=job_name, **kwargs)
        self.slurm_job.launch()
