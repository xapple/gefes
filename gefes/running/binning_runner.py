# Built-in modules #

# Internal modules #
from gefes.running import Runner
from gefes.common.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
class BinningRunner(Runner):
    """Will run stuff on a project"""
    default_time = '7-00:00:00'

    default_steps = [
        {'cluster':     {}},
        {'annotate':   {}},
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.binning = parent, parent
        self.project = self.binning.parent.parent

    def run_slurm(self, steps=None, **kwargs):
        # Make script #
        command = """steps = %s
                     binning = [binning for binning in gefes.projects['%s'].binner if binning.name=='%s'][0]
                     binning.runner(steps)""" % (steps,self.project.name, self.binning.name)
        # Test case #
        if 'test' in self.project.name:
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'

        # Send it #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        job_name = "gefes_%s_%s" % (self.project.name,self.binning.name)
        self.slurm_job = SLURMJob(command, self.binning.p.logs_dir, job_name=job_name, **kwargs)
        self.slurm_job.launch()
