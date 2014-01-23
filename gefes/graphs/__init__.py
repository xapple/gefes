# Built-in modules #
import os, time, getpass, locale

# Third party modules #
import matplotlib
from matplotlib import animation


# No need for an X display #
matplotlib.use('Agg', warn=False)

################################################################################
class Graph(object):

    def __init__(self, parent, base_dir=None, short_name=None):
        # Save parent #
        self.parent = parent
        # Base dir #
        if not base_dir: self.base_dir = self.parent.p.graphs_dir
        else: self.base_dir = base_dir
        # Short name #
        if short_name: self.short_name = short_name
        # Paths #
        self.path = self.base_dir + self.short_name + '.pdf'
        self.avi_path = self.base_dir + self.short_name + '.avi'
        self.csv_path = self.base_dir + self.short_name + '.csv'
        self.json_path = self.base_dir + self.short_name + '.json'
        # Extra #
        self.dev_mode = True

    def save_plot(self, fig, axes, width=12.0, height=7.0, bottom=0.14, top=0.93, left=0.06, right=0.98, sep=()):
        # Adjust #
        fig.set_figwidth(width)
        fig.set_figheight(height)
        fig.subplots_adjust(hspace=0.0, bottom=bottom, top=top, left=left, right=right)
        # Data and source #
        if self.dev_mode:
            fig.text(0.99, 0.98, time.asctime(), horizontalalignment='right')
            job_name = os.environ.get('SLURM_JOB_NAME', 'Unnamed')
            user_msg = 'user: %s, job: %s' % (getpass.getuser(), job_name)
            fig.text(0.01, 0.98, user_msg, horizontalalignment='left')
        # Nice digit grouping #
        if 'x' in sep:
            locale.setlocale(locale.LC_ALL, '')
            seperate = lambda x,pos: locale.format("%d", x, grouping=True)
            axes.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(seperate))
        if 'y' in sep:
            locale.setlocale(locale.LC_ALL, '')
            seperate = lambda x,pos: locale.format("%d", x, grouping=True)
            axes.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(seperate))
        # Save it #
        print(self.path)
        fig.savefig(self.path)

    def save_anim(self, fig, animate, init, width=15, height=15, bitrate=10000,fps=10):

        fig.set_figwidth(height)
        fig.set_figheight(width)

        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=360, interval=20)

        FFMpegWriter = animation.writers['ffmpeg']
        writer = FFMpegWriter( bitrate= bitrate, fps=fps)

        # Save
        anim.save(self.avi_path, writer=writer, codec='x264')#,  extra_args=['--verbose-debug','-vcodec', 'libx264'])
