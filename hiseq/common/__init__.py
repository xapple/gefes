# Built-in modules #
import os, shutil, re, platform, random

# Third party modules #
import sh

# Expositions #
from color import Color
from autopaths import AutoPaths

################################################################################
def isubsample(full_sample, k, full_sample_len=None):
    """Downsample an enumerable list of things"""
    # Determine length #
    if not full_sample_len: full_sample_len = len(full_sample)
    # Check size coherence #
    if not 0 <= k <= full_sample_len:
        raise ValueError('Required that 0 <= k <= full_sample_length')
    # Do it #
    picked = 0
    for i, element in enumerate(full_sample):
        prob = (k-picked) / (full_sample_len-i)
        if random.random() < prob:
            yield element
            picked += 1
    # Did we pick the right amount #
    assert picked == k

################################################################################
def get_git_tag(directory):
    if os.path.exists(directory + '/.git'):
        return sh.git("--git-dir=" + directory + '/.git', "describe", "--tags", "--dirty", "--always").strip('\n')
    else:
        return None

###############################################################################
def reversed_lines(path):
    """Generate the lines of file in reverse order."""
    with open(path, 'r') as handle:
        part = ''
        for block in reversed_blocks(handle):
            for c in reversed(block):
                if c == '\n' and part:
                    yield part[::-1]
                    part = ''
                part += c
        if part: yield part[::-1]

def reversed_blocks(handle, blocksize=4096):
    """Generate blocks of file's contents in reverse order."""
    handle.seek(0, os.SEEK_END)
    here = handle.tell()
    while 0 < here:
        delta = min(blocksize, here)
        here -= delta
        handle.seek(here, os.SEEK_SET)
        yield handle.read(delta)

###############################################################################
def move_with_overwrite(source, dest):
    if os.path.exists(dest):
        if os.path.isdir(dest): shutil.rmtree(dest)
        else: os.remove(dest)
    shutil.move(source, dest)

###############################################################################
def replace_extension(path, new_ext):
    if not new_ext.startswith('.'): new_ext = '.' + new_ext
    base, ext = os.path.splitext(path)
    return base + new_ext

###############################################################################
def find_file_by_name(name, root=os.curdir):
    for dirpath, dirnames, filenames in os.walk(os.path.abspath(root)):
        if name in filenames: return os.path.join(dirpath, name)
    raise Exception("Could not find file '%s' in '%s'") % (name, root)

###############################################################################
def natural_sort(item):
    """
    Sort strings that contain numbers correctly.

    >>> l = ['v1.3.12', 'v1.3.3', 'v1.2.5', 'v1.2.15', 'v1.2.3', 'v1.2.1']
    >>> l.sort(key=natural_sort)
    >>> l.__repr__()
    "['v1.2.1', 'v1.2.3', 'v1.2.5', 'v1.2.15', 'v1.3.3', 'v1.3.12']"
    """
    if item is None: return 0
    def try_int(s):
        try: return int(s)
        except ValueError: return s
    return map(try_int, re.findall(r'(\d+|\D+)', item))

###############################################################################
def which(cmd):
    """https://github.com/jc0n/python-which"""
    def is_executable(path):
        return os.path.exists(path) and os.access(path, os.X_OK)
    if platform.system() == 'Windows':
        def possible_file_extensions(path):
            yield path
            for ext in os.environ.get('PATHEXT', '').split(os.pathsep):
                yield path + ext
    else:
        def possible_file_extensions(path):
            yield path
    path, name = os.path.split(cmd)
    if path:
        if is_executable(cmd):
            return cmd
    else:
        for path in os.environ['PATH'].split(os.pathsep):
            cmd_file = os.path.join(path, cmd)
            for candidate in possible_file_extensions(cmd_file):
                if is_executable(candidate):
                    return candidate
    raise Exception('which failed to locate a proper command path "%s"' % cmd)