# Built-in modules #
import sys, os, time, shutil, re, random

# Third party modules #
import sh

################################################################################
def average(iterator):
    """Iterative mean"""
    count = 0
    total = 0
    for num in iterator:
        count += 1
        total += num
    return float(total)/count

################################################################################
def wait(predicate, interval=1, message=lambda: "Waiting..."):
    ball, next_ball = u"|/-\\", "|"
    sys.stdout.write("    \033[K")
    sys.stdout.flush()
    while not predicate():
        time.sleep(1)
        next_ball = ball[(ball.index(next_ball) + 1) % len(ball)]
        sys.stdout.write("\r " + str(message()) + " " + next_ball + " \033[K")
        sys.stdout.flush()
    print "\r Done. \033[K"
    sys.stdout.flush()

################################################################################
def flatten(L):
    for sublist in L:
        if hasattr(sublist, '__iter__'):
            for item in flatten(sublist): yield item
        else: yield sublist

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
    path, name = os.path.split(cmd)
    if path:
        if is_executable(cmd): return cmd
    else:
        for path in os.environ['PATH'].split(os.pathsep):
            candidate = os.path.join(path, cmd)
            if is_executable(candidate): return candidate
    raise Exception('which failed to locate a proper command path "%s"' % cmd)

################################################################################
def tail(path, window=20):
    with open(path, 'r') as f:
        BUFSIZ = 1024
        f.seek(0, 2)
        num_bytes = f.tell()
        size = window + 1
        block = -1
        data = []
        while size > 0 and num_bytes > 0:
            if num_bytes - BUFSIZ > 0:
                # Seek back one whole BUFSIZ
                f.seek(block * BUFSIZ, 2)
                # Read BUFFER
                data.insert(0, f.read(BUFSIZ))
            else:
                # File too small, start from beginning
                f.seek(0,0)
                # Only read what was not read
                data.insert(0, f.read(num_bytes))
            linesFound = data[0].count('\n')
            size -= linesFound
            num_bytes -= BUFSIZ
            block -= 1
        return '\n'.join(''.join(data).splitlines()[-window:])

################################################################################
def head(path, window=20):
    with open(path, 'r') as handle:
        return ''.join(handle.next() for line in xrange(window))

################################################################################
def is_integer(string):
    try: int(string)
    except ValueError: return False
    return True
