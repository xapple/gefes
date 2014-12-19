# Futures #

# Built-in modules #

# Internal modules #

###############################################################################
class Collection(object):
    """A collection of aggregates."""

    def __repr__(self): return 'Collection: %s' % (self.children)
    def __iter__(self): return iter(self.children)
    def __len__(self): return len(self.children)
    def __add__(self, other): return self.__class__(self.children + other.children)

    def __init__(self, children):
        # Base attribute #
        self.children = children
        # Sort them #
        if all(hasattr(c, "num") for c in self.children): self.children.sort(key = lambda x: x.num)

    @property
    def first(self): return self.children[0]

    def __getitem__(self, key):
        if isinstance(key, basestring):
            return [c for c in self.children if c.name == key][0]
        elif isinstance(key, int):
            if hasattr(self.first, 'num'): return [c for c in self.children if c.num == key][0]
            else: return self.children[key]
        else: raise TypeError('key')
