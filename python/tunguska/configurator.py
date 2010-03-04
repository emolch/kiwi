
class Config:
    def __init__(self, base=None, **kwargs):
        self.base = base
        for k,v in kwargs.iteritems():
            setattr(self,k,v)
        
    def __getitem__(self, k):
        try:
            return getattr(self,k)
        except AttributeError:
            return '%('+k+')s'
    
    def __getattr__(self, k):
        if self.base is not None:
            return getattr(self.base, k)
        else:
            raise AttributeError(k)

    def __hasattr__(self, k):
        if self.base is not None:
            return hasattr(self.base, k)
        else:
            raise AttributeError(k)
    
    def path(self, name, additional=None):
        return self.mkpath(getattr(self, name), additional)
    
    def mkpath(self, template, additional=None):
        last = None
        s = template
        maxrepl = 100
        for i in xrange(maxrepl):
            s = s % self
            if last == s:
                if additional is not None:
                    return s % additional
                else:
                    return s
            last = s
            
        raise Exception('Maximum number of replacements reached (recusive directory naming?)')
