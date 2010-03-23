
import time, calendar

class ConfigAttributeError(AttributeError):
    pass

class Config:
    def __init__(self, base=None, **kwargs):
        self.base = base
        for k,v in kwargs.iteritems():
            setattr(self,k,v)
        
    def __getitem__(self, k):
        try:
            return getattr(self,k)
        except ConfigAttributeError:
            return '%('+k+')s'
    
    def __getattr__(self, k):
        if self.base is not None:
            return getattr(self.base, k)
        else:
            raise ConfigAttributeError(k)

    def __hasattr__(self, k):
        if self.base is not None:
            return hasattr(self.base, k)
        else:
            raise ConfigAttributeError(k)
    
    def has(self, k):
        return hasattr(self,k) and getattr(self,k) is not None
    
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
                    try:
                        return s % additional
                    except KeyError, e:
                        raise ConfigAttributeError(*e.args)
                else:
                    return s
            last = s
            
        raise Exception('Maximum number of replacements reached (recusive directory naming?)')

    def mktime(self, s):
        if isinstance(s, tuple):
            base, offset = s
        else:
            base, offset = s, 0
        
        if base == 'now':
            tbase = time.time()
        else:
            tbase = calendar.timegm(time.strptime(s, "%Y-%m-%d %H:%M:%S"))
        
        return tbase + offset
        
    def timerange(self, name):
        stbeg, stend = getattr(self, name)
        tbeg, tend = self.mktime(stbeg), self.mktime(stend)
        return (tbeg, tend)
        
        
    