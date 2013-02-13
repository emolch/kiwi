
import time, calendar, os

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
    
    def get_or_none(self, k):
        if self.has(k):
            return getattr(self,k)
        else:
            return None
        
    def get(self, k, default=None):
        if self.has(k):
            return getattr(self,k)
        else:
            return default
        
    def get_avail(self, *keys):
        d = {}
        for k in keys:
            if self.has(k):
                d[k] = getattr(self,k)
        return d
            
    def path(self, name, additional=None):
        return self.mkpath(getattr(self, name), additional)

    def path_or_none(self, name, additional=None):
        if not self.has(name):
            return None
        else:
            return self.path(name, additional=additional)
    
    def path_check_file(self, name, additional=None):
        p = self.mkpath(getattr(self, name), additional)
        if not os.path.isfile(p):
            raise Exception('No such file: %s' % p)
        return p
    
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
        
        
    
