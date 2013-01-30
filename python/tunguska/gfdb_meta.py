
from guts import *

class StringID(StringPattern):
    pattern = r'^[A-Za-z][A-Za-z0-9.-]*$'

class ScopeType(StringChoice):
    choices = [
            'global',
            'regional',
            'local',
            'undefined',
        ]


class WaveformType(StringChoice):
    choices = [
            'full waveform',
            'bodywave', 
            'P wave', 
            'S wave', 
            'surface wave',
            'undefined',
        ]


class GFType(StringChoice):
    choices = [
            'Kiwi-HDF',
        ]


class Citation(Object):
    id = StringID.T()
    type = String.T()
    title = Unicode.T()
    journal = Unicode.T(optional=True)
    volume = Unicode.T(optional=True)
    number = Unicode.T(optional=True)
    pages = Unicode.T(optional=True)
    year = String.T()
    note = Unicode.T(optional=True)
    issn = String.T(optional=True)
    doi = String.T(optional=True)
    url = String.T(optional=True)
    eprint = String.T(optional=True)
    authors = List.T(Unicode.T())
    publisher = Unicode.T(optional=True)
    keywords = Unicode.T(optional=True)
    note = Unicode.T(optional=True)
    abstract = Unicode.T(optional=True)
    
    @classmethod
    def from_bibtex(cls, filename=None, stream=None):
        
        from pybtex.database.input import bibtex

        parser = bibtex.Parser()

        if filename is not None:
            bib_data = parser.parse_file(filename)
        elif stream is not None:
            bib_data = parser.parse_stream(stream)

        citations = []

        for id_, entry in bib_data.entries.iteritems():
            d = {} 
            avail = entry.fields.keys()
            for prop in cls.T.properties:
                if prop.name == 'authors' or prop.name not in avail:
                    continue

                d[prop.name] = entry.fields[prop.name]

            if 'author' in entry.persons:
                d['authors'] = []
                for person in entry.persons['author']:
                    d['authors'].append(unicode(person))
            
            c = Citation(id=id_, type=entry.type, **d)
            citations.append(c)

        return citations

class EarthModel(Object):
    id = StringID.T()
    description = String.T(default='', optional=True)
    citation_ids = List.T(StringID.T())


class ModellingCode(Object):
    id = StringID.T()
    name = String.T(optional=True)
    version = String.T(optional=True)
    method = String.T(optional=True)
    author = Unicode.T(optional=True)
    author_email = String.T(optional=True)
    citation_ids = List.T(StringID.T())


class GFSet(Object):
    id = StringID.T()
    derived_from_id = StringID.T(optional=True)
    version = String.T(default='1.0', optional=True)
    author = Unicode.T(optional=True)
    author_email = String.T(optional=True)
    type = GFType.T(default='Kiwi-HDF', optional=True)
    modelling_code_id = StringID.T(optional=True)
    scope_type = ScopeType.T(default='undefined')
    waveform_type = WaveformType.T(default='undefined')
    can_interpolate_source = Bool.T(default=False)
    can_interpolate_receiver = Bool.T(default=False)
    frequency_min = Float.T(optional=True)
    frequency_max = Float.T(optional=True)
    sample_rate = Float.T(optional=True)
    size = Int.T(optional=True)
    citation_ids = List.T(StringID.T())
    description = String.T(default='', optional=True)


class GFSetTypeA(GFSet):
    '''Rotational symmetry, fixed receiver depth'''

    earth_model_id = StringID.T(optional=True)
    distance_min = Float.T()
    distance_max = Float.T()
    distance_delta = Float.T()
    source_depth_min = Float.T()
    source_depth_max = Float.T()
    source_depth_delta = Float.T()
    receiver_depth = Float.T(default=0.0)


class Inventory(Object):
    citations = List.T(Citation.T())
    earth_models = List.T(EarthModel.T())
    modelling_codes = List.T(ModellingCode.T())
    gf_sets = List.T(GFSet.T())



