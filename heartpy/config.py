'''
config file for heartpy
'''


#based on https://venngage.com/blog/color-blind-friendly-palette/#4

__all__ = ['get_colorpalette_poincare']

def init(): # pragma: no cover
    global colorblind 
    colorblind = False
    global colorblind_type 
    colorblind_type = 'deuteranopia'
    global color_style
    color_style = 'default'


def get_colorpalette_poincare():
    '''returns color palettes for poincare plotting

    Function that returns color palettes for poincare plotting.
    Takes arguments from config settings globals.

    Parameters
    ----------
    None

    Returns
    -------
    color_palette : list
        list conntaining color palette for poincare plot, in order
        of scatterplot, SD1 line, SD2 line.

    Examples
    --------
    >>> import heartpy as hp
    >>> hp.config.colorblind = False
    >>> palette = hp.config.get_colorpalette_poincare()
    >>> palette
    ['gray', 'blue', 'red']

    >>> hp.config.colorblind = True
    >>> hp.config.colorblind_type = 'protanopia'
    >>> palette = hp.config.get_colorpalette_poincare()
    >>> palette
    ['#4C4C9B', '#EBAFBE', '#DCDCC7']
    '''

    #poincare color palettes
    #palette color order: ['scatter', 'SD1', 'SD2']
    #identity line is always gray, ellipse always black
    poincare = {'regular': {'default': ['gray', 'blue', 'red'],
                            'retro': ['#63ACBE', '#601A4A', '#EE442F'],
                            'elegant': ['#ABC3C9', '#E0DCD3', '#CCBE9F'],
                            'corporate': ['#44749D', '#BDB8AD', '#EBE7E0'],
                            'zesty': ['#85C0F9', '#F5793A', '#A95AA1']
                            },
                'deuteranopia': {'default': ['#43439C', '#C7C78E', '#787811'],
                                 'retro': ['#9E9CC2', '#383745', '#A17724'],
                                 'elegant': ['#CAB8CB', '#F4D4D4', '#DCB69F'],
                                 'corporate': ['#636D97', '#BDB6AB', '#EDE6DE'],
                                 'zesty': ['#C59434', '#092C48', '#6F7498']
                                 },
                'protanopia': {'default': ['#4C4C9B', '#EBAFBE', '#DCDCC7'],
                               'retro': ['#9C9EB5', '#2A385B', '#8B7F47'],
                               'elegant': ['#BEBCC5', '#E2DAD1', '#C9BD9E'],
                               'corporate': ['#636D97', '#BDB6AB', '#EDE6DE'],
                               'zesty': ['#AE9C45', '#052955', '#6073B1']
                               },
                'tritanopia': {'default': ['#959595', '#46DBFF', '#DE2253'],
                               'retro': ['#6AAECF', '#9E3C50', '#DE2253'],
                               'elegant': ['#E1BE91', '#CD913C', '#78500F'],
                               'corporate': ['#256077', '#F8EAEC', '#E3FAFF'],
                               'zesty': ['#9E3C50', '#CD913C', '#46DBFF']
                               }
                }
    
    if colorblind:
        return poincare[colorblind_type.lower()][color_style.lower()]
    else:
        return poincare['regular'][color_style.lower()]


def get_colorpalette_plotter():
    '''returns color palettes for regular plotting
    
    Function that returns color palettes for regular plotting coloring.
    Takes arguments from config settings globals.
    
    Parameters
    ----------
    None

    Returns
    -------
    color_palette : list
        list conntaining color palette for plotter function, in order
        of line color, accepted peaks color, rejected peaks color.

    Examples
    --------
    >>> import heartpy as hp
    >>> hp.config.colorblind = False
    >>> palette = hp.config.get_colorpalette_plotter()
    >>> palette
    ['#7F7FFF', 'green', 'red']

    >>> hp.config.colorblind = True
    >>> hp.config.colorblind_type = 'protanopia'
    >>> palette = hp.config.get_colorpalette_plotter()
    >>> palette
    ['#4C4C9B', '#EBAFBE', '#DCDCC7']
    '''

    #plotter color palettes
    #color order: ['line color', 'accepted peaks color', 'rejected peaks color']
    plotter_colors = {'regular': {'default': ['#7F7FFF', 'green', 'red'],
                                  'retro': ['#601A4A', '#63ACBE', '#EE442F'],
                                  'elegant': ['#382119', '#70B8CA', '#CCBE9F'],
                                  'corporate': ['#93A7BA', '#44749D', '#CAAB68'],
                                  'zesty': ['#A95AA1', '#0F2080', '#F5793A']
                                  },
                      'deuteranopia': {'default': ['#43439C', '#C7C78E', '#787811'],
                                       'retro': ['#383745', '#9E9CC2', '#A17724'],
                                       'elegant': ['#342A1F', '#CAB8CB', '#DCB69F'],
                                       'corporate': ['#5D6E9E', '#CDB1AD', '#DECBE3'],
                                       'zesty': ['#C59434', '#092C48', '#6F7498']
                                        },
                      'protanopia': {'default': ['#4C4C9B', '#EBAFBE', '#DCDCC7'],
                                     'retro': ['#9C9EB5', '#2A385B', '#8B7F47'],
                                     'elegant': ['#2E2B21', '#C9BD9E', '#BEBCC5'],
                                     'corporate': ['#636D97', '#BDB6AB', '#D1D0DE'],
                                     'zesty': ['#AE9C45', '#052955', '#6073B1']
                                     },
                    'tritanopia': {'default': ['#959595', '#46DBFF', '#DE2253'],
                                   'retro': ['#6AAECF', '#9E3C50', '#DE2253'],
                                   'elegant': ['#E1BE91', '#78500F', '#CD913C'],
                                   'corporate': ['#256077', '#9AEBFD', '#F59AA7'],
                                   'zesty': ['#CD913C', '#46DBFF', '#9E3C50']
                                   }
                    }

    if colorblind:
        return plotter_colors[colorblind_type.lower()][color_style.lower()]
    else:
        return plotter_colors['regular'][color_style.lower()]