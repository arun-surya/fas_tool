# -*- coding: utf-8 -*-

from numpy import arange, sin, pi
import matplotlib
from random import randint
import time
from matplotlib.backends.backend_wx import NavigationToolbar2Wx as NavigationToolbar
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import numpy
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import wx
import copy
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
from matplotlib.patches import Arrow
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.anchored_artists import AnchoredAuxTransformBox
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from mpl_toolkits.axes_grid.anchored_artists import AnchoredDrawingArea
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse
import random
from astropy.time import Time,TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u

def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence 
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circle in data unit.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence 
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)  
        `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
        norm, cmap, transform, etc.

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
    plt.colorbar()

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None
    if 'fc' in kwargs: kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs: kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs: kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs: kwargs.setdefault('linewidth', kwargs.pop('lw'))

    patches = [Circle((x_, y_), s_) for x_, y_, s_ in np.broadcast(x, y, s)]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        collection.set_array(np.asarray(c))
        collection.set_clim(vmin, vmax)
    return collection





class RedirectText(object):
    def __init__(self,aWxTextCtrl):
        self.out=aWxTextCtrl
 
    def write(self,string):
        self.out.WriteText(string)

#Matplotlib Class
class CanvasPanel(wx.Panel):
    def __init__(self, parent):
	global red_x,red_y,blue_x,blue_y
        wx.Panel.__init__(self, parent)
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.toolbar = NavigationToolbar(self.canvas)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
	self.sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Fit()
	self.flexscale=1000
	self.resscale=10000	
	self.pathscale=1000		


    def fibplot(self):
	global ra,dec,coll_x,coll_y,coll_rad,out,coll_stat_coll_x,rest,coll_y_rest,skyxvec,skyyvec,trstar_az,trstar_alt,trstar_ra,trstar_dec,rafld,decfld,alt,az
	#out = self.circles(coll_x, coll_y, coll_rad, alpha=0.7, fc='none',ec='black')
	self.axes.cla()
	starind=np.where(coll_stat=='Star')[0]
	skyind=np.where(coll_stat=='Sky')[0]
	noneind=np.where(coll_stat=='None')[0]	
	coll_x=np.array(coll_x)
	coll_y=np.array(coll_y)	
	if starind.size:
	    fib_star_x=coll_x[starind]
	    fib_star_y=coll_y[starind]
	if  skyind.size:
	    fib_sky_x=coll_x[skyind]
	    fib_sky_y=coll_y[skyind]	
	if  noneind.size:    
	    fib_none_x=coll_x[noneind]
	    fib_none_y=coll_y[noneind]	
	self.axes.scatter(az,alt,c='b')
	if starind.size:
	    self.axes.scatter(fib_star_x,fib_star_y,c='g',s=100,alpha=0.9)
	if skyind.size:
	    self.axes.scatter(fib_sky_x,fib_sky_y,c='b',s=100,alpha=0.5)
	if noneind.size:
	    self.axes.scatter(fib_none_x,fib_none_y,c='r',s=100,alpha=0.5)	
	self.axes.add_collection(out)
	self.axes.scatter(coll_x_rest,coll_y_rest,marker='+')
	box = AnchoredAuxTransformBox(self.axes.transData, loc=2,frameon=False)
        el = Arrow(0, 0, skyxvec, skyyvec, width=.001) 
        box.drawing_area.add_artist(el)
	self.axes.add_artist(box)
	
	self.axes.set_xlim(trstar_az-((rafld/2)+2*coll_rad[0]),trstar_az+((rafld/2)+2*coll_rad[0]))
	self.axes.set_ylim(trstar_alt-((decfld/2)+2*coll_rad[0]),trstar_alt+((decfld/2)+2*coll_rad[0]))
	self.canvas.draw()

class TabFour(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)

	font = wx.Font(15, wx.MODERN, wx.NORMAL, wx.NORMAL, False, u'Consolas')
	self.quote1 = wx.StaticText(self, label="Targets in Field: ", pos=(20, 30))	
	self.tottargets = wx.StaticText(self, label="      ", pos=(20, 30))	
	self.tottargets.SetFont(font)

	self.quote2 = wx.StaticText(self, label="Field Density  (targets/arcmin^2)    ", pos=(20, 30))	
	self.targetdensity = wx.StaticText(self, label="       ", pos=(20, 30))
	self.targetdensity.SetFont(font)

	self.quote8 = wx.StaticText(self, label="Fiber Density  (units/arcmin^2)    ", pos=(20, 30))	
	self.fiberdensity = wx.StaticText(self, label="       ", pos=(20, 30))
	self.fiberdensity.SetFont(font)


	self.quote3 = wx.StaticText(self, label="Fiber Units ", pos=(20, 30))	
	self.fiberunits = wx.StaticText(self, label=" ", pos=(20, 30))
	self.fiberunits.SetFont(font)
	
	self.quote4 = wx.StaticText(self, label="Observed Targets:  ", pos=(20, 30))	
	self.obstargets = wx.StaticText(self, label=" ", pos=(20, 30))	
	self.obstargets.SetFont(font)
	
	self.quote7 = wx.StaticText(self, label="Retention Percentage", pos=(20, 30))	
	self.ret = wx.StaticText(self, label="       ", pos=(20, 30))
	self.ret.SetFont(font)	
	
	self.quote5 = wx.StaticText(self, label="Completeness Percentage", pos=(20, 30))	
	self.comp = wx.StaticText(self, label=" ", pos=(20, 30))	
	font = wx.Font(15, wx.MODERN, wx.NORMAL, wx.NORMAL, False, u'Consolas')
	self.comp.SetFont(font)
	
	self.quote9 = wx.StaticText(self, label="", pos=(20, 30))
	
	self.quote6 = wx.StaticText(self, label="Efficiency percentage", pos=(20, 30))	
	self.eff = wx.StaticText(self, label=" ", pos=(20, 30))		
	self.eff.SetFont(font)

	
	sizerhor1 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor1.Add(self.quote1, 0, wx.ALL|wx.CENTER, 5)
	sizerhor1.Add(self.tottargets, 0, wx.ALL|wx.CENTER, 5)				
				
  
	sizerhor1.Add(self.quote2, 0, wx.ALL|wx.CENTER, 5)
	sizerhor1.Add(self.targetdensity, 0, wx.ALL|wx.CENTER, 5)									
	
	sizerhor3 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor3.Add(self.quote3, 0, wx.ALL|wx.CENTER, 5)
	sizerhor3.Add(self.fiberunits, 0, wx.ALL|wx.CENTER, 5)
	
    
	sizerhor3.Add(self.quote8, 0, wx.ALL|wx.CENTER, 5)
	sizerhor3.Add(self.fiberdensity, 0, wx.ALL|wx.CENTER, 5)	
	
	
	sizerhor4 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor4.Add(self.quote4, 0, wx.ALL|wx.CENTER, 5)
	sizerhor4.Add(self.obstargets, 0, wx.ALL|wx.CENTER, 5)
	
        
	sizerhor4.Add(self.quote7, 0, wx.ALL|wx.CENTER, 5)
	sizerhor4.Add(self.ret, 0, wx.ALL|wx.CENTER, 5)	
	
	sizerhor5 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor5.Add(self.quote5, 0, wx.ALL|wx.CENTER, 5)
	sizerhor5.Add(self.comp, 0, wx.ALL|wx.CENTER, 5)
						
	
	sizerhor6 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor6.Add(self.quote6, 0, wx.ALL|wx.CENTER, 5)
	sizerhor6.Add(self.eff, 0, wx.ALL|wx.CENTER, 5)
							 				                                   
	sizerver = wx.BoxSizer(wx.VERTICAL)	
	sizerver.Add(sizerhor1, 0, wx.ALL|wx.LEFT, 5)							
	#sizerver.Add(sizerhor2, 0, wx.ALL|wx.LEFT, 5)	
	sizerver.Add(sizerhor3, 0, wx.ALL|wx.LEFT, 5)							
	sizerver.Add(sizerhor4, 0, wx.ALL|wx.LEFT, 5)
	sizerver.Add(sizerhor5, 0, wx.ALL|wx.LEFT, 5)							
	sizerver.Add(sizerhor6, 0, wx.ALL|wx.LEFT, 5)	
	
	sizerver.Add(self.quote9, 0, wx.ALL|wx.LEFT, 5)	
	
        self.SetSizer(sizerver)
    def onGetSelection(self, event):
	print 'hello'
    def onShowComp(self, event):
	print 'hello'


class TabThree(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)

	
	self.quote10 = wx.StaticText(self, label="Allocation Mode: ", pos=(20, 30))	
	DropDownList = []
        Options = {"One-to-One","Sky-Nod"}
	for TextNode in Options:
             DropDownList.append( TextNode )
	self.configcmbx=wx.ComboBox(self,value="Select Mode",choices=DropDownList)
	self.quote1 = wx.StaticText(self, label="re-cfg time (min)", pos=(20, 30))
	self.recfg = wx.TextCtrl(self,value='30')
	
	self.selectBtn1 = wx.Button(self, label="Deploy")
        self.selectBtn1.Bind(wx.EVT_BUTTON, self.onGetSelection)
	
	self.selectBtn2 = wx.Button(self, label="Allocate")
        self.selectBtn2.Bind(wx.EVT_BUTTON, self.onGetSelection)
	
	self.selectBtn3 = wx.Button(self, label="Re-Configure")	
	
	
	sizerhor7 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor7.Add(self.selectBtn1, 0, wx.ALL|wx.CENTER, 5)
	sizerhor7.Add(self.selectBtn2, 0, wx.ALL|wx.CENTER, 5)				
	sizerhor7.Add(self.selectBtn3, 0, wx.ALL|wx.CENTER, 5)
				
	sizerhor8 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor8.Add(self.quote10, 0, wx.ALL|wx.CENTER, 5)
	sizerhor8.Add(self.configcmbx, 0, wx.ALL|wx.CENTER, 5)	
	sizerhor8.Add(self.quote1, 0, wx.ALL|wx.CENTER, 5)									
	sizerhor8.Add(self.recfg, 0, wx.ALL|wx.CENTER, 5)					
							 				                                   
	sizerver = wx.BoxSizer(wx.VERTICAL)	
	sizerver.Add(sizerhor8, 0, wx.ALL|wx.LEFT, 5)							
	sizerver.Add(sizerhor7, 0, wx.ALL|wx.CENTER, 5)	
	
        self.SetSizer(sizerver)
    def onGetSelection(self, event):
	print 'hello'
    def onShowComp(self, event):
	print 'hello'

	
class TabTwo(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
	
	lblList = ['Grid', 'Radial']
	self.rbox = wx.RadioBox(self, label = 'Choose Geometry', pos = (80,10), choices = lblList, majorDimension = 1, style = wx.RA_SPECIFY_ROWS) 
        self.rbox.Bind(wx.EVT_RADIOBOX,self.onRadioBox) 
	self.rbox.SetStringSelection('Radial')
	
	self.quote10 = wx.StaticText(self, label="Number of fibers: ", pos=(20, 30))	
	DropDownList = [' 1 ' , ' 7 ' , ' 19 ' , ' 37 ' , ' 61 ' , ' 91 ' , ' 127 ' , ' 169 ' , ' 217 ' , ' 271 ' , ' 331 ' , ' 397 ' , ' 469 ' , ' 547 ' , ' 631 ' , ' 721 ' , ' 817 ' , ' 919 ' , ' 1027 ' , ' 1141 ' , ' 1261 ' , ' 1387 ' , ' 1519 ' , ' 1657 ' , ' 1801 ' , ' 1951 ' , ' 2107 ' , ' 2269 ' , ' 2437 ' , ' 2611 ' , ' 2791 ' , ' 2977 ' , ' 3169 ' , ' 3367 ' , ' 3571 ' , ' 3781 ' , ' 3997 ' , ' 4219 ' , ' 4447 ' , ' 4681 ' , ' 4921 ' , ' 5167 ' , ' 5419 ' , ' 5677 '  ]
	self.fibnum=wx.ComboBox(self,value="721",choices=DropDownList)
	self.selectBtn1 = wx.Button(self, label="Scale to Field")
	
	self.quote6 = wx.StaticText(self, label="Fiber Units:                    Rows", pos=(20, 30))
	self.selectfibrows = wx.TextCtrl(self,value='11')
	self.selectBtn2 = wx.Button(self, label="Scale to Field")
	
	self.quote7 = wx.StaticText(self, label="Coloumns", pos=(20, 30))
	self.selectfibcols = wx.TextCtrl(self,value='11')	
	
	self.quote8 = wx.StaticText(self, label="Distance between units (arcseconds)", pos=(20, 30))
	self.selectfibdist = wx.TextCtrl(self,value='20')	
	
	self.quote9 = wx.StaticText(self, label="Patrol Rius (arcseconds)", pos=(20, 30))
	self.selectfibpat = wx.TextCtrl(self,value='15')	
	

	self.quote11 = wx.StaticText(self, label="Nod  Distance: (arcseconds)", pos=(20, 30))
	self.selectnoder = wx.TextCtrl(self,value='15')	
	self.quote12 = wx.StaticText(self, label="Direction (degree)", pos=(20, 30))	
	self.selectnodetheta = wx.TextCtrl(self,value='10')		
	
	self.quote6.Hide()
	self.selectfibrows.Hide()
	self.quote7.Hide()
	self.selectfibcols.Hide()
	self.fibnum.Show()
	self.quote10.Show()
	self.selectBtn1.Show()
	self.selectBtn2.Hide()

		
	sizerhor3 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor3.Add(self.quote10, 0, wx.ALL|wx.CENTER, 5)
	sizerhor3.Add(self.fibnum, 0, wx.ALL|wx.CENTER, 5)
	sizerhor3.Add(self.selectBtn1, 0, wx.ALL|wx.CENTER, 5)

	sizerhor4 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor4.Add(self.quote6, 0, wx.ALL|wx.CENTER, 5)
	sizerhor4.Add(self.selectfibrows, 0, wx.ALL|wx.CENTER, 5)	
	sizerhor4.Add(self.quote7, 0, wx.ALL|wx.CENTER, 5)
	sizerhor4.Add(self.selectfibcols, 0, wx.ALL|wx.CENTER, 5)
	sizerhor4.Add(self.selectBtn2, 0, wx.ALL|wx.CENTER, 5)
	
	
	
	sizerhor5 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor5.Add(self.quote8, 0, wx.ALL|wx.CENTER, 5)
	sizerhor5.Add(self.selectfibdist, 0, wx.ALL|wx.CENTER, 5)
	
	sizerhor6 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor6.Add(self.quote9, 0, wx.ALL|wx.CENTER, 5)
	sizerhor6.Add(self.selectfibpat, 0, wx.ALL|wx.CENTER, 5)	
									
					
	sizerhor9 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor9.Add(self.quote11, 0, wx.ALL|wx.CENTER, 5)
	sizerhor9.Add(self.selectnoder, 0, wx.ALL|wx.CENTER, 5)	
	sizerhor9.Add(self.quote12, 0, wx.ALL|wx.CENTER, 5)
	sizerhor9.Add(self.selectnodetheta, 0, wx.ALL|wx.CENTER, 5)					

							 				                                   
	sizerver = wx.BoxSizer(wx.VERTICAL)
	sizerver.Add(self.rbox, 0, wx.ALL|wx.LEFT, 5)
	sizerver.Add(sizerhor3, 0, wx.ALL|wx.LEFT, 5)
	sizerver.Add(sizerhor4, 0, wx.ALL|wx.LEFT, 5)
	sizerver.Add(sizerhor5, 0, wx.ALL|wx.LEFT, 5)
	sizerver.Add(sizerhor6, 0, wx.ALL|wx.LEFT, 5)
	sizerver.Add(sizerhor9, 0, wx.ALL|wx.LEFT, 5)						
	
	
        self.SetSizer(sizerver)
    def onGetSelection(self, event):
	print 'hello'
    def onRadioBox(self, event):
	rb = event.GetEventObject() 
        if rb.GetStringSelection() =='Grid':
	    self.quote6.Show()
            self.selectfibrows.Show()
            self.quote7.Show()
            self.selectfibcols.Show()
	    self.fibnum.Hide()
	    self.quote10.Hide()
	    self.selectBtn1.Hide()
	    self.selectBtn2.Show()
	    self.Layout()
	if rb.GetStringSelection() =='Radial':
	    self.quote6.Hide()
            self.selectfibrows.Hide()
            self.quote7.Hide()
            self.selectfibcols.Hide()
	    self.fibnum.Show()
	    self.quote10.Show()
	    self.selectBtn1.Show()
	    self.selectBtn2.Hide()	    
	    self.Layout()    
	
class TabOne(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)

	self.quote1 = wx.StaticText(self, label="Central Field:  RA (degree)", pos=(20, 30))
	self.selectcenra = wx.TextCtrl(self,value='150.20')

	self.quote2 = wx.StaticText(self, label="Dec (degree) ", pos=(20, 30))
	self.selectcendec = wx.TextCtrl(self,value='2.2')
  
	self.selectBtn1 = wx.Button(self, label="Randomize")

	self.quote3 = wx.StaticText(self, label="Field Size (arcmin):                ", pos=(20, 30))
	self.selectfieldsize = wx.TextCtrl(self,value='10')
	
	self.quote4 = wx.StaticText(self, label="Redshift                             low:", pos=(20, 30))
	self.selectlowz = wx.TextCtrl(self,value='2')
	
	self.quote5 = wx.StaticText(self, label="high:", pos=(20, 30))
	self.selecthighz = wx.TextCtrl(self,value='5')	

	self.quote13 = wx.StaticText(self, label="R Mag                                   high:", pos=(20, 30))
	self.selectrmaglow = wx.TextCtrl(self,value='19')
	
	self.quote14 = wx.StaticText(self, label="low:", pos=(20, 30))
	self.selectrmaghigh = wx.TextCtrl(self,value='25')
	
	lblList = ['Square', 'Circular']
	self.rbox = wx.RadioBox(self, label = 'Choose Field Type', pos = (80,10), choices = lblList, majorDimension = 1, style = wx.RA_SPECIFY_ROWS) 
	self.rbox.SetStringSelection('Circular')	
	sizerhor1 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor1.Add(self.quote1, 0, wx.ALL|wx.CENTER|wx.ALIGN_RIGHT, 5)
	sizerhor1.Add(self.selectcenra, 0, wx.ALL|wx.CENTER|wx.ALIGN_RIGHT, 5)	
	sizerhor1.Add(self.quote2, 0, wx.ALL|wx.CENTER|wx.ALIGN_RIGHT, 5)
	sizerhor1.Add(self.selectcendec, 0, wx.ALL|wx.CENTER|wx.ALIGN_RIGHT, 5)
	sizerhor1.Add(self.selectBtn1, 0, wx.ALL|wx.CENTER|wx.ALIGN_RIGHT, 5)

	sizerhor2 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor2.Add(self.quote3, 0, wx.ALL|wx.CENTER, 5)
	sizerhor2.Add(self.selectfieldsize, 0, wx.ALL|wx.CENTER, 5)	
	
	sizerhor3 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor3.Add(self.quote4, 0, wx.ALL|wx.CENTER, 5)
	sizerhor3.Add(self.selectlowz, 0, wx.ALL|wx.CENTER, 5)	
	sizerhor3.Add(self.quote5, 0, wx.ALL|wx.CENTER, 5)
	sizerhor3.Add(self.selecthighz, 0, wx.ALL|wx.CENTER, 5)	
			
	sizerhor10 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor10.Add(self.quote13, 0, wx.ALL|wx.CENTER, 5)
	sizerhor10.Add(self.selectrmaglow, 0, wx.ALL|wx.CENTER, 5)	
	sizerhor10.Add(self.quote14, 0, wx.ALL|wx.CENTER, 5)
	sizerhor10.Add(self.selectrmaghigh, 0, wx.ALL|wx.CENTER,5) 
							 				                                   
	sizerver = wx.BoxSizer(wx.VERTICAL)
	sizerver.Add(sizerhor1, 0, wx.ALL|wx.LEFT, 5)
	sizerver.Add(sizerhor2, 0, wx.ALL|wx.LEFT, 5)
	sizerver.Add(sizerhor3, 0, wx.ALL|wx.LEFT, 5)
	sizerver.Add(sizerhor10, 0, wx.ALL|wx.LEFT, 5)	
	sizerver.Add(self.rbox, 0, wx.ALL|wx.CENTER, 5)
	
        self.SetSizer(sizerver)
    def onGetSelection(self, event):
	print 'hello'
    def onShowComp(self, event):
	print 'hello'




class MyPanel(wx.Panel):
    """"""
 
    #----------------------------------------------------------------------
    def __init__(self, parent):
        """Constructor"""
        wx.Panel.__init__(self, parent)
    
	#log = wx.TextCtrl(self, wx.ID_ANY,style = wx.TE_MULTILINE|wx.TE_READONLY|wx.HSCROLL)
	#log.SetBackgroundColour((240,230,140))
	#log.SetForegroundColour((131,3,0))
	#font = wx.Font(15, wx.MODERN, wx.NORMAL, wx.NORMAL, False, u'Consolas')
	#log.SetFont(font)
	#redir=RedirectText(log)
	#sys.stdout=redir

        self.panell = CanvasPanel(self)
	nb1 = wx.Notebook(self)
	self.tab1 = TabOne(nb1)
	self.tab2 = TabTwo(nb1)
	self.tab3 = TabThree(self)
	self.tab4 = TabFour(self)
	nb1.AddPage(self.tab1, "Target Field")
	nb1.AddPage(self.tab2, "Fiber Units")
        sizerv = wx.BoxSizer(wx.VERTICAL)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
	#self.tab1.myGrid.Bind(gridlib.EVT_GRID_SELECT_CELL, self.onSingleSelect)
        #self.tab2.myGrid.Bind(gridlib.EVT_GRID_RANGE_SELECT, self.onDragSelection)
 	self.tab3.selectBtn1.Bind(wx.EVT_BUTTON, self.onCompButton)
	self.tab3.selectBtn2.Bind(wx.EVT_BUTTON, self.onResButton)
	self.tab3.selectBtn3.Bind(wx.EVT_BUTTON, self.onCfgButton)	
	self.tab1.selectBtn1.Bind(wx.EVT_BUTTON, self.onRandButton)
	self.tab2.selectBtn1.Bind(wx.EVT_BUTTON, self.onscalerad)
	self.tab2.selectBtn2.Bind(wx.EVT_BUTTON, self.onscalegrid)
	sizerv.Add(nb1, 2, wx.EXPAND)
	sizerv.Add(self.tab3, 1, wx.EXPAND)	
	sizerv.Add(self.tab4, 2, wx.EXPAND)	
        sizer.Add(sizerv, 1, wx.EXPAND)
	sizer.Add(self.panell, 1, wx.EXPAND)
        self.SetSizer(sizer)

    def onCfgButton(self, event):
	global trstar_ra,trstar_dec,trstar_az,trstar_alt,rafld,decfld,zlow,zhigh,coll_rows,coll_cols,coll_dist,pat_rad,ra,dec,coll_x,coll_y,coll_rad,circles,data,ind,ra,dec,out,coll_stat,coll_x_rest,coll_y_rest,skyxvec,skyyvec,az,alt,time,ind
	bear_mountain = EarthLocation(lat=19.8206*u.deg, lon=-155.468*u.deg, height=4000*u.m)
	exp=int(self.tab3.recfg.GetValue())*u.minute
	fieldtype=self.tab1.rbox.GetStringSelection()
        time = time + exp
	trstar_ra=float(self.tab1.selectcenra.GetValue())
	trstar_dec=float(self.tab1.selectcendec.GetValue())
	rafld=float(self.tab1.selectfieldsize.GetValue())* 0.0166667
	decfld=float(self.tab1.selectfieldsize.GetValue())* 0.0166667
	zlow=float(self.tab1.selectlowz.GetValue())
	zhigh=float(self.tab1.selecthighz.GetValue())
	rmaglow=float(self.tab1.selectrmaglow.GetValue())	
	rmaghigh=float(self.tab1.selectrmaghigh.GetValue())
	geometry=self.tab2.rbox.GetStringSelection()
	coll_num=int(self.tab2.fibnum.GetValue())	
	coll_rows=int(self.tab2.selectfibrows.GetValue())
	coll_cols=int(self.tab2.selectfibcols.GetValue())
	coll_dist=float(self.tab2.selectfibdist.GetValue())* 0.00027777833333
	pat_rad=float(self.tab2.selectfibpat.GetValue())* 0.00027777833333
	sky_r=float(self.tab2.selectnoder.GetValue())* 0.00027777833333
	sky_theta=float(self.tab2.selectnodetheta.GetValue())
	skyxvec=sky_r* (np.cos(np.deg2rad(sky_theta)))
	skyyvec=sky_r* (np.sin(np.deg2rad(sky_theta)))	
	coll_x=[]
	coll_y=[]
	coll_rad=[]
	
	fc=SkyCoord(ra=trstar_ra*u.degree, dec=trstar_dec*u.degree, frame='icrs')
	fc=fc.transform_to(AltAz(obstime=time,location=bear_mountain))
	trstar_alt=fc.alt.deg
	trstar_az=fc.az.deg
	

	stx=trstar_az-(round(coll_cols/2)*coll_dist)
	sty=trstar_alt-(round(coll_rows/2)*coll_dist)
	
	if geometry=='Grid':
	    
	    for i in range(coll_cols):
		    for j in range(coll_rows):
			    coll_x.append(stx+(i*coll_dist))
			    coll_y.append(sty+(j*coll_dist))
			    coll_rad.append(pat_rad)
	else:
	     #print 'in'
	     dt={1 : 1 , 7 : 2 , 19 : 3 , 37 : 4 , 61 : 5 , 91 : 6 , 127 : 7 , 169 : 8 , 217 : 9 , 271 : 10 , 331 : 11 , 397 : 12 , 469 : 13 , 547 : 14 , 631 : 15 , 721 : 16 , 817 : 17 , 919 : 18 , 1027 : 19 , 1141 : 20 , 1261 : 21 , 1387 : 22 , 1519 : 23 , 1657 : 24 , 1801 : 25 , 1951 : 26 , 2107 : 27 , 2269 : 28 , 2437 : 29 , 2611 : 30 , 2791 : 31 , 2977 : 32 , 3169 : 33 , 3367 : 34 , 3571 : 35 , 3781 : 36 , 3997 : 37 , 4219 : 38 , 4447 : 39 , 4681 : 40 , 4921 : 41 , 5167 : 42 , 5419 : 43 , 5677 : 44}
	     lp=dt[int(coll_num)]
	     coll_x=[trstar_az]
	     coll_y=[trstar_alt]
	     coll_rad=[pat_rad]
    
	     for k in range(1,lp):		    
		    num=int(float(360)/(float(60)/float(k)))
		    for i in range(num):
			    x=trstar_az+ coll_dist*k*np.cos(np.deg2rad(i*(float(60)/float(k))))
			    y=trstar_alt+ coll_dist*k*np.sin(np.deg2rad(i*(float(60)/float(k))))
			    #print x,y
			    coll_x.append(x)
			    coll_y.append(y)
			    coll_rad.append(pat_rad)
	c = SkyCoord(np.array(data.ra),np.array(data.dec), frame="icrs", unit="deg")
   	caltaz = c.transform_to(AltAz(obstime=time,location=bear_mountain))
	data_alt=caltaz.alt.deg 
	data_az=caltaz.az.deg 
	
	cent_dist=np.sqrt(np.square(data_az-trstar_az)+np.square(data_alt-trstar_alt))
	if fieldtype=='Square':
	      ind=np.where(np.logical_and(np.logical_and(np.logical_and(np.logical_and(data_az>(trstar_az-(rafld/2)) ,data_az<(trstar_az+(rafld/2))),np.logical_and(data_alt>(trstar_alt-(decfld/2)) ,data_alt<(trstar_alt+(decfld/2)))),np.logical_and(data.photoz_best>zlow, data.photoz_best<zhigh)),np.logical_and(data.rmag_psf>rmaglow, data.rmag_psf<rmaghigh)))[0] 
	else:
	      ind=np.where(np.logical_and((cent_dist<(rafld/2)),np.logical_and(np.logical_and(data.photoz_best>zlow, data.photoz_best<zhigh),np.logical_and(data.rmag_psf>rmaglow, data.rmag_psf<rmaghigh))))[0] 
	    
	coll_rad=np.array(coll_rad)
	coll_x_rest=np.array(copy.copy(coll_x))
        coll_y_rest=np.array(copy.copy(coll_y))
	out = circles(coll_x, coll_y, coll_rad, alpha=0.7, fc='none',ec='black')
	coll_stat= ["None" for x in range(len(coll_x))]
	coll_stat=np.array(coll_stat)
	ra=np.array(data.ra[ind])
	dec=np.array(data.dec[ind])
	alt=np.array(data_alt[ind])
	az=np.array(data_az[ind])
	self.tab4.tottargets.SetLabel(str(len(ra)))
	self.tab4.targetdensity.SetLabel('%2.2f    ' %(float(len(ra))/(np.pi*np.square(float(self.tab1.selectfieldsize.GetValue())/2))))	    
	self.tab4.fiberunits.SetLabel(str(len(coll_x)))
	out = circles(coll_x_rest, coll_y_rest, coll_rad, alpha=0.7, fc='none',ec='black')
	self.panell.fibplot()

    def onCompButton(self, event):
	global trstar_ra,allobjold,assign,oldstars,newstars,assignold,trstar_dec,trstar_az,trstar_alt,rafld,decfld,zlow,zhigh,coll_rows,coll_cols,coll_dist,pat_rad,ra,dec,coll_x,coll_y,coll_rad,circles,data,ind,ra,dec,out,coll_stat,coll_x_rest,coll_y_rest,skyxvec,skyyvec,az,alt,time,ind
	bear_mountain = EarthLocation(lat=19.8206*u.deg, lon=-155.468*u.deg, height=4000*u.m)
        utcoffset = -4*u.hour
	fieldtype=self.tab1.rbox.GetStringSelection()
        time = Time('2012-7-13 01:00:00') - utcoffset
	trstar_ra=float(self.tab1.selectcenra.GetValue())
	trstar_dec=float(self.tab1.selectcendec.GetValue())
	rafld=float(self.tab1.selectfieldsize.GetValue())* 0.0166667
	decfld=float(self.tab1.selectfieldsize.GetValue())* 0.0166667
	zlow=float(self.tab1.selectlowz.GetValue())
	zhigh=float(self.tab1.selecthighz.GetValue())
	rmaglow=float(self.tab1.selectrmaglow.GetValue())	
	rmaghigh=float(self.tab1.selectrmaghigh.GetValue())
	geometry=self.tab2.rbox.GetStringSelection()
	coll_num=int(self.tab2.fibnum.GetValue())	
	coll_rows=int(self.tab2.selectfibrows.GetValue())
	coll_cols=int(self.tab2.selectfibcols.GetValue())
	coll_dist=float(self.tab2.selectfibdist.GetValue())* 0.00027777833333
	pat_rad=float(self.tab2.selectfibpat.GetValue())* 0.00027777833333
	sky_r=float(self.tab2.selectnoder.GetValue())* 0.00027777833333
	sky_theta=float(self.tab2.selectnodetheta.GetValue())
	skyxvec=sky_r* (np.cos(np.deg2rad(sky_theta)))
	skyyvec=sky_r* (np.sin(np.deg2rad(sky_theta)))	
	coll_x=[]
	coll_y=[]
	coll_rad=[]
	allobjold=float('inf')
	fc=SkyCoord(ra=trstar_ra*u.degree, dec=trstar_dec*u.degree, frame='icrs')
	fc=fc.transform_to(AltAz(obstime=time,location=bear_mountain))
	trstar_alt=fc.alt.deg
	trstar_az=fc.az.deg
	

	stx=trstar_az-(round(coll_cols/2)*coll_dist)
	sty=trstar_alt-(round(coll_rows/2)*coll_dist)
	
	if geometry=='Grid':
	    
	    for i in range(coll_cols):
		    for j in range(coll_rows):
			    coll_x.append(stx+(i*coll_dist))
			    coll_y.append(sty+(j*coll_dist))
			    coll_rad.append(pat_rad)
	else:
	     #print 'in'
	     dt={1 : 1 , 7 : 2 , 19 : 3 , 37 : 4 , 61 : 5 , 91 : 6 , 127 : 7 , 169 : 8 , 217 : 9 , 271 : 10 , 331 : 11 , 397 : 12 , 469 : 13 , 547 : 14 , 631 : 15 , 721 : 16 , 817 : 17 , 919 : 18 , 1027 : 19 , 1141 : 20 , 1261 : 21 , 1387 : 22 , 1519 : 23 , 1657 : 24 , 1801 : 25 , 1951 : 26 , 2107 : 27 , 2269 : 28 , 2437 : 29 , 2611 : 30 , 2791 : 31 , 2977 : 32 , 3169 : 33 , 3367 : 34 , 3571 : 35 , 3781 : 36 , 3997 : 37 , 4219 : 38 , 4447 : 39 , 4681 : 40 , 4921 : 41 , 5167 : 42 , 5419 : 43 , 5677 : 44}
	     lp=dt[int(coll_num)]
	     coll_x=[trstar_az]
	     coll_y=[trstar_alt]
	     coll_rad=[pat_rad]
    
	     for k in range(1,lp):		    
		    num=int(float(360)/(float(60)/float(k)))
		    for i in range(num):
			    x=trstar_az+ coll_dist*k*np.cos(np.deg2rad(i*(float(60)/float(k))))
			    y=trstar_alt+ coll_dist*k*np.sin(np.deg2rad(i*(float(60)/float(k))))
			    #print x,y
			    coll_x.append(x)
			    coll_y.append(y)
			    coll_rad.append(pat_rad)
	c = SkyCoord(np.array(data.ra),np.array(data.dec), frame="icrs", unit="deg")
   	caltaz = c.transform_to(AltAz(obstime=time,location=bear_mountain))
	data_alt=caltaz.alt.deg 
	data_az=caltaz.az.deg 
	
	cent_dist=np.sqrt(np.square(data_az-trstar_az)+np.square(data_alt-trstar_alt))
	if fieldtype=='Square':
	      ind=np.where(np.logical_and(np.logical_and(np.logical_and(np.logical_and(data_az>(trstar_az-(rafld/2)) ,data_az<(trstar_az+(rafld/2))),np.logical_and(data_alt>(trstar_alt-(decfld/2)) ,data_alt<(trstar_alt+(decfld/2)))),np.logical_and(data.photoz_best>zlow, data.photoz_best<zhigh)),np.logical_and(data.rmag_psf>rmaglow, data.rmag_psf<rmaghigh)))[0] 
	else:
	      ind=np.where(np.logical_and((cent_dist<(rafld/2)),np.logical_and(np.logical_and(data.photoz_best>zlow, data.photoz_best<zhigh),np.logical_and(data.rmag_psf>rmaglow, data.rmag_psf<rmaghigh))))[0] 
	    
	coll_rad=np.array(coll_rad)
	coll_x_rest=np.array(copy.copy(coll_x))
        coll_y_rest=np.array(copy.copy(coll_y))
	out = circles(coll_x, coll_y, coll_rad, alpha=0.7, fc='none',ec='black')
	coll_stat= ["None" for x in range(len(coll_x))]
	coll_stat=np.array(coll_stat)
	ra=np.array(data.ra[ind])
	dec=np.array(data.dec[ind])
	alt=np.array(data_alt[ind])
	az=np.array(data_az[ind])
	assign=numpy.full((len(ra)),False,dtype=bool)
	assignold=numpy.full((len(ra)),False,dtype=bool)
	oldstars=[]
	newstars=[]
	if geometry=='Grid':
	    fibden=float(len(coll_x))/float(coll_rows*coll_dist*60*coll_cols*coll_dist*60)
	else:
	    dt={1 : 1 , 7 : 2 , 19 : 3 , 37 : 4 , 61 : 5 , 91 : 6 , 127 : 7 , 169 : 8 , 217 : 9 , 271 : 10 , 331 : 11 , 397 : 12 , 469 : 13 , 547 : 14 , 631 : 15 , 721 : 16 , 817 : 17 , 919 : 18 , 1027 : 19 , 1141 : 20 , 1261 : 21 , 1387 : 22 , 1519 : 23 , 1657 : 24 , 1801 : 25 , 1951 : 26 , 2107 : 27 , 2269 : 28 , 2437 : 29 , 2611 : 30 , 2791 : 31 , 2977 : 32 , 3169 : 33 , 3367 : 34 , 3571 : 35 , 3781 : 36 , 3997 : 37 , 4219 : 38 , 4447 : 39 , 4681 : 40 , 4921 : 41 , 5167 : 42 , 5419 : 43 , 5677 : 44}
	    lp=dt[int(coll_num)]
	    fibden=float(len(coll_x))/float(np.pi*np.square(coll_dist*60*lp))	
        #fibden=float(len(coll_x))/float(np.pi*np.square(float(self.tab1.selectfieldsize.GetValue())/2))
	self.tab4.tottargets.SetLabel(str(len(ra)))
	self.tab4.targetdensity.SetLabel('%2.2f    ' %(float(len(ra))/(np.pi*np.square(float(self.tab1.selectfieldsize.GetValue())/2))))	   
	self.tab4.fiberdensity.SetLabel('%2.2f    ' %(fibden))	 
	self.tab4.fiberunits.SetLabel(str(len(coll_x)))
	out = circles(coll_x_rest, coll_y_rest, coll_rad, alpha=0.7, fc='none',ec='black')
	self.panell.fibplot()
    def onResButton(self, event):
	global trstar_ra,allobjold,newstars,oldstars,ind,assign,assignold,trstar_dec,trstar_az,trstar_alt,rafld,decfld,zlow,zhigh,coll_rows,coll_cols,coll_dist,pat_rad,ra,dec,coll_x,coll_y,coll_rad,circles,data,ind,ra,dec,out,coll_stat,coll_x_rest,coll_y_rest,skyxvec,skyyvec,alt,az
	mode=self.tab3.configcmbx.GetValue()
	sky_r=float(self.tab2.selectnoder.GetValue())* 0.00027777833333
	sky_theta=float(self.tab2.selectnodetheta.GetValue())	
	coll_x=np.array(copy.copy(coll_x_rest))
        coll_y=np.array(copy.copy(coll_y_rest))
	assign=numpy.full((len(ra)),False,dtype=bool)	
	allobj=0
	coll_stat= ["None" for x in range(len(coll_x))]
	coll_stat=np.array(coll_stat)
	if mode=='Sky-Nod':
		for i in range(len(coll_x)):
		    if coll_stat[i] == 'None':
		        dist=np.sqrt(np.square((coll_x[i]-az))+np.square((coll_y[i]-alt)))			
			oldint=np.intersect1d(oldstars,ind)
			oldind=[np.where(ind==star)[0][0] for star in oldint]
			assignind=np.where(assign==True)
			dist[assignind]=float('inf')
			tarinds=np.where(dist<coll_rad[i])[0]	
			dist[oldind]=0	
			tardists=dist[tarinds]
			tarinds=[x for (y,x) in sorted(zip(tardists,tarinds))]
			for k in tarinds:
			    sky_x=az[k]+sky_r* (np.cos(np.deg2rad(sky_theta)))
			    sky_y=alt[k]+sky_r* (np.sin(np.deg2rad(sky_theta)))
			    skyxvec=sky_r* (np.cos(np.deg2rad(sky_theta)))
                            skyyvec=sky_r* (np.sin(np.deg2rad(sky_theta)))
			    dist_sky=np.sqrt(np.square((sky_x-coll_x))+np.square((sky_y-coll_y)))
			    collassignind=np.where(coll_stat != 'None')	
			    dist_sky[collassignind]=float('inf')
			    dist_sky[i]=float('inf')
			    distskyind=dist_sky.argmin()
			    dx_sky= abs(coll_x[distskyind]-sky_x)
			    dy_sky= abs(coll_y[distskyind]-sky_y)
			    if (np.square(dx_sky)+np.square(dy_sky))<np.square(coll_rad[distskyind]) :
				assign[k]=True
				coll_stat[i]='Star'
				coll_stat[distskyind]='Sky'
				coll_x[i]=az[k]				
				coll_y[i]=alt[k]
				coll_x[distskyind]=sky_x				
				coll_y[distskyind]=sky_y
			        allobj=allobj+1	   
				break 
		allfib=2*allobj
	else:	       
		for i in range(len(coll_x)):
		    dist=np.sqrt(np.square((coll_x[i]-az))+np.square((coll_y[i]-alt)))			
		    oldint=np.intersect1d(oldstars,ind)
		    oldind=[np.where(ind==star)[0][0] for star in oldint]
		    assignind=np.where(assign==True)
		    dist[assignind]=float('inf')
		    tarinds=np.where(dist<coll_rad[i])[0]	
		    dist[oldind]=0	
		    tardists=dist[tarinds]
		    tarinds=[x for (y,x) in sorted(zip(tardists,tarinds))]
		    if not not tarinds:
			distind=tarinds[0]
			coll_x[i]=az[distind]
			coll_y[i]=alt[distind]
			coll_stat[i]='Star'			
			assign[distind]=True
			allobj=allobj+1
		allfib=allobj
	newstars=ind[np.where(assign==True)]
        ret=len(np.intersect1d(newstars,oldstars))    
	if allobjold==float('inf'):
	    oldstars=copy.copy(newstars)
	    allobjold=allobj
	else:
	    oldstars=np.intersect1d(newstars,oldstars)    
	self.tab4.obstargets.SetLabel(str(allobj))
	self.tab4.eff.SetLabel('%2.2f %%' %((allfib/float(len(coll_x)))*100))
	self.tab4.comp.SetLabel('%2.2f %%' %((allobj/float(len(ra)))*100))
	self.tab4.ret.SetLabel('%2.2f %%' %((float(ret)/float(allobjold))*100))
	self.tab4.quote9.SetLabel('Time: '+str(time))
	out = circles(coll_x_rest, coll_y_rest, coll_rad, alpha=0.7, fc='none',ec='black')		
	self.panell.fibplot()    

    def onRandButton(self, event):
	global trstar_ra,trstar_dec,rafld,decfld,zlow,zhigh,coll_rows,coll_cols,coll_dist,pat_rad,ra,dec,coll_x,coll_y,coll_rad,circles,data,ind,ra,dec,out,coll_stat,coll_x_rest,coll_y_rest,skyxvec,skyyvec
	self.tab1.selectcenra.SetValue('%2.2f' %(random.uniform(ra.min(),ra.max())))
	self.tab1.selectcendec.SetValue('%2.2f' %(random.uniform(dec.min(),dec.max())))	   
    def onscalerad(self, event):
	global trstar_ra,trstar_dec,rafld,decfld,zlow,zhigh,coll_rows,coll_cols,coll_dist,pat_rad,ra,dec,coll_x,coll_y,coll_rad,circles,data,ind,ra,dec,out,coll_stat,coll_x_rest,coll_y_rest,skyxvec,skyyvec
	fld=float(self.tab1.selectfieldsize.GetValue())* 0.0166667
	d=float(self.tab2.selectfibdist.GetValue())* 0.00027777833333
	dt={1 : 1 , 7 : 2 , 19 : 3 , 37 : 4 , 61 : 5 , 91 : 6 , 127 : 7 , 169 : 8 , 217 : 9 , 271 : 10 , 331 : 11 , 397 : 12 , 469 : 13 , 547 : 14 , 631 : 15 , 721 : 16 , 817 : 17 , 919 : 18 , 1027 : 19 , 1141 : 20 , 1261 : 21 , 1387 : 22 , 1519 : 23 , 1657 : 24 , 1801 : 25 , 1951 : 26 , 2107 : 27 , 2269 : 28 , 2437 : 29 , 2611 : 30 , 2791 : 31 , 2977 : 32 , 3169 : 33 , 3367 : 34 , 3571 : 35 , 3781 : 36 , 3997 : 37 , 4219 : 38 , 4447 : 39 , 4681 : 40 , 4921 : 41 , 5167 : 42 , 5419 : 43 , 5677 : 44}
	dti={v:k for k, v in dt.iteritems()}
	num=int((fld/2)/d)
	self.tab2.fibnum.SetValue(str(dti[num+1]))	  
	
	 
    def onscalegrid(self, event):
	global trstar_ra,trstar_dec,rafld,decfld,zlow,zhigh,coll_rows,coll_cols,coll_dist,pat_rad,ra,dec,coll_x,coll_y,coll_rad,circles,data,ind,ra,dec,out,coll_stat,coll_x_rest,coll_y_rest,skyxvec,skyyvec
	
	self.tab1.selectcenra.SetValue('%2.2f' %(random.uniform(ra.min(),ra.max())))
	self.tab1.selectcendec.SetValue('%2.2f' %(random.uniform(dec.min(),dec.max())))	   

		   	   
	#self.panell.drawredfull()
	#self.panell2.drawbluefull()
       
    #----------------------------------------------------------------------
    
class MyFrame(wx.Frame):
    """"""
 
    #----------------------------------------------------------------------
    def __init__(self):
        """Constructor"""
        wx.Frame.__init__(self, parent=None, title="WFOS Fiber Allocation Simulation Tool", size=(1200, 800))
	
	#p = wx.Panel(self)
        self.panel = MyPanel(self)
        self.Show()
 
if __name__ == "__main__":
    app = wx.App(False)
    frame = MyFrame()
    
    
    data= pd.read_csv('cosmos.csv',sep=',')
    trstar_ra=35.0
    trstar_dec=-4.7
    trstar_az=0
    trstar_alt=0
    rafld=15 * 0.0166667
    decfld=15 * 0.0166667
    zlow=2
    zhigh=5    
    coll_rows=11
    coll_cols=11
    coll_dist=14* 0.00027777833333
    pat_rad= 14* 0.00027777833333
    skyxvec=0
    skyyvec=0
    alt=[]
    az=[]
    coll_x=[]
    coll_y=[]
    coll_rad=[]
    stx=trstar_ra-(round(coll_cols/2)*coll_dist)
    sty=trstar_dec-(round(coll_rows/2)*coll_dist)
    
    for i in range(coll_cols):
	    for j in range(coll_rows):
		    coll_x.append(stx+(i*coll_dist))
		    coll_y.append(sty+(j*coll_dist))
		    coll_rad.append(pat_rad)
    coll_stat= ["None" for x in range(len(coll_x))]
    coll_stat=np.array(coll_stat)
    ind=np.where(np.logical_and(np.logical_and(np.logical_and(data.ra>(trstar_ra-(rafld/2)) ,data.ra<(trstar_ra+(rafld/2))),np.logical_and(data.dec>(trstar_dec-(decfld/2)) ,data.dec<(trstar_dec+(decfld/2)))),np.logical_and(data.photoz_best>zlow, data.photoz_best<zhigh)))[0]
    out = circles(coll_x, coll_y, coll_rad, alpha=0.7, fc='none',ec='black')
    coll_x_rest=np.array(copy.copy(coll_x))
    coll_y_rest=np.array(copy.copy(coll_y))
    ra=np.array(data.ra[ind])
    dec=np.array(data.dec[ind])
    app.MainLoop() 






