"""
A finder chart maker.


To Do:
 - verify that astrometry.net upload works properly
"""

import aplpy
import os
import time
from subprocess import Popen, PIPE
import coord
#from astrometryClient.client import Client as anClient
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import leastsq

# a hackaround for a compass rose in aplpy
def initCompass(self, parent):
    self._ax1 = parent._ax1
    self._wcs = parent._wcs
    self.world2pixel = parent.world2pixel
    self.pixel2world = parent.pixel2world
    self._initialize_compass()
aplpy.overlays.Compass.__init__ = initCompass

# define some pretty colors
red = '#990000'
blue = '#0000FF'
green = '#006600'
orange = '#996600'

class FinderMaker(object):
    """
    A class to make finder charts.

    image: path to input image
    ra,dec: coordinates of object
    name: name of object
    diagnostics: if True, will plot diagnostic plots of the fit to each star's location

    NOTE: Either (ra, dec) OR image must be given (both is ok too).  Will use astrometry.net to find WCS solution
     for any image without one.  If (ra, dec) not given, the object of interest must be clearly discernable
     from the background in the given image.
    """    
    def __init__(self, image=None, ra=None, dec=None, name=None, diagnostics=False):
        self.name = name
        self.diagnostics = diagnostics
        self.annotations = []
        # make sure we have one of the allowed input combinations
        if image != None:
            if image.rsplit('.')[-1].lower() not in ['fit','fits']:
                raise ValueError('Image must be in fits format.')
            else:
                self.image = image
        else:
            self.image = None
        # parse input ra and dec (can be sexagesimal or decimal degrees)
        if (ra != None) & (dec != None):
            try:
                self.ra = float(ra)
                self.dec = float(dec)
            except:
                self.ra, self.dec = coord.s2dec( ra, dec )
        else:
            self.ra = None
            self.dec = None
        if (self.image == None) and any([v==None for v in [self.ra, self.dec]]):
            raise ValueError('Must include either a fits image or the target coordinates.')
        # run it!
        self.go()
    
    def go(self):
        if self.image == None:
            # get an image of the field
            self.get_dss_image()
        else:
            # make sure we have WCS info
            header = pf.open(self.image)[0].header
            if header.get('WCSAXES') == None:
                self.get_astrometry()
        self.build_plot()
        if (self.ra == self.dec == None):
            # have the user choose the object 
            raw_input( '\n\nHit enter and then identify the target.\n\n' )
            while True:
                res = self.get_star()
                if res != None:
                    self.ra, self.dec = res
                    break
                else:
                    print "\nThat didn't work. Try again.\n"
        self.add_object( self.ra, self.dec, wcs=True, marker='a' )
        self.add_offset_stars()
    
    def get_astrometry(self):
        """
        Interact with astrometry.net to get the WCS for our image, using the 
         astrometry client code.
        """
        # connect to astrometry.net
        supernova_key = 'jhvrmcmwgufgmsrw'
        supernova_url = 'http://supernova.astrometry.net/api/'
        nova_key = 'tugzsuwnbcykkeuy'
        nova_url = 'http://nova.astrometry.net/api/'
        new_image = self.image.replace('.fits','.wcs.fits')

        # The short routines below are pulled from astrometryClient/client.py in the __main__
        c = anClient(apiurl=nova_url)
        c.login(nova_key)
        
        # upload the image
        print '\n\nUploading image to astrometry.net\n\n'
        kwargs = {'publicly_visible': 'y', 'allow_modifications': 'd', 'allow_commercial_use': 'd'}
        upres = c.upload(self.image, **kwargs)
        stat = upres['status']
        if stat != 'success':
            raise IOError('Upload failed: status %s\n %s\n' %(str(stat), str(upres)))
        subID = upres['subid']
        print '\n\nUpload successful. Submission id:',subID,'\n\n'
        
        # Wait for the response
        while True:
            stat = c.sub_status(subID, justdict=True)
            jobs = stat.get('jobs', [])
            if len(jobs):
                for j in jobs:
                    if j is not None:
                        break
                if j is not None:
                    print '\n\nReceived job id',j,'\n\n'
                    jobID = j
                    break
            time.sleep(5)

        # wait for the calculation to finish
        success = False
        while True:
            stat = c.job_status(jobID, justdict=True)
            if stat.get('status','') in ['success']:
                success = (stat['status'] == 'success')
                break
            time.sleep(5)
        if not success:
            raise IOError('astrometry.net query failed: status %s'%str(stat))
        
        # download the new image
        print '\n\nGrabbing solved image\n\n'
        url = nova_url.replace('api','new_fits_file/%i' %jobID)
        cmd = 'wget -O %s %s' %(new_image, url)
        os.system( cmd )
        self.image = new_image
    
    def get_dss_image(self, name='finder_field.fits', size=10.0):
        """
        Get a DSS finder chart, if not given input image. Size is width in arcminutes.
        """
        url = "http://archive.stsci.edu/cgi-bin/dss_search?v=3&r=%.8f&d=%.8f$" %(self.ra, self.dec) +\
              "&h=%.2f&w=%.2f&f=fits&c=none&fov=NONE&e=J2000" %(size, size)
        cmd = 'wget "%s" -O %s' %(url, name)
        print 'Downloading image.'
        os.system( cmd )
        self.image = name
    
    def build_plot(self):
        """
        Creates the plot with which you interact.
        """
        self.hdu = pf.open( self.image )
        # assume, from here on out, that the first table is the relevant one
        self.fig = aplpy.FITSFigure( self.image )
        self.fig.show_grayscale(stretch='log', invert=True)
        self.fig.compass = aplpy.overlays.Compass(self.fig)
        self.fig.compass.show_compass(color='black', corner=2, length=0.1)
        self.fig.show_grid()
        self.fig.grid.set_color('k')
        self.fig.grid.set_alpha(0.25)
        if self.name != None:
            plt.title(self.name)
    
    def _gauss2D(self, params, x, y):
        A, x0, y0, sigmaX, sigmaY, C = params
        return A*np.exp(-((x-x0)**2/(2*sigmaX**2) + (y-y0)**2/(2*sigmaY**2))) + C

    def _gauss2D_residuals(self, params, z, x, y):
        return z - self._gauss2D(params, x, y) 
        
    def get_star(self, cutout=10.0, error_plot=None):
        """
        Get a single star from the user interacting with the image.
        Star must be bright enough to fit a Gaussian to it!
        cutout is the region size to consider during fitting (in arcseconds).
        Returns best-fit coordinates (in WCS).
        If error_plot == True, will plot up cutout and the best-fit Gaussian.
        """
        if error_plot == None:
            error_plot = self.diagnostics
        # convert arcseconds into degrees
        cutout = cutout/3600.0
        print "Click on the object."
        [(x,y)] = plt.ginput()  # in pixels
        # map the coutout size (in arcseconds) onto pixel size
        #  Use the fact that delta(Dec) = delta(true angle)
        ra, dec = self.fig.pixel2world(x,y)
        _,y1 = self.fig.world2pixel( ra, dec+cutout )
        _,y2 = self.fig.world2pixel( ra, dec )
        w = abs(y1-y2)
        # get a subarray to fit the Gaussian to
        [xmin,xmax, ymin,ymax] = map( lambda l: int(round(l)), 
                                      [x-w,x+w, y-w,y+w] )
        subarray = self.hdu[0].data[ymin:ymax, xmin:xmax]
        X = np.arange(xmin,xmax) # keeping track of the true pixel
        Y = np.arange(ymin,ymax) #  numbers of the subarray
        XX,YY = np.meshgrid( X,Y )
        # now fit a 2D Gaussian to it
        maximum = np.max(subarray)
        inparams = [maximum, x,y, w/5, w/5, np.median(subarray)] # amp, x, y, sigx, sigy, constant
        # need to make everything 1D
        X = XX.reshape( XX.shape[0]*XX.shape[1] )
        Y = YY.reshape( YY.shape[0]*YY.shape[1] )
        Z = subarray.reshape( subarray.shape[0]*subarray.shape[1] )
        fit_params, success = leastsq(self._gauss2D_residuals, inparams, args=(Z, X, Y))
        if not success:
            return None
        else:
            # I am pretty sure the indexing for fits files starts at 1, so need to
            #  add that in here to match Python index with fits image index
            ra, dec = self.fig.pixel2world( fit_params[1]+1.0, fit_params[2]+1.0 )
            if error_plot:
                curax = plt.gca()  # need to manage current axis variable
                fig = plt.figure()
                ax = Axes3D(fig)
                ax.plot_wireframe(XX,YY,subarray, color='k')
                ax.plot_wireframe(XX,YY, self._gauss2D(fit_params, XX, YY), color='r')
                ax.scatter3D( fit_params[1], fit_params[2], 
                              fit_params[0]+fit_params[5], color='r', s=40 )
                plt.xlabel('RA (px)')
                plt.ylabel('Dec (px)')
                plt.show()
                plt.sca( curax )
            return ra, dec
        
    def add_object(self, x, y, wcs=False, marker='s', s=10, center=None):
        """
        Marks the target (at pixel number x,y) on the plot.
        If WCS == True, x,y must be in WCS.
        Marker can be one of:
         's': a green square
         'c': a blue circle
         'a': two red arrows
         'd': an purple diamond
        s is the size scale for each marker (in pixels)
        if center is True, will place a small cross on the object's center
        """
        if center == None:
            center = self.diagnostics
        if wcs:
            ra,dec = x,y
            x,y = self.fig.world2pixel( ra, dec )
        else:
            ra,dec = self.fig.pixel2world( x, y )
        # find the width and height in RA,Dec coords
        r2,d2 = self.fig.pixel2world( x+s,y+s )
        w = np.abs(r2 - ra)
        h = np.abs(d2 - dec)
        # now make your mark
        if marker == 's':
            self.fig.show_rectangles( [ra], [dec], 2*h, 2*w, color=green, lw=3, layer='s1' )
            if center: self.fig.show_markers( ra, dec, marker='+', c=green, layer='s2' )
        elif marker == 'c':
            self.fig.show_circles( [ra], [dec], w, color=blue, lw=3, layer='c1' )
            if center: self.fig.show_markers( ra, dec, marker='+', c=blue, layer='c2' )
        elif marker == 'd':
            diamond = np.array([ [ra-w, dec], [ra, dec-h], [ra+w, dec], [ra, dec+h] ])
            self.fig.show_polygons( [diamond], color=orange, lw=3, layer='d1' )
            if center: self.fig.show_markers( ra, dec, marker='+', c=orange, layer='d2' )
        elif marker == 'a':
            self.fig.show_arrows( ra+5*w, dec, -4*w, 0.0, color=red, lw=2, layer='t1' )
            self.fig.show_arrows( ra, dec+5*h, 0.0, -4*h, color=red, lw=2, layer='t2' )
            if center: self.fig.show_markers( ra, dec, marker='+', c=red, layer='t3' )
        else:
            raise ValueError('marker must be one of [s,c,a,d]!')
    
    def add_offset_stars(self):
        """
        Have the user choose up to 3 offset stars and add
         them to the plot
        """
        self.remove_offset_stars()  # clear any old ones that may be present
        for i in range(3):
            if i == 0:
                sra,sdec = coord.dec2s(self.ra, self.dec)
                plt.annotate('RA:     %s\nDEC: %s'%(sra,sdec), (0.35, 0.9), xycoords='axes fraction',
                             size='large', weight='bold', ha='left', color=red)
                self.annotations.append( plt.annotate('From star\nto target:', (0.8, 0.9), 
                                         xycoords='axes fraction', weight='bold', color=red) )
                raw_input( '\n\nHit enter, then choose your first offset star.\n' )
            else:
                inn = raw_input('\nHit enter and choose another offset star, or type "q" or "d" to quit.\n')
                if 'q' in inn.lower() or 'd' in inn.lower():
                    return
            ra,dec = self.get_star()
            self.add_object( ra, dec, wcs=True, marker=['c','d','s'][i] )

            osRA, osDec = self.calc_distance(ra, dec)
            # figure out whether star is N,E,S,W
            #  The labels are direction from offset star to target!
            if ra>self.ra:
                rdir = 'W'
            else:
                rdir = 'E'
            if dec>self.dec:
                ddir = 'S'
            else:
                ddir = 'N'
            self.annotations.append( plt.annotate('RA: %.2f"%s \nDec: %.2f"%s' %(osRA, rdir, osDec, ddir),
                         (0.8, 0.8-0.08*i), xycoords='axes fraction', 
                         weight='bold', color=[blue,orange,green][i]) )

    def remove_offset_stars(self):
        """
        Removes all offsets stars from image.
        """
        for ann in self.annotations:
            ann.remove()
        self.annotations = []
        for layer in ['s1','s2','c1','c2','d1','d2']:
            try:
                self.fig.remove_layer( layer )
            except:
                # get an error if that one doesn't exist
                pass

    def calc_distance(self, ra, dec):
        """
        Calculates the distances (in arcseconds) between
         self.ra,self.dec and ra,dec, in each dimension.
        Returns delta(RA), delta(Dec)
        """
        dRA = coord.ang_sep(self.ra, self.dec, ra, self.dec)*3600.0
        dDec = coord.ang_sep(self.ra, self.dec, self.ra, dec)*3600.0
        return dRA, dDec

