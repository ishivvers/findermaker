"""
A finder chart maker.


To Do:
 - verify that astrometry.net upload works properly
 - calculate offsets
 - get it all to look nice and pretty
"""

import aplpy
import os
from astro import coord
from astrometryClient.client import Client as anClient
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import leastsq


class FinderMaker(object):
    
    def __init__(self, image=None, ra=None, dec=None):
        if image != None:
            if image.rsplit('.')[1].lower() != 'fits':
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
        raw_input('ready?')
        if (self.ra == self.dec == None):
            # have the user choose the object 
            print 'You must choose a target.'
            self.ra, self.dec = self.get_star()
        self.add_target( self.ra, self.dec, wcs=True, marker='a' )
        # have the user choose up to 3 offset stars
        # self.offset_stars = []
        # for i in range(3):
        #     if i == 0:
        #         print 'Choose the first offset star.'
        #     else:
        #         inn = raw_input('Hit enter to choose the next offset star, or type "d" to be done.')
        #         if 'd' in inn.lower():
        #             break
        #     ra,dec = self.get_star()
        #     self.add_target( ra, dec, wcs=True, marker=['s','c','d'][i] )
        #     self.offset_stars.append( [ra,dec] )
    
    def get_astrometry(self):
        """
        Interact with astrometry.net to get the WCS for our image.
        """
        # connect to astrometry.net
        supernova_key = 'jhvrmcmwgufgmsrw'
        supernova_url = 'http://supernova.astrometry.net/api/'
        nova_key = 'tugzsuwnbcykkeuy'
        nova_url = 'http://nova.astrometry.net/api/'
        new_image = self.image.rsplit('.')[0] + '.wcs.' + self.image.rsplit('.')[1]
        c = anClient(apiurl=nova_url)
        c.login(nova_key)
        
        # upload the image
        upres = c.upload(self.image)
        stat = upres['status']
        if stat != 'success':
            raise IOError('Upload failed: status %s\n %s\n' %(stat, str(upres)))
        self.anID = upres['subid']
        print 'upload successful. job id:',self.anID
        
        # wait for the response
        print 'waiting for response ...',
        while True:
            stat = c.sub_status(self.anID, justdict=True)
            if stat.get('status','') == 'success':
                break
            print '.',
        
        # download the new image
        url = nova_url.replace('api','new_fits_file/%i' %self.anID)
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
        self.fig = aplpy.FITSFigure( data=self.hdu )
        self.fig.show_grayscale(stretch='log')
        self.fig.show_grid()
    
    def _gauss2D(self, params, x, y):
        A, x0, y0, sigmaX, sigmaY, C = params
        return A*np.exp(-((x-x0)**2/(2*sigmaX**2) + (y-y0)**2/(2*sigmaY**2))) + C

    def _gauss2D_residuals(self, params, z, x, y):
        return z - self._gauss2D(params, x, y) 
        
    def get_star(self, cutout=10.0, error_plot=True):
        """
        Get a single star from the user interacting with the image.
        Star must be bright enough to fit a Gaussian to it!
        cutout is the region size to consider during fitting (in arcseconds).
        Returns best-fit coordinates (in WCS).
        If error_plot == True, will plot up cutout and the best-fit Gaussian.
        """
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
        inparams = [maximum, x,y, w/2, w/2, np.median(subarray)] # amp, x, y, sigx, sigy, constant
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
        
    def add_target(self, x, y, wcs=False, marker='s', w=0.005):
        """
        Marks the target (at pixel number x,y) on the plot.
        If WCS == True, x,y must be in WCS.
        Marker can be one of:
         's': a green square
         'c': a blue circle
         'a': two red arrows
         'd': an orange diamond
        w is the size scale for each marker
        """
        if wcs:
            ra,dec = x,y
        else:
            ra,dec = self.fig.pixel2world( x, y )
        if marker == 's':
            self.fig.show_rectangles( [ra], [dec], 2*w, 2*w, color='green', lw=3 )
            self.fig.show_markers( ra, dec, marker='+', c='green' )
        elif marker == 'c':
            self.fig.show_circles( [ra], [dec], w, color='blue', lw=3 )
            self.fig.show_markers( ra, dec, marker='+', c='blue' )
        elif marker == 'd':
            diamond = np.array([ [ra-w, dec], [ra, dec-w], [ra+w, dec], [ra, dec+w] ])
            self.fig.show_polygons( [diamond], color='orange', lw=3 )
            self.fig.show_markers( ra, dec, marker='+', c='orange' )
        elif marker == 'a':
            self.fig.show_arrows( ra+5*w, dec, -4*w, 0.0, color='red', lw=2 )
            self.fig.show_arrows( ra, dec+5*w, 0.0, -4*w, color='red', lw=2 )
            self.fig.show_markers( ra, dec, marker='+', c='red' )
        else:
            raise ValueError('marker must be one of [s,c,a,d]!')
    
    def add_offsets(self):
        """
        Marks the offset stars on the plot
        """
        # choose offset stars
        # fit for local Gaussian
        # get angular seperation (convert to arcseconds)
        # mark center on plot and add annotation
        pass