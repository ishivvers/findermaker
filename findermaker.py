"""
A finder chart maker.


Plan:

- can point to a fits file with or without WCS info
 - will use astrometry.net to get wcs if needed
- can click to define SN (if clearly detected) or can
 input coords (in any format)
- if don't have local image, must input exact coords, and will
 query online image servers to get background image
- asks you to interact with it to choose offset stars, 
 and will fit for the location of each
"""

import aplpy
from astro import coord
# astrometry.net client script
from client import Client as anClient
import pyfits as pf
import matplotlib.pyplot as plt

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
        c = anClient(api_url=nova_url)
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
    
    def get_dss_image(self):
        """
        Get a DSS finder chart, if not given input image
        """
        pass
    
    def get_object_coords(self):
        """
        Have the user interact with the image to fit for an object's coords
        Used to get target coords if needed, and offset star coords
        """
        pass
        
    def add_target(self):
        """
        Marks the target on the plot
        """
        pass
    
    def add_offsets(self):
        """
        Marks the offset stars on the plot
        """
        # choose offset stars
        # fit for local Gaussian
        # get angular seperation (convert to arcseconds)
        # mark center on plot and add annotation
        pass