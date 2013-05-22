#! /usr/bin/python
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pyfits
from scipy.interpolate import InterpolatedUnivariateSpline

class SpecNormalize():
    """Interactivly fit the continuum of a spectra and normalize it."""
    def __init__(self, rawspec=0, header=0, order=0):
        self.editting_fit = False
        self.pick = 0
        self.order = order
        self.rawspec = rawspec
        self.objectn = header['object'].split('-')[0]
        self.objectd = header['UTSHUT'].split('T')[0]
        self.norm = np.copy(self.rawspec)
        self.fit = np.copy(self.rawspec)
        num_orders = int(np.shape(self.rawspec)[0])
        self.fitted = [False] * num_orders
        self.fitpoints = np.zeros((num_orders,50,2)) # Max # fit points
        self.fig1 = plt.figure(1)
        self.ax = plt.subplot2grid((5,1), (0, 0), rowspan=4)
        self.ax2 = self.fig1.add_subplot(5,1,5)
        self.cid = self.fig1.canvas.mpl_connect('button_press_event',
                                                 self._click_event)
        self.cid2 = self.fig1.canvas.mpl_connect('key_press_event',
                                                 self._key_press)
        self.base_draw()
        plt.show()


    def read_pickle(self, rawspec):
        self.norm = pickle.dump( open(objectn+'-'+objectd+'-norm.p', 'rb'))
        self.fit = pickle.dump( open(objectn+'-'+objectd+'-fit.p', 'rb'))
        self.fitted = pickle.dump( open(objectn+'-'+objectd+'-fitted.p', 'rb'))
        self.fitpoints = pickle.dump( open(objectn+'-'+objectd+'-fitpoints.p',
                    'rb'))
        self.objectn = pickle.dump( open(objectn+'-'+objectd+'-objectn.p',
                    'rb'))
        self.objectd = pickle.dump( open(objectn+'-'+objectd+'-objectd.p',
                    'rb'))


    def write_pickle(self):
        print('pickling')
        pickle.dump( self.norm, open(self.objectn+'-'+self.objectd+'-norm.p', 'wb'))
        pickle.dump( self.fit, open(self.objectn+'-'+self.objectd+'-fit.p', 'wb'))
        pickle.dump( self.fitted, open(self.objectn+'-'+self.objectd+'-fitted.p', 'wb'))
        pickle.dump( self.fitpoints, open(self.objectn+'-'+self.objectd+'-fitpoints.p',
                    'wb'))
        pickle.dump( self.objectn, open(self.objectn+'-'+self.objectd+'-objectn.p',
                    'wb'))
        pickle.dump( self.objectd, open(self.objectn+'-'+self.objectd+'-objectd.p',
                    'wb'))


    def _click_event(self, event):
        """Call functions to handle mouse clicks based on mode"""
        nav = self.ax.get_navigate_mode()
        if ( nav == None and self.editting_fit == False and
             self.fitted[self.order] == False ):
            self.find_points(event)
        elif nav == None and self.editting_fit == True:
            self.edit_fit_points(event)


    def _key_press(self, event):
        """Map keys to various methods/operations"""

        # If you are in edit mode (hitting e after a fit has been made with f)
        if (self.editting_fit == True and event.key in
            ['up','right','down','left']):
            ystep = .005*np.mean(self.rawspec[self.order,:,1])
            xstep = 0.05
            if event.key == 'up':
                self.fitpoints[self.order,self.pick,1] = \
                                 self.fitpoints[self.order,self.pick,1] + ystep
            elif event.key == 'right':
                self.fitpoints[self.order,self.pick,0] = \
                                 self.fitpoints[self.order,self.pick,0] + xstep
            elif event.key == 'down':
                self.fitpoints[self.order,self.pick,1] = \
                                 self.fitpoints[self.order,self.pick,1] - ystep
            elif event.key == 'left':
                self.fitpoints[self.order,self.pick,0] = \
                                 self.fitpoints[self.order,self.pick,0] - xstep
            self.spline_fit_and_plot()
            mark_size = .1*np.amax(self.rawspec[self.order,:,1])
            self.ax.scatter(self.fitpoints[self.order,self.pick,0],
                            self.fitpoints[self.order,self.pick,1],s=mark_size,
                            marker='+', color='black', zorder=3)
            self.ax2.scatter(self.fitpoints[self.order,self.pick,0], 1,
                             s=mark_size, marker='+', color='black', zorder=3)
            self.fig1.canvas.draw()
            return

        # Clear the plot and the selected points.
        if event.key == 'C':
            self.ax.cla()
            self.ax2.cla()
            self.base_draw()
            self.fig1.canvas.draw()
            self.fitpoints[self.order,:,:] = 0
            self.fitted[self.order] = False

        # Safely disconect the canvas and close the figure or, if editting the
        # fit, leave editting mode and go back to browse/select mode.
        if event.key == 'q':
            if self.editting_fit == False:
                self.quit()
            else:
                self.editting_fit = False
                self.base_draw()
                self.spline_fit_and_plot()
                self.fig1.canvas.draw()

        # Fit a spline and redarw plots.
        if event.key == 'f':
            self.fitted[self.order] = True
            self.spline_fit_and_plot()
            self.fig1.canvas.draw()

        # Edit the allready fit data.
        if self.fitted[self.order] == True and event.key == 'e':
            self.editting_fit = True
            self.base_draw()
            self.spline_fit_and_plot()
            self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                          str(self.order)+' Editting Fit')
            self.fig1.canvas.draw()

        # (Next) Advance to the next order. Caps in Next and Previous are to
        # prevent unintended trigering.
        if event.key == 'N':
            if self.editting_fit == False:
                if self.order < ((np.shape(self.rawspec)[0])-1):
                    self.order = self.order + 1
                    self.ipx = []
                    self.ipy = []
                    self.base_draw()
                    plt.draw()

        # (Previous) Go back an order.
        if event.key == 'P':
            if self.editting_fit == False:
                if self.order >= 1:
                    self.ipx = []
                    self.ipy = []
                    self.order = self.order - 1
                    self.base_draw()
                    plt.draw()
                    self.editting_fit = False

        if event.key == 'W':
            print('Pickling')
            self.write_pickle()


    def quit(self):
        """Cleanly quit the normalization processes."""
        self.fig1.canvas.mpl_disconnect(self.cid)
        self.fig1.canvas.mpl_disconnect(self.cid2)
        plt.close()


    def base_draw(self):
        """Sets up the basic plot with the un-normalized spectra and labels."""
        self.ax.cla()
        self.ax2.cla()
        self.redline = np.copy(self.rawspec[self.order,:,:])
        self.redline[:,1] = 1
        self.xmin = self.rawspec[self.order,0,0]
        self.xmax = self.rawspec[self.order,800,0]
        if self.fitted[self.order] == True:
            self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                              str(self.order)+' Mode: Fitted/Browse')
        elif self.editting_fit == True:
            self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                              str(self.order)+' Mode: Editting Fit')
        else:
            self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                              str(self.order)+' Mode: Pick Points/Browse')
        self.ax.set_ylabel('Flux')
        self.ax2.set_ylabel('Normalized Flux')
        self.ax2.set_xlabel(r'Wavelength $(\AA)$')
        self.ax2.plot(self.redline[:,0], self.redline[:,1], color='green')
        self.ax2.set_ylim(0.5, 1.5)
        self.ax2.set_xlim(self.xmin, self.xmax)
        self.ax2.plot(self.norm[self.order,:,0,],
                              self.norm[self.order,:,1], color='blue')
        self.ax.plot(self.rawspec[self.order,:,0],self.rawspec[self.order,:,1])


    def spline_fit_and_plot(self):
        """Fit a cubic-spline to the selected points and then plot the fit."""
        pp = np.where(self.fitpoints[self.order,:,0] == 0)
        pr = int(np.amin(pp))
        a = self.fitpoints[self.order,:,:]
        fit = InterpolatedUnivariateSpline(a[:pr,0], a[:pr,1], k=3)
        self.fit[self.order,:,1] = fit(self.rawspec[self.order,:,0])
        self.base_draw()
        self.norm[self.order,:,1] = (self.rawspec[self.order,:,1] /
                                     self.fit[self.order,:,1])
        self.ax.plot(self.rawspec[self.order,:,0], self.fit[self.order,:,1],
                     color='green')
        self.ax.plot(self.fitpoints[self.order,:pr,0],
                        self.fitpoints[self.order,:pr,1],'o', color='green')
        self.ax2.plot(self.norm[self.order,:,0],self.norm[self.order,:,1],
                      color='blue')
        self.ax2.plot(self.fitpoints[self.order,:pr,0],
                         np.ones((len(self.fitpoints[self.order,:pr,1]))),
                         'o', color='green')


    def find_points(self, event):
        """Convert mouse clicks in the figure to (x,y) coords for fitting."""
        # Populates self.fitpoints with the selected points. It will place each
        # new point one index further down axis 1 as they are added. If you
        # need to fit more than 50 continuum points change the dimensionality
        # of self.fitpoints in __init__ from (num_orders,50,2) to (num_orders,
        # Your # of Points ,2).
        if self.fitpoints[self.order,0,0] == 0:
            self.fitpoints[self.order,0,0] = event.xdata
            self.fitpoints[self.order,0,1] = event.ydata
            pr = 0
        else:
            pp = np.where(self.fitpoints[self.order,:,0] == 0)
            pr = int(np.amin(pp))
            try:
                self.fitpoints[self.order,pr,0] = event.xdata
                self.fitpoints[self.order,pr,1] = event.ydata
            except IndexError:
                print("You've all ready picked the maximum # of points.")
        self.base_draw()
        self.ax.plot(self.fitpoints[self.order,:(pr+1),0],
                     self.fitpoints[self.order,:(pr+1),1],'o', color='red')
        a = self.fitpoints[self.order,:,:]

        # Draws a preview spline fit. The if stament is here because a spline
        # needs at least [order of fit + 1] points to work and I've set the
        # order to 3 for a cubic spline. If you change the fit order, change
        # the conditional acordingly or the fit cause and error before there
        # are enough points for the order of the fit.
        if pr >= 3:
            fit = InterpolatedUnivariateSpline(a[:(pr+1),0], a[:(pr+1),1], k=3)
            fr = np.where(self.rawspec[self.order,:,0]
                          <= int(np.amax(a[:(pr+1),0])))
            fr = int(np.amax(fr))
            fitrange = self.rawspec[self.order,:fr,0]
            pfit = fit(fitrange)
            pfitn = self.rawspec[self.order,:fr,1] / pfit
            self.ax2.plot(self.rawspec[self.order,:fr,0],
                          pfitn, color='red')
            self.ax.plot(self.rawspec[self.order,:fr,0],
                         pfit, '--', color='red')
        event.canvas.draw()

    def edit_fit_points(self, event):
        pp = np.where(self.fitpoints[self.order,:,0] == 0)
        pr = int(np.amin(pp))
        x = event.xdata
        order = self.order
        points = self.fitpoints

        # Only finds points based on x values so that points can be selected in
        # both plots. This also prevents resolution problms between points at
        # similar y values.
        pickx = np.abs(points[order,:pr,0]-x).argmin()
        self.spline_fit_and_plot()
        mark_size = .1*np.amax(self.rawspec[self.order,:,1])
        self.ax.plot(self.rawspec[self.order,:,0], self.fit[self.order,:,1],
                     color='green')
        self.ax.scatter(self.fitpoints[self.order,pickx,0],
                        self.fitpoints[self.order,pickx,1], s=mark_size,
                        marker='+', color='black', zorder=3)
        self.ax2.scatter(self.fitpoints[self.order,pickx,0], 1,
                         s=mark_size, marker='+', color='black' ,zorder=3)
        self.fig1.canvas.draw()
        self.pick = pickx


if __name__ == '__main__':
    rawspec, header = pyfits.getdata('HRCar-2013-02-19-Median.fits', 0,
                                      header=True)
    x = SpecNormalize(rawspec, header)
