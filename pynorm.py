import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pyfits
from scipy.interpolate import InterpolatedUnivariateSpline

class SpecNormalize():
    """Interactivly fit the continuum of a spectra and normalize it."""
    def __init__(self, rawspec=0, order=0):
        self.editting_fit = False 
        self.pick = 0     
        self.order = order
        self.rawspec, header = pyfits.getdata(rawspec, 0, header=True)
        self.objectn = header['object'].split('-')[0]
        self.objectd = header['UTSHUT'].split('T')[0]
        self.norm = np.copy(self.rawspec)
        self.fit = np.copy(self.rawspec)
        num_orders = int(np.shape(self.rawspec)[0])
        self.fitted = [False] * num_orders
        self.fitpoints = np.zeros((num_orders,50,2)) #max # fit points is 50.
        self.fig1 = plt.figure(1)
        self.ax = plt.subplot2grid((5,1), (0, 0), rowspan=4)
        self.ax2 = self.fig1.add_subplot(5,1,5)
        self.cid = self.fig1.canvas.mpl_connect('button_press_event',
                                                 self._click_event)
        self.cid2 = self.fig1.canvas.mpl_connect('key_press_event', 
                                                 self._key_press)
        self.base_draw()
        plt.show()
          
          
    def _click_event(self, event):
        """Call functions to handle mouse clicks based on context"""
        nav = self.ax.get_navigate_mode()
        if nav == None and self.editting_fit == False:
            self.find_points(event)      
        elif nav == None and self.editting_fit == True:
            self.edit_fit_points(event)
          
                    
    def _key_press(self, event):
        """Map key precess to various functions"""
        
        # Clear the plot and the selected points.
        if event.key == 'C':
            self.ax.cla()
            self.ax2.cla() 
            self.base_draw()
            self.fig1.canvas.draw()
            self.fitpoints[self.order,:,:] = 0
            self.fitted[self.order] = False
            
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
            self.spline_fit_and_plot()
            self.fig1.canvas.draw()
            self.fitted[self.order] = True
            
        # Edit the allready fit data.
        if self.fitted[self.order] == True and event.key == 'e':
            self.editting_fit = True
            self.base_draw()
            self.spline_fit_and_plot()
            self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                          str(self.order)+' Editting Fit')
            self.fig1.canvas.draw()
            
        
        # Advance to the next order. (Next) 
        if event.key == 'N':
            if self.editting_fit == False:
                if self.order < ((np.shape(self.rawspec)[0])-1):
                    self.order = self.order + 1
                    self.ipx = []
                    self.ipy = []
                    self.base_draw()
                    plt.draw()
        
        # Go back an order. (Previous)        
        if event.key == 'P':
            if self.editting_fit == False:
                if self.order >= 1:
                    self.ipx = []
                    self.ipy = []
                    self.order = self.order - 1
                    self.base_draw()
                    plt.draw()
                    self.editting_fit = False
                
                
        if (self.editting_fit == True and 
            event.key in ['up','right','down','left']):
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
            mark_size = .05*np.amax(self.rawspec[self.order,:,1])
            self.ax.scatter(self.fitpoints[self.order,self.pick,0],
                            self.fitpoints[self.order,self.pick,1],s=mark_size, 
                            marker='+', color='black' ,zorder = 10)
            self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                              str(self.order)+' Editting Fit')
            self.fig1.canvas.draw()
        
            
        
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
        self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                          str(self.order))
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
        # populates self.fitpoints with the selected points
        if self.fitpoints[self.order,0,0] == 0:
            self.fitpoints[self.order,0,0] = event.xdata
            self.fitpoints[self.order,0,1] = event.ydata
            pr = 0
        else: 
            pp = np.where(self.fitpoints[self.order,:,0] == 0)
            pr = int(np.amin(pp))
            self.fitpoints[self.order,pr,0] = event.xdata
            self.fitpoints[self.order,pr,1] = event.ydata
        self.base_draw()
        self.ax.plot(self.fitpoints[self.order,:(pr+1),0],
                     self.fitpoints[self.order,:(pr+1),1],'o', color='red')
        a = self.fitpoints[self.order,:,:]
        if pr >= 3:  # Spline needs at least order of fit + 1 points to work.
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
        y = event.ydata
        order = self.order
        points = self.fitpoints
        
        # Only find point based on x values so that points can be selected in
        # both plots. This also prevents resolution problms between points at 
        # similar y values.
        pickx = np.abs(points[order,:pr,0]-x).argmin()
        self.spline_fit_and_plot()
        mark_size = .05*np.amax(self.rawspec[self.order,:,1])
        self.ax.plot(self.rawspec[self.order,:,0], self.fit[self.order,:,1],
                     color='green')
        self.ax.scatter(self.fitpoints[self.order,pickx,0],
                        self.fitpoints[self.order,pickx,1], s=mark_size, 
                        marker='+', color='black' ,zorder = 10)
        self.ax2.scatter(self.fitpoints[self.order,pickx,0], 1,
                         s=mark_size, marker='+', color='black' ,zorder = 10)
        self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                          str(self.order)+' Editting Fit')
        self.fig1.canvas.draw()
        self.pick = pickx


if __name__ == '__main__':
    fits = SpecNormalize('HR-Car/HRCar-2013-02-19-Median.fits')

    

