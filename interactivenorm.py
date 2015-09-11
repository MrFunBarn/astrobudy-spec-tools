#! /usr/bin/python
import math
import pickle
import copy
import time

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import numpy.ma as ma
from scipy.interpolate import InterpolatedUnivariateSpline
import pyfits



# InterativeNorm version 1.0
# Written by Brandon Bell, updated circa June 2015.
#
#     SpecNormalize is a class to define a plotting window and interactive
# fitting routines for multi-order flux vs wavelength spectra. The spectra
# are stored as 2D numpy arrays with a seprate array for each order stored
# as an elements in a list.
#
# Initialization:
# SpecNormalize only requires the spectra to initialize, but an object name,
# date of observation, wavelength units, and starting order can all be 
# specified as key workds. The spectra needs to passed to SpecNormalize
# as a list of 2D arrays (wavelength,flux) with each order/band as a seprate
# list element. 
# 
# Once normalized, you can save images that mimich the user interface for
# each and every order fitted, and text exports of both the normalized spectra
# and the the fit function for possible future use/evaluation. The porgram aloso
# creates pickles of current state when closed with Q, finished editing
# something or change the current order.

class SpecNormalize():
    """Interactivly fit the continuum of a spectra and normalize it."""
    def __init__(self,rawspec,objectn='',objectd='',units=r'$\AA$',order=0):
        self.pick = 0
        self.order = order
        self.rawspec = ma.masked_array(rawspec)
        self.objectn = objectn
        self.objectd = objectd
        self.units = units
        self.norm = copy.deepcopy(self.rawspec)
        self.fit = copy.deepcopy(self.rawspec)
        self.trimmed_spec = copy.deepcopy(self.rawspec)
        self.num_orders = len(self.rawspec)
        self.sm = [0] * self.num_orders
        self.fitpoints = [0] * self.num_orders
        self.spec_trim_points = [0] * self.num_orders
        self.fig1 = plt.figure(1)
        self.ax = plt.subplot2grid((3,1), (0, 0), rowspan=2)
        self.ax2 = self.fig1.add_subplot(3,1,3)
        self.cid = self.fig1.canvas.mpl_connect('button_press_event',
                                                 self._click_event)
        self.cid2 = self.fig1.canvas.mpl_connect('key_press_event',
                                                 self._key_press)
        fitted = [False] * self.num_orders
        trimmed  = [False] * self.num_orders
        w_smoothed = [False] * self.num_orders
        self.state = {'fitted': fitted, 'editting_fit': False,
                'w_smoothed': w_smoothed, 'trimmed': trimmed,
                'edditing_trim': False, 'trimming': False,
                'smoothed': False,
                'del_trim': False}


    def start_norm(self):
        """Starts the matplotlib window, nothing happens until start_norm
        is called."""
        self.base_draw()
        plt.show()


    def read_pickle(self):
        self.norm = ma.load( open(
                                self.objectn+'-'+self.objectd+'-norm.p', 'rb'))
        self.fit = ma.load( open(
                               self.objectn+'-'+self.objectd+'-fit.p', 'rb'))
        self.state = pickle.load( open(
                                  self.objectn+'-'+self.objectd+'-state.p',
                                  'rb'))
        self.fitpoints = pickle.load( open(self.objectn+'-'+
                                     self.objectd+'-fitpoints.p', 'rb'))
        if True in self.state['trimmed']:
            self.rawspec = ma.load( open(self.objectn+'-'+
                                     self.objectd+'-trimmed-spec.p', 'rb'))
            self.spec_trim_points = pickle.load( open(self.objectn+'-'+
                                     self.objectd+'-spec-trim-points.p', 'rb'))


    def write_pickle(self):
        # Pickle the current state of the of the data so that progres can be
        # saved.
        print('pickling')
        self.state['trimming'] = False
        self.state['editting_fit'] = False
        ma.dump( self.norm, open(self.objectn+'-'+self.objectd+'-norm.p',
                    'wb'))
        ma.dump( self.fit, open(self.objectn+'-'+self.objectd+'-fit.p',
                    'wb'))
        pickle.dump( self.state,
                    open(self.objectn+'-'+self.objectd+'-state.p', 'wb'))
        pickle.dump( self.fitpoints,
                    open(self.objectn+'-'+self.objectd+'-fitpoints.p',
                    'wb'))
        # The second part of the conditions is to prevent a partial trim
        # section from being saved wich could cause problems.
        if (True in self.state['trimmed'] and
                                self.state['trimming'] == False):
            ma.dump(self.rawspec, open(self.objectn+'-'+
                                     self.objectd+'-trimmed-spec.p', 'wb'))
            pickle.dump(self.spec_trim_points, open(self.objectn+'-'+
                                     self.objectd+'-spec-trim-points.p', 'wb'))



    def _click_event(self, event):
        """Call functions to handle mouse clicks based on mode"""
        nav = self.ax.get_navigate_mode()
        if ( nav == None and self.state['editting_fit'] == False and
             self.state['fitted'][self.order] == False and
             self.state['trimming'] == False):
            self.find_points(event)
        elif (nav == None and self.state['editting_fit'] == True and
                self.state['del_trim'] == False):
            self.edit_fit_points(event)
        elif (nav == None and self.state['editting_fit'] == True and
                self.state['del_trim'] == True):
                self.delete_trim_reg(event)
        elif nav == None and self.state['trimming'] == True:
            self._trim_spec(event)
            if (self.state['trimmed'][self.order] == True and
                    self.state['trimming'] == False):
                self.base_draw()
                self.fig1.canvas.draw()


    def _key_press(self, event):
        """Map keys to various methods/operations"""

        # If you are in edit mode (hitting e after a fit has been made with f)
        if (self.state['editting_fit'] == True and event.key in
            ['up','right','down','left','control','alt','D','A']):
            ystep = 0.05*(self.rawspec[self.order][self.pick,1])
            xstep = ( (self.rawspec[self.order][-1,0]-
                      self.rawspec[self.order][0,0]) /
                      len(self.rawspec[self.order][:,0]) ) # Finds 1px value.
            if event.key == 'up':
                self.fitpoints[self.order][self.pick,1] = \
                                 self.fitpoints[self.order][self.pick,1]+ystep
            elif event.key == 'right':
                self.fitpoints[self.order][self.pick,0] = \
                                 self.fitpoints[self.order][self.pick,0]+xstep
            elif event.key == 'down':
                self.fitpoints[self.order][self.pick,1] = \
                                 self.fitpoints[self.order][self.pick,1]-ystep
            elif event.key == 'left':
                self.fitpoints[self.order][self.pick,0] = \
                                 self.fitpoints[self.order][self.pick,0]-xstep
            elif event.key == 'alt' and self.pick < (len(self.fitpoints) - 1):
                self.pick += 1
            elif event.key == 'control' and self.pick > 0:
                self.pick -= 1
            elif event.key == 'D' and len(self.fitpoints[self.order][:,0]) > 4:
                print('Deleteing Point '+
                      str(self.fitpoints[self.order][self.pick,:]))
                print(self.fitpoints[self.order])
                print(self.pick)
                cut = np.array([self.pick])
                print(cut)
                self.fitpoints[self.order] = \
                           np.delete(self.fitpoints[self.order],
                                   cut, 0)
                self.fitpoints[self.order].view('f8,f8').sort(0, order=['f0'])
                # Reset to self.pick so that it maked sense and lays inside the
                # allowed range.
                if( (self.pick+1) <= np.shape(self.fitpoints[self.order])[0]
                   and self.pick > 1):
                    self.pick = self.pick - 1
                else:
                    self.pick = 0
            elif event.key == 'D' and len(self.fitpoints[self.order][:,0]) <=4:
                print("""Deletion of point would not leave enough points to
                      calculate fit. Add at least one more point before
                      deleting this one.""")
            elif event.key == 'A':
                self.fitpoints[self.order] = \
                                         np.vstack((self.fitpoints[self.order],
                                         [event.xdata, event.ydata]))
                self.fitpoints[self.order].view('f8,f8').sort(0,order=['f0'])
                # Make the Added Point the selected point for editting.
                self.pick = np.abs(self.fitpoints[self.order][:,0]-
                                   event.xdata).argmin()
                print('Adding Point '+
                      str(self.fitpoints[self.order][self.pick,:]))

            self.spline_fit_and_plot()
            self.ax.plot(self.fitpoints[self.order][self.pick,0],
                         self.fitpoints[self.order][self.pick,1], 'ro')
            self.ax2.plot(self.fitpoints[self.order][self.pick,0], 1, 'ro')
            self.fig1.canvas.draw()
            return

        # (Trim) Define a section of data to trim.
        if event.key == 'T':
            if self.state['editting_fit'] == False:
                self.state['trimming'] = True
                print("Select trim Region by picking the blue side first and "
                      "then the red side of the region.")
                self.base_draw()
                self.fig1.canvas.draw()
            elif self.state['editting_fit'] == True:
                self.state['del_trim'] = True
                self.base_draw()
                self.fig1.canvas.draw()

        # (Clear) the plot and the selected points.
        if event.key == 'C':
            self.fitpoints[self.order] = 0
            self.state['editting_fit'] = False
            self.state['fitted'][self.order] = False
            self.state['trimming'] = False
            self.state['trimmed'][self.order] = False
            self.spec_trim_points[self.order] = 0
            self.rawspec[self.order].mask = ma.nomask
            self.ax.cla()
            self.ax2.cla()
            self.base_draw()
            self.fig1.canvas.draw()

        # (Quit,quit) Safely disconect the canvas and close the figure or, if
        # editting the fit, leave editting mode and go back to browse/select
        # mode.
        if event.key == 'Q':
            if self.state['editting_fit'] == False:
                self.quit()
        if event.key ==  'q':
            self.state['editting_fit'] = False
            self.write_pickle()
            self.base_draw()
            self.spline_fit_and_plot()
            self.fig1.canvas.draw()

        # (Fit) a spline and redarw plots.
        if event.key == 'l' and self.state['editting_fit'] == False:
            self.state['fitted'][self.order] = True
            self.spline_fit_and_plot()
            self.fig1.canvas.draw()

        # (Edit) the allready fit data.
        if self.state['fitted'][self.order] == True and event.key == 'e':
            self.state['editting_fit'] = True
            self.base_draw()
            self.spline_fit_and_plot()
            self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                          str(self.order)+' Editting Fit')
            self.fig1.canvas.draw()

        # (Next) Advance to the next order.
        if event.key in ['n','up','right']:
            if self.state['editting_fit'] == False:
                if self.order < (self.num_orders-1):
                    self.order = self.order + 1
                    if self.state['w_smoothed'][self.order] == False:
                        self.smooth3
                    self.base_draw()
                    plt.draw()

        # (Previous) Go back an order.
        if event.key in ['p','down','left']:
            if self.state['editting_fit'] == False:
                if self.order >= 1:
                    self.order = self.order - 1
                    if self.state['w_smoothed'][self.order] == False:
                        self.smooth3
                    self.base_draw()
                    plt.draw()

        # Pickle work so far. write pickle is also called when an order is
        # fitted, so this is just to save progress in finding points etc...
        if event.key == 'w':
            self.write_pickle()

        # Smooth the spectra with a 3px boxcar. Doesn't actually change the
        # data, just makes a smoothed version an plots it. The norm preview
        # plot always shows un-smoothed data.
        if event.key == 'b':
            self.smooth3(self.rawspec[self.order])
            self.state['smoothed'] = True
            self.base_draw()
            if (self.state['editting_fit'] == True or
                self.state['fitted'][self.order] == True):
                self.spline_fit_and_plot()
            self.fig1.canvas.draw()

        # Un-Smooth the spectra. Just changes the state so that the un-smoothed
        # spectra is plotted again.
        if event.key == 'B' and self.state['smoothed'] == True:
            self.state['smoothed'] = False
            self.base_draw()
            if (self.state['editting_fit'] == True or
                self.state['fitted'][self.order] == True):
                self.spline_fit_and_plot()
            self.fig1.canvas.draw()

    def quit(self):
        """Cleanly quit the normalization processes."""
        self.write_pickle()
        self.fig1.canvas.mpl_disconnect(self.cid)
        self.fig1.canvas.mpl_disconnect(self.cid2)
        plt.close()


    def base_draw(self):
        """Sets up the basic plot with the un-normalized spectra and labels."""
        self.ax.cla()
        self.ax2.cla()
        self.redline = np.copy(self.rawspec[self.order])
        self.redline[:,1] = 1
        self.xmin = self.rawspec[self.order][0,0]
        self.xmax = np.max(self.rawspec[self.order][:,0])
        self.ax.set_ylabel('Flux')
        self.ax2.set_ylabel('Normalized Flux')
        self.ax2.set_xlabel('Wavelength ('+self.units+')')
        self.ax2.plot(self.redline[:,0], self.redline[:,1], color='green')
        self.ax2.set_ylim(0.5, 1.5)
        self.ax2.set_xlim(self.xmin, self.xmax)
        if self.state['smoothed'] == True:
            self.ax.plot(self.sm[self.order][:,0],
                         self.sm[self.order][:,1])
        elif ((self.state['trimming'] == True or
              self.state['editting_fit'] == True) and
              self.state['trimmed'][self.order] == True):
            self.ax.plot(self.rawspec.data[self.order][:,0],
                         self.rawspec.data[self.order][:,1])
            for r in range(np.shape(self.spec_trim_points[self.order])[0]):
                self.ax.axvline(self.spec_trim_points[self.order][r,0],
                                color='red')
                self.ax.axvline(self.spec_trim_points[self.order][r,1],
                                color='red')
                self.ax.axvspan(self.spec_trim_points[self.order][r,0],
                        self.spec_trim_points[self.order][r,1],
                                color='black',alpha=0.5)
        else:
            self.ax.plot(self.rawspec[self.order][:,0],
                         self.rawspec[self.order][:,1])
        if (self.state['fitted'][self.order] == True
            and self.state['editting_fit'] == False):
            self.ax.set_title(self.objectn+' '+self.objectd+' Order/Band-'+
                              str(self.order+1)+' Mode: Fitted/Browse')
            self.ax2.plot(self.norm[self.order][:,0],
                                  self.norm[self.order][:,1], color='blue')
            if type(self.fitpoints[self.order]) != int:
                self.ax2.plot(self.fitpoints[self.order][:,0],
                                 np.ones((len(self.fitpoints[self.order][:,1])
                                         )),'o', color='green')
                self.ax.plot(self.fitpoints[self.order][:,0],
                                self.fitpoints[self.order][:,1],'o',
                                color='green')
                self.ax.plot(self.rawspec[self.order][:,0],
                             self.fit[self.order][:,1],
                             color='green')
        elif self.state['editting_fit'] == True:
            self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                              str(self.order+1)+' Mode: Editting Fit Points')
            self.ax2.plot(self.norm[self.order][:,0],
                                  self.norm[self.order][:,1], color='blue')
            if self.state['del_trim'] == True:
                self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                              str(self.order+1)+' Mode: Editting Fit Points'
                                              '-> Del-Trim')
        else:
            self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                              str(self.order+1)+' Mode: Pick Points/Browse')
            if type(self.fitpoints[self.order]) != int:
                self.ax.plot(self.fitpoints[self.order][:,0],
                                self.fitpoints[self.order][:,1],'o',
                                color='green')
                self.ax.plot(self.fitpoints[self.order][:,0],
                             self.fitpoints[self.order][:,1],'o', color='red')
        if self.state['trimming'] == True:
            self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                              str(self.order+1)+' Mode: Trimming')


    def spline_fit_and_plot(self):
        """Fit a cubic-spline to the selected points and then plot the fit."""
        self.fitpoints[self.order].view('f8,f8').sort(0,order=['f0'])
        a = self.fitpoints[self.order][:,:]
        fit = InterpolatedUnivariateSpline(a[:,0],a[:,1],k=3)
        self.fit[self.order][:,1] = fit(self.rawspec[self.order][:,0])
        self.norm[self.order][:,1] = (self.rawspec[self.order][:,1] /
                                     self.fit[self.order][:,1])
        self.base_draw()
        self.ax.plot(self.rawspec[self.order][:,0], self.fit[self.order][:,1],
                     color='green')
        self.ax.plot(self.fitpoints[self.order][:,0],
                        self.fitpoints[self.order][:,1],'o', color='green')
        self.ax2.plot(self.fitpoints[self.order][:,0],
                         np.ones((len(self.fitpoints[self.order][:,1]))),
                         'o', color='green')


    def find_points(self, event):
        """Convert mouse clicks in the figure to (x,y) coords for fitting."""
        # Populates self.fitpoints with the selected points. Builds a list with
        # an array of x,y values at each order. This to accomidata the need
        # for differnt numbers of fit points for different orders/bands.
        if type(self.fitpoints[self.order]) == int:
            self.fitpoints[self.order] = np.array([[event.xdata, event.ydata]])
        else:
            self.fitpoints[self.order] = np.vstack((self.fitpoints[self.order],
                                                  [event.xdata, event.ydata]))
            self.fitpoints[self.order].view('f8,f8').sort(0,order=['f0'])
        self.base_draw()
        self.ax.plot(self.fitpoints[self.order][:,0],
                     self.fitpoints[self.order][:,1],'o', color='red')
        a = self.fitpoints[self.order][:,:]

        # Draws a preview spline fit. The if stament is here because a spline
        # needs at least [order of fit + 1] points to work and I've set the
        # order to 3 for a cubic spline. If you change the fit order, change
        # the conditional acordingly or the fit cause and error before there
        # are enough points for the order of the fit.
        if len(self.fitpoints[self.order]) > 3:
            fit = InterpolatedUnivariateSpline(a[:,0], a[:,1], k=3)
            am = a[-1,0]
            fr = np.where(self.rawspec[self.order][:,0] <= am)
            fr = fr[0][-1] + 1
            fitrange = self.rawspec[self.order][:fr,0]
            pfit = fit(fitrange)
            pfitn = self.rawspec[self.order][:fr,1] / pfit
            self.ax2.plot(self.rawspec[self.order][:fr,0],
                          pfitn, color='red')
            self.ax.plot(self.rawspec[self.order][:fr,0],
                         pfit, '--', color='red')
        event.canvas.draw()


    def edit_fit_points(self, event):
        x = event.xdata
        order = copy.deepcopy(self.order)
        points = copy.deepcopy(self.fitpoints)

        # Only finds points based on x values so that points can be selected in
        # either plot. This also prevents resolution problms between points at
        # similar y values.
        pickx = np.abs(points[order][:,0]-x).argmin()
        self.pick = pickx
        self.spline_fit_and_plot()
        mark_size = .1*np.amax(self.rawspec[self.order][:,1])
        self.ax.plot(self.rawspec[self.order][:,0], self.fit[self.order][:,1],
                     color='green')
        self.ax.plot(self.fitpoints[self.order][self.pick,0],
                     self.fitpoints[self.order][self.pick,1], 'ro')
        self.ax2.plot(self.fitpoints[self.order][self.pick,0], 1, 'ro')
        self.fig1.canvas.draw()


    def delete_trim_reg(self, event):
        """Delete the whole trim region nearest the click."""
        x = event.xdata
        pick = np.abs(self.spec_trim_points[self.order]-x).argmin()
        points = self.spec_trim_points[self.order]
        row = (pick-(pick%2))/2
        print('Deleting trim range '+str(
            self.spec_trim_points[self.order][row,0])+' to '+ str(
            self.spec_trim_points[self.order][row,1]))
        if np.shape(self.spec_trim_points[self.order])[0] >= 1:
            self.spec_trim_points[self.order] = \
                    np.delete(self.spec_trim_points[self.order],row,axis=0)
            print(self.spec_trim_points[self.order])

            self.rawspec[self.order] = ma.masked_array(
                                      self.rawspec[self.order].data, mask=False)
            self.norm[self.order] = ma.masked_array(
                                         self.norm[self.order].data, mask=False)
            self.fit[self.order] = ma.masked_array(
                                          self.fit[self.order].data, mask=False)
        else:
            self.spec_trim_points[self.order] = 0
        if np.shape(self.spec_trim_points[self.order])[0] == 0:
            self.spec_trim_points[self.order] = 0
            self.state['trimmed'][self.order] = False
        else:
            for i in range((np.shape(points)[0]) -1 ):
                start = self.spec_trim_points[self.order][i,0]
                end = self.spec_trim_points[self.order][i,1]
                self.rawspec[self.order] = ma.masked_inside(
                                            self.rawspec[self.order],start,end)
                self.fit[self.order] = ma.masked_inside(
                                                self.fit[self.order],start,end)
                self.norm[self.order] = ma.masked_inside(
                                               self.norm[self.order],start,end)
                self.sm[self.order] = ma.masked_inside(
                                                 self.sm[self.order],start,end)
        self.state['del_trim'] = False
        self.base_draw()
        self.spline_fit_and_plot()
        self.fig1.canvas.draw()
        return


    def _trim_spec(self, event):
        """Trim out the section of data between a set of x values."""
        x = event.xdata
        trimline = np.abs(self.rawspec[self.order][:,0]-x).argmin()
        # This section handels the case where there is currently no trim
        # sections defined for the data. This is detected because the list
        # that holds an array of trim points for each order/band is initally
        # set 0 at for a band/order. This Trimming algorithm is not especially
        # robust if something unexpected happens between defining the first and
        # second boundry of the region..
        if type(self.spec_trim_points[self.order]) == int:
            self.spec_trim_points[self.order] = np.array(
                                [[self.rawspec[self.order][trimline,0]]])
            self.base_draw
            self.ax.axvline(self.spec_trim_points[self.order][0,0],color='red')
            self.fig1.canvas.draw()
        elif (np.shape(self.spec_trim_points[self.order])[1] == 1
              and self.state['trimmed'][self.order] != True):
            self.spec_trim_points[self.order]= \
                      np.append(self.spec_trim_points[self.order],
                      [[self.rawspec[self.order][trimline,0]]], axis=1)
            self.base_draw
            self.ax.axvline(self.spec_trim_points[self.order][0,0],color='red')
            self.ax.axvline(self.spec_trim_points[self.order][0,1],color='red')
            self.ax.axvspan(self.spec_trim_points[self.order][0,0],
                            self.spec_trim_points[self.order][0,1],
                            color='black',alpha=0.5)
            self.fig1.canvas.draw()
            start = int(np.where(self.spec_trim_points[self.order][0,0] ==
                                self.rawspec[self.order][:,0])[0])
            start = self.rawspec[self.order][start,0]
            end = int(np.where(self.spec_trim_points[self.order][0,1] ==
                                self.rawspec[self.order][:,0])[0])
            end = self.rawspec[self.order][end,0]
            del_splice = np.array([start,end])
            self.rawspec[self.order] = ma.masked_inside(
                                            self.rawspec[self.order],start,end)
            self.fit[self.order] = ma.masked_inside(
                                                self.fit[self.order],start,end)
            self.norm[self.order] = ma.masked_inside(
                                               self.norm[self.order],start,end)
            self.sm[self.order] = ma.masked_inside(
                                               self.sm[self.order],start,end)

            T = str(np.shape(self.spec_trim_points[self.order])[0])
            print(T+' Sections Are Trimmed From Data. \n'
                    'Trimming is now finished. To trim another section press'
                    ' T again. \n')
            self.state['trimming'] = False
            self.state['trimmed'][self.order] = True
            return

        # This sections Handles the case where there is allready at least one
        # trim section in the data. The new section is built and append to the
        # list of regions.
        elif (self.state['trimmed'][self.order] == True and
            self.spec_trim_points[self.order][-1,1] != 0):
            new = np.array([self.rawspec[self.order][trimline,0],0])
            self.spec_trim_points[self.order] = np.vstack(
                                        (self.spec_trim_points[self.order],
                                        new))
            self.base_draw
            self.ax.axvline(self.spec_trim_points[self.order][-1,0],
                            color='red')
            self.fig1.canvas.draw()
        elif self.spec_trim_points[self.order][-1,1] == 0:
            self.spec_trim_points[self.order][-1,1] = \
                                           self.rawspec[self.order][trimline,0]
            self.base_draw
            self.ax.axvline(self.spec_trim_points[self.order][-1,0],
                            color='red')
            self.ax.axvline(self.spec_trim_points[self.order][-1,1],
                            color='red')
            self.ax.axvspan(self.spec_trim_points[self.order][-1,0],
                            self.spec_trim_points[self.order][-1,1],
                            color='black',alpha=0.5)
            self.fig1.canvas.draw()
            start = int(np.where(self.spec_trim_points[self.order][-1,0] ==
                                self.rawspec[self.order][:,0])[0])
            start = self.rawspec[self.order][start,0]
            end = int(np.where(self.spec_trim_points[self.order][-1,1] ==
                                self.rawspec[self.order][:,0])[0])
            end = self.rawspec[self.order][end,0]
            self.rawspec[self.order] = ma.masked_inside(
                                            self.rawspec[self.order],start,end)
            self.fit[self.order] = ma.masked_inside(
                                                self.fit[self.order],start,end)
            self.norm[self.order] = ma.masked_inside(
                                               self.norm[self.order],start,end)
            self.sm[self.order] = ma.masked_inside(
                                               self.sm[self.order],start,end)
            T = str(np.shape(self.spec_trim_points[self.order])[0])
            print(T+' Sections Are Trimmed From Data. \n'
                    'Trimming is now finished. To trim another section press'
                    ' T again. \n')
            self.state['trimming'] = False
            return

    def smooth3(self, x, smorder=False):
        # the order keyword and acompanying if statement allow the function be
        # called on a specific order (as in plotting) or on the current order
        # specified by state (as during interacvie normalization).
        if smorder == False: order=self.order
        else: order = smorder
        x = ma.compress_rows(x)
        nelm = np.shape(x)[0]
        self.sm[order] = np.zeros( (int(math.floor(nelm/3)),2) )
        sm = self.sm[order]
        for i in range(int(math.floor(nelm/3))):
          n=3*i
          sm[i,0] = ( x[n,0] + x[n+1,0] + x[n+2,0] ) / 3
          sm[i,1] = ( x[n,1] + x[n+1,1] + x[n+2,1] ) / 3
        #self.sm[order] = sm
        self.state['w_smoothed'][order] = True
        if self.state['trimmed'][order] == True:
            masked = np.empty(np.shape(self.sm[order]))
            masked.fill(False)
            for r in range(np.shape(self.spec_trim_points[order])[0]):
                start = self.spec_trim_points[order][r,0]
                end = self.spec_trim_points[order][r,1]
                mask = np.empty(np.shape(self.sm[order]))
                mask[:,0] = np.logical_and(self.sm[order][:,0] >= start,
                                           self.sm[order][:,0] <= end)
                mask[:,1] = np.logical_and(self.sm[order][:,0] >= start,
                                           self.sm[order][:,0] <= end)
                maskednew1 = np.where(mask[:,0] == True)
                maskednew2 = np.where(mask[:,1] == True)
                masked[maskednew1,0] = True
                masked[maskednew2,1] = True
            self.sm[order] = ma.masked_array(self.sm[order],
                                                      mask=masked)

    
    def norm_smooth3(self, x, smorder=False):
        # A version of smooth3x that is designed to simply return a smoothed
        # version of the passed array without, altering any execution state
        # variables or objects. Was made to solve a problem in whole_spec_plot
        # where an array seprate from self.sm was needed.
        if smorder == False: order=self.order
        else: order = smorder
        x = ma.compress_rows(x)
        nelm = np.shape(x)[0]
        sm = np.zeros( (int(math.floor(nelm/3)),2) )
        for i in range(int(math.floor(nelm/3))):
          n=3*i
          sm[i,0] = ( x[n,0] + x[n+1,0] + x[n+2,0] ) / 3
          sm[i,1] = ( x[n,1] + x[n+1,1] + x[n+2,1] ) / 3
        return sm
        if self.state['trimmed'][order] == True:
            masked = np.empty(np.shape(self.sm[order]))
            masked.fill(False)
            for r in range(np.shape(self.spec_trim_points[order])[0]):
                start = self.spec_trim_points[order][r,0]
                end = self.spec_trim_points[order][r,1]
                mask = np.empty(np.shape(self.sm[order]))
                mask[:,0] = np.logical_and(self.sm[order][:,0] >= start,
                                           self.sm[order][:,0] <= end)
                mask[:,1] = np.logical_and(self.sm[order][:,0] >= start,
                                           self.sm[order][:,0] <= end)
                maskednew1 = np.where(mask[:,0] == True)
                maskednew2 = np.where(mask[:,1] == True)
                masked[maskednew1,0] = True
                masked[maskednew2,1] = True
            #self.sm[order] = ma.masked_array(self.sm[order],
                            #                          mask=masked)
            smoo = ma.masked_array(self.sm[order],mask=masked)
            return smoo


    # Make plots of the rawdata fit and resulting normalized spectra for future
    # reference. Will only make a plot for orders that have been fit.
    def plot_order_fits(self):
        for order in range(self.num_orders):
            if type(self.fitpoints[order]) != int:
                normspecname = str(self.objectn+'-'+self.objectd+'-'+'Order'+str(order+1)+
                                '-Normalized-Data.dat')
                fitfilename = str(self.objectn+'-'+
                                self.objectd+'-'+'Order'+str(order+1)+'-Fit-Function.dat')
                
                # Make a line array for plotting on the normaliztion fit plots.
                redline = np.copy(self.norm[order,:,:])
                redline[:,1] = 1
        
                # Make a plot of the fit and un-normalized data to check the fit at
                # latter times if something seems fishy in the data.
                xmin = self.norm[order][0,0]
                xmax = self.norm[order][800,0]
                ymax = np.amax(self.rawspec[order][:,1])
            
                plt.clf()
                fig=plt.figure(1)
                fig.suptitle(self.objectn+' '+self.objectd+' Order '+str(order+1))
                plt.subplot(211)
                plt.ylabel('Flux')
                plt.title('Data and Fit')
                plt.plot(self.rawspec[order][:,0],self.rawspec[order][:,1])
                plt.plot(self.fit[order][:,0],self.fit[order][:,1])
                plt.plot(self.fitpoints[order][:,0], self.fitpoints[order][:,1], 'r+')
                plt.ylim(0,(1.25*ymax))
                plt.xlim(xmin, xmax)
                plt.subplot(212)
                plt.title('Normalized Spectra')
                plt.ylabel('Normalized Flux')
                plt.xlabel(r'Wavelength '+self.units)
                plt.plot(self.norm[order][:,0],self.norm[order][:,1])
                plt.plot(redline[:,0], redline[:,1])
                plt.plot(self.fitpoints[order][:,0],np.ones(len(self.fitpoints[order][:,0])), 'go')
                plt.ylim(0.5, 1.5)
                plt.xlim(xmin, xmax)
                pltfilename = (self.objectn+'-'+self.objectd+
                            'Reference-Data-Fit-'+'Order'+str(order+1)+'.png')
                plt.savefig(pltfilename, dpi=300)
                print('Saving Fit-Plot and data files for order '+str(order+1))
        return
    
    def whole_spec_plot(self,minor,major):
        print('Generating Plot of all orders/bands')
        minorLocator=MultipleLocator(minor)
        majorLocator=MultipleLocator(major)
        xmin = self.norm[0][0,0]
        xmax = np.amax(self.norm[:][-1,0])
        ymax = 1.5
        fig = plt.figure(figsize=(16, 2))
        ax = fig.add_subplot(1,1,1,autoscale_on=False,
                             ylim=(0,ymax),xlim=(xmin, xmax))
        plt.title(self.objectn+' '+self.objectd, fontsize=9)
        plt.xlabel(r'Wavelength '+self.units,fontsize=8)
        plt.ylabel('Flux', fontsize=8)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        for order in range(self.num_orders):
            if self.state['fitted'][order] == False: continue
            y = self.norm_smooth3(self.norm[order], smorder=order)
            middle_index = int(len(y[:,0]) // 2)
            ro = str(order + 1)
            if order%2==0:
                ax.plot(y[:,0],y[:,1],
                       color='green', linewidth=.07)
                ax.text(y[middle_index,0], .9, ro,
                        fontsize=4, color='green',
                        horizontalalignment='center',
                        va='center')
            else:
                ax.plot(y[:,0],y[:,1],
                        color='blue', linewidth=.07)
                ax.text(y[middle_index,0], 1.1, ro,
                        fontsize=4, color='blue',
                        horizontalalignment='center',
                        va='center')
        pltfilename = self.objectn+'-'+self.objectd+'-All-Data.png'
        plt.savefig(pltfilename, bbox_inches='tight', dpi=800)
        return

    def dump_text_files(self, header=False): 
    # write the normalized spectra and the fit function to a (wave flux) 
    # columb text file (seprate for each band/order.
        for order in range(self.num_orders):
            if type(self.fitpoints[order]) != int:
                normspecname = str(self.objectn+'-'+self.objectd+'-'+'Order'
                                   +str(order+1)+'-Normalized-Data.dat')
                fitfilename = str(self.objectn+'-'+
                                  self.objectd+'-'+'Order'+str(order+1)+
                                  '-Fit-Function.dat')
                np.savetxt(normspecname, self.norm[order][:,:])
                np.savetxt(fitfilename, self.fit[order][:,:])
