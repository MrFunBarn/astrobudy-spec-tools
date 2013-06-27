#! /usr/bin/python
import math
import pickle
import copy

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pyfits
from scipy.interpolate import InterpolatedUnivariateSpline

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
                'smoothed': False}


    def start_norm(self):
        self.base_draw()
        plt.show()


    def read_pickle(self):
        self.norm = pickle.load( open(
                                self.objectn+'-'+self.objectd+'-norm.p', 'rb'))
        self.fit = pickle.load( open(
                               self.objectn+'-'+self.objectd+'-fit.p', 'rb'))
        self.state = pickle.load( open(
                                  self.objectn+'-'+self.objectd+'-state.p',
                                  'rb'))
        self.fitpoints = pickle.load( open(self.objectn+'-'+
                                     self.objectd+'-fitpoints.p', 'rb'))
        if True in self.state['trimmed']:
            self.trimmed_spec = pickle.load( open(self.objectn+'-'+
                                     self.objectd+'-trimmed-spec.p', 'rb'))
            self.spec_trim_points = pickle.load( open(self.objectn+'-'+
                                     self.objectd+'-spec-trim-points.p', 'rb'))


    def write_pickle(self):
        print('pickling')
        pickle.dump( self.norm, open(self.objectn+'-'+self.objectd+'-norm.p',
                    'wb'))
        pickle.dump( self.fit, open(self.objectn+'-'+self.objectd+'-fit.p',
                    'wb'))
        pickle.dump( self.state,
                    open(self.objectn+'-'+self.objectd+'-state.p', 'wb'))
        pickle.dump( self.fitpoints,
                    open(self.objectn+'-'+self.objectd+'-fitpoints.p',
                    'wb'))
        # The second part of the conditions is to prevent a partial trim
        # section from being saved wich could case problems.
        if (True in self.state['trimmed'] and
                                self.state['trimming'] == False):
            pickle.dump(self.trimmed_spec, open(self.objectn+'-'+
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
        elif nav == None and self.state['editting_fit'] == True:
            self.edit_fit_points(event)
        elif nav == None and self.state['trimming'] == True:
            self._trim_spec(event)


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

        if event.key == 'T':
            self.state['trimming'] = True
            print('trimming')

        # Clear the plot and the selected points.
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

        # Safely disconect the canvas and close the figure or, if editting the
        # fit, leave editting mode and go back to browse/select mode.
        if event.key == 'Q':
            if self.state['editting_fit'] == False:
                self.quit()
        if event.key ==  'q':
            self.state['editting_fit'] = False
            self.write_pickle()
            self.base_draw()
            self.spline_fit_and_plot()
            self.fig1.canvas.draw()

        # Fit a spline and redarw plots.
        if event.key == 'f' and self.state['editting_fit'] == False:
            self.state['fitted'][self.order] = True
            self.spline_fit_and_plot()
            self.fig1.canvas.draw()

        # Edit the allready fit data.
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
        if (self.state['smoothed'] == True ):
            self.ax.plot(self.sm[self.order][:,0],
                         self.sm[self.order][:,1])
        else:
            self.ax.plot(self.rawspec[self.order][:,0],
                         self.rawspec[self.order][:,1])
        if (self.state['fitted'][self.order] == True
            and self.state['editting_fit'] == False):
            self.ax.set_title(self.objectn+' '+self.objectd+' Order/Band-'+
                              str(self.order)+' Mode: Fitted/Browse')
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
                              str(self.order)+' Mode: Editting Fit Points')
            self.ax2.plot(self.norm[self.order][:,0],
                                  self.norm[self.order][:,1], color='blue')
        else:
            self.ax.set_title(self.objectn+' '+self.objectd+' Order-'+
                              str(self.order)+' Mode: Pick Points/Browse')
            if type(self.fitpoints[self.order]) != int:
                self.ax.plot(self.fitpoints[self.order][:,0],
                                self.fitpoints[self.order][:,1],'o',
                                color='green')
                self.ax.plot(self.fitpoints[self.order][:,0],
                             self.fitpoints[self.order][:,1],'o', color='red')


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
        # and array of x,y values at each order. This to accomidata the need
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
        # both plots. This also prevents resolution problms between points at
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


    def _trim_spec(self, event):
        x = event.xdata
        trimline = np.abs(self.rawspec[self.order][:,0]-x).argmin()
        if type(self.spec_trim_points[self.order]) == int:
            self.spec_trim_points[self.order] = np.array(
                                [[self.rawspec[self.order][trimline,0]]])
            self.base_draw
            self.ax.axvline(self.spec_trim_points[self.order][0,0],color='red')
            self.fig1.canvas.draw()
        elif np.shape(self.spec_trim_points[self.order])[0] == 1:
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
            print('Just Trimmed')
            self.state['trimming'] = False
            self.state['trimmed'][self.order] = True
        if self.state['trimmed'][self.order] == True:
            pass

    def smooth3(self, x):
        print('called')
        x = ma.compress_rows(x)
        nelm = np.shape(x)[0]
        self.sm[self.order] = np.zeros( (int(math.floor(nelm/3)),2) )
        sm = self.sm
        o = self.order
        for i in range(int(math.floor(nelm/3))):
          n=3*i
          sm[o][i,0] = ma.mean( x[n,0] + x[n+1,0] + x[n+2,0] )
          sm[o][i,1] = ma.mean( x[n,1] + x[n+1,1] + x[n+2,1] )
        self.sm = sm
        self.state['w_smoothed'][self.order] = True
        start = int(np.where(self.spec_trim_points[self.order][0,0] <=
                            self.rawspec[self.order][:,0])[0].argmin())
        start = self.rawspec[self.order][start,0]
        end = int(np.where(self.spec_trim_points[self.order][0,1] >=
                            self.rawspec[self.order][:,0])[0])
        end = self.rawspec[self.order][end,0]
        del_splice = np.array([start,end])
        print(del_splice)
        self.sm[self.order] = ma.masked_inside(
                                        self.sm[self.order],start,end)
        return
