# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 14:10:41 2016

@author: Visscher
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 11:46:25 2016

@author: Visscher
"""

import wx
import numpy as np
import GUI_baseclass as guib
from helixmc.pose import HelixPose
from helixmc.random_step import RandomStepSimple
import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from helixmc import util
from helixmc import random, kBT
import FiberMC as fiber
import NucleosomeMC as nuc

#############################
#   DNA creation            #
#############################

def create(bpturn = 10.356, amp = 0.239, phase = 1.199, freq = 1.11, mod = 0.959,
           n_nuc = 2, NRL = 197, unwrap  = -20, n_ribs = 1, diameter = 330, \
           rise = 100, nld = 17):
    """
    Create a chromatin fiber (2 nucleosomes for now) with curved linker DNA
    defined by the parameters in
    
    Parameters
    ----------
    bpturn : float
        basepairs per full turn 
    amp : float
        amplitude of roll modulation
    phase : float
    freq: float
    mod : float
    
    returns
    -------
    DNA : HelixPose
    """
    a = np.array([bpturn,amp,phase,freq,mod])
    #create DNA and nucleosomes    
    n_bp = 147 + NRL
    DNA, dyads, nucl = fiber.create_dna(n_bp, n_nuc, NRL, unwrap = unwrap)
    start = []
    end = []
    link_len = []
    #curve linker dna
    pars1 = DNA.params
    for i in range(0, len(dyads)-1):
        end.append(dyads[1]-nucl.d_index)
        start.append(dyads[0]+(147-nucl.d_index))
        link_len.append(end[i] - start[i])
        pars1[start[i]+unwrap:end[i]-unwrap,:] = fiber.create_curved_linker(link_len[i]-2*unwrap,a)
    DNA.set_params(pars1)
    
    #Set fiber at right spot
    fiber_params = fiber.set_fiber_params(n_nuc, diameter = diameter, rise = rise, nld = nld, ribs = n_ribs)
    fiber_dyad_frames = fiber.get_fiber_dyad_frames(fiber_params, nucl.d_index, nucl.nuc_type + '.3DNA')
    pars = fiber.cast_nucs_fiber(DNA, dyads, fiber_dyad_frames, cut = start[0]-dyads[0]+unwrap)
    DNA.set_params(pars)
    return DNA, dyads, nucl, start, end
    

        
###########################################
#########   DYNAMIC UPDATE CLASS   ########
###########################################

plt.ion()
    
    
    
    
    
    
########################################   
######    CONTROL FRAME CLASS     ######
########################################

class ControlFrame(guib.ControlFrame):
    #constructor
    def __init__(self,parent):
        #initialize parent class
        guib.ControlFrame.__init__(self,parent)
        self.apars = np.array([10.4,0.0,0.0,1.0,1.0])
        self.n_nuc = 2
        self.NRL = 197
        self.unwrap = -20
        self.zoom = 20
        self.fiber_params =  fiber.set_fiber_params(self.n_nuc, diameter = 330.0, rise = 100.0, nld = 17)
        self.fig = plt.figure(0)
        self.ax = p3.Axes3D(self.fig)
        self.DNA, self.dyads, self.nucl, self.start, self.end = create(unwrap = self.unwrap)
        self.coord = self.DNA.coord/10.0
        self.cutcoords = self.DNA.coord[self.start[0]+self.unwrap]/10.0
        self.fiber_dyad_frames = fiber.get_fiber_dyad_frames(self.fiber_params, self.nucl.d_index, self.nucl.nuc_type + '.3DNA')
        self.CurveDNA()
        self.ChangePlot()    
    
    
    
    
    def GetCut(self):
        
        
        fiber.cast_nucs_fiber(self.DNA, self.dyads, self.fiber_dyad_frames, cut =self.start[0]-self.dyads[0]+self.unwrap)
    
    def Boltzmann_Energy(self, DNApars = "DNA_gau.npy"):
        '''
        get the bending energy of a dna basepair using the gaussian approximation
        
        Parameters
        ----------
        bp_step : int 
            number of the basepair
        DNApars : string
            name of a file containing a 7 by 6 matrix (upper row average DNA base-
            pair step parameters, lower 6 rows correlation matrix)
        
        Returns
        -------
        energy : float
            bending energy of the basepair
        '''
        file_working = util.locate_data_file(DNApars)
        params = np.load(file_working)
        avg = params[0]
        cov = params[1:]
        prec = np.linalg.inv(cov)
        return 0.5*kBT*np.dot((self.DNA.params[self.start[0]+self.unwrap]-avg),np.dot(prec,(self.DNA.params[self.start[0]+self.unwrap]-avg)))
            
    
    #Get distance of the cut
    def GetDist(self):
        coord = self.DNA.coord/10.0
        dx = coord[self.start[0]+1+self.unwrap,0]-self.coord[self.start[0]+self.unwrap,0]
        dy = coord[self.start[0]+1+self.unwrap,1]-self.coord[self.start[0]+self.unwrap,1]
        dz = coord[self.start[0]+1+self.unwrap,2]-self.coord[self.start[0]+self.unwrap,2]
        return np.sqrt(dx**2+dy**2+dz**2)
    
    #Get the inproduct of the two basepair normal vector of the cut(should be close to one for a good cut)
    def GetDot(self):
        A1 = self.DNA.frames[self.start[0]+self.unwrap]
        A2 = self.DNA.frames[1+self.start[0]+self.unwrap]
        z = np.array([0,0,1])
        return np.dot(np.dot(A1,z), np.dot(A2,z))
        

    def ChangePlot(self):
        """
        Update the matplot window
        """
        coord = self.DNA.coord/10.0
        
        self.ax.clear()
        #self.ax.plot3D(coord[self.dyads[0]:self.dyads[1], 0], coord[self.dyads[0]:self.dyads[1], 1], coord[self.dyads[0]:self.dyads[1], 2], 'k-')
        self.ax.plot3D(coord[:,0], coord[:, 1], coord[:, 2], 'k-')
        self.ax.scatter3D(coord[self.start[0]+self.unwrap:self.start[0]+2+self.unwrap, 0], coord[self.start[0]+self.unwrap:self.start[0]+2+self.unwrap, 1], coord[self.start[0]+self.unwrap:self.start[0]+2+self.unwrap, 2])
        self.ax.set_xlim3d(-self.zoom+self.cutcoords[0],self.zoom+self.cutcoords[0])
        self.ax.set_ylim3d(-self.zoom+self.cutcoords[1],self.zoom+self.cutcoords[1])
        self.ax.set_zlim3d(-self.zoom+self.cutcoords[2],self.zoom+self.cutcoords[2])
        self.fig.canvas.draw_idle()  
        
    def CurveDNA(self):
        """
        Curve the linker DNA
        """
        
        newpars = self.DNA.params
        for i in np.arange(0,len(self.dyads)-1):
            newpars[self.start[i]+self.unwrap:self.end[i]-self.unwrap,:] = fiber.create_curved_linker(self.end[i] - self.start[i]-2*self.unwrap,self.apars)
            self.DNA.set_params(newpars)
        self.GetCut()
        self.m_textDot.SetValue('%.3f' %self.GetDot())
        self.m_textDist.SetValue('%.2f' %self.GetDist())
        self.m_textShift.SetValue('%.3f' %self.DNA.params[self.start[0]+self.unwrap,0])
        self.m_textSlide.SetValue('%.3f' %self.DNA.params[self.start[0]+self.unwrap,1])
        self.m_textRise.SetValue('%.3f' %self.DNA.params[self.start[0]+self.unwrap,2])
        self.m_textEnergy.SetValue('%.3f' %self.Boltzmann_Energy())
        self.m_textTwist2.SetValue('%.3f' %self.DNA.params[self.start[0]+self.unwrap,5])

    def Scroll(self,parnumb):
        """
        Change variable upon scrolling event
        """
        sliders = [self.m_sliderTwist,self.m_sliderAmp,self.m_sliderPhase
                    ,self.m_sliderFreq,self.m_sliderMod]
        texts = [self.m_textBptwist, self.m_textAmp,self.m_textPhase,
                 self.m_textFreq,self.m_textMod]
        self.apars[parnumb] = sliders[parnumb].GetValue()/1000.0
        texts[parnumb].SetValue('%.3f' %self.apars[parnumb])

        self.CurveDNA()
        
        #updateplot
        self.ChangePlot()
    
    def TextChange(self,parnumb):
        """
        Change variables upon changing in the textboxes event
        """
        sliders = [self.m_sliderTwist,self.m_sliderAmp,self.m_sliderPhase
                    ,self.m_sliderFreq,self.m_sliderMod]
        texts = [self.m_textBptwist, self.m_textAmp,self.m_textPhase,
                 self.m_textFreq,self.m_textMod]
        
        self.apars[parnumb] = float(texts[parnumb].GetValue())
        sliders[parnumb].SetValue(int(self.apars[parnumb]*1000))
        
        self.CurveDNA()
        
        
        #updateplot
        self.ChangePlot()
    
    def PlusPress(self,parnumb):
        """
        Change variables upon pressing plus button events
        """
        sliders = [self.m_sliderTwist,self.m_sliderAmp,self.m_sliderPhase
                    ,self.m_sliderFreq,self.m_sliderMod]
        texts = [self.m_textBptwist, self.m_textAmp,self.m_textPhase,
                 self.m_textFreq,self.m_textMod]
                 
        if self.apars[parnumb]*1000 < sliders[parnumb].GetMax():
            self.apars[parnumb]+=.001
            sliders[parnumb].SetValue(self.apars[parnumb]*1000)
            texts[parnumb].SetValue('%.3f' %(self.apars[parnumb]))
        
            self.CurveDNA()
        
        
            #updateplot
            self.ChangePlot()

    
    def MinPress(self,parnumb):
        """
        Change variables upon pressing plus button events
        """
        sliders = [self.m_sliderTwist,self.m_sliderAmp,self.m_sliderPhase
                    ,self.m_sliderFreq,self.m_sliderMod]
        texts = [self.m_textBptwist, self.m_textAmp,self.m_textPhase,
                 self.m_textFreq,self.m_textMod]
                 
        if self.apars[parnumb]*1000 > sliders[parnumb].GetMin():
            self.apars[parnumb]-=.001
            sliders[parnumb].SetValue(self.apars[parnumb]*1000)
            texts[parnumb].SetValue('%.3f' %(self.apars[parnumb]))
            
            self.CurveDNA()
        
        
            #updateplot
            self.ChangePlot()
    
        
    #What to do when scrolled
    def TwistScroll(self,event):
        self.Scroll(0)

        
    def AmpScroll(self,event):
        self.Scroll(1)

    def PhaseScroll(self,event):
        self.Scroll(2)

    def FreqScroll(self,event):
        self.Scroll(3)
        
    def ModScroll(self,event):
        self.Scroll(4)
    
    #what to do when text is changed
    def ChangeTextTwist(self,event):
        self.TextChange(0)

    def ChangeTextAmp(self,event):
        self.TextChange(1)

    def ChangeTextPhase(self,event):
        self.TextChange(2)

    def ChangeTextFreq(self,event):
        self.TextChange(3)
        
    def ChangeTextMod(self,event):
        self.TextChange(4)
        
    #what to do when + or - is pressed
    
    def PlusTwist(self,event):
        self.PlusPress(0)
    
    def PlusAmp(self,event):
        self.PlusPress(1)

    def PlusPhase(self,event):
        self.PlusPress(2)
        
    def PlusFreq(self,event):
        self.PlusPress(3)
        
    def PlusMod(self,event):
        self.PlusPress(4)
        
    def MinTwist(self,event):
        self.MinPress(0)
        
    def MinAmp(self,event):
        self.MinPress(1)        
        
    def MinPhase(self,event):
        self.MinPress(2)

    def MinFreq(self,event):
        self.MinPress(3)
        
    def MinMod(self,event):
        self.MinPress(4)
        
    #what to do when create is clicked
    def Create(self,event):
        self.n_ribs = int(self.m_textRibs.GetValue())
        self.rise = float(self.m_textFiberRise.GetValue())
        self.nld = float(self.m_textNLD.GetValue())
        self.diameter = float(self.m_textDiameter.GetValue())
        
        self.NRL = int(self.m_textNRL.GetValue())
        self.n_nuc = int(self.m_textNucs.GetValue())
        self.unwrap = int(self.m_textUnwrap.GetValue())
        self.DNA, self.dyads, self.nucl, self.start, self.end = create(n_nuc = self.n_nuc, NRL = self.NRL, unwrap = self.unwrap, rise = self.rise, nld = self.nld, n_ribs = self.n_ribs)
        self.cutcoords = self.DNA.coord[self.start[0]+self.unwrap]/10.0
        self.fiber_params =  fiber.set_fiber_params(self.n_nuc, diameter = self.diameter, rise = self.rise, nld = self.nld, ribs = self.n_ribs)
        if self.m_checkStack.GetValue():
            self.fiber_params[0,5] = float(self.m_textFiberTwist.GetValue())
        self.fiber_dyad_frames = fiber.get_fiber_dyad_frames(self.fiber_params, self.nucl.d_index, self.nucl.nuc_type + '.3DNA')

    def ZoomScroll(self,event):
        scrollval = self.m_sliderZoom.GetValue()
        self.zoom = 20-.15*scrollval
        self.ax.set_xlim3d(-self.zoom+self.cutcoords[0],self.zoom+self.cutcoords[0])
        self.ax.set_ylim3d(-self.zoom+self.cutcoords[1],self.zoom+self.cutcoords[1])
        self.ax.set_zlim3d(-self.zoom+self.cutcoords[2],self.zoom+self.cutcoords[2])
        self.fig.canvas.draw_idle()
    
    def Reset(self,event):
        self.apars[0] = 10.4
        self.apars[1] = 0.0
        self.apars[2] = 0.0
        self.apars[3] = 1.0
        self.apars[4] = 1.0

        #update value in textbox
        self.m_textBptwist.SetValue('%.1f' %self.apars[0])
        self.m_textAmp.SetValue('%.1f' %self.apars[1])
        self.m_textFreq.SetValue('%.1f' %self.apars[2])
        self.m_textPhase.SetValue('%.1f' %self.apars[3])
        self.m_textMod.SetValue('%.1f' %self.apars[4])
        
        self.m_sliderTwist.SetValue(self.apars[0]*1000)
        self.m_sliderFreq.SetValue(self.apars[1]*1000)
        self.m_sliderMod.SetValue(self.apars[2]*1000)
        self.m_sliderPhase.SetValue(self.apars[3]*1000)
        self.m_sliderAmp.SetValue(self.apars[4]*1000)
        
        
        #curve DNA
        self.CurveDNA()
        
        #updateplot
        self.ChangePlot()
        
        
#mandatory in wx, create an app, False stands for not deteriction stdin/stdout
#refer manual for details
app = wx.App(False)
 
#create an object of ControlFrame
frame0 = ControlFrame(None)
#show the frame
frame0.Show(True)



#start the applications
app.MainLoop()