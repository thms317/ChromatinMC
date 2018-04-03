# -*- coding: utf-8 -*- 

###########################################################################
## Python code generated with wxFormBuilder (version Jun 17 2015)
## http://www.wxformbuilder.org/
##
## PLEASE DO "NOT" EDIT THIS FILE!
###########################################################################

import wx
import wx.xrc

###########################################################################
## Class ControlFrame
###########################################################################

class ControlFrame ( wx.Frame ):
	
	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = wx.EmptyString, pos = wx.DefaultPosition, size = wx.Size( 1167,697 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		
		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		
		bSizer6 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer5 = wx.BoxSizer( wx.HORIZONTAL )
		
		fgSizer2 = wx.FlexGridSizer( 0, 5, 0, 0 )
		fgSizer2.SetFlexibleDirection( wx.BOTH )
		fgSizer2.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		self.m_staticText1 = wx.StaticText( self, wx.ID_ANY, u"Basepairs/Turn", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText1.Wrap( -1 )
		self.m_staticText1.SetForegroundColour( wx.SystemSettings.GetColour( wx.SYS_COLOUR_BTNTEXT ) )
		self.m_staticText1.SetBackgroundColour( wx.SystemSettings.GetColour( wx.SYS_COLOUR_APPWORKSPACE ) )
		
		fgSizer2.Add( self.m_staticText1, 0, wx.ALL, 5 )
		
		self.m_button3 = wx.Button( self, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_button3, 0, wx.ALL, 5 )
		
		self.m_sliderTwist = wx.Slider( self, wx.ID_ANY, 10400, 9400, 11400, wx.DefaultPosition, wx.DefaultSize, wx.SL_HORIZONTAL )
		self.m_sliderTwist.SetBackgroundColour( wx.SystemSettings.GetColour( wx.SYS_COLOUR_APPWORKSPACE ) )
		
		fgSizer2.Add( self.m_sliderTwist, 0, wx.ALL, 5 )
		
		self.m_button4 = wx.Button( self, wx.ID_ANY, u"+", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_button4, 0, wx.ALL, 5 )
		
		self.m_textBptwist = wx.TextCtrl( self, wx.ID_ANY, u"10.4", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_textBptwist, 0, wx.ALL, 5 )
		
		self.m_staticText4 = wx.StaticText( self, wx.ID_ANY, u"Amplitude", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText4.Wrap( -1 )
		fgSizer2.Add( self.m_staticText4, 0, wx.ALL, 5 )
		
		self.m_button5 = wx.Button( self, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_button5, 0, wx.ALL, 5 )
		
		self.m_sliderAmp = wx.Slider( self, wx.ID_ANY, 0, 0, 1000, wx.DefaultPosition, wx.DefaultSize, wx.SL_HORIZONTAL )
		fgSizer2.Add( self.m_sliderAmp, 0, wx.ALL, 5 )
		
		self.m_button6 = wx.Button( self, wx.ID_ANY, u"+", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_button6, 0, wx.ALL, 5 )
		
		self.m_textAmp = wx.TextCtrl( self, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_textAmp, 0, wx.ALL, 5 )
		
		self.m_staticText5 = wx.StaticText( self, wx.ID_ANY, u"Phase", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText5.Wrap( -1 )
		fgSizer2.Add( self.m_staticText5, 0, wx.ALL, 5 )
		
		self.m_button7 = wx.Button( self, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_button7, 0, wx.ALL, 5 )
		
		self.m_sliderPhase = wx.Slider( self, wx.ID_ANY, 0, 0, 6280, wx.DefaultPosition, wx.DefaultSize, wx.SL_HORIZONTAL )
		fgSizer2.Add( self.m_sliderPhase, 0, wx.ALL, 5 )
		
		self.m_button8 = wx.Button( self, wx.ID_ANY, u"+", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_button8, 0, wx.ALL, 5 )
		
		self.m_textPhase = wx.TextCtrl( self, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_textPhase, 0, wx.ALL, 5 )
		
		self.m_staticText6 = wx.StaticText( self, wx.ID_ANY, u"Frequency", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText6.Wrap( -1 )
		fgSizer2.Add( self.m_staticText6, 0, wx.ALL, 5 )
		
		self.m_button9 = wx.Button( self, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_button9, 0, wx.ALL, 5 )
		
		self.m_sliderFreq = wx.Slider( self, wx.ID_ANY, 1000, 500, 1500, wx.DefaultPosition, wx.DefaultSize, wx.SL_HORIZONTAL )
		fgSizer2.Add( self.m_sliderFreq, 0, wx.ALL, 5 )
		
		self.m_button10 = wx.Button( self, wx.ID_ANY, u"+", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_button10, 0, wx.ALL, 5 )
		
		self.m_textFreq = wx.TextCtrl( self, wx.ID_ANY, u"1.0", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_textFreq, 0, wx.ALL, 5 )
		
		self.m_staticText7 = wx.StaticText( self, wx.ID_ANY, u"Modulation period", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText7.Wrap( -1 )
		fgSizer2.Add( self.m_staticText7, 0, wx.ALL, 5 )
		
		self.m_button11 = wx.Button( self, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_button11, 0, wx.ALL, 5 )
		
		self.m_sliderMod = wx.Slider( self, wx.ID_ANY, 1000, 500, 1500, wx.DefaultPosition, wx.DefaultSize, wx.SL_HORIZONTAL )
		fgSizer2.Add( self.m_sliderMod, 0, wx.ALL, 5 )
		
		self.m_button12 = wx.Button( self, wx.ID_ANY, u"+", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_button12, 0, wx.ALL, 5 )
		
		self.m_textMod = wx.TextCtrl( self, wx.ID_ANY, u"1.0", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_textMod, 0, wx.ALL, 5 )
		
		self.m_staticText9 = wx.StaticText( self, wx.ID_ANY, u"Zoom", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText9.Wrap( -1 )
		fgSizer2.Add( self.m_staticText9, 0, wx.ALL, 5 )
		
		
		fgSizer2.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )
		
		self.m_sliderZoom = wx.Slider( self, wx.ID_ANY, 0, 0, 100, wx.DefaultPosition, wx.DefaultSize, wx.SL_HORIZONTAL )
		fgSizer2.Add( self.m_sliderZoom, 0, wx.ALL, 5 )
		
		
		fgSizer2.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )
		
		self.m_textZoom = wx.TextCtrl( self, wx.ID_ANY, u"0", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_textZoom, 0, wx.ALL, 5 )
		
		self.m_buttonReset = wx.Button( self, wx.ID_ANY, u"Reset", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.m_buttonReset, 0, wx.ALL, 5 )
		
		
		bSizer5.Add( fgSizer2, 1, wx.EXPAND, 5 )
		
		gSizer1 = wx.GridSizer( 0, 2, 0, 0 )
		
		self.m_staticText12 = wx.StaticText( self, wx.ID_ANY, u"Distance between cut Basepairs", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText12.Wrap( -1 )
		gSizer1.Add( self.m_staticText12, 0, wx.ALL, 5 )
		
		self.m_textDist = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer1.Add( self.m_textDist, 0, wx.ALL, 5 )
		
		self.m_staticText13 = wx.StaticText( self, wx.ID_ANY, u"Cut dot product", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText13.Wrap( -1 )
		gSizer1.Add( self.m_staticText13, 0, wx.ALL, 5 )
		
		self.m_textDot = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer1.Add( self.m_textDot, 0, wx.ALL, 5 )
		
		self.m_staticText121 = wx.StaticText( self, wx.ID_ANY, u"Shift", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText121.Wrap( -1 )
		gSizer1.Add( self.m_staticText121, 0, wx.ALL, 5 )
		
		self.m_textShift = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer1.Add( self.m_textShift, 0, wx.ALL, 5 )
		
		self.m_staticText131 = wx.StaticText( self, wx.ID_ANY, u"Slide", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText131.Wrap( -1 )
		gSizer1.Add( self.m_staticText131, 0, wx.ALL, 5 )
		
		self.m_textSlide = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer1.Add( self.m_textSlide, 0, wx.ALL, 5 )
		
		self.m_staticText141 = wx.StaticText( self, wx.ID_ANY, u"Rise", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText141.Wrap( -1 )
		gSizer1.Add( self.m_staticText141, 0, wx.ALL, 5 )
		
		self.m_textRise = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer1.Add( self.m_textRise, 0, wx.ALL, 5 )
		
		self.m_staticText17 = wx.StaticText( self, wx.ID_ANY, u"deviation from average twist", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText17.Wrap( -1 )
		gSizer1.Add( self.m_staticText17, 0, wx.ALL, 5 )
		
		self.m_textTwist2 = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer1.Add( self.m_textTwist2, 0, wx.ALL, 5 )
		
		self.m_staticText15 = wx.StaticText( self, wx.ID_ANY, u"Energy", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText15.Wrap( -1 )
		gSizer1.Add( self.m_staticText15, 0, wx.ALL, 5 )
		
		self.m_textEnergy = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer1.Add( self.m_textEnergy, 0, wx.ALL, 5 )
		
		
		bSizer5.Add( gSizer1, 1, wx.EXPAND, 5 )
		
		
		bSizer6.Add( bSizer5, 1, wx.EXPAND, 5 )
		
		fgSizer4 = wx.FlexGridSizer( 0, 2, 0, 0 )
		fgSizer4.SetFlexibleDirection( wx.BOTH )
		fgSizer4.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		self.m_staticText8 = wx.StaticText( self, wx.ID_ANY, u"Number of nucleosomes", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText8.Wrap( -1 )
		fgSizer4.Add( self.m_staticText8, 0, wx.ALL, 5 )
		
		self.m_textNucs = wx.TextCtrl( self, wx.ID_ANY, u"2", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer4.Add( self.m_textNucs, 0, wx.ALL, 5 )
		
		self.m_staticText10 = wx.StaticText( self, wx.ID_ANY, u"Nucleosome repeat length", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText10.Wrap( -1 )
		fgSizer4.Add( self.m_staticText10, 0, wx.ALL, 5 )
		
		self.m_textNRL = wx.TextCtrl( self, wx.ID_ANY, u"197", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer4.Add( self.m_textNRL, 0, wx.ALL, 5 )
		
		self.m_staticText14 = wx.StaticText( self, wx.ID_ANY, u"Unwrap", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText14.Wrap( -1 )
		fgSizer4.Add( self.m_staticText14, 0, wx.ALL, 5 )
		
		self.m_textUnwrap = wx.TextCtrl( self, wx.ID_ANY, u"-20", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer4.Add( self.m_textUnwrap, 0, wx.ALL, 5 )
		
		
		fgSizer4.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )
		
		self.m_buttonCreate = wx.Button( self, wx.ID_ANY, u"Create DNA!", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer4.Add( self.m_buttonCreate, 0, wx.ALL, 5 )
		
		
		bSizer6.Add( fgSizer4, 1, wx.EXPAND, 5 )
		
		bSizer3 = wx.BoxSizer( wx.VERTICAL )
		
		self.m_buttonMinimize = wx.Button( self, wx.ID_ANY, u"Minimize!", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer3.Add( self.m_buttonMinimize, 0, wx.ALL, 5 )
		
		self.m_staticText171 = wx.StaticText( self, wx.ID_ANY, u"Progress", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText171.Wrap( -1 )
		bSizer3.Add( self.m_staticText171, 0, wx.ALL, 5 )
		
		self.m_gaugeProgress = wx.Gauge( self, wx.ID_ANY, 100, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL )
		self.m_gaugeProgress.SetValue( 0 ) 
		bSizer3.Add( self.m_gaugeProgress, 0, wx.ALL, 5 )
		
		self.m_staticEnergy2 = wx.StaticText( self, wx.ID_ANY, u"Energy (kT)", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticEnergy2.Wrap( -1 )
		bSizer3.Add( self.m_staticEnergy2, 0, wx.ALL, 5 )
		
		self.m_textEnergy2 = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer3.Add( self.m_textEnergy2, 0, wx.ALL, 5 )
		
		
		bSizer6.Add( bSizer3, 1, wx.EXPAND, 5 )
		
		
		self.SetSizer( bSizer6 )
		self.Layout()
		
		self.Centre( wx.BOTH )
		
		# Connect Events
		self.m_button3.Bind( wx.EVT_BUTTON, self.MinTwist )
		self.m_sliderTwist.Bind( wx.EVT_SCROLL, self.TwistScroll )
		self.m_button4.Bind( wx.EVT_BUTTON, self.PlusTwist )
		self.m_textBptwist.Bind( wx.EVT_TEXT_ENTER, self.ChangeTextTwist )
		self.m_button5.Bind( wx.EVT_BUTTON, self.MinAmp )
		self.m_sliderAmp.Bind( wx.EVT_SCROLL, self.AmpScroll )
		self.m_button6.Bind( wx.EVT_BUTTON, self.PlusAmp )
		self.m_textAmp.Bind( wx.EVT_TEXT_ENTER, self.ChangeTextAmp )
		self.m_button7.Bind( wx.EVT_BUTTON, self.MinPhase )
		self.m_sliderPhase.Bind( wx.EVT_SCROLL, self.PhaseScroll )
		self.m_button8.Bind( wx.EVT_BUTTON, self.PlusPhase )
		self.m_textPhase.Bind( wx.EVT_TEXT_ENTER, self.ChangeTextPhase )
		self.m_button9.Bind( wx.EVT_BUTTON, self.MinFreq )
		self.m_sliderFreq.Bind( wx.EVT_SCROLL, self.FreqScroll )
		self.m_button10.Bind( wx.EVT_BUTTON, self.PlusFreq )
		self.m_textFreq.Bind( wx.EVT_TEXT_ENTER, self.ChangeTextFreq )
		self.m_button11.Bind( wx.EVT_BUTTON, self.MinMod )
		self.m_sliderMod.Bind( wx.EVT_SCROLL, self.ModScroll )
		self.m_button12.Bind( wx.EVT_BUTTON, self.PlusMod )
		self.m_textMod.Bind( wx.EVT_TEXT_ENTER, self.ChangeTextMod )
		self.m_sliderZoom.Bind( wx.EVT_SCROLL, self.ZoomScroll )
		self.m_buttonReset.Bind( wx.EVT_BUTTON, self.Reset )
		self.m_buttonCreate.Bind( wx.EVT_BUTTON, self.Create )
		self.m_buttonMinimize.Bind( wx.EVT_BUTTON, self.Minimize )
	
	def __del__( self ):
		pass
	
	
	# Virtual event handlers, overide them in your derived class
	def MinTwist( self, event ):
		event.Skip()
	
	def TwistScroll( self, event ):
		event.Skip()
	
	def PlusTwist( self, event ):
		event.Skip()
	
	def ChangeTextTwist( self, event ):
		event.Skip()
	
	def MinAmp( self, event ):
		event.Skip()
	
	def AmpScroll( self, event ):
		event.Skip()
	
	def PlusAmp( self, event ):
		event.Skip()
	
	def ChangeTextAmp( self, event ):
		event.Skip()
	
	def MinPhase( self, event ):
		event.Skip()
	
	def PhaseScroll( self, event ):
		event.Skip()
	
	def PlusPhase( self, event ):
		event.Skip()
	
	def ChangeTextPhase( self, event ):
		event.Skip()
	
	def MinFreq( self, event ):
		event.Skip()
	
	def FreqScroll( self, event ):
		event.Skip()
	
	def PlusFreq( self, event ):
		event.Skip()
	
	def ChangeTextFreq( self, event ):
		event.Skip()
	
	def MinMod( self, event ):
		event.Skip()
	
	def ModScroll( self, event ):
		event.Skip()
	
	def PlusMod( self, event ):
		event.Skip()
	
	def ChangeTextMod( self, event ):
		event.Skip()
	
	def ZoomScroll( self, event ):
		event.Skip()
	
	def Reset( self, event ):
		event.Skip()
	
	def Create( self, event ):
		event.Skip()
	
	def Minimize( self, event ):
		event.Skip()
	

