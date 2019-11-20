#@ Integer (label="SCA1 C:", value=1, persist=true) channel_sca1
#@ Float (label="SCA1 Min:", value=0, persist=true) min_sca1
#@ Integer (label="Hoescht C:", value=2, persist=true) channel_hoescht
#@ Float (label="Hoescht Min:", value=0, persist=true) min_hoescht
#@ Integer (label="OPN C:", value=3, persist=true) channel_opn
#@ Float (label="OPN Min:", value=0, persist=true) min_opn
#@ Integer (label="GFP C:", value=4, persist=true) channel_gfp
#@ Float (label="GFP Min:", value=0, persist=true) min_gfp

# Measure distances of OPN positive cells to the nearest SCA1/GFP double positive region in z-projected stacks
#
# 	-by Richard Butler, Gurdon Institute Imaging Facility


import sys
import math as maths

from java.awt import Color, Font

from ij import IJ, WindowManager, Prefs, ImagePlus, ImageStack
from ij.plugin import ImageCalculator, Duplicator, ZProjector, RoiEnlarger, Straightener, Selection
from ij.plugin.filter import GaussianBlur, MaximumFinder, ThresholdToSelection, Binary, EDM
from ij.process import ImageStatistics, Blitter, ImageProcessor, ShortProcessor, ByteProcessor, AutoThresholder, FloodFiller
from ij.measure import ResultsTable, Measurements
from ij.gui import Roi, ShapeRoi, TextRoi, Overlay

from fiji.process3d import SEDT


MINA = 5	#cell area range, µm²
MAXA = 100

if((channel_sca1 + channel_hoescht + channel_opn + channel_gfp != 10) or
	channel_sca1<1 or channel_hoescht<1 or channel_opn<1 or channel_gfp<1 or
	channel_sca1>4 or channel_hoescht>4 or channel_opn>4 or channel_gfp>4) :
	IJ.error("Invalid channels")
	exit(0)

def maxProject(image, chan):
	proj = ShortProcessor(image.getWidth(), image.getHeight())
	stack = image.getStack()
	for z in range(1,imp.getNSlices()+1):
		ip = stack.getProcessor(image.getStackIndex(chan,z,1))
		proj.copyBits(ip, 0,0, Blitter.MAX)
	return proj

def fillHoles(mask):
	width = mask.getWidth()
	height = mask.getHeight()
	ff = FloodFiller(mask)
	mask.setColor(127)
	foreground = 127
	background = 0
	for y in range(height):
	    if mask.getPixel(0,y)==background:
	    	ff.fill(0, y)
	    if mask.getPixel(width-1,y)==background:
	    	ff.fill(width-1, y)
	for x in range(width):
	    if mask.getPixel(x,0)==background:
	    	ff.fill(x, 0)
	    if mask.getPixel(x,height-1)==background:
	    	ff.fill(x, height-1)
	n = width*height
	for i in range(n):
		if mask.get(i)==127:
		    mask.set(i, 0)
		else:
		    mask.set(i, 255)

def mask2D(ip, sigmaPx, k, method, minimum, doFillHoles, doWatershed):
	mask = ip.duplicate()
	sub = mask.duplicate()
	mask.blurGaussian(sigmaPx)
	if k > 0:
		sub.blurGaussian(k*sigmaPx)
		mask.copyBits(sub, 0,0, Blitter.SUBTRACT)
	
	stats = mask.getStatistics()
	hist = mask.getStatistics().histogram
	thresh = AutoThresholder().getThreshold(method, hist)
	thresh = (thresh/float(255)) * (stats.max-stats.min) + stats.min

	mask.threshold( int(thresh) )
	mask = mask.convertToByte(False)

	if doFillHoles:
		fillHoles(mask)

	if doWatershed:
		floatEdm = EDM().makeFloatEDM(mask, 0, False)
		maxIp = MaximumFinder().findMaxima(floatEdm, 0.5, ImageProcessor.NO_THRESHOLD, MaximumFinder.SEGMENTED, False, True)
		if (maxIp != None):
			mask.copyBits(maxIp, 0, 0, Blitter.AND)

	mask.dilate()
	mask.erode()

	mask.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
	roi = ThresholdToSelection().convert(mask)
	ip.setRoi(roi)
	mean = ip.getStatistics().mean
	
	if mean < minimum:	#if the mask area intensity mean in the original image is less than the minimum required
		mask = ByteProcessor(ip.getWidth(), ip.getHeight())	#return empty mask
	
	return mask

def getRois(mask):
	mask.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
	composite = ThresholdToSelection().convert(mask)
	rois = ShapeRoi(composite).getRois()
	return rois

def getRoi(mask):
	mask.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
	roi = ThresholdToSelection().convert(mask)
	return roi


imp = IJ.getImage()
cal = imp.getCalibration()
ol = Overlay()

proj_sca1 = maxProject(imp, channel_sca1)
proj_hoescht = maxProject(imp, channel_hoescht)
proj_opn = maxProject(imp, channel_opn)
proj_gfp = maxProject(imp, channel_gfp)

mask_sca1 = mask2D(proj_sca1, 1.0/cal.pixelWidth, 6, AutoThresholder.Method.Otsu, min_sca1, False, False)	#surrounding lumen
mask_hoescht = mask2D(proj_hoescht, 1.5/cal.pixelWidth, 5, AutoThresholder.Method.Otsu, min_hoescht, True, True)		#Hoescht
mask_opn = mask2D(proj_opn, 1.0/cal.pixelWidth, 0, AutoThresholder.Method.MaxEntropy, min_opn, False, False)		#surrounding lumen fragments
mask_gfp = mask2D(proj_gfp, 0.5/cal.pixelWidth, 5, AutoThresholder.Method.Otsu, min_gfp, False, True)	#cells

roi_sca1 = getRoi(mask_sca1)
roi_sca1.setStrokeColor(Color.YELLOW)
ol.add(roi_sca1)

#roi_hoescht = getRoi(mask_hoescht)
#roi_hoescht.setStrokeColor(Color.GREEN)
#ol.add(roi_hoescht)

roi_opn = getRoi(mask_opn)
roi_opn.setStrokeColor(Color.CYAN)
ol.add(roi_opn)

roi_gfp = getRoi(mask_gfp)
roi_gfp.setStrokeColor(Color.BLUE)
ol.add(roi_gfp)

mask_sca1_gfp = mask_sca1.duplicate()
mask_sca1_gfp.copyBits(mask_gfp, 0,0, Blitter.AND)
fillHoles(mask_sca1_gfp)
cells_sca1_gfp = getRois(mask_sca1_gfp)

mask_hoescht_opn = mask_hoescht.duplicate()
mask_hoescht_opn.copyBits(mask_opn, 0,0, Blitter.AND)
cells_hoescht_opn = getRois(mask_hoescht_opn)

sedt = SEDT()
edm_sca1_gfp = sedt.compute(ImagePlus("", mask_sca1_gfp).getStack())

rt = ResultsTable.getResultsTable()
celli = 0
FONT = Font(Font.SANS_SERIF, Font.BOLD, 16)
for cell in cells_hoescht_opn:
	area = cell.getStatistics().area * cal.pixelWidth * cal.pixelHeight
	if area>=MINA and area<=MAXA:	
		cell.setStrokeColor(Color.MAGENTA)
		ol.add(cell)

		edm_sca1_gfp.setRoi(cell)
		dist = edm_sca1_gfp.getStatistics().mean * cal.pixelWidth

		row = rt.getCounter()
		rect = cell.getBounds()
		rt.setValue("Image", row, imp.getTitle())
		rt.setValue("OPN+ Cell", row, celli)
		rt.setValue("X", row, (rect.x+rect.width/2.0)*cal.pixelWidth)
		rt.setValue("Y", row, (rect.y+rect.height/2.0)*cal.pixelHeight)
		rt.setValue("OPN+ cell distance from SCA1+ GFP+ ("+u"\u00b5"+"m"+")", row, dist)

		label = TextRoi(rect.x, rect.y-16, str(celli), FONT)
		label.setStrokeColor(Color.MAGENTA)
		ol.add(label)
		
		celli += 1

rt.show("Results")
imp.setOverlay(ol)


