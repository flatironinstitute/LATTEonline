import os
import re
import ast
import csv
import sys
import time
import json
#import js2py
import random
import sqlite3
import matplotlib
import numpy as np
import pandas as pd
import configparser

import seaborn as sb
matplotlib.use('Agg')
#import sqlite3 as sql
from glob import glob
import lightkurve as lk
from matplotlib import cm
from sqlite3 import Error
import pandas.io.sql as psql
import astropy.io.fits as pf
from functools import partial
import matplotlib.pyplot as plt
from astropy.table import Table
from os.path import dirname, join
from argparse import ArgumentParser
from os.path import basename, exists
from scipy.optimize import curve_fit
from astroquery.mast import Catalogs

# bokeh tools - - - - - - - - - - -
from jinja2 import Template
from bokeh.io import curdoc
from bokeh.resources import CDN
from bokeh.embed import file_html

from bokeh.models.widgets import Select
from bokeh.themes import built_in_themes
from bokeh.plotting import figure, curdoc
from bokeh.events import DoubleTap, ButtonClick, DocumentReady
from bokeh.layouts import column, row, layout, Spacer
from bokeh.models import TextInput, DataTable, Dropdown, DateFormatter, TableColumn, Paragraph, Button, Slider, CustomJS, Range1d, Select, Toggle, ColumnDataSource, Panel, Tabs, Span, Div, Select, Slider, Div

# custom functions - - - - - - - - -
import LATTEutils, LATTEbrew

global target_name_input
global datapath
global syspath
global TIC_url
global button_loaded
global twitter_logo_text1

# scp -r /Users/neisner/Documents/code/LATTEonline ccalin029:/mnt/home/neisner/projects/
# scp -r /mnt/home/neisner/projects/LATTEonline popeye:/mnt/home/neisner/projects/
# scp -r /Users/neisner/Documents/code/LATTEonline popeye:/mnt/home/neisner/projects/
# scp -r /Users/neisner/Documents/code/LATTEonline/app/*.py popeye:/mnt/home/neisner/projects/LATTEonline/app/
# open the confidguration file so that this can be used on different machines as long as there is a cofnigureation file pointing to the LATTE output folder
#rsync -r /Users/neisner/Documents/code/LATTEonline/ popeye:/mnt/home/neisner/projects/LATTEonline

syspath = str(os.path.abspath(LATTEutils.__file__))[0:-14]

# open the argument parser thing to determine what paths to use

config = configparser.ConfigParser()
config.read('_config.txt')
config.sections()

# configuration files (see whether this is running locally on Nora's computer or not)
if "/Users/neisner/Documents/code/" in syspath:
	datapath = config['local']['datapath'] # where the data is
	outpath = config['local']['outpath'] # where to store to the temp files ( in '/static')
	codepath = config['local']['codepath'] # where the codes are to access the data (might not need anymore now)

	logo_text = """<div id="logo">
				<img alt="LATTE Logo" title="Client Logo" src="applatte/static/LATTE_imgs/LATTE_logo_small.png" style="width: 90%;max-height: 100%"/>
				</div>
	"""
	PHT_logo_text = """<div id="PHT logo">
				<img alt="PHT logo" title="PHT Logo" src="applatte/static/LATTE_imgs/pht_logo_small.png" style="width: 83%;max-height: 83%"/>
				</div>
	"""

else:
	datapath = config['workstation']['datapath']
	outpath = config['workstation']['outpath']
	codepath = config['workstation']['codepath']

	logo_text = """<div id="logo">
				<img alt="LATTE Logo" title="Client Logo" src="/app/static/LATTE_imgs/LATTE_logo_small.png" style="width: 90%;max-height: 100%"/>
				</div>
	"""
	PHT_logo_text = """<div id="PHT logo">
				<img alt="PHT logo" title="PHT Logo" src="/app/static/LATTE_imgs/pht_logo_small.png" style="width: 83%;max-height: 83%"/>
				</div>
	"""



#https://css-tricks.com/simple-social-sharing-links/

# - - - - - - - - - - - - needed functions - - - - - - - - - - - - 
def rebin(arr, new_shape):
	'''
	bin the data.
	'''
	shape = (
		new_shape[0],
		arr.shape[0] // new_shape[0],
		new_shape[1],
		arr.shape[1] // new_shape[1],
	)
	return arr.reshape(shape).mean(-1).mean(1)

# - - - - - - - - - - 
# - - - - - - - - - -
# START the document

doc = curdoc()


with open('./templates/index.html') as f:
    index_template = Template(f.read())
    doc.template = index_template


# ------- LOGOS -------
# ---------------------


# add the LATTE logo
latte_logo = Div(text=logo_text)


# add the PHT and PHCC logos 
pht_latte_logo = Div(text=PHT_logo_text)

# ------------------------------------------------------------------------------
# - - - - - - - - - - -	LOAD SQL TABLE FOR DATA ACCESS- - - - - - - - - - - - - 
# ------------------------------------------------------------------------------


DB_FILE = datapath + '/tess_database.db'

# database to get the data
class Database:
	def __init__(self, db_file):
		self.connection = None
		self.db_file = db_file
		self.re_lc = re.compile('data/s000(\d+)/tess\d{13}-s(\d{4})-(\d+)-\d+-s_lc.fits')
		self.open_connection()

	def __del__(self):
		self.close_connection()

	def open_connection(self):

		try:
			self.connection = sqlite3.connect(self.db_file)
			print("Database file opened, SQLite version:", sqlite3.version)
		except Error as e:
			print("Database connect failed with error", e)

	def close_connection(self):
		if self.connection:
			self.connection.close()

	# tid needs to be a int. Just cast it this way it takes care of leading 0's
	def search(self, tic_id):

		SUBFOLDER_SPLIT = 1e7

		#print("start db search")
		sql_select_by_tic = """ SELECT sector_id, tp_filename FROM fits
												WHERE tic_id = ? """

		try:
			c = self.connection.cursor()
			c.execute(sql_select_by_tic, (tic_id,))
			rows = c.fetchall()
			
			if len(rows) == 0:
				print('Could not find tic', tic_id, "in data")
				return [0,0,0]

			sector = ([i[0] for i in rows])
			tp_filename_list = []
			lc_filename_list = []

			for row in rows:
					tp_format = 'data/s{sector:04d}/{subfolder:09d}/{filename}'
					tp_filename = tp_format.format(sector = row[0], subfolder = int(tic_id / SUBFOLDER_SPLIT), filename = row[1])
					lc_filename = tp_filename.replace('tp', 'lc')
					tp_filename_list.append(tp_filename)
					lc_filename_list.append(lc_filename)

			#lc_filename_list = ([tp.replace('tp', 'lc') for tp in tp_filename])
		
		except Error as e:
			print("INSERT failed with error", e)

		#print("end db search")
		return (sector, lc_filename_list, tp_filename_list)


	def search_random(self):

		#sql_select_by_tic = """ SELECT RANDOM() FROM fits LIMIT 1"""

		sql_select_by_tic = """SELECT * FROM fits ORDER BY RANDOM() LIMIT 1;"""

		c = self.connection.cursor()
		c.execute(sql_select_by_tic)
		rows = c.fetchall() # change this to make it iterate over multiple 

		sector = ([i[0] for i in rows])
		tp_filename_list = []
		lc_filename_list = []

		SUBFOLDER_SPLIT = 1e7

		for row in rows:
				tp_format = 'data/s{sector:04d}/{subfolder:09d}/{filename}'
				tp_filename = tp_format.format(sector = row[0], subfolder = int(row[1] / SUBFOLDER_SPLIT), filename = row[2])
				lc_filename = tp_filename.replace('tp', 'lc')
				tp_filename_list.append(tp_filename)
				lc_filename_list.append(lc_filename)

		return (sector, lc_filename_list, tp_filename_list)

# load the database - this is a quick way to search for all of the needed urls
db = Database(DB_FILE)

# add text for the information of the target! Empty at the begingin but as you load data this will update

# parameters used throughout 
width1 = 150
width2 = 300

color1 = 'dimgray'
color2 = '#293040'
fontsize = '110%'
text_height = 15
width_gadgets = 250


target_targetinfo_name_text  = " <!DOCTYPE html> \
<html> \
<style> \
.tooltip { \
  position: relative; \
  display: inline-block; \
  border-bottom: 1px dotted black; \
} \
 \
.tooltip .tooltiptext { \
  visibility: hidden; \
  width: 300px; \
  background-color: #293040; \
  color: #fff; \
  text-align: center; \
  border-radius: 6px; \
  font-size: 14px; /* Set a font-size */ \
  padding: 5px 0; \
   \
  /* Position the tooltip */ \
  position: absolute; \
  z-index: 1; \
  top: 100%; \
  left: 50%; \
  margin-left: -60px; \
} \
 \
.tooltip:hover .tooltiptext { \
  visibility: visible; \
} \
</style> \
<body style='text-align:center;'>  "
 
text_end = " <b> Target Information <b> <b> <br> \
<div class='tooltip'> <br> <b> TIC ID <b> <br> \
  <span class='tooltiptext'>{}</span> \
</div> \
<br> <b> RA <br> Dec <br> TESS mag <br> Teff <br> Stellar Radius <br> <br> ExoFOP <br> SIMBAD <b>\
 \
</body> \
</html> \
".format("other names")


# ------------------------------------------------------------------------------------
# - - - - - - - - - - -	DEFINE THE DATA AND PLOTTING REGION - - - - - - - - - - - - - 
# ------------------------------------------------------------------------------------

source = ColumnDataSource(data=dict(x=[], y=[])) # unbinned data
source_binned = ColumnDataSource(data=dict(x_binned=[], y_binned=[])) # binned data
source_md = ColumnDataSource(data=dict(md_x1 = [], md_y1 =[], md_x2 = [], md_y2 =[])) # momentun dumps
source_dummy = ColumnDataSource(data=dict(x_dummy=[1e2], y_dummy=[1e-10])) # dummy data
# - - - - - - - 
# check whether the URL has any information in it that we care about
# - - - - - - -

# - - - - make the figures - - - - -

p = figure(height=350, 
	width=700, 
	title="TESS LC", 
	toolbar_location="below", 
	sizing_mode="scale_both",
	tools='box_zoom,wheel_zoom,pan,reset, tap',
	)


main_points = p.circle(x="x", 
	y="y", 
	source=source, 
	size=3, 
	color="darkorange", 
	line_color=None, 
	alpha=0.8,
	)


binned_points = p.circle(x="x_binned", 
	y="y_binned", 
	source=source_binned, 
	size=3, 
	color="black", 
	line_color=None, 
	alpha=0.5)


main_segment = p.segment(x0 = "md_x1",
	y0 = "md_y1",
	x1 = "md_x2",
	y1 = "md_y2",
	level = 'underlay',
	source = source_md,
	color="blue", 
	line_width=1)

p.xaxis.axis_label = r"$$\text{Time (BJD - 2457000)}$$"
p.yaxis.axis_label = r"$$\text{Normalised Flux}$$"

# when LATTE is first launched, we want to prompt the user to enter a TIC ID on the left 
# to do that we will put some text on the empty plot 

#welcome_text_source = ColumnDataSource(
#			pd.DataFrame.from_records([dict(
#			x=[100],
#			y=[50],
#			text="1. Enter a TIC ID to start. \n 2. Double click on figure to select events.",
#			color="black")]))
#
#welcome_text = p.text(
#		text="text",
#		text_align="center",
#		text_font_size="15px",
#		source= welcome_text_source)

# - - - - -	ZOOM TAB - - - - - - 

source_zoom = ColumnDataSource(data=dict(x_zoom=[], y_zoom=[]))
source_zoom_binned = ColumnDataSource(data=dict(x_zoom_binned=[], y_zoom_binned=[]))

plot_zoom = figure(
	plot_height=300,
	min_width=600,
	min_height=300,
	title="",
	sizing_mode="stretch_both",
	#x_range=(0,5),
	toolbar_location="below", 
	tools='box_zoom,wheel_zoom,pan,reset'
)

zoom_points	= plot_zoom.circle(x="x_zoom", 
	y="y_zoom", 
	source=source_zoom, 
	size=5, 
	color="darkorange", 
	line_color=None, 
	alpha=0.8,
	)

zoom_points_binned = plot_zoom.circle(x="x_zoom_binned", 
	y="y_zoom_binned", 
	source=source_zoom_binned, 
	size=5, 
	color="black", 
	line_color=None, 
	alpha=0.5)

#zoom_segment = plot_zoom.segment(x0 = "md_x1", 
#	y0 = "md_y1",
#	x1 = "md_x2",
#	y1 = "md_y2",
#	source = source_md,
#	level = 'underlay',
#	color="blue", 
#	line_width=1)

# add x and y labels
plot_zoom.xaxis.axis_label = r"$$\text{Time (BJD - 2457000)}$$"
plot_zoom.yaxis.axis_label = r"$$\text{Normalised Flux}$$"

# - - - - -	BACKGROUND TAB - - - - - - 

source_bkg = ColumnDataSource(data=dict(x_bkg=[], y_bkg=[]))

plot_bkg = figure( 
	plot_height=300,
	min_width=600,
	min_height=300,
	title="",
	sizing_mode="stretch_both",
	x_range=plot_zoom.x_range, 
	toolbar_location="below", 
	tools='box_zoom,wheel_zoom,pan,reset'
)

bkg_points = plot_bkg.circle(
	x="x_bkg",
	y="y_bkg",
	source=source_bkg,
	line_color=None,
	color="blue",
	alpha=0.5,
	size=5,
)

# add x and y labels
plot_bkg.xaxis.axis_label = r"$$\text{Time (BJD - 2457000)}$$"
plot_bkg.yaxis.axis_label = r"$$\text{Background Flux}$$"


# - - - - -	CENTROID TAB - - - - - - 

source_xcen1 = ColumnDataSource(data=dict(x_xcen1=[], y_xcen1=[]))
source_xcen2 = ColumnDataSource(data=dict(x_xcen2=[], y_xcen2=[]))
source_ycen1 = ColumnDataSource(data=dict(x_ycen1=[], y_ycen1=[]))
source_ycen2 = ColumnDataSource(data=dict(x_ycen2=[], y_ycen2=[]))

plot_xcen = figure( 
	plot_height=150,
	min_width=600,
	min_height=150,
	title="",
	sizing_mode="stretch_both",
	x_range=plot_zoom.x_range,
	toolbar_location="below", 
	tools='box_zoom,wheel_zoom,pan,reset'
)

plot_ycen = figure(
	plot_height=150,
	min_width=600,
	min_height=150,
	title="",
	sizing_mode="stretch_both",
	x_range=plot_zoom.x_range,
	toolbar_location="below", 
	tools='box_zoom,wheel_zoom,pan,reset'
)

xcen1_points = plot_xcen.circle(
	x="x_xcen1",
	y="y_xcen1",
	source=source_xcen1,
	line_color=None,
	color="red",
	alpha=0.4,
	size=5,
)

xcen2_points = plot_xcen.circle(
	x="x_xcen2",
	y="y_xcen2",
	source=source_xcen2,
	line_color=None,
	color="black",
	alpha=0.3,
	size=5,
)

ycen1_points = plot_ycen.circle(
	x="x_ycen1",
	y="y_ycen1",
	source=source_ycen1,
	line_color=None,
	color="red",
	alpha=0.4,
	size=5,
)

ycen2_points = plot_ycen.circle(
	x="x_ycen2",
	y="y_ycen2",
	source=source_ycen2,
	line_color=None,
	color="black",
	alpha=0.3,
	size=5,
)

# add x and y labels
plot_xcen.yaxis.axis_label = r"$$\text{x-centroid}$$"
plot_ycen.xaxis.axis_label = r"$$\text{Time (BJD - 2457000)}$$"
plot_ycen.yaxis.axis_label = r"$$\text{y-centroid}$$"


# - - - - - PERIODOGRAM TAB - - - - - - 

source_periodgrm = ColumnDataSource(
	data=dict(x_periodgrm=[], y_periodgrm=[])
)

#source_periodgrm_smooth = ColumnDataSource(
#	data=dict(x_periodgrm_smooth=[], y_periodgrm_smooth=[])
#)

# add an extra figure to plot an additional parameter - such as the background or the centroid shifts.
plot_periodgrm = figure( 
	plot_height=300,
	min_width=600,
	min_height=300,
	title="",
	sizing_mode="stretch_both",
	x_axis_type="log",
	y_axis_type="log",
	toolbar_location="below", 
	tools='box_zoom,wheel_zoom,pan,reset'
)

pg_points = plot_periodgrm.line(
	x="x_periodgrm",
	y="y_periodgrm",
	source=source_periodgrm,
	line_color="black",
	color="black",
	alpha=1,
)

dummy_points = plot_periodgrm.circle(x="x_dummy", 
	y="y_dummy", 
	source=source_dummy, 
	size=30, 
	color="white", 
	line_color=None, 
	alpha=0,
	)



# - - - -	EVOLUTIONARY TRACKS TAB - - - - 

# Ddownload the files for the evolutionary tracks
phase0 = pd.read_csv('{}/LATTE_eep_data/eep_phase0.csv'.format(syspath)) # these are the main-sequence tracks
phase2 = pd.read_csv('{}/LATTE_eep_data/eep_phase2.csv'.format(syspath)) # these are the post main-sequence tracks

# we need the TOI file which has the stellar parameters (unfortunately this is a different file to the other TOI file - both are installed with --new-data)
infile = "{}/data/TOI_list_star_params.txt".format(datapath)
exoplanets = pd.read_csv(infile, comment = '#')

# if either the radius or the temperature is unknown, ust ignore it because it can't be plotted.
# also uknown radii default to 1, so get rid of these as they are wrong.
exoplanets_r = exoplanets[~(np.isnan(exoplanets['Stellar Radius (R_Sun) err']) & (exoplanets['Stellar Radius (R_Sun)'] == 1))]

source_evol_points = ColumnDataSource(
	data=dict(x_evol_points=list(exoplanets_r['Stellar Eff Temp (K)']), y_evol_points=list(exoplanets_r['Stellar Radius (R_Sun)']))
)

source_evol_target = ColumnDataSource(
	data=dict(x_evol_target=[], y_evol_target=[])
)

# add an extra figure to plot an additional parameter - such as the background or the centroid shifts.
plot_evol = figure( 
	plot_height=300,
	min_width=600,
	min_height=300,
	title="",
	sizing_mode="stretch_both",
	#x_axis_type="log",
	#y_axis_type="log",
	toolbar_location="below", 
	tools='box_zoom,wheel_zoom,pan,reset'
)


for i in range(300,1700, 100): # we have evolutionary tracks from 0.3 M_sun to 1.6 M_sun in steps of 0.1
	colname_phase0_T = (str(i).rjust(6, "0")) + '0M_phase0_T'
	colname_phase0_R = (str(i).rjust(6, "0")) + '0M_phase0_R'

	colname_phase2_T = (str(i).rjust(6, "0")) + '0M_phase2_T'
	colname_phase2_R = (str(i).rjust(6, "0")) + '0M_phase2_R'

	phase0_T = phase0[colname_phase0_T]
	phase0_R = phase0[colname_phase0_R]

	phase2_T = phase2[colname_phase2_T]
	phase2_R = phase2[colname_phase2_R]

	# the -99 values are just there as padding to make the lists the same length -- we dont' need them to mask out.
	mask0 = (phase0_T != -99) & (phase0_R != -99)
	mask2 = (phase2_T != -99) & (phase2_R != -99)

	source_evol_MS = ColumnDataSource(
		data=dict(x_evol_MS=list(phase0_T[mask0]), y_evol_MS=list(phase0_R[mask0]))
	)

	source_evol_pMS = ColumnDataSource(
		data=dict(x_evol_pMS=list(phase2_T[mask2]), y_evol_pMS=list(phase2_R[mask2]))
	)

	# plot the 1 Solar Mass track in a different colour to make it stand out - easier for reference.
	if i == 1000:
		color_MS = 'red'
		color_pMS = 'red'
	else:
		color_MS = 'black'
		color_pMS = 'grey'

	evol_targets_MS = plot_evol.line(
		x="x_evol_MS",
		y="y_evol_MS",
		source=source_evol_MS,
		line_color=color_MS,
		color=color_MS,
		alpha=1,
	)

	evol_targets_pMS = plot_evol.line(
		x="x_evol_pMS",
		y="y_evol_pMS",
		source=source_evol_pMS,
		line_color=color_pMS,
		color=color_pMS,
		line_dash = 'dashed',
		alpha=1,
	)

evol_points = plot_evol.circle(x="x_evol_points", 
	y="y_evol_points", 
	source=source_evol_points, 
	size=2, 
	color="navy", 
	line_color=None, 
	alpha=0.3,
	)

evol_targets = plot_evol.scatter(x="x_evol_target", 
	y="y_evol_target", 
	source=source_evol_target, 
	size=15, 
	color="magenta", 
	line_color=None, 
	alpha=1,
	marker='star'
	)


# add x and y labels

plot_evol.x_range.start = 8000
plot_evol.x_range.end = 3000
plot_evol.y_range.start = 0.2
plot_evol.y_range.end = 4

plot_evol.xaxis.axis_label = r"$$\text{T}_{\text{eff}} (K)$$"
plot_evol.yaxis.axis_label = r"$$\text{Radius (R}_{\odot})$$"

div_error = Div(text=""" this is not a TESS target """, visible=False, width=150, height=15, align = 'start', style={'font-size': '90%', 'color': 'red'})
div_error_astroquery = Div(text=""" can't connect to astroquery """, visible=False, width=150, height=15, align = 'start', style={'font-size': '80%', 'color': 'red'})


# - - - -	PHASE FOLD TAB - - - - 

source_phased = ColumnDataSource(data=dict(x_phase=[], y_phase=[]))
source_phased_fit = ColumnDataSource(data=dict(x_phase_fit=[], y_phase_fit=[]))

#source_phased_binned = ColumnDataSource(data=dict(x_phase_binned=[], y_phase_binned=[]))


plot_phase = figure(
	plot_height=300,
	min_width=650,
	min_height=200,
	title="",
	sizing_mode="stretch_both",
	toolbar_location="below", 
	tools='box_zoom,wheel_zoom,pan,reset'
)

phase_points = plot_phase.circle(x="x_phase", 
	y="y_phase", 
	source=source_phased, 
	size=5, 
	color="darkslategray", 
	line_color=None, 
	alpha=0.8,
	)

phase_fit = plot_phase.line(x="x_phase_fit", 
	y="y_phase_fit", 
	source=source_phased_fit, 
	color="red", 
	line_color="red", 
	alpha=1,
	line_width = 3
	)

plot_phase.xaxis.axis_label = r"\text{Phase (days)}"
plot_phase.yaxis.axis_label = r"\text{Flux}"


#phase_points_binned = plot_phase.circle(x="x_phase_binned", 
#	y="y_phase_binned", 
#	source=source_phased_binned, 
#	size=5, 
#	color="black", 
#	line_color=None, 
#	alpha=0.5)


# - - - - - - - - - - - - - - - - 

# function to select the data	- only run when the enter button is pressed (takes a while at the moment to download data)
# to do: data locally or need s pinnign wheel, progress bar etc


color_theme_menu = Select(title="Color theme", value="light",
				 options=['light', 'dark'], width=width_gadgets)


plotting_options_menu = Select(title="Report plot options", value="split axis by hemispheres",
				 options=["split axis by hemispheres", "split axis by sectors", "don't split axis", "plot sectors with events"], width=width_gadgets)


def color_theme():

	if color_theme_menu.value == 'light':

		for plot in [p, plot_zoom, plot_bkg, plot_xcen, plot_ycen, plot_periodgrm, plot_evol, plot_phase]:
			plot.background_fill_color = 'white'
			plot.outline_line_color = 'black'
			plot.ygrid.grid_line_color = 'gainsboro'
			plot.xgrid.grid_line_color = 'gainsboro'

		binned_points.glyph.fill_color = 'black'
		main_points.glyph.fill_color = 'darkorange'
		main_segment.glyph.line_color = 'blue'

		zoom_points.glyph.fill_color = 'darkorange' 
		zoom_points_binned.glyph.fill_color = 'black'
		#zoom_segment.glyph.line_color = 'blue'

		bkg_points.glyph.fill_color = 'blue'

		xcen1_points.glyph.fill_color = 'red'
		ycen1_points.glyph.fill_color = 'red'
		xcen2_points.glyph.fill_color = 'black'
		ycen2_points.glyph.fill_color = 'black'

		pg_points.glyph.line_color = 'black'

		#evol_targets_MS.glyph.line_color = 'black'
		evol_targets_pMS.glyph.line_color = 'grey'
		evol_points.glyph.fill_color = 'navy'
		evol_targets.glyph.fill_color = 'maroon'

		phase_points.glyph.fill_color = 'darkslategray'

	elif color_theme_menu.value == 'dark':

		for plot in [p, plot_zoom, plot_bkg, plot_xcen, plot_ycen, plot_periodgrm, plot_evol, plot_phase]:
			plot.background_fill_color = '#293040'
			plot.outline_line_color = 'white'
			plot.ygrid.grid_line_color = 'dimgray'
			plot.xgrid.grid_line_color = 'dimgray'

			main_points.glyph.fill_color = 'white'
			binned_points.glyph.fill_color = 'dimgray'

		binned_points.glyph.fill_color = 'white'
		main_points.glyph.fill_color = 'dimgray'
		main_segment.glyph.line_color = 'blue'

		zoom_points.glyph.fill_color = 'dimgray' 
		zoom_points_binned.glyph.fill_color = 'white'
		#zoom_segment.glyph.line_color = 'paleturquoise'

		bkg_points.glyph.fill_color = 'paleturquoise'

		xcen1_points.glyph.fill_color = 'salmon'
		ycen1_points.glyph.fill_color = 'salmon'
		xcen2_points.glyph.fill_color = 'white'
		ycen2_points.glyph.fill_color = 'white'
		
		pg_points.glyph.line_color = 'white'

		#evol_targets_MS.glyph.line_color = 'white'
		evol_targets_pMS.glyph.line_color = 'dimgray'
		evol_points.glyph.fill_color = 'paleturquoise'
		evol_targets.glyph.fill_color = 'pink'

		phase_points.glyph.fill_color = 'white'

# plot_zoom, plot_evol

color_theme_menu.on_change('value', lambda attr, old, new: color_theme())


# -----------------------------------
# - - - - - - FUNCTIONS - - - - - - -
# -----------------------------------

def select_data(lucky = False, single = True, url = False):
	
	'''
	Function to load the data needed to plot the initial LCs displayed in the interface.
	'''
	# take note that the data is being dowloaded 

	#print ("start select data")

	global transit_time_list

	# delete the markings
	transit_time_list = []
	tr_text.text = "{}".format(str(np.around(transit_time_list,2))[1:-1])

	# -----------------
	# grey out the DV report button if no events
	if len(transit_time_list) > 0:
		button_done.disabled=False
	else:
		button_done.disabled=True


	if (url == False):
		input_name = str(target_name_input.value)
	else:
		input_name = str(dummy_name_input.value)

	# if lucky is false then see what the tic is (otherwise there is no tic so don't do this)
	if lucky == False:
		try:
			int(input_name)
			tic_num_only = True
		except:
			if len(input_name) > 0:
				tic_num_only = False
			else:
				all_data = [-99]
				div_spinner_loading_data.visible = False
				return all_data

		if tic_num_only == False:
			#input_name = str(target_name_input.value)
	
			try: # try with simbad 
				#TIC = 'TIC {}'.format(TICID)
				
				#print ("start query")
				catalog_data = Catalogs.query_object(input_name, catalog="TIC", radius = 0.001)
				
				TICID	= str(catalog_data[0]['ID'])
				TIC = "TIC " + str(TICID) # numerical only
				
				#print ("start db")
	
				_ , lc_paths0, tp_paths0 = db.search(int(TICID))
	
				#print ("end db")
	
				if single == True:
					lc_paths0 = [lc_paths0[-1]]
					tp_paths0 = [tp_paths0[-1]]
	
				df_other_names = (pd.DataFrame(np.array(catalog_data)))
				
				other_names_ref = ['HIP','TYC','UCAC','TWOMASS','SDSS','ALLWISE','GAIA','APASS','KIC']
				
				names_str = ''
				
				for name in other_names_ref:
					
					if len(df_other_names[name][0]) > 0:
						names_str = names_str +  (str(name) + ' ' + str(df_other_names[name][0]) + "<br>")
					else:
						names_str = names_str + (str(name) + ' ' + "--- <br>")
	
			except: # if the TIC thing didn't work either, then check whether SIMBAD is currently working
				try:
					catalog_data = Catalogs.query_object('55525572', catalog="TIC", radius = 0.001)
					all_data = [-99]
					div_error_astroquery.visible=True
					div_spinner_loading_data.visible = False
					return all_data
				except:
					#p.title.text_color="red"
					#p.title.text = "THIS TARGET HASN'T BEEN OSBERVED BY TESS"
					all_data = [-99]
					div_error.visible=True
					div_spinner_loading_data.visible = False
					return all_data
	
		elif tic_num_only == True:
	
			#if it fails with SIMABD then the target can't be resolved so check whether its just the TIC UD
			try:
	
				#print ("start query - tic only")
	
				TIC = "TIC " + str(input_name) # numerical only
				_ , lc_paths0, tp_paths0 = db.search(int(input_name))
				
				if single == True:
					lc_paths0 = [lc_paths0[-1]]
					tp_paths0 = [tp_paths0[-1]]
				
				names_str = 'unavailable'
				#print ("end db")
	
			except: # if the TIC thing didn't work either, then check whether SIMBAD is currently working
				try:
					catalog_data = Catalogs.query_object('55525572', catalog="TIC", radius = 0.001)
					all_data = [-99]
					div_error_astroquery.visible=True
					div_spinner_loading_data.visible = False
					return all_data
				except:
					#p.title.text_color="red"
					#p.title.text = "THIS TARGET HASN'T BEEN OSBERVED BY TESS"
					all_data = [-99]
					div_error.visible=True
					div_spinner_loading_data.visible = False
					return all_data

		else:
			all_data = [-99]
			div_error_astroquery.visible=True
			div_spinner_loading_data.visible = False
	
			return all_data

	else: # if lucky is true

		_ , lc_paths0, tp_paths0 = db.search_random()

		names_str = 'unavailable'
	

	div_error.visible = False
	div_error_astroquery.visible = False
	lc_paths0 = lc_paths0
	tp_paths0 = tp_paths0

	# if no new data - delete old data and print that this is not a valid TIC ID
	
	if lc_paths0 == 0:

		source.data = dict(
			x=[0],
			y=[0]
			)

		source_binned.data = dict(
			x_binned=[0],
			y_binned=[0]
			)

		source_md.data = dict(
			md_x1=[0],
			md_y1=[0],
			md_x2=[0],
			md_y2=[0]
			)

		source_zoom.data = dict(
			x_zoom=[0],
			y_zoom=[0]
			)

		zoom_points_binned.data = dict(
			x_zoom_binned=[0],
			y_zoom_binned=[0]
			)

		# give some indication that this is not a valid TIC ID
		#p.title.text_color="red"
		#p.title.text = "THIS TARGET HASN'T BEEN OSBERVED BY TESS"
		all_data = [-99]

		global twitter_logo_text1
		twitter_logo_text2 = "</style> \
		<script> \
		 \
		</script> \
		</head> \
		<body> \
		 \
		<!-- Add font awesome icons --> \
		<a href='https://twitter.com/intent/tweet?url=http://latte-online.flatironinstitute.org&text=Looking+at+TESS+light+curves+using+LATTE%21&hashtags=TESS,NASA,LATTE' target=_blank' class='fa fa-twitter'></a>\
		<a href='https://www.facebook.com/sharer.php?u=http://latte-online.flatironinstitute.org' target=_blank' class='fa fa-facebook'></a> \
		</body> \
		</html> \
		"
		
		#		<a href='https://www.linkedin.com/shareArticle?url=http://latte-online.flatironinstitute.org&title=LATTE&summary=Looking+at+TESS+light+curves+using+LATTE%21&source=LATTE' target=_blank' class='fa fa-linkedin'></a> \

		twitter_logo_text = twitter_logo_text1 + twitter_logo_text2

		twitter_logo.text = twitter_logo_text

		return all_data

	else:

		lc_paths = [datapath + '/' + lcp for lcp in lc_paths0]
		
		#print("download data")
		# extract the data from the fits files

		alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, ra, dec, TICID = LATTEutils.download_data(lc_paths,binfac = 10)

		#print("end download data")
		# change things into arrays because they're better

		# we can't handle nans so get rid of those
		finite_mask = np.isfinite(alltime) & np.isfinite(allflux) & np.isfinite(allflux_err)
		finite_binned = np.isfinite(alltimebinned) & np.isfinite(allfluxbinned) 

		alltime = alltime[finite_mask]
		allflux = allflux[finite_mask]
		allfbkg = allfbkg[finite_mask]
		allflux_err = allflux_err[finite_mask]
		
		# do some sigma clipping on the background flux so that it is easier to plot in the quick zoom tabs
		sig_clipping = (allfbkg < (np.median(allfbkg) + (5 * np.std(allfbkg))))
		bkg_clipped = allfbkg[sig_clipping]
		time_clipped = alltime[sig_clipping]

		#freq_smooth = smooth.frequency
		#power_smooth = smooth.power

		# add a sensical title to the plot about the parameters of the target we're looking at
		p.title.text_color="black"
		p.title.text = "TIC %s" %(TICID)
		
		#p.title.text = "%s		mag = %.2f		Teff = %.0f		Rstar = %.2f" %(TIC, tessmag, teff, srad)

		minf = np.nanmin(allflux)
		maxf = np.nanmax(allflux)
		height = maxf - minf

		md_x = all_md
		md_y1 = [minf]
		md_y2 = [minf + height*0.3]

		# -	 -	 -	 -	 -	 -
		# change/update the plotting data	text (i.e. delete the text)
		
		# welcome_text_source.data = pd.DataFrame.from_records([dict(text="", x = )])

		# -	 -	 -	 -	 -	 - UPDATE PLOT DATA -	 -	 -	 -	 -	 -
		
		# main data (binned + unbinned)
		source.data = dict(
			x=alltime,
			y=allflux
			)

		source_binned.data = dict(
			x_binned=alltimebinned[finite_binned],
			y_binned=allfluxbinned[finite_binned]
			)

		# update the plotting region axis limits and reset the 'sector selector'
		# start each not target by showing the whole data set and not an indiciual sector
		select_sector.value = 'all'
		p.x_range.start = np.nanmin(alltime)
		p.x_range.end = np.nanmax(alltime)

		# update red line to the start of the sector - this is where the zoom in will happen
		start_line = np.nanmin(alltime) + 5 # note that this is an entirely arbitrary offset from the start of the sector

		transit_time.location = start_line
		#plot_zoom.x_range.start = start_line - 1
		#plot_zoom.x_range.end = start_line + 1

		transit_mask = (alltime > (start_line - 1)) & (alltime < (start_line + 1))
		transit_mask_binned = (alltimebinned > (start_line - 1)) & (alltimebinned < (start_line + 1))
		alltime12_mask = (alltime12 > (start_line - 1)) & (alltime12 < (start_line + 1))

		# change the y axis limits if the centroid and backrgound plots - otherwise they're a bit crazy

		time_clipped_time_mask = (time_clipped > (start_line - 1)) & (time_clipped < (start_line + 1))
		
		#plot_bkg.y_range.start = np.nanmin(bkg_clipped[transit_time_mask])
		#plot_bkg.y_range.end.  = np.nanmax(bkg_clipped[transit_time_mask])

		source_zoom.data = dict(
			x_zoom=alltime[transit_mask],
			y_zoom=allflux[transit_mask]
			)

		source_zoom_binned.data = dict(
			x_zoom_binned=alltimebinned[transit_mask_binned],
			y_zoom_binned=allfluxbinned[transit_mask_binned]
			)

		# momentum dumps -	 -	 -	 -	 -	 -
		source_md.data = dict(
			md_x1=all_md,
			md_y1=[md_y1]*len(all_md), 
			md_x2 = all_md,
			md_y2 = [md_y2]*len(all_md)
			)

		# background data	-	 -	 -	 -	 -	 -
		source_bkg.data = dict(
			x_bkg=time_clipped[time_clipped_time_mask],
			y_bkg=bkg_clipped[time_clipped_time_mask]
			)

		# centroids
		source_xcen1.data = dict(
			x_xcen1=alltime12[alltime12_mask], 
			y_xcen1=allx1[alltime12_mask])

		source_xcen2.data = dict(
			x_xcen2=alltime12[alltime12_mask], 
			y_xcen2=allx2[alltime12_mask])


		source_ycen1.data = dict(
			x_ycen1=alltime12[alltime12_mask], 
			y_ycen1=ally1[alltime12_mask])

		source_ycen2.data = dict(
			x_ycen2=alltime12[alltime12_mask], 
			y_ycen2=ally2[alltime12_mask])


		# turn the time and flux into a lightkurve object to make the periodogram
		lc = lk.lightcurve.LightCurve(time = alltime, flux = allflux)
		ls = lc.to_periodogram(normalization="psd")
		
		best_period = (1/ls.frequency_at_max_power).to('d').value
		t0 = 0

		phased_time = np.array([-0.5+( ( t - t0-0.5*best_period) % best_period) / best_period for t in alltime]) * best_period

		#smooth = ls.smooth(method="boxkernel", filter_width=20.0) # remove the smoothing for now
		
		freq = ls.frequency
		power = ls.power

		# perioddogram -	 -	 -	 -	 -	 -
		source_periodgrm.data = dict(
			x_periodgrm=freq.value, 
			y_periodgrm=power.value
		)

		# evolutionary track -	 -	 -	 -	 -	 -
		source_periodgrm.data = dict(
			x_periodgrm=freq.value, 
			y_periodgrm=power.value
		)

		source_evol_target.data = dict(
			x_evol_target=[float(teff)],
			y_evol_target=[float(srad)]
			)

		# phase plot -	 -	 -	 -	 -	 -	 - 
		source_phased.data = dict(
			x_phase=phased_time,
			y_phase=allflux
			)

		source_phased_fit.data = dict(x_phase_fit=[], y_phase_fit=[])

		plot_zoom.x_range.start = start_line - 1
		plot_zoom.x_range.end = start_line + 1


		# delete anly selected transit times 
		tr_text.text = "{}".format(str(np.around(list(),2))[1:-1])

		# print the information about the target to the right of the target (instead of in the plot title?)

		exofop_url = "https://exofop.ipac.caltech.edu/tess/target.php?id={}".format(TICID)
		exofop_link = '<a href="%s"target="_blank">ExoFOP</a>' % (exofop_url)

		simbad_url = 'https://simbad.u-strasbg.fr/simbad/sim-coo?Coord={}+{}&Radius=2&Radius.unit=arcmin'.format(ra, dec)
		simbad_text = '<a href="%s"target="_blank">SIMBAD</a>' % (simbad_url)


		target_tic_text = r"$$%s$$" %(TICID)
		target_ra_text= r"$${}~^\circ$$".format(round(ra,5))
		target_dec_text = r"$${}~^\circ$$".format(round(dec,5))
		target_mag_text = r"$$%.2f~$$"%(tessmag)
		target_teff_text = r"$$%.0f~\text{K}$$"%(teff)
		target_radius_text = r"$$%.2f~\text{R}_{\odot}$$"%(srad)
		target_exofop_text = "{}".format(exofop_link)
		target_simbad_text = "{}".format(simbad_text)
		estimated_planet_rad = " "
		estimated_planet_rad2 = " "


		target_targetinfo_data.text = """ <br> <br> {} <br> {} <br> {} <br> {} <br> {} <br> {} <br> <br> {} <br> {} """.format(target_tic_text, target_ra_text, target_dec_text, target_mag_text, target_teff_text, target_radius_text, target_exofop_text, target_simbad_text)


		target_targetinfo_name_text  = " <!DOCTYPE html> \
		<html> \
		<style> \
		.tooltip { \
		  position: relative; \
		  display: inline-block; \
		  border-bottom: 1px dotted black; \
		} \
		 \
		.tooltip .tooltiptext { \
		  visibility: hidden; \
		  width: 300px; \
		  background-color: #293040; \
		  color: #fff; \
		  text-align: center; \
		  font-size: 14px; /* Set a font-size */ \
		  border-radius: 6px; \
		  padding: 5px 0; \
		   \
		  /* Position the tooltip */ \
		  position: absolute; \
		  z-index: 1; \
		  top: 100%; \
		  left: 50%; \
		  margin-left: -60px; \
		} \
		 \
		.tooltip:hover .tooltiptext { \
		  visibility: visible; \
		} \
		</style> \
		<body style='text-align:center;'>  "
		 
		text_end = " <b> Target Information <b> <b> <br> \
		<div class='tooltip'> <br> <b> TIC ID <b> <br> \
		  <span class='tooltiptext'>{}</span> \
		</div> \
		<br> <b> RA <br> Dec <br> TESS mag <br> Teff <br> Stellar Radius <br> <br> ExoFOP <br> SIMBAD <b>\
		 \
		</body> \
		</html> \
		".format(names_str)

		target_targetinfo_name.text = target_targetinfo_name_text + text_end

		
		twitter_logo_text2 = "</style> \
		<script> \
		 \
		</script> \
		</head> \
		<body> \
		 \
		<!-- Add font awesome icons --> \
		<a href='https://twitter.com/intent/tweet?url=http://latte-online.flatironinstitute.org%23ID%3D{}&text=Looking+at+the+light+curve+of+TIC+{}+using+LATTE%21&hashtags=TESS,NASA,LATTE' target=_blank' class='fa fa-twitter'></a>\
		<a href='https://www.facebook.com/sharer.php?u=http://latte-online.flatironinstitute.org%23ID%3D{}' target=_blank' class='fa fa-facebook'></a> \
		</body> \
		</html> \
		".format(p.title.text.split(' ')[1], p.title.text.split(' ')[1], p.title.text.split(' ')[1])
		
		#<a href='https://www.linkedin.com/shareArticle?url=http://latte-online.flatironinstitute.org%23ID%3D{}&title=LATTE&summary=Looking+at+the+light+curves+of+TIC+{}+using+LATTE%21&source=LATTE' target=_blank' class='fa fa-linkedin'></a> \

		twitter_logo_text = twitter_logo_text1 + twitter_logo_text2

		twitter_logo.text = twitter_logo_text


	# bundle up the data to pass on
	all_data = [alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, datapath, lc_paths,in_sec, TICID, ra, dec, phased_time]
		
	# take note that the data has finished being dowloaded 
	div_spinner_loading_data.visible = False

	return all_data


# - - - - - - - - - MAKE A SPINNER- - - - - - - - - - -

# run a spinner when the report is being made - to show the user that something is happening (progress bar for the future?)
# the spinner has been coded in java scipt and was added to bokeh using CustomJS


spinner_text_report = """
<!DOCTYPE html>
<html>
<head>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">
<meta name="viewport" content="width=device-width, initial-scale=0">
<style>
.buttonload {
  background-color: grey; /* grey background */
  border: none; /* Remove borders */
  border-radius: 20px; 
  color: white; /* White text */
  padding: 0px 20px; /* Some padding */
  font-size: 14px; /* Set a font-size */
  text-align: center;
}

/* Add a right margin to each icon */
.fa {
  margin-left: 0px;
  margin-right: 0px;
}
</style>
</head>
<body>

<button class="buttonload">
  <i class="fa fa-spinner fa-spin"></i>generating report...
</button>

</body>
</html>
"""

spinner_text_loading = """
<!DOCTYPE html>
<html>
<head>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">
<meta name="viewport" content="width=device-width, initial-scale=0">
<style>
.buttonload {
  background-color: grey; /* grey background */
  border: none; /* Remove borders */
  border-radius: 20px; 
  color: white; /* White text */
  padding: 0px 20px; /* Some padding */
  font-size: 14px; /* Set a font-size */
  text-align: center;
}

/* Add a right margin to each icon */
.fa {
  margin-left: 0px;
  margin-right: 0px;
}
</style>
</head>
<body>

<button class="buttonload">
  <i class="fa fa-spinner fa-spin"></i>loading data...
</button>

</body>
</html>
"""

div_spinner = Div(text=spinner_text_report,width=width_gadgets,height=40,visible=False, align = 'center')
div_spinner_loading_data = Div(text=spinner_text_loading,width=width_gadgets,height=40,visible=False, align = 'center')

#download_text = """
#<a href="bokeh_latte/static/LATTE_logo_small.png" download="w3logo">
#"""

cb = CustomJS(args=dict(div_spinner=div_spinner)
				,code='''
				console.log('cb triggered!')
				div_spinner.change.emit()''')


cb_data = CustomJS(args=dict(div_spinner=div_spinner_loading_data)
				,code='''
				console.log('cb_data triggered!')
				div_spinner_loading_data.change.emit()''')

# the title of the image is not needed but the code doesn't work without it so I'm goign to move on with my life and leave it there before I go insane trygin to figure out why
# important - the temperoary file needs to be stored in the static folder!

#tic=target_name_input, random_number=random_number

JScode_fetch = """
var filename = 'DV_report_';
filename = filename.concat(p.title.text.split(' ')[1])
filename = filename.concat('.pdf')
alert(filename);

var path = '/app/static/online_output_latte/'
path = [path + p.title.text.split(' ')[1] + '_rn' + rn_text.text + '/DV_report_' + p.title.text.split(' ')[1] + '.pdf'].join()

fetch(path, {cache: "no-store"}).then(response => response.blob())
					.then(blob => {
						//addresses IE
						if (navigator.msSaveBlob) {
							navigator.msSaveBlob(blob, filename);
						}

						else {
							var link = document.createElement("a");
							link = document.createElement('a')
							link.href = URL.createObjectURL(blob);
							window.open(link.href, '_blank');
							link.download = filename
							link.target = "_blank";
							link.style.visibility = 'hidden';
							link.dispatchEvent(new MouseEvent('click'))
						}
						return response.text();
					});
"""

# - - - - - - - - - - - - - 

# add the capability that if the size of a 'dummy' point is changed the report is downloaded 
# crucial step to download the report
#dummy_points.glyph.js_on_change('size', CustomJS(args=dict(p=p),code=JScode_fetch))


# ---------------------------------------
# ---------- CALLBACK FUNCTIONS ---------
# ---------------------------------------

#filename = filename.concat(target_name_input); filename = filename.concat('.pdf');
def make_report(all_data):

	'''
	This function triggers the full report being generated. (this triggers the 'brew' LATTE scipt)
	'''

	# only trigger the report making if at least one transit time has been added to the list (the button shouldn't be 'pressable' if this is not the case)

	if len(transit_time_list) > 0 :

		# start the spinner to show that things are happening
		div_spinner.visible = True
		
		def brew():
			# if the length of the data is == 1, there is no data and the report cannot be generated.

			if len(all_data) != 1: # add arguments (this is necessary because the LATTE legacy code had arguments so we will just simulated that)
				class Args:
					model = False
					noshow = True
					o = False
					nickname = False
					FFI = False
					save = True
					north = False
					mpi = False
					axis_split = True

				args=Args()
			
				simple = False	# we don't want to run the simple version - that's option is simply to do quick test
				save = True	# always save the files - no point running this if you're not going to save them
				DV = True	 # we'll make a DV report for all of them
				args.noshow = True	# don't show these, just save them
				args.auto = True # the code will automatically choose it's own apertures in this mode
				BLS = False
				model = False # modelling currently isn't availanbe - should we change that? 


				if plotting_options_menu.value == "split axis by hemispheres":
					args.axis_split = 150
				elif plotting_options_menu.value == "split axis by sectors":
					args.axis_split = 15
				elif plotting_options_menu.value == "plot sectors with events":
					args.axis_split = 0
				else:
					args.axis_split = -999

				#if plotting_options_menu.value == "split":

				#	args.axis_split = True
				#else:
				#	args.axis_split = False

				alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, datapath, lc_paths,sectors, tic, ra, dec, phased_time = all_data
				#	0		1			2			 3		 4			5				6		7		8	 9		 10			11		 12		 13		14		15		 16	17	 18		19		20	 21	22
				
				# run brew to make hte report
				global random_number
				#generate a random number that will be used to keep track of what the file is called (for storing)
				# this is in case multiple people run the same target at the same time - files need to be stored uniquely. 
				random_number = random.randint(1000000, 9999999) # 7 digit random number
				
				#print ("the random number is {}".format(random_number))
				LATTEbrew.brew_LATTE(tic, random_number, db, lc_paths, datapath, outpath, syspath, transit_time_list, simple, BLS, model, save, DV, sectors, sectors, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, np.nanmax(allflux), np.nanmin(allflux), tessmag, teff, srad, ra, dec, [-111], [-111], [-111], args)
				
				rn_text.text = str(random_number)
				# once the report has been made, remove the spinner again! 

				# then delete everything in the folder that isn't the pdf file 
				# be VERY careful when editing this line of code
				os.system('rm {}/{}_rn{}/{}_*.png'.format(outpath, tic, random_number,tic)) # THIS REMOVES GENERATED PNG FILES
				div_spinner.visible=False

		doc.add_next_tick_callback(brew)

		def change_size():
			# change size to trigger the download of the report 

			dummy_points.glyph.size = dummy_points.glyph.size + 0.1

		doc.add_next_tick_callback(change_size)

	else:
		pass


# when the figure data needs to be re-loaded
def update_data(lc_collection):
	'''
	Function to change the binning of the data. 
	# take in the unbinned data and rebin it to the desired cadence!
	'''

	time = lc_collection[0]
	flux = lc_collection[1]

	binfac = int(bin_factor.value)

	N = len(time)
	n = int(np.floor(N / binfac) * binfac)
	X = np.zeros((2, n))
	X[0, :] = time[:n]
	X[1, :] = flux[:n]
	Xb = rebin(X, (2, int(n / binfac)))

	time_binned = Xb[0]
	flux_binned = Xb[1]

	finite_binned = np.isfinite(time_binned) & np.isfinite(flux_binned) 

	source.data = dict(
		x=time,
		y=flux)

	source_binned.data = dict(
		x_binned=time_binned[finite_binned],
		y_binned=flux_binned[finite_binned]
		)

	# the axes go back to the full range, so make sure that the drop down menu reflects that		
	select_sector.value = 'all'

	p.x_range.start = np.nanmin(lc_collection[12])
	p.x_range.end = np.nanmax(lc_collection[13])


# button to change the TIC ID (data)

button_loaded = False
# main enter button to load multiple sectors
def button_callback():
	'''
	button to change the sectors that are available. 
	'''
	
	#print ("start_button")
	global source_data
	global button_loaded 

	button_loaded = True
	#print ("start spinner")

	div_spinner_loading_data.visible = True
	#print ("start data")

	source_data = select_data(single = False)
	
	#print ("end data")

	if len(source_data) != 1: # this means that the TIC ID was valid

		sectors = source_data[14]
	
		strings = ["{} - Sector {}".format(i, str(x)) for i,x in enumerate(sectors)]
	
		select_sector.options = ["all"] + list(strings)

# main enter button to load one sector
def button_callback_single():
	'''
	button to change the sectors that are available. 
	'''
	#print ("single")

	global button_loaded 
	global source_data


	button_loaded = True
	div_spinner_loading_data.visible = True

	source_data = select_data(single = True)
	
	if len(source_data) != 1: # this means that the TIC ID was valid

		sectors = source_data[14]
	
		strings = ["{} - Sector {}".format(i, str(x)) for i,x in enumerate(sectors)]
	
		select_sector.options = ["all"] + list(strings)

# callback to read the url
def url_callback():


	if (button_loaded == False):
		div_spinner_loading_data.visible = True
	
		global source_data
		source_data = select_data(single = False, url = True)
		
		if len(source_data) != 1: # this means that the TIC ID was valid
	
			sectors = source_data[14]
		
			strings = ["{} - Sector {}".format(i, str(x)) for i,x in enumerate(sectors)]
		
			select_sector.options = ["all"] + list(strings)

	else:
		#print ("nope")
		button_loaded == False

# choose a random tic id to plot
def button_lucky_callback():
	'''
	button to change the sectors that are available. 
	'''
	global button_loaded # don't ask why this needs to b repeated but the code won't work without this 
	global source_data
	
	button_loaded = True

	source_data = select_data(lucky = True)
	
	if len(source_data) != 1: # this means that the TIC ID was valid

		sectors = source_data[14]
	
		strings = ["{} - Sector {}".format(i, str(x)) for i,x in enumerate(sectors)]
	
		select_sector.options = ["all"] + list(strings)

global current_transit_time

# draw the line and move the zoom in window
def line_callback(event, source_data):
	
	'''
	Function to change the location of the line when the figure is double clicked
	'''
	
	global current_transit_time

	transit_time.location = event.x
	current_transit_time = event.x

	#plot_zoom.x_range.start = np.nanmin(current_transit_time - 1)
	#plot_zoom.x_range.end = np.nanmax(current_transit_time + 1)

	alltime = source_data[0]
	allflux = source_data[1]

	alltimebinned = source_data[4]
	allfluxbinned = source_data[5]

	allfbkg = source_data[11]

	sig_clipping = (allfbkg < (np.median(allfbkg) + (5 * np.std(allfbkg))))
	bkg_clipped = allfbkg[sig_clipping]
	time_clipped = alltime[sig_clipping]

	alltime12 = source_data[10]
	allx1 = source_data[6]
	allx2 = source_data[7]
	ally1 = source_data[8]
	ally2 = source_data[9]

	transit_mask = (alltime > (current_transit_time - 1)) & (alltime < (current_transit_time + 1))
	transit_mask_binned = (alltimebinned > (current_transit_time - 1)) & (alltimebinned < (current_transit_time + 1))
	alltime12_mask = (alltime12 > (current_transit_time - 1)) & (alltime12 < (current_transit_time + 1))

	# change the y axis limits if the centroid and backrgound plots - otherwise they're a bit crazy

	time_clipped_time_mask = (time_clipped > (current_transit_time - 1)) & (time_clipped < (current_transit_time + 1))
	
	#plot_bkg.y_range.start = np.nanmin(bkg_clipped[transit_time_mask])
	#plot_bkg.y_range.end.  = np.nanmax(bkg_clipped[transit_time_mask])

	
	source_zoom.data = dict(
		x_zoom=alltime[transit_mask],
		y_zoom=allflux[transit_mask]
		)

	source_zoom_binned.data = dict(
		x_zoom_binned=alltimebinned[transit_mask_binned],
		y_zoom_binned=allfluxbinned[transit_mask_binned]
		)

	# background data	-	 -	 -	 -	 -	 -
	source_bkg.data = dict(
		x_bkg=time_clipped[time_clipped_time_mask],
		y_bkg=bkg_clipped[time_clipped_time_mask]
		)

	# centroids
	source_xcen1.data = dict(
		x_xcen1=alltime12[alltime12_mask], 
		y_xcen1=allx1[alltime12_mask])

	source_xcen2.data = dict(
		x_xcen2=alltime12[alltime12_mask], 
		y_xcen2=allx2[alltime12_mask])


	source_ycen1.data = dict(
		x_ycen1=alltime12[alltime12_mask], 
		y_ycen1=ally1[alltime12_mask])

	source_ycen2.data = dict(
		x_ycen2=alltime12[alltime12_mask], 
		y_ycen2=ally2[alltime12_mask])


	plot_zoom.x_range.start = current_transit_time - 1
	plot_zoom.x_range.end = current_transit_time + 1

	#transit_time_mask = (source_data[0] > (current_transit_time - 1)) & (source_data[0] < (current_transit_time + 1))

	#time = source_data[0]
	#bkg = source_data[11]
	
	# change the y range of the background plot (because the background flux sometimes has strange scales)
	
	#if len(source_data[11][transit_time_mask] >0): # you can only do this if the transit time that has been selected is close to data...
	#	plot_bkg.y_range.start = np.nanmin(source_data[11][transit_time_mask])
	#	plot_bkg.y_range.end	 = np.nanmax(source_data[11][transit_time_mask])



global transit_time_list

transit_time_list = []

# add delete the event times
def callback_times(foo):
	'''
	Function to add/remove transit times from the printed list
	'''
	global current_transit_time
	global transit_time_list

	if type(current_transit_time) != list:

		if foo == 'add_time_button':
	
			transit_time_list.append(current_transit_time)
	
		else:
			if len(transit_time_list)>0:
				transit_time_list.pop()
	
		tr_text.text = "\n {}".format(str(np.around(transit_time_list,2))[1:-1])
		
		if len(transit_time_list) > 0:
			button_done.disabled=False
		else:
			button_done.disabled=True

# move through the sectorcs
def callback_sector_skip(foo):
	'''
	Function to click through the sectors (or choose to view them all)
	'''

	if foo == 'back_button':

		if select_sector.value == 'all':

			index = len(source_data[12]) - 1

			start_sec = float(source_data[12][index][0])
			end_sec = float(source_data[13][index][0])
			sector = int(source_data[14][index])
		
			p.x_range.start = start_sec
			p.x_range.end = end_sec
		
			select_sector.value = "{} - Sector {}".format(index, str(sector))

		else:
			index = int(select_sector.value.split(" ")[0])

			if index == 0:

				select_sector.value = 'all'
				p.x_range.start = np.nanmin(source_data[12])
				p.x_range.end = np.nanmax(source_data[13])

			else:
				index = index - 1

				start_sec = float(source_data[12][index][0])
				end_sec = float(source_data[13][index][0])
				sector = int(source_data[14][index])
				
				p.x_range.start = start_sec
				p.x_range.end = end_sec
				
				select_sector.value = "{} - Sector {}".format(index, str(sector))

	else:

		if select_sector.value == 'all':

			index = 0

			start_sec = float(source_data[12][index][0])
			end_sec = float(source_data[13][index][0])
			sector = int(source_data[14][index])
		
			p.x_range.start = start_sec
			p.x_range.end = end_sec
		
			select_sector.value = "{} - Sector {}".format(index, str(sector))


		else:
			index = int(select_sector.value.split(" ")[0])

			if index == len(source_data[12]) - 1:

				select_sector.value = 'all'
				p.x_range.start = np.nanmin(source_data[12])
				p.x_range.end = np.nanmax(source_data[13])

			else:
				index = index + 1

				start_sec = float(source_data[12][index][0])
				end_sec = float(source_data[13][index][0])
				sector = int(source_data[14][index])
		
				p.x_range.start = start_sec
				p.x_range.end = end_sec
		
				select_sector.value = "{} - Sector {}".format(index, str(sector))

# what sectors are available
def sector_function(source_data):
	

	if select_sector.value == 'all':
		p.x_range.start = np.nanmin(source_data[12])
		p.x_range.end = np.nanmax(source_data[13])

	else:	
		index = int(select_sector.value.split(" ")[0])
		sector = int(select_sector.value.split(" ")[-1])
	
		start_sec = float(source_data[12][index][0])
		end_sec = float(source_data[13][index][0])

		p.x_range.start = start_sec
		p.x_range.end = end_sec

# make the phase plot
def update_phase_data(source_data):

	# Get the current slider values


	try: # if no values have been entered don't do anything
		period = float(period_input.value)
		t0 = float(t0_input.value)
	
		phased_time = np.array([-0.5+( ( t - t0-0.5*period) % period) / period for t in np.array(source_data[0])]) * period
	
		source_phased.data = dict(x_phase=phased_time, y_phase=source_data[1])
		source_phased_fit.data = dict(x_phase_fit=[], y_phase_fit=[])
		estimated_planet_rad.text = " "
		estimated_planet_rad2.text = " "
	except:
		pass

# function for fitting 
def cosh_gauss_fit(x, *args):
		
		c0,c1,t0,d,gamma = args
		
		fit =c0 + c1 * ( 1 - pow( 1 - np.exp( 1 - np.cosh( (x-t0) / d ) ) ,gamma) ) 
		return fit

# update the transirt fit 
def update_transit_fit(source_data):

	# Get the current slider values
	try: # if no values have been entered don't do anything

		period = float(period_input.value)
		t0 = float(t0_input.value)
	
		# initial guesses...
		#c0 = 0 # offset in the y 
		#c1 = -1 # depth
		#t0 = 0
		#d = 0.5 # transit duration at the bottom of the transit
		#gamma = 0.7
	
		c0 = 0
		t0 = 0
	
		#c0 = 0 # offset in the y 
		#c1 = -1 # depth
		#t0 = 0
		#d = 0.5 # transit duration at the bottom of the transit
		#gamma = 0.7
	
		if menu_phase_fit.value == "Transit":
			p0 = [c0,-0.1,t0,0.3,3]
			transit = True
	
		elif	menu_phase_fit.value == "Short transit":
			p0 = [c0,-0.1,t0,0.1,2]
			transit = True
	
		elif	menu_phase_fit.value == "Long transit":
			p0 = [c0,-0.1,t0,0.5,3]
			transit = True
	
		elif	menu_phase_fit.value == "Eclipse":
			p0 = [c0,-0.1,t0,0.1,2]
			transit = False
	
		elif	menu_phase_fit.value == "Deep eclipse":
			p0 = [c0,-0.5,t0,0.5,2]
			transit = False
	
		time = source_phased.data['x_phase']
		flux = source_phased.data['y_phase']
		
		try:
			coeff, cm = curve_fit(cosh_gauss_fit, time, flux, p0=p0)
			coeff_error = np.sqrt(np.diag(cm))
			
			fine_time_axis = np.linspace(np.min(time),np.max(time),10000)
			flux_fit = cosh_gauss_fit(fine_time_axis, *coeff)
		
			source_phased_fit.data = dict(x_phase_fit=fine_time_axis, y_phase_fit=flux_fit)
		
			# given the transit depth, calculate the radius of the object
		
			srad = source_data[17] * 109.076 # convert to earth radii
	
			if transit == True:
		
				pl_rad = np.sqrt(coeff[1]*-1) * srad
		
				estimated_planet_rad.text = "Estimated radius:"
				estimated_planet_rad2.text = r"$$%.2f~\text{R}_{\odot}$$" % float(pl_rad)
			
			if transit == False:
				#target_simbad_text.text = r"$$\text{Estimated radius of eclipser}: {} \text{R}\_{odot}$$".format(simbad_text)
				estimated_planet_rad.text = " "
				estimated_planet_rad2.text = " "
	
		except:
			estimated_planet_rad.text = "Fit failed."
			estimated_planet_rad2.text = " "

	except: # if no values have been entered don't do anything
		pass

# make the PDF report downloadable


# ---------------------------------------
# -------- END CALLBACK FUNCTIONS -------
# ---------------------------------------



# ------------------------------------------------------------------------------
# - - - - - - - - - - -	MAKE BUTTONS + DROP DOWN MENUS- - - - - - - - - - - - - 
# ------------------------------------------------------------------------------


target_targetinfo_name = Div(text= target_targetinfo_name_text + text_end, width=width1, height=120, style={'font-size': fontsize, 'color': 'black'})
target_targetinfo_data = Div(text= """ """, width=width1, height=120, style={'font-size': fontsize, 'color': 'black'})

estimated_planet_rad  = Div(text=""" """, width=150, height=15, style={'font-size': fontsize, 'color': color2})
estimated_planet_rad2 = Div(text=""" """, width=150, height=15, style={'font-size': fontsize, 'color': color2})

# make a div for the text undenear the image
# logo links (PHT and PHCC)
pht_url = "https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess"
phcc_url = "https://www.planethunters.coffee"
pht_link = '<a href="%s"target="_blank">Coffee Chat</a>' % (phcc_url)
phcc_link = '<a href="%s"target="_blank">Planet Hunters TESS</a>'% (pht_url)

pht_text  = Div(text=pht_link,	width=125, height=20, style={'font-size': '90%', 'color': 'grey'})
phcc_text = Div(text=phcc_link, width=125, height=20, style={'font-size': '90%', 'color': 'grey'})

target_text_sep_name = Div(text="""<b>- - - - - - - - - - - - - - -<b>""", width=width1, height=30, style={'font-size': fontsize, 'color': 'darkorange'})
target_text_sep = Div(text=""" """, width=width1, height=30, style={'font-size': fontsize, 'color': 'darkorange'})

# ----

# buttons to add and delete time
add_time_button = Button(label = "Add time" ,width=width_gadgets)
add_time_button.on_click(partial(callback_times, foo="add_time_button"))

del_time_button = Button(label = "Delete time", width=width_gadgets)
del_time_button.on_click(partial(callback_times, foo="del_time_button"))

# ----

# button to move forward and backwards 
# flip through the sector options
next_button = Button(label = "Next", width=120)
next_button.on_click(partial(callback_sector_skip, foo="next_button"))

back_button = Button(label = "Back", width=120)
back_button.on_click(partial(callback_sector_skip, foo="back_button"))


# ----

joss_url = 'https://joss.theoj.org/papers/10.21105/joss.02101'
reference_text = "Copyright: Nora Eisner. If you use this tool for work that is being published please cite <a href='%s'target='_blank'>Eisner et al. 2020</a>." % joss_url
joss_text = Div(text=reference_text, style={'font-size': '75%', 'color': 'grey'})

# ----

# enter the TIC ID button 
button_single = Button(label = "Enter (latest sector only)", button_type = 'primary', width = width_gadgets - 42, name = 'SingleButton')
button_single.on_click(button_callback_single)

# ----

target_name_input = TextInput(title='Enter Simbad or TIC ID of a star:', width=width_gadgets) #, id = 'IDInput', name = 'IDInput')
#target_name_input.on_change('value', lambda attr, old, new: button_callback_single()) # if you press enter after entering an ID it will load one sector

# ----

button = Button(label = "All sectors (takes longer)", id = 'button_all_sectors', name = 'button_all_sectors', button_type = 'primary', width = width_gadgets)
button.on_click(button_callback)

# ----

#I'm feeling lucky button 
button_lucky = Button(label = "I'm feeling lucky", button_type = 'warning', width = width_gadgets)
button_lucky.on_click(button_lucky_callback)

# ----

# help button to primpt the user to add a TIC ID
tic_prompt_button = Button(label="?", width = 5)
tic_prompt_text = CustomJS(args={}, code='alert("(1) Enter a TIC ID or any ID recognised by SIMBAD. If the target was observed by TESS the data will be shown. (If you\'re feelig lucky, click the orange button below to select a random target!) (2) Double click on the figure to select an event time. (3) Press the Add time button to add the time to the list of event times to analyse. (4) Press the Make report button to generate a pdf summary.");')
tic_prompt_button.js_on_click(tic_prompt_text)

# ----

# final button to trigger report making
button_done = Button(label = "Make Report", button_type = 'primary', width=width_gadgets - 42)
button_done.on_event(ButtonClick, lambda event: make_report(source_data)) 

# ----

# help button 
report_prompt_button = Button(label="?", width = 5)
report_prompt_text = CustomJS(args={}, code='alert("Select transit times by clicking on the Add time button. Added times appear to the right of the light curve plots. Once at least one event time has beeen added a pdf report can be generated. Reports are automatically downloaded.");')
report_prompt_button.js_on_click(report_prompt_text)

# ----

# button to change the bin factor
bin_factor = Select(title="Bin factor", value="10",
				 options=["2", "5", '10', '20'], width=width_gadgets)
bin_factor.on_change('value', lambda attr, old, new: update_data(source_data))

# ----

# make a drop bown menu to select which sector to look at 
select_sector = Select(title="Select sector", value="all", options=["all"], width=width_gadgets)
select_sector.on_change('value', lambda attr, old, new: sector_function(source_data))

# ----

# the text to add/remove transit events
tr_text_name = Div(text="""<b>Event times<b> \n""", width=120, height=100, style={'font-size': '120%', 'color': 'black'})
tr_text = Div(text=""" """, width=100, height=100, style={'font-size': '120%', 'color': 'black'})
rn_text = Div(text=""" - """, width=1, height=1, style={'font-size': '0%', 'color': '0'}) # this is dummy text for the naming of the planets

# ----

# DUMMY INPUT FILE FOR THE URL 
# if the dummy changes - which changes when the url changes then reload the data
dummy_name_input = TextInput(value = '0', id = 'IDInput', name = 'IDInput', visible = False)
dummy_name_input.on_change('value', lambda attr, old, new: url_callback())

# ----

button_phase = Button(label = "Enter", button_type = 'primary', height = 31, width = 160)
button_phase.on_event(ButtonClick, lambda event: update_phase_data(source_data))

# ----

button_phase_fit = Button(label = "Fit transit", button_type = 'primary', height = 31, width = 110)
button_phase_fit.on_event(ButtonClick, lambda event: update_transit_fit(source_data))

# ----

# select what kind of fit - drop down menu
menu_phase_fit = Select(title="Signal type", value="Transit",
				 options=["Transit", "Short transit", "Long transit", 'Eclipse', 'Deep eclipse'], height = 31, width = 160)

# ----

# help button to prompt the user to add a TIC ID
fit_question_button = Button(label="?", width = 5)
fit_question_text = CustomJS(args={}, code='alert("Once you have phase folded the data, you can generate a very BASIC model for the fit. Note that the transits/eclipses have to be at a phase of 0, and that the model is NOT a full transit model and only an initial guess. If the model is sucessfull, it will provide an estimate of the planet radius. Note that this is only an approximation!");')
fit_question_button.js_on_click(fit_question_text)

# ----

#dummy_input = TextInput(value = '0', height = 0, width = 0)

t0_input = TextInput(title='T0', height = 35, width = 160)
period_input = TextInput(title='Period', height = 35, width = 160)


#target_name_input.on_change('value', lambda attr, old, new: url_callback())


# draw on the red line to indicate the transit if the figrue is double clicked
transit_time = Span(location=0, dimension='height', line_color='red', line_dash='dashed', line_width=3)

p.add_layout(transit_time)
plot_zoom.add_layout(transit_time)
plot_bkg.add_layout(transit_time)
plot_xcen.add_layout(transit_time)
plot_ycen.add_layout(transit_time)

p.on_event(DoubleTap, lambda event: line_callback(event, source_data))
plot_zoom.on_event(DoubleTap, lambda event: line_callback(event, source_data))

#quick_zoom_toggle.on_click(quick_zoom)

# initialize some parameters
# source_data = np.array(0)
current_transit_time = []

# slider for the marker size - lots of these as it needs to change all of the marker sizes (on all of the subplots)
marker_size_callback = CustomJS(args=dict(renderer=zoom_points), code="""
	renderer.glyph.size = cb_obj.value;
""")

marker_size_binned_callback = CustomJS(args=dict(renderer=zoom_points_binned), code="""
	renderer.glyph.size = cb_obj.value;
""")

# background
marker_size_callback_cen = CustomJS(args=dict(renderer=bkg_points), code="""
	renderer.glyph.size = cb_obj.value;
""")

# centroid marker size
marker_size_xcen1 = CustomJS(args=dict(renderer=xcen1_points), code="""
	renderer.glyph.size = cb_obj.value;
""")
marker_size_xcen2 = CustomJS(args=dict(renderer=xcen2_points), code="""
	renderer.glyph.size = cb_obj.value;
""")
marker_size_ycen1= CustomJS(args=dict(renderer=ycen1_points), code="""
	renderer.glyph.size = cb_obj.value;
""")
marker_size_ycen2 = CustomJS(args=dict(renderer=ycen2_points), code="""
	renderer.glyph.size = cb_obj.value;
""")

marker_size_phase = CustomJS(args=dict(renderer=phase_points), code="""
	renderer.glyph.size = cb_obj.value;
""")

# one slider to rule them all
slider = Slider(start=1, end=10, value=5, step=.1, title="Marker size", width=width_gadgets)
slider.js_on_change('value', marker_size_callback)
slider.js_on_change('value', marker_size_binned_callback)
slider.js_on_change('value', marker_size_callback_cen)
slider.js_on_change('value', marker_size_xcen1)
slider.js_on_change('value', marker_size_xcen2)
slider.js_on_change('value', marker_size_ycen1)
slider.js_on_change('value', marker_size_ycen2)
slider.js_on_change('value', marker_size_phase)


dummy_points.glyph.js_on_change('size', CustomJS(args=dict(p=p, rn_text=rn_text),code=JScode_fetch))

#------ SOCIAL MEDIA ---------



twitter_logo_text1 = """<!DOCTYPE html>
<html>
<head>
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">
<style>
.fa {
	padding: 5px;
	font-size: 30px;
	width: 30px;
	text-align: center;
	text-decoration: none;
	margin: 2px 2px;
	border-radius: 80%;
}
.fa:hover {
		opacity: 0.7;
}

.fa-twitter {
	background: #55ACEE;
	color: white;
}

.fa-facebook {
	background: #3B5998;
	color: white;
}
.fa-linkedin {
	background: #007bb5;
	color: white;
}
"""

twitter_logo_text2 = "</style> \
<script> \
 \
</script> \
</head> \
<body> \
 \
<!-- Add font awesome icons --> \
<a href='https://twitter.com/intent/tweet?url=http://latte-online.flatironinstitute.org&text=Looking+at+TESS+light+curves+using+LATTE%21&hashtags=TESS,NASA,LATTE' target=_blank' class='fa fa-twitter'></a>\
<a href='https://www.facebook.com/sharer.php?u=http://latte-online.flatironinstitute.org' target=_blank' class='fa fa-facebook'></a> \
</body> \
</html> \
"

#<a href='https://www.linkedin.com/shareArticle?url=http://latte-online.flatironinstitute.org&title=LATTE&summary=Looking+at+TESS+light+curves+using+LATTE%21&source=LATTE' target=_blank' class='fa fa-linkedin'></a> \

twitter_logo_text = twitter_logo_text1 + twitter_logo_text2

twitter_logo = Div(text=twitter_logo_text)



#  -----------------  change the url  ----------------- 
JScode_url = """
var url_name = '#ID=' + p.title.text.split(' ')[1];
window.location.hash = url_name;
"""

change_url = CustomJS(args=dict(p=p), code=JScode_url)
p.title.js_on_change('text', change_url)


# -----------------

# grey out the DV report button if the length is 1

if len(transit_time_list) > 0:
	button_done.disabled=False
else:
	button_done.disabled=True

# ----------------- LAYOUT ----------------- 

# make the 'other' plots used as diagnostics
# these will be shown in panels that can easily be flicked through for quick diagnostics
# once we have a quicker way to download data we want ALL of the diagnostics to be shown this way

panels = [None, None, None, None, None, None]

# Main panel: data
panels[0] = Panel(child=plot_zoom, title="Zoom Plot")
panels[1] = Panel(child=plot_bkg, title="Background Flux")
panels[2] = Panel(child=column( plot_xcen, plot_ycen, sizing_mode="stretch_width"),title="Centroid Position",)
panels[3] = Panel(child=plot_periodgrm, title="Periodogram")
panels[4] = Panel(child=plot_evol, title="Evolutionary tracks")

spacer1 = Spacer(width = 0, height = 18)

panels[5] = Panel(child=column(plot_phase, row(column(t0_input,spacer1,estimated_planet_rad), column(period_input, spacer1,estimated_planet_rad2), column(spacer1, button_phase), column(menu_phase_fit, spacer1, row(button_phase_fit, fit_question_button)))), title="Phase fold") # phase folde 

tabs = Tabs(tabs=panels)

#return column(plot, tabs, sizing_mode="stretch_width")
# to align in the middle align = "center"
instructions_txt = Div(text="""Double click on figure to select times.""", width=width_gadgets, height=10, align = 'start', style={'font-size': '90%', 'color': 'darkgrey'})

settings_txt = Div(text="""Settings - - - - - - - - - - - - - - - - - -""", width=width_gadgets, height=20, align = 'start', style={'font-size': '130%', 'color': 'darkorange'})

transit_times_txt = Div(text="""Select event times - - - - - - - - - - - -""", width=width_gadgets, height=20, align = 'start', style={'font-size': '130%', 'color': 'darkorange'})

#target_info_name = column(target_targetinfo_name,target_text_sep_name, target_tic_name,target_ra_name, target_dec_name,target_mag_name,target_teff_name,target_radius_name,target_exofop_name,target_simbad_name, tr_text_name)
#target_info_text = column(target_targetinfo_text,target_text_sep,target_tic_text,target_ra_text,target_dec_text, target_mag_text,target_teff_text,target_radius_text,target_exofop_text,target_simbad_text,tr_text)

space = Spacer(width = 15, height = 0)
h_space = Spacer(width = 0, height = 20)

target_info_name = column(h_space, target_targetinfo_name, Spacer(width = 0, height = 120), tr_text_name)
target_info_text = column(h_space, target_targetinfo_data, Spacer(width = 0, height = 120), tr_text)


l = column(row(column(latte_logo, dummy_name_input, target_name_input, div_error, div_error_astroquery, row(button_single, tic_prompt_button), button, button_lucky, instructions_txt, div_spinner, div_spinner_loading_data, settings_txt, bin_factor, slider, select_sector, row(back_button, next_button), transit_times_txt, add_time_button, del_time_button, row(button_done, report_prompt_button), rn_text, color_theme_menu, plotting_options_menu, pht_latte_logo, row(pht_text, phcc_text)), column(p,tabs, twitter_logo, joss_text), space, target_info_name, target_info_text), sizing_mode="scale_both")

#select_data()	# initial load of the data
curdoc().add_root(l)



