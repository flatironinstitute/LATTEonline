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
import sqlite3 as sql
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
from astroquery.mast import Catalogs
from astroquery.simbad import Simbad

# bokeh tools - - - - - - - - - - -
from bokeh.io import curdoc
from bokeh.resources import CDN
from bokeh.embed import file_html

from bokeh.models.widgets import Select
from bokeh.themes import built_in_themes
from bokeh.plotting import figure, curdoc
from bokeh.events import DoubleTap, ButtonClick
from bokeh.layouts import column, row, layout, Spacer
from bokeh.models import TextInput, DataTable, Dropdown, DateFormatter, TableColumn, Paragraph, Button, Slider, CustomJS, Range1d, Select, Toggle, ColumnDataSource, Panel, Tabs, Span, Div, Select, Slider, Div

# custom functions - - - - - - - - -
import LATTEutils, LATTEbrew

global target_name_input
global datapath
global syspath


# scp -r /Users/neisner/Documents/code/LATTEonline ccalin029:/mnt/home/neisner/projects/

# open the confidguration file so that this can be used on different machines as long as there is a cofnigureation file pointing to the LATTE output folder

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
				<img alt="LATTE Logo" title="Client Logo" src="/applatte/static/LATTE_imgs/LATTE_logo_small.png" style="width: 90%;max-height: 100%"/>
				</div>
	"""
	PHT_logo_text = """<div id="PHT logo">
				<img alt="PHT logo" title="PHT Logo" src="/applatte/static/LATTE_imgs/pht_logo_small.png" style="width: 83%;max-height: 83%"/>
				</div>
	"""


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

# START the document

doc = curdoc()

#curdoc().theme = 'dark_minimal'

# make it look pretty !!

# add the LATTE logo
latte_logo = Div(text=logo_text)

#tooltips = [
#			('Player', '3443'),
#			('Three-Pointers Made', '3424'),
#			('Three-Pointers Attempted', '3424'),
#			('Three-Point Percentage', '25243')   
#		   ]

#latte_logo.add_tools(HoverTool(tooltips=tooltips))

# add the PHT and PHCC logos 
pht_latte_logo = Div(text=PHT_logo_text)


# ------------------------------------------------------------------------------
# - - - - - - - - - - -  LOAD SQL TABLE FOR DATA ACCESS- - - - - - - - - - - - - 
# ------------------------------------------------------------------------------

DB_FILE = datapath + '/my_database.db'

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
		self.connection.close()

	# tid needs to be a int. Just cast it this way it takes care of leading 0's
	def search(self, tic_id):
		sql_select_by_tic = """ SELECT tic_id,sector_id,lc_filename,tp_filename FROM fits
							WHERE tic_id = ? """
		try:
			c = self.connection.cursor()
			c.execute(sql_select_by_tic, (tic_id,))
			row = c.fetchall() # change this to make it itterate over multiple 

			if len(row) == 0:
				print('Could not find tic', tic_id, "in data")
				return [0,0,0]
			
			sector = ([i[1] for i in row])
			lc_path = ([i[2] for i in row])
			tp_path = ([i[2] for i in row])

			return (sector, lc_path, tp_path)

			#return (row[0][1], row[0][2], row[0][3])

		except Error as e:
			print("INSERT failed with error", e)


# load the database - this is a quick way to search for all of the needed urls
db = Database(DB_FILE)

# ------------------------------------------------------------------------------
# - - - - - - - - - - -  MAKE BUTTONS + DROP DOWN MENUS- - - - - - - - - - - - - 
# ------------------------------------------------------------------------------

width_gadgets = 250

# enter the TIC ID button 
button = Button(label = "Enter", button_type = 'primary', width = width_gadgets - 42)

target_name_input = TextInput(title='Enter Simbad or TIC ID of a star:', width=width_gadgets)

# help button to primpt the user to add a TIC ID

tic_prompt_button = Button(label="?", width = 5)
tic_prompt_text = CustomJS(args={}, code='alert("(1) Enter a TIC ID or any ID recognised by SIMBAD. If the target was observed by TESS the data will be shown. (2) Double click on the figure to select a transit time. (3) Press the Add time button to add the time to the list of transit times to analyse. (4) Press the Make report button to generate a pdf summary.");')
tic_prompt_button.js_on_click(tic_prompt_text)

# button to enter 
button_done = Button(label = "Make Report", button_type = 'success', width=width_gadgets - 42)

report_prompt_button = Button(label="?", width = 5)
report_prompt_text = CustomJS(args={}, code='alert("Select transit times by clicking on the Add time button. Added times appear to the right of the light curve plots. Once at least one transit time has beeen added a pdf report can be generated. Reports are automatically downloaded.");')
report_prompt_button.js_on_click(report_prompt_text)


# button to change the bin factor
bin_factor = Select(title="Bin factor", value="10",
			   options=["2", "5", '10', '20'], width=width_gadgets)


# make a drop bown menu to select which sector to look at 
select_sector = Select(title="Select sector", value="all", options=["all"], width=width_gadgets)

# the text to add/remove transit events
tr_text_name = Div(text="""<b>Transit list<b> \n""", width=100, height=100, style={'font-size': '120%', 'color': 'black'})
tr_text = Div(text=""" """, width=100, height=100, style={'font-size': '120%', 'color': 'black'})
rn_text = Div(text=""" - """, width=1, height=1, style={'font-size': '0%', 'color': '0'}) # this is dummy text for the naming of the planets

# add text for the information of the target! Empty at the begingin but as you load data this will update

width1 = 150
width2 = 300

color1 = 'dimgray'
color2 = 'darkslategrey'
fontsize = '110%'
text_height = 15

target_targetinfo_name = Div(text=""" <b>Target Information <b>""", width=width1, height=10, style={'font-size': fontsize, 'color': 'black'})
target_tic_name = Div(text=""" TIC ID """, width=width1, height=text_height, style={'font-size': fontsize, 'color': color1})
target_ra_name = Div(text=""" RA  """, width=width1, height=text_height, style={'font-size': fontsize, 'color': color1})
target_dec_name = Div(text=""" Dec """, width=width1, height=text_height, style={'font-size': fontsize, 'color': color1})
target_mag_name = Div(text=""" TESS Mag """, width=width1, height=text_height, style={'font-size': fontsize, 'color': color1})
target_teff_name = Div(text=""" Teff """, width=width1, height=text_height, style={'font-size': fontsize, 'color': color1})
target_radius_name = Div(text=""" Stellar Radius """, width=width1, height=50, style={'font-size': fontsize, 'color': color1})
target_exofop_name = Div(text=""" ExoFOP """, width=width1, height=text_height, style={'font-size': fontsize, 'color': color1})
target_simbad_name = Div(text=""" SIMBAD """, width=width1, height=60, style={'font-size': fontsize, 'color': color1})

target_targetinfo_text = Div(text=""" """, width=width1, height=10, style={'font-size': fontsize, 'color': 'black'})
target_tic_text = Div(text=""" """, width=width2, height=text_height, style={'font-size': fontsize, 'color': color2})
target_ra_text = Div(text=""" """, width=width2, height=text_height, style={'font-size': fontsize, 'color': color2})
target_dec_text = Div(text=""" """, width=width2, height=text_height, style={'font-size': fontsize, 'color': color2})
target_mag_text = Div(text=""" """, width=width2, height=text_height, style={'font-size': fontsize, 'color': color2})
target_teff_text = Div(text=""" """, width=width2, height=text_height, style={'font-size': fontsize, 'color': color2})
target_radius_text = Div(text=""" """, width=width2, height=50, style={'font-size': fontsize, 'color': color2})
target_exofop_text = Div(text=""" """, width=width2, height=text_height, style={'font-size': fontsize, 'color': color2})
target_simbad_text = Div(text=""" """, width=width2, height=60, style={'font-size': fontsize, 'color': color2})

# make a div for the text undenear the image
# logo links (PHT and PHCC)
# make a div for the text undenear the image
# logo links (PHT and PHCC)
pht_url = "https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess"
phcc_url = "https://www.planethunters.coffee"
pht_link = '<a href="%s">Coffee Chat</a>' % (phcc_url)
phcc_link = '<a href="%s">Planet Hunters TESS</a>'% (pht_url)

pht_text = Div(text=pht_link,  width=125, height=20, style={'font-size': '90%', 'color': 'grey'})
phcc_text = Div(text=phcc_link, width=125, height=20, style={'font-size': '90%', 'color': 'grey'})

target_text_sep_name = Div(text="""<b>- - - - - - - - - - - - - - -<b>""", width=width1, height=30, style={'font-size': fontsize, 'color': 'darkorange'})
target_text_sep = Div(text=""" """, width=width1, height=30, style={'font-size': fontsize, 'color': 'darkorange'})

# buttons to add and delete time
add_time_button = Button(label = "Add time" ,width=width_gadgets)
del_time_button = Button(label = "Delete time", width=width_gadgets)

# button to move forward and backwards 
next_button = Button(label = "Next", width=120)
back_button = Button(label = "Back", width=120)

#quick_zoom_toggle = Toggle(label = "zoom", button_type = 'warning', width=100)

joss_url = 'https://joss.theoj.org/papers/10.21105/joss.02101'
reference_text = "Copyright: Nora Eisner. If you use this tool for work that is being published please cite <a href='%s'>Eisner et al. 2020.</a>." % joss_url
joss_text = Div(text=reference_text, style={'font-size': '75%', 'color': 'grey'})

# ------------------------------------------------------------------------------------
# - - - - - - - - - - -  DEFINE THE DATA AND PLOTTING REGION - - - - - - - - - - - - - 
# ------------------------------------------------------------------------------------

source = ColumnDataSource(data=dict(x=[], y=[])) # unbinned data
source_binned = ColumnDataSource(data=dict(x_binned=[], y_binned=[])) # binned data
source_md = ColumnDataSource(data=dict(md_x1 = [], md_y1 =[], md_x2 = [], md_y2 =[])) # momentun dumps
source_dummy = ColumnDataSource(data=dict(x_dummy=[1e2], y_dummy=[1e-10])) # dummy data


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
	source = source_md,
	color="blue", 
	line_width=1)



p.xaxis.axis_label = "Time (BJD - 2457000)"
p.yaxis.axis_label = "Normalised Flux"


# when LATTE is first launched, we want to prompt the user to enter a TIC ID on the left 
# to do that we will put some text on the empty plot 

welcome_text_source = ColumnDataSource(
			pd.DataFrame.from_records([dict(
			x=100,
			y=50,
			text="Enter a TIC ID to start.",
			color="black")]))


welcome_text = p.text(
		x="x",
		y="y",
		text="text",
		text_align="center",
		text_font_size="20px",
		text_font="helvetica",
		text_font_style="normal",
		source= welcome_text_source)



# - - - - -  ZOOM TAB - - - - - - 

source_zoom = ColumnDataSource(data=dict(x_zoom=[], y_zoom=[]))


plot_zoom = figure(
	plot_height=300,
	min_width=600,
	min_height=300,
	title="",
	sizing_mode="stretch_both",
	x_range=(0,5),
	toolbar_location="below", 
	tools='box_zoom,wheel_zoom,pan,reset'
)


zoom_points  = plot_zoom.circle(x="x", 
	y="y", 
	source=source, 
	size=5, 
	color="darkorange", 
	line_color=None, 
	alpha=0.8,
	)


zoom_points_binned = plot_zoom.circle(x="x_binned", 
	y="y_binned", 
	source=source_binned, 
	size=5, 
	color="black", 
	line_color=None, 
	alpha=0.5)


zoom_segment = plot_zoom.segment(x0 = "md_x1", 
	y0 = "md_y1",
	x1 = "md_x2",
	y1 = "md_y2",
	source = source_md,
	color="blue", 
	line_width=1)


# add x and y labels
plot_zoom.xaxis.axis_label = "Time (BJD - 2457000)"
plot_zoom.yaxis.axis_label = "Normalised Flux"


# - - - - -  BACKGROUND TAB - - - - - - 

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
plot_bkg.xaxis.axis_label = "Time (BJD - 2457000)"
plot_bkg.yaxis.axis_label = "Background Flux"



# - - - - -  CENTROID TAB - - - - - - 

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
plot_xcen.yaxis.axis_label = "x-centroid"
plot_ycen.xaxis.axis_label = "Time (BJD - 2457000)"
plot_ycen.yaxis.axis_label = "y-centroid"


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



# add the smoothes line to the periodogram (not needed most of the time so delete for now)
#plot_periodgrm.line(
#	x="x_periodgrm_smooth",
#	y="y_periodgrm_smooth",
#	source=source_periodgrm_smooth,
#	line_color="red",
#	color="red",
#	alpha=1,
#)

# add x and y labels
plot_periodgrm.xaxis.axis_label = r"$$\text{Frequency}~(\mu Hz)$$"
plot_periodgrm.yaxis.axis_label = (r"$$\text{Power Density}~(A{^2} \mu Hz^{-1})$$")


# - - - -  EVOLUTIONARY TRACKS TAB - - - - 

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
	color="maroon", 
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

# - - - - - - - - - - - - - - - - 

# function to select the data  - only run when the enter button is pressed (takes a while at the moment to download data)
# to do: data locally or need s pinnign wheel, progress bar etc


color_theme_menu = Select(title="Color theme", value="light",
			   options=['light', 'dark'], width=width_gadgets)


def color_theme():

	if color_theme_menu.value == 'light':

		for plot in [p, plot_zoom, plot_bkg, plot_xcen, plot_ycen, plot_periodgrm, plot_evol]:
			plot.background_fill_color = 'white'
			plot.outline_line_color = 'black'
			plot.ygrid.grid_line_color = 'gainsboro'
			plot.xgrid.grid_line_color = 'gainsboro'

		binned_points.glyph.fill_color = 'black'
		main_points.glyph.fill_color = 'darkorange'
		main_segment.glyph.line_color = 'blue'

		zoom_points.glyph.fill_color = 'darkorange' 
		zoom_points_binned.glyph.fill_color = 'black'
		zoom_segment.glyph.line_color = 'blue'

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


	elif color_theme_menu.value == 'dark':

		for plot in [p, plot_zoom, plot_bkg, plot_xcen, plot_ycen, plot_periodgrm, plot_evol]:
			plot.background_fill_color = 'darkslategray'
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
		zoom_segment.glyph.line_color = 'paleturquoise'

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

# plot_zoom, plot_evol

color_theme_menu.on_change('value', lambda attr, old, new: color_theme())


# -----------------------------------
# - - - - - - FUNCTIONS - - - - - - -
# -----------------------------------

def select_data():

	'''
	Function to load the data needed to plot the initial LCs displayed in the interface.
	'''

	global transit_time_list

	# delete the markings
	transit_time_list = []
	tr_text.text = "{}".format(str(np.around(transit_time_list,2))[1:-1])

	# -----------------
	
	# grey out the DV report button if no  (i.e. )
	if len(transit_time_list) > 0:
		button_done.disabled=False
	else:
		button_done.disabled=True

	target_name = str(target_name_input.value)
	#TIC = 'TIC {}'.format(TICID)

	if target_name[0:3] == 'TIC' or target_name[0:3] == 'tic':

		TICID = str(int(target_name[3:]))
		TIC = "TIC " + str(TICID) # numerical only

		_ , lc_paths0, tp_paths0 = db.search(int(TICID))

	else:
		target_names_table = Simbad.query_objectids(target_name)
		
		try: # see if the name is a simbad name
			results = list(target_names_table['ID'])
			
			for r in results: 
				if 'TIC' in r:
					TIC = r # with the words
		
			TICID = TIC[4:] # numerical only

			_ , lc_paths0, tp_paths0 = db.search(int(TICID))

		except:
			
			TICID = int(target_name)
			TIC = "TIC " + str(TICID) # numerical only
			_ , lc_paths0, tp_paths0 = db.search(int(TICID))
			

			if lc_paths0 == 0: # if the length of the path is still zero then the target has not been osberved by tess
				#p.title.text_color="red"
				#p.title.text = "THIS TARGET HASN'T BEEN OSBERVED BY TESS"
				all_data = [-99]
				div_error.visible=True

				return all_data

	div_error.visible=False

	lc_paths0 = lc_paths0
	tp_paths0 = tp_paths0

	# if no new data - delete old data and print that this is not a valid TIC ID
	if len(lc_paths0) == 0:

		source.data = dict(
			x=[0],
			y=[0]
			)

		source_binned.data = dict(
			x_binned=[0],
			y_binned=[0]
			)

		source_md.data = dict(
			x=[0],
			y=[0], 
			xm01 = [0],
			ym01 = [0]
			)

		# give some indication that this is not a valid TIC ID
		#p.title.text_color="red"
		#p.title.text = "THIS TARGET HASN'T BEEN OSBERVED BY TESS"
		all_data = [-99]

		return all_data

	else:

		lc_paths = [datapath + '/' + lcp for lcp in lc_paths0]
		
		# extract the data from the fits files
		alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, ra, dec = LATTEutils.download_data(lc_paths)

		# change things into arrays because they're better
		
		allx1 = np.array(allx1)
		allx2 = np.array(allx2)
		ally1 = np.array(ally1)
		ally2 = np.array(ally2)
		allfbkg = np.array(allfbkg)
		alltime = np.array(alltime)
		allflux = np.array(allflux)
		alltime12 = np.array(alltime12)
		allflux_err = np.array(allflux_err)
		alltimebinned = np.array(alltimebinned)
		allfluxbinned = np.array(allfluxbinned)
		
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


		# turn the time and flux into a lightkurve object to make the periodogram
		lc = lk.lightcurve.LightCurve(time = alltime, flux = allflux)
		ls = lc.to_periodogram(normalization="psd")
		
		#smooth = ls.smooth(method="boxkernel", filter_width=20.0) # remove the smoothing for now
		
		freq = ls.frequency
		power = ls.power

		#freq_smooth = smooth.frequency
		#power_smooth = smooth.power

		# add a sensical title to the plot about the parameters of the target we're looking at
		p.title.text_color="black"
		p.title.text = "%s" %(TIC)
		#p.title.text = "%s		mag = %.2f		Teff = %.0f		Rstar = %.2f" %(TIC, tessmag, teff, srad)


		minf = np.nanmin(allflux)
		maxf = np.nanmax(allflux)
		height = maxf - minf

		md_x = all_md
		md_y1 = [minf]
		md_y2 = [minf + height*0.3]

		# -   -   -   -   -   -

		# change/update the plotting data  text (i.e. delete the text)
		welcome_text_source.data = pd.DataFrame.from_records([dict(text="")])


		# -   -   -   -   -   - UPDATE PLOT DATA -   -   -   -   -   -
		
		# main data (binned + unbinned)
		source.data = dict(
			x=alltime,
			y=allflux
			)

		source_binned.data = dict(
			x_binned=alltimebinned[finite_binned],
			y_binned=allfluxbinned[finite_binned]
			)

		# momentum dumps -   -   -   -   -   -
		source_md.data = dict(
			x=all_md,
			y=[md_y1]*len(all_md), 
			xm01 = all_md,
			ym01 = [md_y2]*len(all_md)
			)

		# background data  -   -   -   -   -   -
		source_bkg.data = dict(
			x_bkg=time_clipped,
			y_bkg=bkg_clipped
			)

		# centroids
		source_xcen1.data = dict(
			x_xcen1=alltime12, 
			y_xcen1=allx1)

		source_xcen2.data = dict(
			x_xcen2=alltime12, 
			y_xcen2=allx2)


		source_ycen1.data = dict(
			x_ycen1=alltime12, 
			y_ycen1=ally1)

		source_ycen2.data = dict(
			x_ycen2=alltime12, 
			y_ycen2=ally2)

		# perioddogram -   -   -   -   -   -
		source_periodgrm.data = dict(
			x_periodgrm=freq.value, 
			y_periodgrm=power.value
		)

		# evolutionary track -   -   -   -   -   -
		source_periodgrm.data = dict(
			x_periodgrm=freq.value, 
			y_periodgrm=power.value
		)


		source_evol_target.data = dict(
			x_evol_target=[float(teff)],
			y_evol_target=[float(srad)]
			)

		#source_periodgrm_smooth.data = dict(
		#	x_periodgrm_smooth=freq_smooth.value, 
		#	y_periodgrm_smooth=power_smooth.value
		#)


		# update the plotting region axis limits and reset the 'sector selector'
		# start each not target by showing the whole data set and not an indiciual sector
		select_sector.value = 'all'
		p.x_range.start = np.nanmin(start_sec)
		p.x_range.end = np.nanmax(end_sec)

		
		# update red line to the start of the sector - this is where the zoom in will happen
		start_line = np.nanmin(alltime) + 5 # note that this is an entirely arbitrary offset from the start of the sector

		transit_time.location = start_line
		plot_zoom.x_range.start = start_line - 1
		plot_zoom.x_range.end = start_line + 1


		# change the y axis limits if the centroid and backrgound plots - otherwise they're a bit crazy

		transit_time_mask = (time_clipped > (start_line - 1)) & (time_clipped < (start_line + 1))
		plot_bkg.y_range.start = np.nanmin(bkg_clipped[transit_time_mask])
		plot_bkg.y_range.end   = np.nanmax(bkg_clipped[transit_time_mask])

		# delete anly selected transit times 
		tr_text.text = "{}".format(str(np.around(list(),2))[1:-1])

		# print the information about the target to the right of the target (instead of in the plot title?)


		exofop_url = "https://exofop.ipac.caltech.edu/tess/target.php?id={}".format(TICID)
		exofop_link = '<a href="%s">ExoFOP</a>' % (exofop_url)

		simbad_url = 'https://simbad.u-strasbg.fr/simbad/sim-coo?Coord={}+{}&Radius=2&Radius.unit=arcmin'.format(ra, dec)
		simbad_text = '<a href="%s">SIMBAD</a>' % (simbad_url)


		target_tic_text.text = r"$$%s$$" %(TICID)
		target_ra_text.text = r"$${}~^\circ$$".format(round(ra,5))
		target_dec_text.text = r"$${}~^\circ$$".format(round(dec,5))
		target_mag_text.text = r"$$%.2f  $$"%(tessmag)
		target_teff_text.text = r"$$%.0f \text{K}$$"%(teff)
		target_radius_text.text = r"$$%.2f~\text{R}_{\odot}$$"%(srad)
		target_exofop_text.text = "{}".format(exofop_link)
		target_simbad_text.text = "{}".format(simbad_text)


	# bundle up the data to pass on
	all_data = [alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, datapath, lc_paths,in_sec, TICID, ra, dec]
	

	return all_data



# - - - - - - - - - MAKE A SPINNER- - - - - - - - - - -

# run a spinner when the report is being made - to show the user that something is happening (progress bar for the future?)
# the spinner has been coded in java scipt and was added to bokeh using CustomJS

spinner_text = """
<!-- https://www.w3schools.com/howto/howto_css_loader.asp -->
<div class="loader">
<style scoped>
.loader {
	border: 16px solid #f3f3f3; /* Light grey */
	border-top: 16px solid orange; /* orange */
	border-radius: 50%;
	width: 100px;
	height: 100px;
	animation: spin 5s linear infinite;
}

@keyframes spin {
	0% { transform: rotate(0deg); }
	100% { transform: rotate(360deg); }
} 
</style>
</div>
"""

div_spinner = Div(text=spinner_text,width=120,height=120,visible=False, align = 'center')


#download_text = """
#<a href="bokeh_latte/static/LATTE_logo_small.png" download="w3logo">
#"""


cb = CustomJS(args=dict(div_spinner=div_spinner)
			  ,code='''
			  console.log('cb triggered!')
			  div_spinner.change.emit()''')


# the title of the image is not needed but the code doesn't work without it so I'm goign to move on with my life and leave it there before I go insane trygin to figure out why
# important - the temperoary file needs to be stored in the static folder!

#tic=target_name_input, random_number=random_number

JScode_fetch = """
var filename = 'DV_report_';
filename = filename.concat(p.title.text.split(' ')[1])
filename = filename.concat('.pdf')
alert(filename);

var path = '/applatte/static/online_output_latte/'
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
							URL.revokeObjectURL(url);
						}
						return response.text();
					});

"""


# - - - - - - - - - - - - - 

# add the capability that if the size of a 'dummy' point is changed the report is downloaded 
# crucial step to download the report
#dummy_points.glyph.js_on_change('size', CustomJS(args=dict(p=p),code=JScode_fetch))


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
			
				args=Args()
			
				simple = False  # we don't want to run the simple version - that's option is simply to do quick test
				save = True  # always save the files - no point running this if you're not going to save them
				DV = True   # we'll make a DV report for all of them
				args.noshow = True  # don't show these, just save them
				args.auto = True # the code will automatically choose it's own apertures in this mode
				BLS = False
				model = False # modelling currently isn't availanbe - should we change that? 
				
	
				alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, tessmag, teff, srad, datapath, lc_paths,sectors, tic, ra, dec = all_data
				#  0	  1		  2			 3		 4		  5			  6	  7	  8	 9	   10		  11		 12	   13	  14	  15	   16	17   18	  19	  20   21  22
	
				# run brew to make hte report
				global random_number
				#generate a random number that will be used to keep track of what the file is called (for storing)
				# this is in case multiple people run the same target at the same time - files need to be stored uniquely. 
				random_number = random.randint(1000000, 9999999) # 7 digit random number
				
				#print ("the random number is {}".format(random_number))

				LATTEbrew.brew_LATTE(tic, random_number, db, lc_paths, datapath, outpath, syspath, transit_time_list, simple, BLS, model, save, DV, sectors, sectors, alltime, allflux, allflux_err, all_md, alltimebinned, allfluxbinned, allx1, allx2, ally1, ally2, alltime12, allfbkg, start_sec, end_sec, in_sec, np.nanmax(allflux), np.nanmin(allflux), tessmag, teff, srad, ra, dec, [-111], [-111], [-111], args)
				
				rn_text.text = str(random_number)
				# once the report has been made, remove the spinner again! 
				div_spinner.visible=False


		
		doc.add_next_tick_callback(brew)

		def change_size():
			# change size to trigger the download of the report 

			dummy_points.glyph.size = dummy_points.glyph.size + 0.1

		doc.add_next_tick_callback(change_size)

	else:
		print ("need at least one transit in the list")


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
def button_callback():
	'''
	button to change the sectors that are available. 
	'''
	global source_data
	source_data = select_data()
	
	if len(source_data) != 1: # this means that the TIC ID was valid

		sectors = source_data[14]
	
		strings = ["{} - Sector {}".format(i, str(x)) for i,x in enumerate(sectors)]
	
		select_sector.options = ["all"] + list(strings)


global current_transit_time


def line_callback(event, source_data):
	'''
	Fucntion to change the location of the line when the figure is double clicked
	'''
	global current_transit_time

	transit_time.location = event.x
	current_transit_time = event.x

	plot_zoom.x_range.start = np.nanmin(current_transit_time - 1)
	plot_zoom.x_range.end = np.nanmax(current_transit_time + 1)

	
	transit_time_mask = (source_data[0] > (current_transit_time - 1)) & (source_data[0] < (current_transit_time + 1))

	#time = source_data[0]
	#bkg = source_data[11]
	
	# change the y range of the background plot (because the background flux sometimes has strange scales)
	
	if len(source_data[11][transit_time_mask] >0): # you can only do this if the transit time that has been selected is close to data...
		plot_bkg.y_range.start = np.nanmin(source_data[11][transit_time_mask])
		plot_bkg.y_range.end   = np.nanmax(source_data[11][transit_time_mask])


# - - - - - - - - - - 

global transit_time_list

transit_time_list = []

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


# make the PDF report downloadable


# - - - - - - - - LINK BUTTONS TO FUNCTIODNS - - - - - - - - 

# get the data for the given LC!
button.on_click(button_callback)

# change the bin factor
bin_factor.on_change('value', lambda attr, old, new: update_data(source_data))

# add.remove times from list
del_time_button.on_click(partial(callback_times, foo="del_time_button"))
add_time_button.on_click(partial(callback_times, foo="add_time_button"))

# flip through the sector options
back_button.on_click(partial(callback_sector_skip, foo="back_button"))
next_button.on_click(partial(callback_sector_skip, foo="next_button"))

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


# one slider to rule them all
slider = Slider(start=1, end=10, value=5, step=.1, title="Marker size", width=width_gadgets)
slider.js_on_change('value', marker_size_callback)
slider.js_on_change('value', marker_size_binned_callback)
slider.js_on_change('value', marker_size_callback_cen)
slider.js_on_change('value', marker_size_xcen1)
slider.js_on_change('value', marker_size_xcen2)
slider.js_on_change('value', marker_size_ycen1)
slider.js_on_change('value', marker_size_ycen2)


# slider to change the zoom window size 
# slider_zoom_window = Slider(start=0.5, end=5, value=1, step=.1, title="Zoom window size")


# change bin factor
bin_factor.on_change('value', lambda attr, old, new: update_data(source_data))

# change sector selection
select_sector.on_change('value', lambda attr, old, new: sector_function(source_data))

# final button to trigger report makign
button_done.on_event(ButtonClick, lambda event: make_report(source_data)) 


dummy_points.glyph.js_on_change('size', CustomJS(args=dict(p=p, rn_text=rn_text),code=JScode_fetch))


# -----------------

# grey out the DV report button if the length is 1

if len(transit_time_list) > 0:
	button_done.disabled=False
else:
	button_done.disabled=True

# -----------------

# make the 'other' plots used as diagnostics
# these will be shown in panels that can easily be flicked through for quick diagnostics
# once we have a quicker way to download data we want ALL of the diagnostics to be shown this way

panels = [None, None, None, None, None]

# Main panel: data
panels[0] = Panel(child=plot_zoom, title="Zoom Plot")
panels[1] = Panel(child=plot_bkg, title="Background Flux")
panels[2] = Panel(child=column( plot_xcen, plot_ycen, sizing_mode="stretch_width"),title="Centroid Position",)
panels[3] = Panel(child=plot_periodgrm, title="Periodogram")
panels[4] = Panel(child=plot_evol, title="Evolutionary tracks")




tabs = Tabs(tabs=panels)

#return column(plot, tabs, sizing_mode="stretch_width")
# to align in the middle align = "center"
settings_txt = Div(text="""Settings - - - - - - - - - - - - - - - - - -""", width=width_gadgets, height=20, align = 'start', style={'font-size': '130%', 'color': 'darkorange'})

transit_times_txt = Div(text="""Select transit times - - - - - - - - - - - -""", width=width_gadgets, height=20, align = 'start', style={'font-size': '130%', 'color': 'darkorange'})


target_info_name = column(target_targetinfo_name,target_text_sep_name, target_tic_name,target_ra_name, target_dec_name,target_mag_name,target_teff_name,target_radius_name,target_exofop_name,target_simbad_name, tr_text_name)
target_info_text = column(target_targetinfo_text,target_text_sep,target_tic_text,target_ra_text,target_dec_text, target_mag_text,target_teff_text,target_radius_text,target_exofop_text,target_simbad_text,tr_text)

#quick_zoom_toggle -- removed quick zoom toggle for now (can always add back in)

space = Spacer(width = 15, height = 0)

l = column(row(column(latte_logo, target_name_input, div_error, row(button, tic_prompt_button), div_spinner, settings_txt, bin_factor, slider, select_sector, row(back_button, next_button), transit_times_txt, add_time_button, del_time_button, row(button_done, report_prompt_button), rn_text, color_theme_menu, pht_latte_logo, row(pht_text, phcc_text)), column(p,tabs, joss_text), space, target_info_name, target_info_text), sizing_mode="scale_both")

#select_data()  # initial load of the data
curdoc().add_root(l)


html = file_html(l, CDN, "latte_html")

with open('latte_html.txt', 'w') as f:
	f.write(html)

#
from bokeh.embed import server_document
script = server_document("https://demo.bokeh.org/sliders")

with open('latte_html_script.txt', 'w') as f:
	f.write(script)

















# -----------


# NOT USED AT THE MOMENT BUT KEEP INCASE WE WANT IT BACK
def quick_zoom(event):

	'''CURRENTLY NOT USED'''

	try:	
		if quick_zoom_toggle.active == True:
			
			line_loc = transit_time.location
			p.x_range.start = line_loc - 1
			p.x_range.end =line_loc + 1
		
			#main_points.glyph.size = 5
			#binned_points.glyph.size = 5

		if quick_zoom_toggle.active == False:
			

			line_loc = transit_time.location
			p.x_range.start = line_loc - 1
			p.x_range.end =line_loc + 1
	
			if select_sector.value == 'all':
	
				select_sector.value = 'all'
				p.x_range.start = np.nanmin(source_data[12])
				p.x_range.end = np.nanmax(source_data[13])
	
	
			else:
				index = int(select_sector.value.split(" ")[0])
	
				start_sec = float(source_data[12][index][0])
				end_sec = float(source_data[13][index][0])
				sector = int(source_data[14][index])
			
				p.x_range.start = start_sec
				p.x_range.end = end_sec
			
				select_sector.value = "{} - Sector {}".format(index, str(sector))
			
			#main_points.glyph.size = 3
			#binned_points.glyph.size = 3

	except:

		print ("must select a point in the LC first")



JScode_fetch_2 = """
var filename = p.title.text; 
filename = 'test.pdf';
alert(filename);

fetch('/applatte/static/DV_report_470710327.pdf', {cache: "no-store"}).then(response => response.blob())
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
							URL.revokeObjectURL(url);
						}
						return response.text();
					});

"""