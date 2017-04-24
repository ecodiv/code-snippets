#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
DESCRIPTION:  A summary of the code used in the tutorial Plotting GRASS data
              in Python (part 1) to plot distance to schools in south-west
              Wake County in North Carolina of the different municipals.
              The code below can be used to create Figure 4 and 10 in the
              tutorial.
              Please see https://tutorials.ecodiv.earth/grass-python-plots.html
              for details.
NOTES:        The code uses the NC sample data set (you can get this data from
              the download section on https://grass.osgeo.org/). To run the
              code, first start GRASS GIS in the NC location (preferably in
              a new mapset) and then start python from the command line. Code
              should work on GRASS GIS 7.2 +

@author: pvbreugel add ecodiv dot earth
"""

# Load libraries
import os
import grass.script as gs
import pandas as pd
from pandas.io import sql
import matplotlib.pyplot as plt
import numpy as np
import sqlite3
from tabulate import tabulate
os.chdir("/home/paulo/")

#------------------------------------------------------------------------------
# Prepare the vector layers
#------------------------------------------------------------------------------
gs.run_command("v.extract", input="boundary_county", where="NAME='WAKE'",
               output="WakeCounty")
gs.run_command("v.overlay", ainput="boundary_municp", binput="WakeCounty",
               operator="and", output="WakeMunicp", olayer="0,1,0")
sqlstat = "MB_NAME='Durham' OR MB_NAME IS NULL"
gs.run_command("v.db.droprow", input="WakeMunicp", where=sqlstat,
            output="WakeMunicp2")
gs.run_command("g.rename", vector="WakeMunicp2,WakeMunicp", overwrite=True)

#------------------------------------------------------------------------------
# Create the map (Figure 4 in tutorial)
#------------------------------------------------------------------------------

# Get the Paired-12 color palet
from palettable.colorbrewer.qualitative import Paired_12
cls = Paired_12.colors
csl = [":".join(map(str, x)) for x in cls]

# Write color codes to a new column  in the WakeMunicp attribute table.
sqlpath = gs.read_command("db.databases", driver="sqlite").replace('\n', '')
con = sqlite3.connect(sqlpath, isolation_level=None)
sqlstat = "ALTER TABLE WakeMunicp ADD RGBcolor varchar;"
sql.execute(sqlstat, con)
for i in xrange(len(df)):
    sql.execute("UPDATE WakeMunicp SET RGBcolor='{}' WHERE TOWN_CODE='{}';".
                format(csl[i], df['ID'].iloc[i]), con)
con.close()

# Write legend file
sqlstat = ("SELECT TOWN_CODE as ID, MB_NAME as Municipality, sum(AREA) as "
           "'Area (ha)' FROM WakeMunicp GROUP BY MB_NAME")
df = pd.read_sql_query(sqlstat, con).dropna()
lgt = "legend/area_curved"
with open("tmplegendfile", "w") as text_file:
    for i in xrange(len(df)):
        text_file.write("{}|{}|5|lf|none|{}|1|area||\n".
                        format(df.iloc[i, 1], lgt, csl[i]))
        
# Print out map
gs.run_command("g.region", vector="WakeCounty", res=100)
gs.run_command("d.mon", start="cairo", width=600, height=450,
			   output="mymap.png")
gs.run_command("d.frame", flags="c", frame="map1", at="0,100,30,100")
gs.run_command("d.vect", map="WakeCounty", fill_color="#fffef1",
			   color="53:53:53", flags="s")
gs.run_command("d.vect", map="WakeMunicp", rgb_column="RGBcolor", font="Sans",
			   label_size="10")

gs.run_command("d.frame", flags="c", frame="map2", at="5,100,0,29")
gs.run_command("d.legend.vect", flags="b", at="5,95", font="Roboto:Light",
			   title="Municipals", input="tmplegendfile", border_color="none",
			   title_fontsize=20, fontsize=16)
gs.run_command("d.mon", stop="cairo")


#------------------------------------------------------------------------------
# Compute distances
#------------------------------------------------------------------------------
# Set region
gs.run_command("g.region", flags="a", vector="WakeCounty", res="50")

# Convert school point layer to raster layer and compute distances to schools
gs.run_command("v.to.rast", input="schools_wake", output="schools", use="attr",
               attribute_column="cat")
gs.run_command("r.grow.distance", flags="m", input="schools",
               distance="dist2schools")

# Convert the WakeMunicp map to raster layer.
gs.run_command("v.to.rast", input="WakeMunicp", output="WakeMunicp",
               use="attr", attribute_column="OBJECTID", type="area",
               label_column="MB_NAME")

#------------------------------------------------------------------------------
# Create boxplots (Figure 11 in tutorial)
#------------------------------------------------------------------------------

# Get the data in a dataframe
gs.run_command("r.stats", input=["WakeMunicp", "dist2schools"],
               output="tmp1.csv", separator="comma", flags="al1gn")
df3 = pd.read_csv("tmp1.csv", sep=",", header=None, usecols=[0,1,2,3,4])
df3.columns = ["x", "y", "ID", "municipals", "distance2school"]
dfw = pd.DataFrame({col:vals['distance2school'] for col, vals in df3.groupby('municipals')})
meds = dfw.median().sort_values()
dfw = dfw[meds.index]

# Create the figure
fig, ax = plt.subplots(figsize=(6.5, 6))
bp = dfw[meds.index].boxplot(rot=90, patch_artist=True, notch=True,
        return_type='dict')

# Set labels and ticks properties
ax.set_ylabel("distance to schools (km)", fontsize=12)
ax.set_xlabel("Municipals", fontsize=12)
ax.tick_params(axis='y', labelsize=11)
ax.tick_params(axis='x', labelsize=11)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
        
# Create series with colors ordered according to median distance of municipals
from palettable.colorbrewer.qualitative import Paired_12
cls = Paired_12.colors 
municp = pd.Series(df3['municipals'].unique()).sort_values()
coltab = pd.DataFrame(
		 data=[[float(x)/255 for x in item] for item in cls],
		 index=municp)        
coltab = coltab.reindex(meds.index)
colseries = coltab.apply(lambda x: (x[0], x[1], x[2]), axis=1)
      
# Set properties of boxplot elements
plt.setp(bp['medians'], color=(0.2,0.2,0.2), linewidth=2, linestyle='-')
plt.setp(bp['fliers'], color=(0.8,0.8,0.8), marker='o', markersize=1)
plt.setp(bp['whiskers'], color=(0.6,0.6,0.6), linestyle='--', linewidth=2)
plt.setp(bp['caps'], color=(0.6,0.6,0.6))
plt.setp(bp['boxes'], color=(0.3,0.3,0.3))
for i, patch in enumerate(bp['boxes']):
	patch.set_facecolor(colseries[i])
plt.tight_layout()
plt.savefig(filename='d2s_municp_ordered.png')
plt.close('all') 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
