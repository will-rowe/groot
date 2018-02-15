# Get this figure: fig = py.get_figure("https://plot.ly/~wpmr/14/")
# Get this figure's data: data = py.get_figure("https://plot.ly/~wpmr/14/").get_data()
# Add data to this figure: py.plot(Data([Scatter(x=[1, 2], y=[2, 3])]), filename ="groot-performance", fileopt="extend")
# Get y data of first trace: y1 = py.get_figure("https://plot.ly/~wpmr/14/").get_data()[0]["y"]

# Get figure documentation: https://plot.ly/python/get-requests/
# Add data documentation: https://plot.ly/python/file-options/

# If you're using unicode in your file, you may need to specify the encoding.
# You can reproduce this figure in Python with the following code!

# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

import plotly.plotly as py
from plotly.graph_objs import *
py.sign_in('username', 'api_key')
trace1 = {
  "x": ["100", "1000", "10000", "100000", "1000000", "10000000"], 
  "y": ["5.77", "5.9", "12.97", "75.18", "829.71", "14504.82"], 
  "mode": "lines", 
  "name": "1 CPU", 
  "type": "scatter", 
  "uid": "6a819c", 
  "xsrc": "wpmr:13:e63bb5", 
  "ysrc": "wpmr:13:17faf6"
}
trace2 = {
  "x": ["100", "1000", "10000", "100000", "1000000", "10000000"], 
  "y": ["5.04", "5.75", "8.06", "41.53", "386.12", "4539.8"], 
  "mode": "lines", 
  "name": "4 CPU", 
  "type": "scatter", 
  "uid": "cd99e0", 
  "xsrc": "wpmr:13:e63bb5", 
  "ysrc": "wpmr:13:04d7ae"
}
trace3 = {
  "x": ["100", "1000", "10000", "100000", "1000000", "10000000"], 
  "y": ["7.21", "6.21", "8.74", "30.88", "299.53", "2769"], 
  "mode": "lines", 
  "name": "8 CPU", 
  "type": "scatter", 
  "uid": "77fdd6", 
  "xsrc": "wpmr:13:e63bb5", 
  "ysrc": "wpmr:13:bf5a84"
}
data = Data([trace1, trace2, trace3])
layout = {
  "autosize": True, 
  "hovermode": "closest", 
  "legend": {
    "x": 0, 
    "y": -0.1, 
    "font": {"family": "Roboto"}, 
    "orientation": "h"
  }, 
  "showlegend": True, 
  "title": "Click to enter Plot title", 
  "xaxis": {
    "anchor": "y", 
    "autorange": True, 
    "dtick": "1", 
    "exponentformat": "SI", 
    "range": [2, 7], 
    "showline": True, 
    "side": "bottom", 
    "tick0": 0, 
    "tickfont": {"family": "Roboto"}, 
    "tickmode": "linear", 
    "ticks": "outside", 
    "title": "reads processed <sub>log(number of reads)</sub>", 
    "titlefont": {
      "family": "Roboto", 
      "size": 16
    }, 
    "type": "log", 
    "zeroline": True
  }, 
  "yaxis": {
    "autorange": True, 
    "dtick": "1", 
    "exponentformat": "none", 
    "linewidth": 1, 
    "range": [0.51025932493, 4.35368355524], 
    "showline": True, 
    "tickfont": {"family": "Roboto"}, 
    "tickmode": "linear", 
    "ticks": "outside", 
    "title": "runtime <sub>log(wall clock seconds)</sub>", 
    "titlefont": {
      "family": "Roboto", 
      "size": 16
    }, 
    "type": "log", 
    "zeroline": True
  }
}
fig = Figure(data=data, layout=layout)
plot_url = py.plot(fig)