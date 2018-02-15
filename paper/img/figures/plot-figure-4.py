# Get this figure: fig = py.get_figure("https://plot.ly/~wpmr/17/")
# Get this figure's data: data = py.get_figure("https://plot.ly/~wpmr/17/").get_data()
# Add data to this figure: py.plot(Data([Scatter(x=[1, 2], y=[2, 3])]), filename ="groot-benchmark", fileopt="extend")
# Get y data of first trace: y1 = py.get_figure("https://plot.ly/~wpmr/17/").get_data()[0]["y"]

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
  "x": ["1", "5", "10", "15", "20"], 
  "y": ["146.2", "781.45", "1541.36", "2460.63", "3373.83"], 
  "line": {"color": "rgb(23, 190, 207)"}, 
  "mode": "lines", 
  "name": "groot", 
  "type": "scatter", 
  "uid": "ca8b75", 
  "xaxis": "x", 
  "xsrc": "wpmr:16:4898c9", 
  "yaxis": "y2", 
  "ysrc": "wpmr:16:29c829"
}
trace2 = {
  "x": ["1", "5", "10", "15", "20"], 
  "y": ["955.5", "5337.92", "10014.87", "14385.57", "19243.29"], 
  "line": {"color": "rgb(148, 103, 189)"}, 
  "mode": "lines", 
  "name": "args-oap", 
  "type": "scatter", 
  "uid": "a5d4aa", 
  "xaxis": "x", 
  "xsrc": "wpmr:16:4898c9", 
  "yaxis": "y2", 
  "ysrc": "wpmr:16:a8830c"
}
trace3 = {
  "x": ["1", "5", "10", "15", "20"], 
  "y": ["1", "1", "2", "2", "2"], 
  "line": {
    "color": "rgb(23, 190, 207)", 
    "dash": "longdashdot"
  }, 
  "marker": {
    "color": "rgb(23, 190, 207)", 
    "line": {
      "color": "rgb(23, 190, 207)", 
      "width": 0
    }, 
    "symbol": "cross"
  }, 
  "mode": "lines+markers", 
  "name": "groot false positives", 
  "type": "scatter", 
  "uid": "7533ed", 
  "xaxis": "x", 
  "xsrc": "wpmr:16:4898c9", 
  "yaxis": "y", 
  "ysrc": "wpmr:16:31f64a"
}
trace4 = {
  "x": ["1", "5", "10", "15", "20"], 
  "y": ["1", "0", "0", "0", "0"], 
  "line": {
    "color": "rgb(23, 190, 207)", 
    "dash": "dot", 
    "width": 2
  }, 
  "marker": {
    "color": "rgba(23, 190, 207, 0)", 
    "line": {
      "color": "rgb(23, 190, 207)", 
      "width": 2
    }, 
    "size": 8, 
    "symbol": "circle"
  }, 
  "mode": "lines+markers", 
  "name": "groot false negatives", 
  "type": "scatter", 
  "uid": "8633dc", 
  "xaxis": "x", 
  "xsrc": "wpmr:16:4898c9", 
  "yaxis": "y", 
  "ysrc": "wpmr:16:00e926"
}
trace5 = {
  "x": ["1", "5", "10", "15", "20"], 
  "y": ["3", "3", "3", "3", "3"], 
  "line": {"dash": "longdashdot"}, 
  "marker": {"symbol": "cross"}, 
  "mode": "lines+markers", 
  "name": "args-oap false positives", 
  "type": "scatter", 
  "uid": "0e3954", 
  "xaxis": "x", 
  "xsrc": "wpmr:16:4898c9", 
  "yaxis": "y", 
  "ysrc": "wpmr:16:abcb7b"
}
trace6 = {
  "x": ["1", "5", "10", "15", "20"], 
  "y": ["73", "94", "104", "111", "116"], 
  "line": {
    "color": "rgb(148, 103, 189)", 
    "dash": "dot"
  }, 
  "marker": {
    "color": "rgba(148, 103, 189, 0)", 
    "line": {
      "color": "rgb(148, 103, 189)", 
      "width": 2
    }, 
    "size": 8
  }, 
  "mode": "lines+markers", 
  "name": "args-oap false negatives", 
  "type": "scatter", 
  "uid": "de1995", 
  "xaxis": "x", 
  "xsrc": "wpmr:16:4898c9", 
  "yaxis": "y", 
  "ysrc": "wpmr:16:621f98"
}
data = Data([trace1, trace2, trace3, trace4, trace5, trace6])
layout = {
  "autosize": True, 
  "bargap": 0.73, 
  "barmode": "group", 
  "hovermode": "closest", 
  "legend": {
    "font": {"family": "Roboto"}, 
    "orientation": "v", 
    "traceorder": "normal"
  }, 
  "margin": {"l": 80}, 
  "showlegend": True, 
  "xaxis": {
    "anchor": "y", 
    "autorange": True, 
    "domain": [0, 1], 
    "mirror": False, 
    "range": [-0.205049261084, 21.2050492611], 
    "showgrid": True, 
    "showline": True, 
    "showticklabels": True, 
    "side": "bottom", 
    "tickfont": {"family": "Roboto"}, 
    "tickmode": "auto", 
    "ticks": "outside", 
    "title": "metagenome coverage<sub> (X times coverage)</sub>", 
    "titlefont": {
      "family": "Roboto", 
      "size": 16
    }, 
    "type": "linear", 
    "zeroline": True
  }, 
  "yaxis": {
    "autorange": True, 
    "domain": [0, 0.45], 
    "dtick": 10, 
    "range": [-8.05611964187, 124.056119642], 
    "showgrid": False, 
    "showline": False, 
    "showticklabels": True, 
    "tick0": 100, 
    "tickfont": {"family": "Roboto"}, 
    "tickmode": "linear", 
    "ticks": "outside", 
    "title": "errors <sub>(number of miss-annotated genes)</sub>", 
    "titlefont": {
      "family": "Roboto", 
      "size": 16
    }, 
    "type": "linear", 
    "zeroline": False
  }, 
  "yaxis2": {
    "anchor": "x", 
    "autorange": True, 
    "domain": [0.55, 1], 
    "dtick": 5000, 
    "exponentformat": "none", 
    "range": [-914.749444444, 20304.2394444], 
    "showgrid": False, 
    "showline": False, 
    "showticklabels": True, 
    "tick0": 10000, 
    "tickfont": {"family": "Roboto"}, 
    "tickmode": "linear", 
    "ticks": "outside", 
    "title": "runtime <sub>(wall clock seconds)</sub>", 
    "titlefont": {
      "family": "Roboto", 
      "size": 16
    }, 
    "type": "linear", 
    "zeroline": False
  }
}
fig = Figure(data=data, layout=layout)
plot_url = py.plot(fig)