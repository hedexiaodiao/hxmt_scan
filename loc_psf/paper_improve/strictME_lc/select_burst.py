# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:55:19 2020

@author: admin
"""
import numpy as np
from pyecharts import options as opts
from pyecharts.charts import Kline
from pyecharts.charts import Scatter
from pyecharts.charts import Line
from pyecharts.charts import Grid
from astropy.io import fits
import re
import os,sys
from astropy.io import fits as pf

ObsID = sys.argv[1]

def plot_lc(ndata,groupname,picname,gtix,gtiy):
    ndata = np.array(ndata)
    x_data = ndata[0,:]
    y_blind = ndata[1,:]
    y_big = ndata[2,:]
    y_small = ndata[3,:]
    lc = (
        Scatter()#init_opts=opts.InitOpts(width="1600px", height="400px"))
        #Line(init_opts=opts.InitOpts(width="1680px", height="800px"))
        .add_xaxis(xaxis_data=x_data)
        .add_yaxis(
            series_name="box0",
            y_axis=y_blind,
            #areastyle_opts=opts.AreaStyleOpts(opacity=0.5),
            #linestyle_opts=opts.LineStyleOpts(),
            label_opts=opts.LabelOpts(is_show=False),
        )
        .add_yaxis(
            series_name="box1",
            y_axis=y_big,
            #yaxis_index=1,
            #areastyle_opts=opts.AreaStyleOpts(opacity=0.5),
            #linestyle_opts=opts.LineStyleOpts(),
            label_opts=opts.LabelOpts(is_show=False),
        )
        .add_yaxis(
            series_name="box2",
            y_axis=y_small,
            #yaxis_index=1,
            #areastyle_opts=opts.AreaStyleOpts(opacity=0.5),
            #linestyle_opts=opts.LineStyleOpts(),
            label_opts=opts.LabelOpts(is_show=False),
        )
        .set_global_opts(
            title_opts=opts.TitleOpts(
                title="光变曲线图",
                subtitle=groupname,
                pos_left="center",
                pos_top="top",
            ),
            tooltip_opts=opts.TooltipOpts(trigger="axis", axis_pointer_type="cross"),
            legend_opts=opts.LegendOpts(pos_left="left"),
            datazoom_opts=[
                opts.DataZoomOpts(range_start=0, range_end=100),
                opts.DataZoomOpts(type_="inside", range_start=0, range_end=100),
            ],
            xaxis_opts=opts.AxisOpts(type_="category", boundary_gap=False),
            yaxis_opts=opts.AxisOpts(name="Counts Rate", type_="value"),
        )
        .set_series_opts(
            markarea_opts=opts.MarkAreaOpts(
                is_silent=False,

            ),
            axisline_opts=opts.AxisLineOpts(),
        )

        #.render(picname+'.html')
    )
    '''
    gti = (
        Scatter()
        .add_xaxis(
        xaxis_data=gtix,
    )
        .add_yaxis(
        series_name="gti",
        y_axis=gtiy,
        xaxis_index=0,
        yaxis_index=0,
    )

    )
    '''
    overlap_gti_lc = lc
    grid_chart = Grid(init_opts=opts.InitOpts(width="1600px", height="700px"))
    grid_chart.add(
        overlap_gti_lc,
        grid_opts=opts.GridOpts(pos_left="3%", pos_right="1%", height="80%"),
                   )
    grid_chart.render(picname+'.html')

dir_head = '/sharefs/hbkg/user/luoqi/psfl/genlc'
Erange_list = ['1_6','2_6','7_12','7_20']
instr_list = ['le','le','me','me']
detname_list = [['0-30','32-62','64-94'],['0-30','32-62','64-94'],['0-17','18-35','36-53'],['0-17','18-35','36-53']]

for i in range(len(Erange_list)):
    Erange = Erange_list[i]
    fitsname = dir_head + '/' + Erange + '/lc/%s%s_g%s_%s.lc' % (ObsID, instr_list[i], str(int(0)), detname_list[i][0])
    if os.path.exists(fitsname) == 0:
        print("No input files:", fitsname)
        sys.exit(0)
    else:
        hdul = pf.open(fitsname)
        time_start = hdul[1].data['TIME'][0]
        x_data = hdul[1].data['TIME'] - time_start
        hdul.close()
        with open('./time_start.txt', 'a+') as f:
            print(time_start, file=f)
        del hdul
    ndata = np.zeros((4,x_data.shape[0]))
    ndata[0,:] = x_data
    for j in range(3):
        fitsname = dir_head+'/'+ Erange + '/lc/%s%s_g%s_%s.lc'%(ObsID,instr_list[i],str(int(j)),detname_list[i][j])
        if os.path.exists(fitsname)==0:
            print("No input files:",fitsname)
            sys.exit(0)
        else:
            hdul = pf.open(fitsname)
            ndata[j+1,:] = hdul[1].data['COUNTS']
            hdul.close()
            del hdul
    lcname = '%s%s_%s'%(ObsID,instr_list[i],Erange_list[i])
    picoutputname = dir_head+'/lc_html/%s%s_%s'%(ObsID,instr_list[i],Erange)
    plot_lc(ndata,lcname,picoutputname,[],[])

