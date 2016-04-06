from bokeh.models import HoverTool
from bokeh.plotting import ColumnDataSource, figure, show, output_file
import bokeh
import argparse
import sys
import os
import pandas
from math import pi
import collections
import numpy as np


def run(clip_deseq_fname, peaks_fname, lib):
    clip_deseq = pandas.read_csv(clip_deseq_fname, sep='\t')
    gl_deseq_df = pandas.read_csv(lib['gl_expression'], sep='\t')
    gl_deseq_l = gl_deseq_df.to_dict('records')
    gl_deseq = {}
    for row in gl_deseq_l:
        gl_deseq[row['WormBase ID (WS240)']] = row
    wbid_to_name = dict(zip(gl_deseq_df['WormBase ID (WS240)'].tolist(),
                            gl_deseq_df['Gene name']))
    oogenic, spermatogenic = read_programs(lib)
    clip_deseq['has_ortiz'] = [id_to_gl_deseq_val(x, gl_deseq)!=0 for x in\
                        clip_deseq['gene_id']]
    clip_deseq = clip_deseq[clip_deseq['has_ortiz']]
    cs = [id_to_color(x, oogenic, spermatogenic)\
          for x in clip_deseq['gene_id'].tolist()]
    from bokeh.sampledata.unemployment1948 import data
    data['Year'] = [str(x) for x in data['Year']]
    years = list(data['Year'])
    months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    data = data.set_index('Year')
    # This is the colormap from the original NYTimes plot
    colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce",
              "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
    month = []
    x_vals = clip_deseq['log2FoldChange'].tolist()
    y_vals = y_vals = [-np.log10(x) for x in clip_deseq['padj'].tolist()]
    labels= [str((x, wbid_to_name[x])) for x in clip_deseq['gene_id'].tolist()]
    source = ColumnDataSource(
        data=dict(x_vals=x_vals, y_vals=y_vals,
                  labels=labels,
                  color=cs)#month=month, year=year, color=color, rate=rate, label=label)
    )
    TOOLS = "resize,crosshair,hover,save,pan,box_zoom,wheel_zoom"
    p = figure(title="Fig 1A. FBF-RNA interactions are influenced by cell fate.",
               #x_range=years,
               y_range=[0, 20],
               title_text_font_size='12pt',
               #x_axis_location="above", plot_width=900, plot_height=400,
               #toolbar_location="left",
               tools=TOOLS)
    p.scatter(x_vals, y_vals, source=source, color="color",
            )#, 1, 1, #source=source,
    p.select_one(HoverTool).tooltips = [
        ('label', '@labels'),
        ('x_vals', '@x_vals'),
    ]
    p.xaxis.axis_label = ' FBF spermatogenic/oogenic (log2)'
    p.yaxis.axis_label = '-log10(p value) for differential binding'
    p.xaxis.axis_label_text_font_size = '12pt'
    p.yaxis.axis_label_text_font_size = '12pt'
    output_file(
        'Fig 1A. FBF-RNA interactions are influenced by cell fate.html',
        title="Fig 1A. FBF-RNA interactions are influenced by cell fate.html")
    bokeh.plotting.save(p)      # show the plot
    ###################
    # Peaks.
    peaks_df = pandas.read_csv(peaks_fname, sep='\t')
    peaks_l = peaks_df.to_dict('records')
    peaks_d = collections.defaultdict(list)
    for row in peaks_l: peaks_d[row['gene_id']].append(row)
    size = [id_to_size(x, peaks_d) for x in clip_deseq['gene_id'].tolist()]
    # Colors.
    cs = [id_to_color(x, oogenic, spermatogenic)\
          for x in clip_deseq['gene_id'].tolist()]
    for row in peaks_l: peaks_d[row['gene_id']].append(row)
    size = [id_to_size(x, peaks_d) for x in clip_deseq['gene_id'].tolist()]
    gl_deseq_vals = [id_to_gl_deseq_val(x, gl_deseq) \
          for x in clip_deseq['gene_id'].tolist()]
    # CUT:
#    plt.scatter(
#        clip_deseq['log2FoldChange'], gl_deseq_vals, linewidth=0,
#        color=cs, s=size, alpha=0.2)
#    plt.ylim([-10, 15])
#    plt.xlim([-4,8])
#    plt.axhline(y=0, linestyle='--', c='k', alpha=0.5)
#    plt.xlabel('FBF spermatogenic/oogenic (log2)')
#    plt.ylabel('RNA-seq spermatogenic/oogenic (log2)')
#    plt.savefig('volcano_sized_colored_by_magnitude.pdf', format='pdf')
#    plt.clf()
#    plt.close()
    # end cut
    x_vals = clip_deseq['log2FoldChange'].tolist()
    y_vals = [id_to_gl_deseq_val(x, gl_deseq) \
          for x in clip_deseq['gene_id'].tolist()]
    source = ColumnDataSource(
        data=dict(x_vals=x_vals, y_vals=y_vals,
                  labels=labels,
                  color=cs)#month=month, year=year, color=color, rate=rate, label=label)
    )
    TOOLS = "resize,crosshair,hover,save,pan,box_zoom,wheel_zoom"
    p = figure(title="Fig 2A. FBF-RNA interactions are influenced by cell fate.",
               y_range=[-10, 15],
               x_range=[-4, 8],
               title_text_font_size='12pt',
               tools=TOOLS)
    p.scatter(x_vals, y_vals, source=source, color="color",
              size=size)
    #, 1, 1, #source=source,
    p.select_one(HoverTool).tooltips = [
        ('label', '@labels'),
        ('x_vals', '@x_vals'),
    ]
    p.xaxis.axis_label = 'FBF spermatogenic/oogenic (log2)'
    p.yaxis.axis_label = 'RNA-seq spermatogenic/oogenic (log2)'
    p.xaxis.axis_label_text_font_size = '12pt'
    p.yaxis.axis_label_text_font_size = '12pt'

    output_file(
        'Fig 2A. FBF-RNA interactions are influenced by cell fate.html',
        title="Fig 2A. FBF-RNA interactions are influenced by cell fate.html")
    bokeh.plotting.save(p)      # show the plot

def id_to_size(wbid, peaks_d):
    #print wbid
    if wbid not in peaks_d:
        return 2
    return 10
    #print "Found size"
    #print min([1e3, max([x['exp'] for x in peaks_d[wbid]])/1.])
#    return min([1e2, max([x['exp'] for x in peaks_d[wbid]])/1.])

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--peaks', help="""Input peaks file (Required).""")
    parser.add_argument('-c', '--config', help="""Folder of config.py.""")
    parser.add_argument(
        '-v', '--no_ui',
        help="""Do not wait for user input (default: false).""",
        action='store_true', default=False)
    args = parser.parse_args()
    return args


def id_to_color(wbid, oogenic, spermatogenic):
    if (wbid in oogenic) and (wbid not in spermatogenic):
        
        return "#933b41"# NYT color
    elif (wbid not in oogenic) and (wbid in spermatogenic):
        return "#6600ff"
        #return "#75968f"
    elif (wbid not in oogenic) and (wbid not in spermatogenic):
        return '#75968f' # NYT color
    else:
        return '#75968f'


def read_programs(lib):
    oogenic = {}
    spermatogenic = {}
    with open(lib['spermatogenic_program'], 'r') as f:
        for li in f:
            s = li.rstrip('\t').split('\t')
            spermatogenic[s[0]] = s[5]
    with open(lib['oogenic_program'], 'r') as f:
        for li in f:
            s = li.rstrip('\t').split('\t')
            oogenic[s[0]] = s[5]
    return oogenic, spermatogenic


def id_to_gl_deseq_val(wbid, gl_deseq):
    if wbid not in gl_deseq:
        return 0
    val = gl_deseq[wbid]['log2 fold change (of normalized reads)']
    try:
        val=float(val)
        return val
    except:
        return 0

if __name__ == '__main__':
    args = get_args()
#    sys.path.insert(0, args.config)
#    import config
#    lib = config.config()
#    del sys.path[0]
    if not os.path.exists('figs/'): os.system('mkdir figs/')
#    rename('counts/')
    fog3_lib = '/groups/Kimble/Common/fog_iCLIP/calls/lib/'
    lib = {
        'oogenic_program': fog3_lib + 'TableS1_oogenic.txt',
        'spermatogenic_program': fog3_lib + 'TableS2_spermatogenic.txt',
        'gl_expression': fog3_lib + '/ortiz/DESeq_genes_in_gonad.txt',
    }
    run(
        '/groups/Kimble/Common/fbf_celltype/oogenic/peak_stats.txt',
        args.peaks,
        lib)
