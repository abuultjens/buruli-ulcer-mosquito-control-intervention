import os
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import contextily as ctx
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
import numpy as np

treatment_circle_colour = 'red'
control_circle_colour = 'blue'
circle_alpha = 0.5

buffer_fill_color = 'grey'
buffer_alpha = 0.5

heatmap_low_colour = 'yellow'
middle_heatmap_colour = 'orange'
heatmap_high_colour = 'darkred'

custom_cmap = LinearSegmentedColormap.from_list("custom_heatmap", [heatmap_low_colour, middle_heatmap_colour, heatmap_high_colour])

# basemap
basemap_source = ctx.providers.OpenStreetMap.Mapnik

# metric CRS for buffering
metric_crs = "EPSG:28355"

# scale bar
scale_bar_length = 2000 
scale_bar_thickness_factor = 0.1

# legend placement
#legend_bbox = (0.72, 1.015)
legend_bbox = (0.045, 0.2)
legend_loc = 'upper left'
legend_fontsize = 13

# scale bar placement
scale_bar_bbox = (0.905, 0.0105)

# north arrow placement
north_arrow_tip_coords = (0.95, 0.15)
north_arrow_text_coords = (0.95, 0.03) 
north_arrow_fontsize = 16
north_arrow_arrow_width = 3
north_arrow_headwidth = 12

# heading font size
heading_fontsize = 20

# axis label and ticks
x_axis_label_fontsize = 14
y_axis_label_fontsize = 14
x_tick_fontsize = 12
y_tick_fontsize = 12 
num_xticks = 6 
num_yticks = 6 

# meshblock colours for cases
meshblock_colour_1_case = 'yellow'
meshblock_colour_2_cases = 'orange'
meshblock_colour_3_cases = 'darkred'

# paths to polygon KML files to overlay
polygon_kml_files = [
     'Treatment_1.kml',
     'Treatment_2.kml',
     'Treatment_3.kml',
     'Treatment_4.kml',
     'Treatment_5.kml',
     'Treatment_6.kml',
     'Control_1.kml',
     'Control_2.kml',
     'Control_3.kml',
     'Control_4.kml',
     'Control_5.kml',
     'Control_6.kml',
]
polygon_edgecolor = 'black'
polygon_facecolor = 'none'
polygon_linewidth = 1
polygon_zorder = 3

# polygon labels
polygon_label_fontsize = 12
polygon_label_color = 'grey'
polygon_label_offset_factor = 0.02 

# file paths
meshblock_shp_path = '/Users/abuultjens/Google Drive/OneDrive/PhD/Bioinformatics/Mozzie_surveillance/QGIS/1270055001_mb_2011_vic_shape/MB_2011_VIC.shp'
cases_csv_path = '66_Inner_northwest_2024_cases_symptom_MEDIAN.csv'

# read the meshblock shapefile
meshblocks = gpd.read_file(meshblock_shp_path)

# read csv of cases and create a GeoDataFrame
cases = pd.read_csv(cases_csv_path)
geometry = [Point(xy) for xy in zip(cases['lon'], cases['lat'])]
cases_gdf = gpd.GeoDataFrame(cases, geometry=geometry, crs="EPSG:4326")

# make both GeoDataFrames use the same CRS
if meshblocks.crs != cases_gdf.crs:
    cases_gdf = cases_gdf.to_crs(meshblocks.crs)

# determine which meshblock each case falls within.
cases_with_mesh = gpd.sjoin(cases_gdf, meshblocks, how='left', predicate='within')

# count the number of cases per meshblock
if 'mesh_id' in meshblocks.columns:
    counts = cases_with_mesh.groupby('mesh_id').size().reset_index(name='case_count')
    meshblocks = meshblocks.merge(counts, on='mesh_id', how='left')
else:
    counts = cases_with_mesh.groupby('index_right').size().reset_index(name='case_count')
    meshblocks = meshblocks.merge(counts, left_index=True, right_on='index_right', how='left')

meshblocks['case_count'] = meshblocks['case_count'].fillna(0).astype(int)

# filter only meshblocks that have at least one case
meshblocks_cases = meshblocks[meshblocks['case_count'] > 0]

# define the mapping area (bounding box) and clip the data accordingly.
min_lon, max_lon = 144.86, 144.986887
min_lat, max_lat = -37.785, -37.714593
meshblocks_cases_clipped = meshblocks_cases.cx[min_lon:max_lon, min_lat:max_lat]

# Classify meshblocks into discrete case categories
def classify_cases(count):
    if count == 1:
        return '1 case'
    elif count == 2:
        return '2 cases'
    else:
        return '3+ cases'

meshblocks_cases_clipped['case_category'] = meshblocks_cases_clipped['case_count'].apply(classify_cases)

# define discrete colours for each category
category_colors = {
    '1 case': meshblock_colour_1_case,
    '2 cases': meshblock_colour_2_cases,
    '3+ cases': meshblock_colour_3_cases
}

# define control and treatment site coordinates
# Control sites
control_sites_data = [
    {"site": "C1", "lat": -37.761268, "lon": 144.892909},
    {"site": "C2", "lat": -37.753892, "lon": 144.906185},
    {"site": "C3", "lat": -37.741625, "lon": 144.959313},
    {"site": "C4", "lat": -37.738252, "lon": 144.908115},
    {"site": "C5", "lat": -37.754925, "lon": 144.928309},
    {"site": "C6", "lat": -37.725209, "lon": 144.882608}
]
# Treatment sites
treatment_sites_data = [
    {"site": "T1", "lat": -37.76749, "lon": 144.915136},
    {"site": "T2", "lat": -37.724524, "lon": 144.938588},
    {"site": "T3", "lat": -37.73955, "lon": 144.939485},
    {"site": "T4", "lat": -37.760762, "lon": 144.943505},
    {"site": "T5", "lat": -37.750055, "lon": 144.891333},
    {"site": "T6", "lat": -37.773462, "lon": 144.942397}
]

# make GeoDataFrames for the sites.
control_sites_df = pd.DataFrame(control_sites_data)
treatment_sites_df = pd.DataFrame(treatment_sites_data)

control_sites_df['geometry'] = control_sites_df.apply(lambda row: Point(row['lon'], row['lat']), axis=1)
treatment_sites_df['geometry'] = treatment_sites_df.apply(lambda row: Point(row['lon'], row['lat']), axis=1)

control_sites_gdf = gpd.GeoDataFrame(control_sites_df, geometry='geometry', crs="EPSG:4326")
treatment_sites_gdf = gpd.GeoDataFrame(treatment_sites_df, geometry='geometry', crs="EPSG:4326")

# create 800m radius buffers for the sites.
control_sites_metric = control_sites_gdf.to_crs(metric_crs)
treatment_sites_metric = treatment_sites_gdf.to_crs(metric_crs)

control_sites_metric['buffer'] = control_sites_metric.geometry.buffer(650)
treatment_sites_metric['buffer'] = treatment_sites_metric.geometry.buffer(650)

control_buffers = gpd.GeoDataFrame(control_sites_metric[['site']], geometry=control_sites_metric['buffer'], crs=metric_crs).to_crs("EPSG:4326")
treatment_buffers = gpd.GeoDataFrame(treatment_sites_metric[['site']], geometry=treatment_sites_metric['buffer'], crs=metric_crs).to_crs("EPSG:4326")

# plot the map with a basemap.
fig, ax = plt.subplots(figsize=(10, 10))

# set plot limits.
ax.set_xlim(min_lon, max_lon)
ax.set_ylim(min_lat, max_lat)

# add basemap.
ctx.add_basemap(ax, source=basemap_source, crs="EPSG:4326", zorder=0)

# plot site buffers
#control_buffers.plot(ax=ax, color=buffer_fill_color, alpha=buffer_alpha,
#                     edgecolor=buffer_fill_color, zorder=1)
#treatment_buffers.plot(ax=ax, color=buffer_fill_color, alpha=buffer_alpha,
#                       edgecolor=buffer_fill_color, zorder=1)

# plot meshblocks with cases
#for cat, color in category_colors.items():
#    subset = meshblocks_cases_clipped[meshblocks_cases_clipped['case_category'] == cat]
#    subset.plot(ax=ax, color=color, edgecolor='black', linewidth=0.8, zorder=2)

# load and plot additional polygon KML files 
if polygon_kml_files:
    for kml_file in polygon_kml_files:
        try:
            polygons = gpd.read_file(kml_file, driver='KML')
            polygons = polygons.to_crs("EPSG:4326")
            # make label from filename
            filename = os.path.basename(kml_file)
            label = os.path.splitext(filename)[0].replace("_", " ")
            # get fill from filename
            if "Treatment" in label:
                poly_fill = treatment_circle_colour
            elif "Control" in label:
                poly_fill = control_circle_colour
            else:
                poly_fill = polygon_facecolor
            polygons.plot(ax=ax, edgecolor=polygon_edgecolor, facecolor=poly_fill,
                          linewidth=polygon_linewidth, zorder=polygon_zorder)
            # Get total bounds
            bounds = polygons.total_bounds
            x_center = (bounds[0] + bounds[2]) / 2
            y_top = bounds[3]
            # Compute offset
            offset = (bounds[3] - bounds[1]) * polygon_label_offset_factor
            ax.text(x_center, y_top + offset, label,
                    ha='center', va='bottom',
                    fontsize=polygon_label_fontsize,
                    color=polygon_label_color,
                    path_effects=[pe.withStroke(linewidth=1, foreground='black')])
        except Exception as e:
            print(f"Error loading {kml_file}: {e}")

# Add combined legend for the site buffers, the 800m radius zone, and the meshblock categories
legend_handles = [
    Patch(facecolor=treatment_circle_colour, edgecolor=treatment_circle_colour, label='Treatment Sites'),
    Patch(facecolor=control_circle_colour, edgecolor=control_circle_colour, label='Control Sites'),
#    Patch(facecolor=buffer_fill_color, edgecolor=buffer_fill_color, label='650m radius'),

#    Patch(facecolor=meshblock_colour_1_case, edgecolor='black', label='1 case'),
#    Patch(facecolor=meshblock_colour_2_cases, edgecolor='black', label='2 cases'),
#    Patch(facecolor=meshblock_colour_3_cases, edgecolor='black', label='3 cases')
]
ax.legend(handles=legend_handles, loc=legend_loc, bbox_to_anchor=legend_bbox,
          prop={'size': legend_fontsize})

# add scale bar
lat_avg = (min_lat + max_lat) / 2.0
meters_per_degree = 111320 * np.cos(np.radians(lat_avg))
scale_bar_length_degrees = scale_bar_length / meters_per_degree
scale_bar_thickness = scale_bar_length_degrees * scale_bar_thickness_factor

fontprops = fm.FontProperties(size=8)
scalebar = AnchoredSizeBar(ax.transData,
                           scale_bar_length_degrees,
                           f'{scale_bar_length} m',
                           loc='lower right',
                           pad=0.1,
                           color='black',
                           frameon=True,
                           size_vertical=scale_bar_thickness,
                           fontproperties=fontprops,
                           bbox_to_anchor=scale_bar_bbox,
                           bbox_transform=ax.transAxes)
ax.add_artist(scalebar)

# add north arrow
ax.annotate('N',
            xy=north_arrow_tip_coords,
            xytext=north_arrow_text_coords,
            arrowprops=dict(facecolor='black',
                            width=north_arrow_arrow_width,
                            headwidth=north_arrow_headwidth),
            ha='center', va='center',
            fontsize=north_arrow_fontsize,
            xycoords=ax.transAxes)

# set axis ticks and labels
ax.set_xticks(np.linspace(min_lon, max_lon, num_xticks))
ax.set_yticks(np.linspace(min_lat, max_lat, num_yticks))
ax.tick_params(axis='x', labelsize=x_tick_fontsize)
ax.tick_params(axis='y', labelsize=y_tick_fontsize)

# add title and axis labels
ax.set_title('Treatment and control sites',
             fontsize=heading_fontsize)
ax.set_xlabel('Longitude', fontsize=x_axis_label_fontsize)
ax.set_ylabel('Latitude', fontsize=y_axis_label_fontsize)

# save figure
output_svg = 'Fig_1_v7_650m_66-MEDIAN_NO-ZONES.svg'
output_png = 'Fig_1_v7_650m_66-MEDIAN_NO-ZONES.png'
plt.savefig(output_svg, format='svg')
plt.savefig(output_png, format='png', dpi=300)

plt.show()
