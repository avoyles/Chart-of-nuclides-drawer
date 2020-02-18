#!/usr/bin/env python3
"""

This script produces chart of nuclides in svg format from data
provided in xml document.

    usage: ChartDrawer.py [-h] [--names] [--halflives] [--magic] [--numbers]
                        [--z Z Z] [--n N N]
                        datafile outfile

    positional arguments:
    datafile     Input data base XML file (required)
    outfile      Output data SVG file (required)

    optional arguments:
    -h, --help   show this help message and exit
    --names      Disable names
    --halflives  Disable half-lives
    --magic      Disable magic numbers
    --numbers    Disable numbers along axis
    --unknown    Disable isotopes with unknown decay mode
    --isomers    Disable detialed isomer data
    --z Z Z      Atomic number Z range (int), default: [0, 120]
    --n N N      Neutron number N range (int), default: [0, 180]
    --p '...'    List of product isotopes to draw
    --t '...'    List of target isotopes to draw (mass 0 draws all stable isotopes)

   The Nubase2xml.py script will generate the input xml file from Nubase format

   The NuBase ascii file can be downloaded from:
   http://amdc.in2p3.fr/nubase/nubtab03.asc

   Please note that the database is old (2003) so many values are outdated.
   As soon as expected NuBase2013 is released the Nubase2xml script will 
   be updated to parse new data format.

   The generated chart of nuclides is following (but not precisely) the format
   of Karlsruher Chart of Nuclides. The nuclides are coded with colors according to
   primary decay mode:
   Black - stable
   Yellow - Alpha
   Red - B+/EC
   Blue - B-
   Green - fission
   Orange - proton / two-proton emission
   Violet - cluster emission
   Light blue - neutron / two-neutron emission

   Secondary decay mode is indicated by triangle. A large triangle is used if
   secondary decay mode branching is larger then 5%. Otherwise small triangle is
   used (likewise for tertiary decay mode).

   A Nuclide class provides parser to read data from NuBase ascii file or xml
   document. Note that some information is not being used on the chart 
   (e.g. mass of nuclides, isomeric states), but is present in the data base.
   The Nuclide class can easyli used to write scripts using these informations 
   as well.
"""

import sys
import argparse
import re
import xml.dom.minidom
from Nuclide import *
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
import cairosvg
import subprocess
import numpy as np


# Definition of colors used for decay modes
COLORS = { 'is': '#000000',
           'b-': '#62aeff',
           '2b-': '#62aeff',
           'b+': '#ff7e75',
           'ec': '#ff7e75',
           '2ec': '#ff7e75',
           'a' : '#fffe49',
           'sf': '#5cbc57',
           'p' : '#ffa425',
           '2p': '#ffa425', 
           'n' : '#9fd7ff',
           '2n': '#9fd7ff',
           'it': '#ffffff',
           'cluster': '#a564cc',
           '?': '#cccccc' }

FONT_COLOR_DARK = '#000000'
FONT_COLOR_BRIGHT = '#aaaaaa'

# Size of rectangle in pixels
SIZE_SHAPE = 30
# Size of margin between rectangles in pixels
SIZE_GAP = 2 
# Total size of one nuclid with margin
SIZE_FIELD = SIZE_SHAPE + SIZE_GAP
# Font size used for element name
SIZE_FONT = 7 
# Font size used for half-life
SIZE_FONT_HL = 5 

MAGIC_NUMBERS = [2, 8, 20, 28, 50, 82, 126]

def load_xml_nuclear_table(datafile, n_range, z_range,
                           n_limits = [None, None], z_limits = [None, None]):
    """Loads data from nuclear table in xml format. Returns list of
    Nuclide objects
    """
    # Make high and low limit oposite
    # Later each point is checked against:
    # n_limits[0] = N if N < n_limits[0]
    # n_limits[1] = N if N > n_limits[1]
    # (Z likewise)
    # So oposite limit here forces first point to set 
    # reasonable limits without loosing any data point
    n_limits[0] = n_range[1]
    n_limits[1] = n_range[0]

    z_limits[0] = z_range[1]
    z_limits[1] = z_range[0]

    try:
        dom = xml.dom.minidom.parse(datafile)
    except (EnvironmentError, xml.parsers.expat.ExpatError) as err:
        print("{0}: import error: {1}".format(datafile, err))
        return None

    data = []
    for nuclide in dom.getElementsByTagName("nuclide"):
        try:
            A = int(nuclide.getAttribute('A'))
            Z = int(nuclide.getAttribute('Z'))
            N = A - Z

            if not(n_range[0] <= N <= n_range[1] and 
                z_range[0] <= Z <= z_range[1]):
                continue
            elif N > n_range[1] and Z > z_range[1]:
                break

            if N < n_limits[0]:
                n_limits[0] = N
            if N > n_limits[1]:
                n_limits[1] = N
            if Z < z_limits[0]:
                z_limits[0] = Z
            if Z > z_limits[1]:
                z_limits[1] = Z

            isotope = NuclideXml(Z, A, nuclide)
            data.append(isotope)
        except (ValueError, LookupError) as err:
            print("{0}: import error: {1}".format(datafile, err))
            return False
    return data

def _draw_isomer_rectangle(layer, position, color, name):
    """Draws rectangle (basic nuclide on map) position is
    given for left top corner """
    x = position[0]
    y = position[1]
    rectangle = svg.createElement("rect")
    rectangle.setAttribute("id", '{}'.format(name))
    rectangle.setAttribute("width", str((SIZE_SHAPE/2.0)-0.5))
    rectangle.setAttribute("height", str((SIZE_SHAPE*0.7)-0.5))
    # rectangle.setAttribute("stroke", "#000000")
    # rectangle.setAttribute("stroke-width", "0.5")
    rectangle.setAttribute("fill", color)
    rectangle.setAttribute("x", str(position[0]+0.3))
    rectangle.setAttribute("y", str(position[1] + (SIZE_SHAPE*0.3) + 0.2))
    layer.appendChild(rectangle)


def _draw_rectangle(layer, position, color, name):
    """Draws rectangle (basic nuclide on map) position is
    given for left top corner """
    x = position[0]
    y = position[1]
    rectangle = svg.createElement("rect")
    rectangle.setAttribute("id", '{}'.format(name))
    rectangle.setAttribute("width", str(SIZE_SHAPE))
    rectangle.setAttribute("height", str(SIZE_SHAPE))
    rectangle.setAttribute("stroke", "#000000")
    rectangle.setAttribute("stroke-width", "0.5")
    rectangle.setAttribute("fill", color)
    rectangle.setAttribute("x", str(position[0]))
    rectangle.setAttribute("y", str(position[1]))
    layer.appendChild(rectangle)

def _draw_triangle(layer, position, color, name, corner = 'rb'):
    """Draws triangle (half-rectangle), position is given
    for left top corner of rectangle, triangle is drawn in
    right bottom corner"""
    x = position[0]
    y = position[1]
    triangle = svg.createElement("polygon")
    triangle.setAttribute("id", '{}'.format(name))
    triangle.setAttribute("stroke", "#000000")
    triangle.setAttribute("stroke-width", "0.0")
    triangle.setAttribute("stroke-linejoin", "bevel")
    triangle.setAttribute("fill", color)
    triangle.setAttribute("x", str(position[0]))
    triangle.setAttribute("y", str(position[1]))
    if corner == 'lt':
        x1 = x + 0.25
        y1 = y + SIZE_SHAPE - 0.25
        x2 = x1 
        y2 = y + 0.25
        x3 = x + SIZE_SHAPE - 0.25
        y3 = y2 
    else:
        #default right bottom corner
        x1 = x + 0.25
        y1 = y + SIZE_SHAPE - 0.25
        x2 = x + SIZE_SHAPE - 0.25
        y2 = y + 0.25
        x3 = x2 
        y3 = y1 
    triangle.setAttribute( "points", 
                           "{},{} {},{} {},{}".format(x1, y1, x2, y2, x3, y3))
    layer.appendChild(triangle)

def _draw_small_triangle(layer, position, color, name, corner = 'rb'):
    """Draws small triangle in the corner of rectangle,
    position is left top corner of rectangle"""
    x = position[0]
    y = position[1]
    small_triangle = svg.createElement("polygon")
    small_triangle.setAttribute("id", '{}'.format(name))
    small_triangle.setAttribute("stroke", "#000000")
    small_triangle.setAttribute("stroke-width", "0.0")
    small_triangle.setAttribute("stroke-linejoin", "bevel")
    small_triangle.setAttribute("fill", color)
    small_triangle.setAttribute("x", str(position[0]))
    small_triangle.setAttribute("y", str(position[1]))
    if corner == 'lt':
        x1 = x + 0.25
        y1 = y +  1 / 3 * SIZE_SHAPE
        x2 = x1 
        y2 = y + 0.25
        x3 = x + 1 / 3 * SIZE_SHAPE
        y3 = y2 
    elif corner == 'rt':
        x1 = x + 2 / 3 * SIZE_SHAPE
        y1 = y + 0.25
        x2 = x + SIZE_SHAPE - 0.25
        y2 = y1
        x3 = x2
        y3 = y + 1 / 3 * SIZE_SHAPE 
    else:
        # default right bottom case
        x1 = x + SIZE_SHAPE * 2 / 3
        y1 = y + SIZE_SHAPE - 0.25
        x2 = x + SIZE_SHAPE - 0.25
        y2 = y1
        x3 = x2
        y3 = y + SIZE_SHAPE * 2 / 3
    small_triangle.setAttribute(
            "points",
            "{},{} {},{} {},{}".format(x1, y1, x2, y2, x3, y3))
    layer.appendChild(small_triangle)

def _draw_small_isomer_triangle(layer, position, color, name, corner = 'rb'):
    """Draws small triangle in the corner of rectangle,
    position is left top corner of rectangle"""
    x = position[0]- (SIZE_SHAPE/2.0) - 0.05
    y = position[1] - 0.05

    # print('x,y: ',x,y)
    small_triangle = svg.createElement("polygon")
    small_triangle.setAttribute("id", '{}'.format(name))
    small_triangle.setAttribute("stroke", "#000000")
    small_triangle.setAttribute("stroke-width", "0.0")
    small_triangle.setAttribute("stroke-linejoin", "bevel")
    small_triangle.setAttribute("fill", color)
    small_triangle.setAttribute("x", str(position[0] - (SIZE_SHAPE/2.0)))
    small_triangle.setAttribute("y", str(position[1]))
    if corner == 'lt':
        x1 = x + 0.25
        y1 = y +  1 / 3 * SIZE_SHAPE
        x2 = x1 
        y2 = y + 0.25
        x3 = x + 1 / 3 * SIZE_SHAPE
        y3 = y2 
    elif corner == 'rt':
        x1 = x + 2 / 3 * SIZE_SHAPE
        y1 = y + 0.25
        x2 = x + SIZE_SHAPE - 0.25
        y2 = y1
        x3 = x2
        y3 = y + 1 / 3 * SIZE_SHAPE 
    else:
        # default right bottom case
        x1 = x + SIZE_SHAPE * 2 / 3
        y1 = y + SIZE_SHAPE - 0.25
        x2 = x + SIZE_SHAPE - 0.25
        y2 = y1
        x3 = x2
        y3 = y + SIZE_SHAPE * 2 / 3
    small_triangle.setAttribute(
            "points",
            "{},{} {},{} {},{}".format(x1, y1, x2, y2, x3, y3))
    layer.appendChild(small_triangle)

    # rectangle.setAttribute("width", str(SIZE_SHAPE/2.0))
    # rectangle.setAttribute("height", str(SIZE_SHAPE*0.7))
    # # rectangle.setAttribute("stroke", "#000000")
    # # rectangle.setAttribute("stroke-width", "0.5")
    # rectangle.setAttribute("fill", color)
    # rectangle.setAttribute("x", str(position[0]+0.5))
    # rectangle.setAttribute("y", str(position[1] + (SIZE_SHAPE*0.3) - 0.5))

    # loc_1 = [position[0]+(SIZE_SHAPE/2.0)-0.5,  position[1] + (SIZE_SHAPE*0.3) +0.5]
    # print('x_loc',loc_1)
    # loc_2 = [position[0]+(SIZE_SHAPE/2.0)-0.5, position[1] + (SIZE_SHAPE) -0.5 ]
    # # loc_2 =[42.5, 34.0]
    # print('y_loc',loc_2)

    # _draw_isomer_line(layers[0], loc_1, loc_2, "jbjjb")

def _draw_text(layer, position, font_color, font_size, text):
    """Draws text"""
    x = position[0]
    y = position[1]
    text_node = svg.createTextNode(text)

    text_el = svg.createElement("text")
    text_el.appendChild(text_node)
    text_el.setAttribute("text-anchor", "middle")
    text_el.setAttribute("font-family", "sans")
    text_el.setAttribute(
                "style", 
                "font-size:{}px; fill:{}".format(font_size, font_color))
    text_el.setAttribute("x", '{0:.2f}'.format(x))
    text_el.setAttribute("y", '{}'.format(y))
    layer.appendChild(text_el)

def _draw_text_superscript(layer, position, font_color, font_size, text,text2):
    """Draws text"""
    x = position[0]
    y = position[1]
    text_node = svg.createTextNode(text)
    sup_node = svg.createTextNode(text2)

    text_el = svg.createElement("text")
    text_el.appendChild(text_node)
    text_el.setAttribute("text-anchor", "left")
    text_el.setAttribute("font-family", "sans")
    text_el.setAttribute(
                "style", 
                "font-size:{}px; fill:{}".format(font_size, font_color))
    text_el.setAttribute("x", '{0:.2f}'.format(x))
    text_el.setAttribute("y", '{}'.format(y))
    layer.appendChild(text_el)    

    sup_el = svg.createElement("text")
    sup_el.appendChild(sup_node)
    # sup_el.setAttribute("text", 'Supscrit Text')
    sup_el.setAttribute("text-anchor", "middle")
    sup_el.setAttribute("font-family", "sans")
    sup_el.setAttribute(
                "style", 
                "font-size:{}px; fill:{}".format(font_size-2, font_color))
    sup_el.setAttribute("x", '{0:.2f}'.format(x))
    sup_el.setAttribute("y", '{}'.format(y))
    # sup_el.setAttribute("dx", '1em')
    # sup_el.setAttribute("dy", '.9em')
    layer.appendChild(sup_el)  

def _draw_line(layer, begin, end, name):
    """Draws line, begin and end should be a lists of [x,y]
    coordinates of line"""
    x1 = begin[0]
    y1 = begin[1]
    x2 = end[0]
    y2 = end[1]
    line = svg.createElement("line")
    line.setAttribute("id", str(name))
    line.setAttribute("stroke", "#000000")
    line.setAttribute("stroke-width", "1.0")
    line.setAttribute("x1", str(x1))
    line.setAttribute("y1", str(y1))
    line.setAttribute("x2", str(x2))
    line.setAttribute("y2", str(y2))
    layer.appendChild(line)

def _draw_isomer_line(layer, begin, end, name):
    """Draws line, begin and end should be a lists of [x,y]
    coordinates of line"""
    # print('isomer begin: ',begin)
    # print('isomer end: ',end)
    x1 = begin[0]
    y1 = begin[1]
    x2 = end[0]
    y2 = end[1]
    line = svg.createElement("line")
    line.setAttribute("id", str(name))
    line.setAttribute("stroke", "#000000")
    line.setAttribute("stroke-width", "0.5")
    line.setAttribute("x1", str(x1))
    line.setAttribute("y1", str(y1))
    line.setAttribute("x2", str(x2))
    line.setAttribute("y2", str(y2))
    layer.appendChild(line)    

def _draw_border_line(layer, begin, end, name):
    """Draws line, begin and end should be a lists of [x,y]
    coordinates of line"""
    # print('begin: ',begin)
    # print('end: ',end)
    x1 = begin[0]
    y1 = begin[1]
    x2 = end[0]
    y2 = end[1]
    line = svg.createElement("line")
    line.setAttribute("id", str(name))
    # line.setAttribute("stroke", "#ff0000")
    line.setAttribute("stroke", "#ffff00")
    line.setAttribute("stroke-width", "2.0")
    line.setAttribute("x1", str(x1))
    line.setAttribute("y1", str(y1))
    line.setAttribute("x2", str(x2))
    line.setAttribute("y2", str(y2))
    layer.appendChild(line)    

def draw_nuclide(nuclide, layers, position,  args, is_isomer=False):
    """ Draws nuclide data, including primary and secondary decay modes,
        and name of nuclide """

    # List of accepted basic decay modes, primary color is chosen on 
    # that basis. A '?' mode is for placeholders.
    basic_decay_modes = ['is', 'a', 'b-', 'b+',
                           'ec', 'p', '2p', 'sf', 
                           'n', '2n', '2ec', '2b-']
    # This reg ex. matches cluster emission marked by isotopic name
    # it matches names starting by at excatly two digits and
    # ending with letter(s). Remember that all decay modes are lower cased!
    # Cluster decays are only secondary or tertiary
    cluster_re = r'[0-9]{2}([a-z]+)$'

    # First decay mode should be largest and should match one
    # of basic decay modes
    if nuclide.decay_modes[0]['mode'] in basic_decay_modes:
        primary_color = COLORS[nuclide.decay_modes[0]['mode']]
    elif nuclide.decay_modes[0]['mode'] == '?' and args.unknown:
        primary_color = COLORS['?']
    else:
        # Order of basic and secondary decay modes is not kept well in NWC data
        # eg. sometimes B+p is given before B+
        # We will swap two element so any basic mode comes first
        if len(nuclide.decay_modes) > 1:
            i = 0
            for i in range(1, len(nuclide.decay_modes)):
                if nuclide.decay_modes[i]['mode'] in basic_decay_modes:
                    nuclide.decay_modes[0], nuclide.decay_modes[i] = nuclide.decay_modes[i], nuclide.decay_modes[0]
                    primary_color = COLORS[nuclide.decay_modes[0]['mode']]
                    break
            else:
                return
        else:
            return

    
    # Ommit p-unstable and n-unstable (this information
    # is in half-life)
    if nuclide.half_life['value'].find('unstable') >= 0:
        return

    # If there is more decay modes, and if at least one matches
    # basic decay modes, a secondary color will be used
    # Large triangle is used for modes with branching > 5%
    # but not for long lived nuclides (quasi-stable)
    # small triangle for other
    # If large triangle is used, a tertiary decay mode might
    # be indicated with small triangle
    secondary_size = None
    tertiary_size  = None
    if len(nuclide.decay_modes) > 1:
        for i in range(1, len(nuclide.decay_modes)):
            try:
                if nuclide.decay_modes[i]['mode'] in basic_decay_modes:
                    v = float(nuclide.decay_modes[i]['value'].strip('#'))
                    secondary_color = COLORS[nuclide.decay_modes[i]['mode']]
                    if v > 5.0 and primary_color != COLORS['is']:
                        secondary_size = 'large'
                    elif v > 0.0:
                        secondary_size = 'small'
                elif (re.search(cluster_re, nuclide.decay_modes[i]['mode']) 
                        is not None):
                    secondary_size = 'small'
                    secondary_color = COLORS['cluster']
                    break
            except ValueError:
                continue

        if ( len(nuclide.decay_modes) > 2 and
             ( secondary_size == 'large' or
               (secondary_size == 'small' and primary_color == COLORS['is'])) ):
            for i in range(2, len(nuclide.decay_modes)):
                try:
                    if nuclide.decay_modes[i]['mode'] in basic_decay_modes:
                        v = float(nuclide.decay_modes[i]['value'].strip('#'))
                        if v > 0.0:
                            tertiary_color = (
                                    COLORS[nuclide.decay_modes[i]['mode']])
                            tertiary_size  = 'small'
                        break
                    elif re.search(cluster_re, 
                                nuclide.decay_modes[i]['mode']) is not None:
                        tertiary_size = 'small'
                        tertiary_color = COLORS['cluster']
                        break
                except ValueError:
                    continue

    _draw_rectangle(layers[0], position,
                    primary_color, '{}0'.format(nuclide))

    if secondary_size == 'large':
        if secondary_color == COLORS['a']:
            corner = 'lt'
        elif ( secondary_color == COLORS['b+'] and 
               primary_color != COLORS['a'] and 
               primary_color != COLORS['p'] ) :
            corner = 'lt'
        else:
            corner = 'rb'
        _draw_triangle(layers[1], position, secondary_color,
                       '{}1'.format(nuclide), corner)
    elif secondary_size == 'small':
        if secondary_color == COLORS['a']:
            corner = 'lt'
        elif ( secondary_color == COLORS['b+'] and 
               primary_color != COLORS['a'] and 
               primary_color != COLORS['p'] ) :
            corner = 'lt'
        elif secondary_color == COLORS['cluster']:
            corner = 'rt'
        else:
            corner = 'rb'
        _draw_small_triangle(layers[1], position, secondary_color,
                             '{}1'.format(nuclide), corner)

    if tertiary_size == 'small':
        if tertiary_color == COLORS['a']:
            corner = 'lt'
        elif ( tertiary_color == COLORS['b+'] and 
               primary_color != COLORS['a'] and 
               secondary_color != COLORS['a']) :
            corner = 'lt'
        elif tertiary_color == COLORS['cluster']:
            corner = 'rt'
        else:
            corner = 'rb'
        _draw_small_triangle(layers[1], position, tertiary_color,
                             '{}2'.format(nuclide), corner)

    font_color = FONT_COLOR_BRIGHT if primary_color == COLORS['is'] else FONT_COLOR_DARK
    if (args.isomers and is_isomer):
    # if is_isomer:
        # print('Plotting isomer...', str(nuclide.A) ,nuclide.element)
        if args.names:
            # print('Plotting isomer...', str(nuclide.A) ,nuclide.element)
            element_name = "$*" + str(nuclide.A) + "*$" + nuclide.element 

            tx = position[0] + SIZE_SHAPE / 2 
            ty = position[1] + SIZE_GAP + 0.9 * SIZE_FONT

            _draw_text(layers[3], [tx, ty], font_color, SIZE_FONT-2, element_name)

            # print('tx: ',tx)
            # print('ty: ',ty)

            isomer_primary_color = primary_color



            target_decay_mode = nuclide.decay_modes
            # [0]
            # print('Nuclide: ',str(Z+N)+target_element)
            # print(nuclide.decay_modes)
            found_isomers = nuclide.isomers
            # found_isomers['half_life']
            # print('Isomeric states: \n', found_isomers)
            # print(len(found_isomers))
            # print(type(found_isomers))
            for possible_isomers in found_isomers:
                # print(type(possible_isomers))
                isomer_key = possible_isomers['half_life']
                # print(isomer_key)
                isomer_relation = isomer_key['relation']
                # print(isomer_relation)
                if isomer_relation=='=':
                    # print('Isomer found! ',str(Z+N)+target_element)
                    # print(isomer_key)
                    isomer_half_life = isomer_key['value']
                    isomer_half_life_units = isomer_key['unit']
                    # print('half-life: ', isomer_half_life, ' ', isomer_half_life_units)
                    # isomer_decay_modes = isomer_key['decay_modes']
                    if isomer_half_life_units in {'d', 'h', 'm','s','Yy','Zy','Ey','Py','Ty','Gy','My','ky','y'}:
                        # is_isomer=True
                        # print(type(possible_isomers['decay_modes']))
                        isomer_decay_modes = possible_isomers['decay_modes']
                        # for dict_index in range(len(isomer_decay_modes)):
                        # print('number of isomers: ',len(isomer_decay_modes))
                        # print('isomer decay modes:\n', isomer_decay_modes)
                        # print('isomer decay modes:', isomer_decay_modes[dict_index])

                        
                        # print("True isomer found for:", str(Z+N)+target_element,"!")

        

                        isomer_secondary_size = None
                        isomer_tertiary_size  = None
                        if len(isomer_decay_modes) > 0:
                            for i in range(0, 1):
                                pri_mode_str = isomer_decay_modes[i]['mode']
                                # print('Mode: ',pri_mode_str)
                                try:
                                    if (pri_mode_str in basic_decay_modes) or (pri_mode_str == 'it'):
                                        v = float(isomer_decay_modes[i]['value'].strip('#'))
                                        isomer_primary_color = COLORS[isomer_decay_modes[i]['mode']]
                                        # print('loop 0 color: ',isomer_primary_color)
                                except ValueError:
                                    continue                    

                        if len(isomer_decay_modes) > 1:
                            for i in range(1, len(isomer_decay_modes)):
                                # print('Secondary Mode: ',isomer_decay_modes[i]['mode'])
                                sec_mode_str = isomer_decay_modes[i]['mode']
                                # print(type(sec_mode_str))
                                # print(type('it'))
                                # print(sec_mode_str == 'it')
                                # print('Testing at line 575')
                                try:
                                    # print(nuclide.decay_modes[i]['mode'] == 'it')
                                    # print('Testing at line 580')
                                    if ((isomer_decay_modes[i]['mode'] in basic_decay_modes) or (sec_mode_str == 'it')):
                                        # print('Testing at line 582')
                                        v = float(isomer_decay_modes[i]['value'].strip('#'))
                                        isomer_secondary_color = COLORS[isomer_decay_modes[i]['mode']]
                                        # print('loop 1 secondary color: ',isomer_secondary_color)
                                        if v > 5.0 and isomer_primary_color != COLORS['is']:
                                            isomer_secondary_size = 'large'
                                        elif v > 0.0:
                                            isomer_secondary_size = 'small'
                                    elif (re.search(cluster_re, isomer_decay_modes[i]['mode']) 
                                            is not None):
                                        isomer_secondary_size = 'small'
                                        isomer_secondary_color = COLORS['cluster']
                                        break
                                except ValueError:
                                    # print('Warning at line 588')
                                    continue

                            if ( len(isomer_decay_modes) > 2 and
                                 ( isomer_secondary_size == 'large' or
                                   (isomer_secondary_size == 'small' and isomer_primary_color == COLORS['is'])) ):
                                for i in range(2, len(isomer_decay_modes)):
                                    try:
                                        if isomer_decay_modes[i]['mode'] in basic_decay_modes:
                                            v = float(isomer_decay_modes[i]['value'].strip('#'))
                                            if v > 0.0:
                                                isomer_tertiary_color = (
                                                        COLORS[isomer_decay_modes[i]['mode']])
                                                isomer_tertiary_size  = 'small'
                                            break
                                        elif re.search(cluster_re, 
                                                    isomer_decay_modes[i]['mode']) is not None:
                                            isomer_tertiary_size = 'small'
                                            isomer_tertiary_color = COLORS['cluster']
                                            break
                                    except ValueError:
                                        continue
                            # print('isomer secondary color: ',isomer_secondary_color)

                            # isomer_secondary_color = isomer_primary_color
                            
                        # _draw_rectangle(layers[0], position,
                        #                 primary_color, '{}0'.format(nuclide))

                        if isomer_secondary_size == 'large':
                            if isomer_secondary_color == COLORS['a']:
                                corner = 'lt'
                            elif ( isomer_secondary_color == COLORS['b+'] and 
                                   isomer_primary_color != COLORS['a'] and 
                                   isomer_primary_color != COLORS['p'] ) :
                                corner = 'lt'
                            else:
                                corner = 'rb'
                            _draw_triangle(layers[1], position, isomer_secondary_color,
                                           '{}1'.format(nuclide), corner)
                        elif isomer_secondary_size == 'small':
                            if isomer_secondary_color == COLORS['a']:
                                corner = 'lt'
                            elif ( isomer_secondary_color == COLORS['b+'] and 
                                   isomer_primary_color != COLORS['a'] and 
                                   isomer_primary_color != COLORS['p'] ) :
                                corner = 'lt'
                            elif isomer_secondary_color == COLORS['cluster']:
                                corner = 'rt'
                            else:
                                corner = 'rb'
                            _draw_small_isomer_triangle(layers[1], position, isomer_secondary_color,
                                                 '{}1'.format(nuclide), corner)

                        if isomer_tertiary_size == 'small':
                            if isomer_tertiary_color == COLORS['a']:
                                corner = 'lt'
                            elif ( isomer_tertiary_color == COLORS['b+'] and 
                                   isomer_primary_color != COLORS['a'] and 
                                   isomer_secondary_color != COLORS['a']) :
                                corner = 'lt'
                            elif isomer_tertiary_color == COLORS['cluster']:
                                corner = 'rt'
                            else:
                                corner = 'rb'
                            _draw_small_isomer_triangle(layers[1], position, isomer_tertiary_color,
                                                 '{}2'.format(nuclide), corner)

                        _draw_isomer_rectangle(layers[0], position, isomer_primary_color, "{}{}ib".format(nuclide.A,nuclide.Z))
                        loc_1 = [position[0]+(SIZE_SHAPE/2.0),  position[1] + (SIZE_SHAPE*0.3) +0.2]
                        # print('x_loc',loc_1)
                        loc_2 = [position[0]+(SIZE_SHAPE/2.0), position[1] + (SIZE_SHAPE) -0.2 ]
                        # loc_2 =[42.5, 34.0]
                        # print('y_loc',loc_2)

                        _draw_isomer_line(layers[1], loc_1, loc_2, "{}{}iline".format(nuclide.A,nuclide.Z))

                        # isomer_text_y = (loc_1[1] + loc_2[1]) * 0.5
                        isomer_text_y1 = position[1] + SIZE_SHAPE * 0.5
                        isomer_text_y2 = isomer_text_y1 + SIZE_FONT_HL*0.75
                        isomer_text_x1 = position[0] + 0.25*SIZE_SHAPE
                        isomer_text_x2 = position[0] + 0.75*SIZE_SHAPE

                        half_life_string = nuclide.half_life['value'] 
                        sci_re = r'^[-+]?[0-9]*\.?[0-9]+([eE]+[-+]?[0-9]+)$'
                        if re.search(sci_re, half_life_string) is not None:
                            try:
                                hl = float(half_life_string)
                                half_life_string = '{0:.1e}'.format(hl)
                            except TypeError:
                                pass
                        if nuclide.half_life['value'] != '?':
                            # half_life_string += 
                            half_life_unit = nuclide.half_life['unit']
                            if nuclide.half_life['relation'] != '=':
                                half_life_string = nuclide.half_life['relation'] 
                                half_life_unit = nuclide.half_life['value'] 
                        else:
                            half_life_string = nuclide.decay_modes[0]['value']+'%'


                        _draw_text(layers[3], [isomer_text_x2, isomer_text_y1], font_color, SIZE_FONT_HL-2, half_life_string)
                        _draw_text(layers[3], [isomer_text_x2, isomer_text_y2], font_color, SIZE_FONT_HL-2, half_life_unit)
                        _draw_text(layers[3], [isomer_text_x1, isomer_text_y1], font_color, SIZE_FONT_HL-2, isomer_half_life)
                        _draw_text(layers[3], [isomer_text_x1, isomer_text_y2], font_color, SIZE_FONT_HL-2, isomer_half_life_units)


            # _draw_isomer_line(layers[0], [145.0, 45.0], [161.0, 45.0], "{}{}il".format(nuclide.A,nuclide.Z))

            # isomer_color = COLORS[nuclide.decay_modes[0]['mode']]
            # print('color: ', isomer_color)
            # print('isomer decay mode: ', )

    # Normal plotting
    else:
        if args.names:
            # Swapping from <element symbol> AAA  to:   ^{AAA} <element symbol>
            # element_name = nuclide.element + " " + str(nuclide.A) 
            # element_name = nuclide.element 
            element_name = "$*" + str(nuclide.A) + "*$" + nuclide.element 
            # element_name = "<sup>" + str(nuclide.A) + "</sup>" + nuclide.element 
            # element_name = "heat".sub() + str(nuclide.A) + nuclide.element 
            # element_name2 = str(nuclide.A) 

            tx = position[0] + SIZE_SHAPE / 2 
            ty = position[1] + SIZE_GAP + 1.25 * SIZE_FONT

            _draw_text(layers[3], [tx, ty], font_color, SIZE_FONT, element_name)
            # _draw_text_superscript(layers[3], [tx, ty], font_color, SIZE_FONT, element_name, element_name2)
            # _draw_text(layers[3], [tx-7, ty-3], font_color, SIZE_FONT-2, element_name2)

        if (args.halflives and not(nuclide.half_life['extrapolated'] == 'True')):
            # For stable and quasi-stable nuclide print isotopic abundance
            # for unstable - half live
            if primary_color != COLORS['is']:
                half_life_string = nuclide.half_life['value'] 
                sci_re = r'^[-+]?[0-9]*\.?[0-9]+([eE]+[-+]?[0-9]+)$'
                if re.search(sci_re, half_life_string) is not None:
                    try:
                        hl = float(half_life_string)
                        half_life_string = '{0:.1e}'.format(hl)
                    except TypeError:
                        pass
                if nuclide.half_life['value'] != '?':
                    half_life_string +=  ' ' + nuclide.half_life['unit']
                    if nuclide.half_life['relation'] != '=':
                        half_life_string = nuclide.half_life['relation'] + ' ' + half_life_string
            else:
                half_life_string = nuclide.decay_modes[0]['value']+'%'

            # text position center
            tx = position[0] + SIZE_SHAPE / 2
            ty = position[1] + SIZE_SHAPE - 1.5 * SIZE_FONT_HL

            _draw_text(layers[3], [tx, ty], font_color, SIZE_FONT_HL,
                       half_life_string)

def draw_target_border(layers, n_magic, z_magic,
                             n_limits, z_limits, size):
    for N, limits in n_magic.items():
        x1 = (N - n_limits[0] + 1) * SIZE_FIELD + SIZE_GAP / 2
        x2 = x1
        y1 = size[1] - (limits[1] - z_limits[0] + 2) * SIZE_FIELD - SIZE_GAP / 2
        y2 = size[1] - (limits[0] - z_limits[0] + 1) * SIZE_FIELD - SIZE_GAP / 2
        _draw_border_line(layers[4], [x1, y1], [x2, y2], "{}n0".format(N))
        x1 = (N - n_limits[0] + 2) * SIZE_FIELD + SIZE_GAP / 2
        x2 = x1
        _draw_border_line(layers[4], [x1, y1], [x2, y2], "{}n1".format(N))

    Z_counter=0
    for Z, limits in z_magic.items():
        x1 = (limits[0] - n_limits[0] + 1) * SIZE_FIELD + SIZE_GAP / 2 
        x2 = (limits[1] - n_limits[0] + 2) * SIZE_FIELD + SIZE_GAP / 2
        y1 = size[1] - (Z - Z_counter - z_limits[0] + 2) * SIZE_FIELD - SIZE_GAP / 2
        y2 = y1
        _draw_border_line(layers[4], [x1, y1], [x2, y2], "{}z0".format(Z))
        y1 = size[1] - (Z - Z_counter - z_limits[0] + 1) * SIZE_FIELD - SIZE_GAP / 2
        y2 = y1
        _draw_border_line(layers[4], [x1, y1], [x2, y2], "{}z1".format(Z))  
        Z_counter+=1  

def draw_magic_lines(layers, n_magic, z_magic,
                             n_limits, z_limits, size):
    for N, limits in n_magic.items():
        x1 = (N - n_limits[0] + 1) * SIZE_FIELD + SIZE_GAP / 2
        x2 = x1
        y1 = size[1] - (limits[1] - z_limits[0] + 2) * SIZE_FIELD - SIZE_GAP / 2
        y2 = size[1] - (limits[0] - z_limits[0] + 1) * SIZE_FIELD - SIZE_GAP / 2
        _draw_line(layers[2], [x1, y1], [x2, y2], "{}n0".format(N))
        x1 = (N - n_limits[0] + 2) * SIZE_FIELD + SIZE_GAP / 2
        x2 = x1
        _draw_line(layers[2], [x1, y1], [x2, y2], "{}n1".format(N))

    for Z, limits in z_magic.items():
        x1 = (limits[0] - n_limits[0] + 1) * SIZE_FIELD + SIZE_GAP / 2 
        x2 = (limits[1] - n_limits[0] + 2) * SIZE_FIELD + SIZE_GAP / 2
        y1 = size[1] - (Z - z_limits[0] + 2) * SIZE_FIELD - SIZE_GAP / 2
        y2 = y1
        _draw_line(layers[2], [x1, y1], [x2, y2], "{}z0".format(Z))
        y1 = size[1] - (Z - z_limits[0] + 1) * SIZE_FIELD - SIZE_GAP / 2
        y2 = y1
        _draw_line(layers[2], [x1, y1], [x2, y2], "{}z1".format(Z))

def draw_numbers(layers, shape, n_limits, z_limits, size):
    for n in range(0, n_limits[1] - n_limits[0] + 1):
        if (n + n_limits[0]) % 2 == 0 and (n + n_limits[0] > 0):
            z_first = 0
            while ( not(shape[n][z_first]) and 
                    z_first < z_limits[1] - z_limits[0] + 1 ):
                z_first += 1
            x = (n + 1) * SIZE_FIELD + SIZE_GAP + SIZE_SHAPE / 2 
            y = size[1] - (z_first + 1) * SIZE_FIELD + SIZE_GAP + 1.25 * SIZE_FONT
            _draw_text(layers[3], [x, y], '#000000', 
                       SIZE_FONT * 1.5, str(n + n_limits[0]))
    for z in range(0, z_limits[1] - z_limits[0] + 1):
        if (z + z_limits[0] % 2) % 2 == 0 and (z + z_limits[0] > 0):
            n_first = 0
            while ( not(shape[n_first][z]) and 
                    n_first < n_limits[1] - n_limits[0] + 1 ):
                n_first += 1
            x = (n_first) * SIZE_FIELD + SIZE_SHAPE / 2 + 3 * SIZE_GAP
            y = size[1] - (z + 2) * SIZE_FIELD + SIZE_SHAPE / 2 +2 * SIZE_GAP
            _draw_text(layers[3], [x, y], '#000000', 
                       SIZE_FONT * 1.5, str(z + z_limits[0]))

class FooAction(argparse.Action):
        def __init__(self, option_strings, dest, nargs=None, **kwargs):
            if nargs is not None:
                raise ValueError("nargs not allowed")
            super(FooAction, self).__init__(option_strings, dest, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            print('%r %r %r' % (namespace, values, option_string))
            setattr(namespace, self.dest, values)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create SVG format chart of nuclides')

    parser.add_argument('datafile', type=argparse.FileType('r'), 
                         help='Input data base XML file (required)')
    parser.add_argument('outfile', type=argparse.FileType('w'), 
                         help='Output data SVG file (required)')
    parser.add_argument('--names', action='store_false', 
                        help='Disable names')
    parser.add_argument('--halflives', action='store_false', 
                        help='Disable half-lives')
    parser.add_argument('--isomers', action='store_false', 
                        help='Disable separate isomer plotting')
    parser.add_argument('--magic', action='store_false', 
                        help='Disable magic numbers')
    parser.add_argument('--numbers', action='store_false', 
                        help='Disable numbers along axis')
    parser.add_argument('--unknown', action='store_false', 
                        help='Disable isotopes with unknown decay mode')
    # parser.add_argument('--t', action='store_true', 
    #                     dest='list_of_targets', type=str,  help='Enable names for target isotope(s)')
    parser.add_argument('--p', 
        # action='FooAction',
            # action='store', 
                        dest="list_of_products", 
                        type=str,  help='Enable names for product isotope(s)')
    parser.add_argument('--t', 
        # action='FooAction',
            # action='store', 
                        dest="list_of_targets", 
                        type=str,  help='Enable names for target isotope(s)')
    parser.add_argument('--z', nargs=2, default=[0,120],
                        dest='Z', type=int, help='Atomic number Z range (%(type)s), default: %(default)s')
    parser.add_argument('--n', nargs=2, default=[0,180],
                        dest='N', type=int, help='Neutron number N range (%(type)s), default: %(default)s')
    args = parser.parse_args()
    # print(parser.parse_args())

    # print(args.N, args.Z)

    if args.list_of_products is not None:
        # Only print target names
        args.names=False
        args.halflives=False
        # print('Drawing products: ',args.list_of_products) 
        # p_string = args.list_of_products.strip()
        # print('*',p_string,'*')
        individual_products = list(set(args.list_of_products.strip().split(' ')))
        print('Drawing products: ',individual_products) 

    if args.list_of_targets is not None:
        # Only print target names
        args.names=False
        args.halflives=False
        # print('Drawing targets : ',args.list_of_targets) 
        individual_targets = list(set(args.list_of_targets.strip().split(' ')))
        print('Drawing targets : ',individual_targets) 


    if args.N[0] > args.N[1]:
        print('Wrong N range {}, {}'.format(args.N[0], args.N[1]))
        print('Try {} -h for more information'.format(sys.argv[0]))
        exit()

    if args.Z[0] > args.Z[1]:
        print('Wrong Z range {}, {}'.format(args.Z[0], args.Z[1]))
        print('Try {} -h for more information'.format(sys.argv[0]))
        exit()


    # Create document type SVG and document itself with proper headers
    dom = xml.dom.minidom.getDOMImplementation()
    doctype_svg = dom.createDocumentType("svg", "-//W3C//DTD SVG 1.1//EN", "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd")
    svg = dom.createDocument("http://www.w3.org/2000/svg", "svg", doctype_svg)
    root = svg.documentElement
    root.setAttribute("xmlns", "http://www.w3.org/2000/svg")
    root.setAttribute("version", "1.1")
    
    # XML parser does not guarantee to preserve order of elements.
    # In fact minidom sorts elements alphabetically.
    # In SVG file order is important (elements added later appear at
    # the top of previous). To keep good order of elements we introduce
    # 4 layers (groups)
    #
    # layer0 is intended for squares (primary decay mode)
    # layer1 for triangles (secondary decay mode)
    # layer2 for magic number lines etc.
    # layer3 for text

    layers = []
    for l in range(5):
        layer = svg.createElement("g")
        layer.setAttribute("id", "layer{}".format(l))
        layer.setAttribute("fill", "none")
        root.appendChild(layer)
        layers.append(layer)

    n_limits = [None, None]
    z_limits = [None, None]
    # print(args.datafile, args.N, args.Z, n_limits, z_limits)
    data = load_xml_nuclear_table(args.datafile, args.N, args.Z, n_limits, z_limits)
    # print(data)
    # lookup_data = load_xml_nuclear_table(args.datafile, args.N, args.Z, n_limits, z_limits) 
    # print(lookup_data)
    # Size of picture is now calculated, and proper attributes
    # are assigned to root element
    # Additional margins are added to provide space for numbers on the 
    # left and bottom
    size = [(n_limits[1] - n_limits[0] + 2) * SIZE_FIELD + SIZE_GAP,
            (z_limits[1] - z_limits[0] + 2) * SIZE_FIELD + SIZE_GAP]
    root.setAttribute("width", str(size[0]))
    root.setAttribute("height", str(size[1]))

    # This variable is used to draw numbers next to last 
    # first element in Z rows, and below first element in N column
    #
    # Could be used instead of n_magic and z_magic?
    shape = []
    for n in range(n_limits[0], n_limits[1] + 1):
        n_list = []
        for z in range(z_limits[0], z_limits[1] + 1):
            n_list.append(False)
        shape.append(n_list)

    n_magic = {}
    z_magic = {}
    n_targets = {}
    z_targets = {}
    Z_counter=0
    for nuclide in data:
        is_isomer=False
        if args.list_of_targets is not None:
            args.names=False
            args.halflives=False
            my_flag=False
        if args.list_of_products is not None:
            args.names=False
            args.halflives=False   
            my_flag=False 
        N = nuclide.N
        Z = nuclide.Z
        target_element = nuclide.element
        target_decay_mode = nuclide.decay_modes[0]
        # print('Nuclide: ',str(Z+N)+target_element)
        # print(nuclide.decay_modes)
        found_isomers = nuclide.isomers
        # found_isomers['half_life']
        # print('Isomeric states: \n', found_isomers)
        # print(len(found_isomers))
        # print(type(found_isomers))
        for possible_isomers in found_isomers:
            # print(type(possible_isomers))
            isomer_key = possible_isomers['half_life']
            # print(isomer_key)
            isomer_relation = isomer_key['relation']
            # print(isomer_relation)
            if isomer_relation=='=':
                # print('Isomer found! ',str(Z+N)+target_element)
                # print(isomer_key)
                isomer_half_life = isomer_key['value']
                isomer_half_life_units = isomer_key['unit']
                # print('half-life: ', isomer_half_life, ' ', isomer_half_life_units)
                if isomer_half_life_units in {'d', 'h', 'm','s','Yy','Zy','Ey','Py','Ty','Gy','My','ky','y'}:
                    is_isomer=True
                    # print("True isomer found for:", str(Z+N)+target_element,"!")
        # print('Isomer half-life:',found_isomers['half_life'])
        if N in MAGIC_NUMBERS:
            if n_magic.get(N) is not None:
                if n_magic[N][1] < Z:
                    n_magic[N][1] = Z
            else:
                n_magic[N] = [Z, Z]
        if Z in MAGIC_NUMBERS:
            if z_magic.get(Z) is not None:
                if z_magic[Z][1] < N:
                    z_magic[Z][1] = N
            else:
                z_magic[Z] = [N, N]

        # Draw names for any specified product isotopes        
        if args.list_of_products is not None:
            # args.names=False
            # args.halflives=False
            for target in individual_products:
                # print('Nuclide: ',str(Z+N)+target_element)
                cleaned_target = re.sub(r'\d+', '', target)
                # print('Cleaned Target: ',cleaned_target)
                try:
                    cleaned_AAA = re.search('(\d+?)',  target).group(1)
                except:
                    cleaned_AAA = '-1'
                # print('Cleaned AAA: ',cleaned_AAA)    
                # print("Cleaned target: ", cleaned_target)
                # print(type(args.list_of_targets))
                # print('Current nuclide: ',str(Z+N),target_element)
                if target_element == cleaned_target: 
                    # print('Cleaned Target: ',cleaned_target)
                    # print('Cleaned AAA: ',cleaned_AAA)
                    if (str(Z+N)+target_element) == target:
                        # in basic_decay_modes:
                        # and target_decay_mode =='is':
                        # print("Product Success!: ", str(Z+N),target_element)
                        try:
                            # my_flag=True
                            args.names=True
                            args.halflives=True
                        # draw_nuclide(nuclide, layers, [x, y], args)
                        except IndexError:
                            print('IndexError: nuclide {}'.format(nuclide))
                        # args.names=False
                        # args.halflives=False
                    # else:
                    #     args.names=False
                    #     args.halflives=False


        # Draw names for any natural abundance target isotopes or specified target isotopes        
        if args.list_of_targets is not None:
            # args.names=False
            # args.halflives=False
            # print("Individual targets: ",individual_targets)
            for target in individual_targets:
                # print('Nuclide: ',str(Z+N)+target_element)
                cleaned_target = re.sub(r'\d+', '', target)
                # print('Cleaned Target: ',cleaned_target)
                try:
                    cleaned_AAA = int(re.findall(r'\d+',  target)[0])
                except:
                    cleaned_AAA = '-1'
                # print('Cleaned AAA: ',cleaned_AAA)    
                # print("Cleaned target: ", cleaned_target)
                # print(type(args.list_of_targets))
                # print('Current nuclide: ',str(Z+N),target_element)
                if target_element == cleaned_target: 
                    # print('Cleaned Target: ',cleaned_target)
                    # print('Cleaned AAA: ',cleaned_AAA)
                    if  int(cleaned_AAA)==0 and nuclide.decay_modes[0]['mode'] =='is':
                        # in basic_decay_modes:
                        # and target_decay_mode =='is':
                        # print("nat Target Success!: ", str(Z+N),target_element)
                        try:
                            my_flag=True
                            args.names=True
                            args.halflives=True


                            # print('Z: ', Z)
                            # print('n expansion check: ',n_targets.get(N))
                            # print('z expansion check: ',z_targets.get(Z))

                            n_targets[N] = [Z, Z]
                            # print('n_targets: ',n_targets)
                            # Z_counter=0
                            # z_targets[Z] = [N, N]
                            # print('z_targets: ',z_targets)

                            # if n_targets.get(N) is not None:
                            #     print('')
                            #     # print('n_targets: ',n_targets)
                            # #     if n_targets[N][1] < Z:
                            # #         n_targets[N][1] = Z
                            # else:
                            #     n_targets[N] = [Z, Z]
                            #     print('n_targets: ',n_targets)
                            if z_targets.get(Z) is not None:
                                Z_counter+=1
                                z_targets[Z+Z_counter] = [N, N]
                                # Z_counter+=1
                                # print('')
                                # print('z_targets: ',z_targets)
                            #     if z_targets[Z][1] < N:
                            #         z_targets[Z][1] = N
                            else:
                                z_targets[Z] = [N, N]
                                # print('z_targets: ',z_targets)

                        # draw_nuclide(nuclide, layers, [x, y], args)
                        except IndexError:
                            print('IndexError: nuclide {}'.format(nuclide))
                        # args.names=False
                        # args.halflives=False
                    elif int(cleaned_AAA)!=0 and (str(Z+N)+target_element) == target and nuclide.decay_modes[0]['mode'] =='is':
                        # in basic_decay_modes:
                        # and target_decay_mode =='is':
                        # print("Stable Target Success!: ", str(Z+N),target_element)
                        try:
                            my_flag=True
                            args.names=True
                            args.halflives=True


                            # print('Z: ', Z)
                            # print('n expansion check: ',n_targets.get(N))
                            # print('z expansion check: ',z_targets.get(Z))

                            n_targets[N] = [Z, Z]
                            # print('n_targets: ',n_targets)
                            # Z_counter=0
                            # z_targets[Z] = [N, N]
                            # print('z_targets: ',z_targets)

                            # if n_targets.get(N) is not None:
                            #     print('')
                            #     # print('n_targets: ',n_targets)
                            # #     if n_targets[N][1] < Z:
                            # #         n_targets[N][1] = Z
                            # else:
                            #     n_targets[N] = [Z, Z]
                            #     print('n_targets: ',n_targets)
                            if z_targets.get(Z) is not None:
                                Z_counter+=1
                                z_targets[Z+Z_counter] = [N, N]
                                # Z_counter+=1
                                # print('')
                                # print('z_targets: ',z_targets)
                            #     if z_targets[Z][1] < N:
                            #         z_targets[Z][1] = N
                            else:
                                z_targets[Z] = [N, N]
                                # print('z_targets: ',z_targets)

                        # draw_nuclide(nuclide, layers, [x, y], args)
                        except IndexError:
                            print('IndexError: nuclide {}'.format(nuclide))
                        # args.names=False
                        # args.halflives=False
                    elif (str(Z+N)+target_element) == target:
                        # in basic_decay_modes:
                        # and target_decay_mode =='is':
                        # print("Matched Success!: ", str(Z+N),target_element)
                        try:
                            my_flag=True
                            args.names=True
                            args.halflives=True
                        # draw_nuclide(nuclide, layers, [x, y], args)
                        except IndexError:
                            print('IndexError: nuclide {}'.format(nuclide))
                        # args.names=False
                        # args.halflives=False
                    # else:
                    #     args.names=False
                    #     args.halflives=False

        shape[N - n_limits[0]][Z - z_limits[0]] = True

        # Position is passed for upper left corner of square
        x = (N - n_limits[0] + 1) * SIZE_FIELD + SIZE_GAP
        y = size[1] - (Z - z_limits[0] + 2) * SIZE_FIELD 
        try:
            draw_nuclide(nuclide, layers, [x, y], args,is_isomer)
        except IndexError:
            print('IndexError: nuclide {}'.format(nuclide))
    # if True:
    # Draw border around target isotopes
    try:
        # print('n_targets: ',n_targets)
        # print('z_targets: ',z_targets)
        draw_target_border(layers, n_targets, z_targets, n_limits, z_limits, size)
        # draw_nuclide(nuclide, layers, [x, y], args)
    except IndexError:
        print('IndexError: nuclide {}'.format(nuclide))
    if args.magic:
        # print('nmagic: ',n_magic)
        # print(type(n_magic))
        draw_magic_lines(layers, n_magic, z_magic, n_limits, z_limits, size)
    if args.numbers:
        draw_numbers(layers, shape, n_limits, z_limits, size)
    # if t:
    #   print(t)    
    #   print(list_of_targets) 
    # print(_element)
    # if args.list_of_targets is not None:
    #   print(args.list_of_targets) 
    #   individual_targets = args.list_of_targets.split(' ')
    #   for target in individual_targets:
    #       cleaned_target = re.sub(r'\d+', '', target)
    #       print("Cleaned target: ", cleaned_target)
    #   # print(type(args.list_of_targets))
       #    for nuclide in data:
       #        target_N = nuclide.N
       #        target_Z = nuclide.Z
       #        target_element = nuclide.element
       #        target_decay_mode = nuclide.decay_modes[0]
       #        # print(target_decay_mode)
       #        # print(str(target_Z+target_N),target_element)
       #        # print(nuclide.decay_modes)

       #        if target_element == cleaned_target: 
       #            if nuclide.decay_modes[0]['mode'] =='is':
       #            # in basic_decay_modes:
       #                # and target_decay_mode =='is':
       #                print("Success!: ", str(target_Z+target_N),target_element)
       #                try:
       #                    args.names=True
       #                    args.halflives=True
       #                    draw_nuclide(nuclide, layers, [x, y], args)
       #                except IndexError:
       #                    print('IndexError: nuclide {}'.format(nuclide))
       #                args.names=False
       #                args.halflives=False

        # print(data)
        # print(args.datafile, args.N, args.Z, n_limits, z_limits)
        # lookup_data = load_xml_nuclear_table(args.datafile, args.N, args.Z, n_limits, z_limits) 


        # lookup_data = load_xml_nuclear_table(args.datafile, args.N, args.Z,
     #                              n_limits, z_limits) 
        # print(lookup_data)

    args.outfile.write(svg.toprettyxml(indent="  ", encoding="utf-8").decode("utf-8"))


    # print(args.outfile)
    outfile_string = str(args.outfile).split("name='")[1].split("' mode='w'")[0]

    tempfile = open(outfile_string).read()
    # print(type(tempfile))
    # print(tempfile)
    tempfile = tempfile.replace('$*','<tspan baseline-shift = "super">')
    tempfile = tempfile.replace('*$','</tspan>')
    open(outfile_string, 'w').write(tempfile)

    # # Replace placeholders for superscript
    # # print(args.outfile)
    # # open('newfile.svg', 'w').write(open('outfile.svg').read().replace('$*','<tspan dx="-1ex" dy="+1ex"><tspan class="subscript">'))
    # # open('newfile.svg', 'w').write(open('outfile.svg').read().replace('$*','<tspan class="subscript" dy="+5px">'))
    # open('newfile.svg', 'w').write(open('outfile.svg').read().replace('$*','<tspan baseline-shift = "super">'))
    # # open('newfile.svg', 'w').write(open('outfile.svg').read().replace('$*','<tspan dy="-7" font-size=".7em">'))
    # # open('outfile.svg', 'w').write(open('newfile.svg').read().replace('*$','</tspan></tspan>'))
    # open('outfile.svg', 'w').write(open('newfile.svg').read().replace('*$','</tspan>'))

    # drawing = svg2rlg("outfile.svg")
    # renderPDF.drawToFile(drawing, "outfile.pdf")
    # cairosvg.svg2pdf(url='outfile.svg', write_to='image.pdf')

    # Save to PDF via inkscape
    bashCommand = "inkscape " + outfile_string + " -A " + outfile_string[:-3] + "pdf"
    # print(bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    # output, error = process.communicate()
