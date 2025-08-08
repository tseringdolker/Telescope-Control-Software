import sys
import random
import numpy as np
from datetime import datetime, timedelta
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QTabWidget, QGridLayout, QLineEdit, QGroupBox,
    QComboBox, QCheckBox, QFrame, QSizePolicy, QScrollArea, QMessageBox,
    QProgressBar, QSlider, QTableWidget, QTableWidgetItem, QHeaderView, QFileDialog, QSplitter
)
from PyQt5.QtCore import Qt, QTimer, pyqtSignal, QThread
from PyQt5.QtGui import QFont, QDoubleValidator, QColor, QPixmap
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from opcua import Client, ua
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_constellation
from astropy.time import Time
import astropy.units as u
from astropy.utils.iers import conf
import logging
import os
from threading import Event, Thread
import time
import shutil
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import requests
import json
from matplotlib.image import imread
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image
from io import BytesIO

# Color theme constants
COLOR_BACKGROUND = "#06141B"  # Dark background
COLOR_DARK = "#11212D"       # Darker seafoam
COLOR_MEDIUM = "#458B74"     # Dark Seafoam
COLOR_LIGHT = "#67CDC9"      # Primary Seafoam Green
COLOR_TEXT_LIGHT = "#E0FFFF" # Seafoam White
COLOR_TEXT_MEDIUM = "#8FBC8F" # Seafoam Grey
COLOR_ACCENT1 = "#67CDC9"    # Seafoam Green
COLOR_ACCENT2 = "#B2EBF2"    # Light Seafoam
COLOR_ACCENT3 = "#458B74"    # Dark Seafoam
COLOR_ACCENT4 = "#8FBC8F"    # Seafoam Grey

class SkyViewer(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent_window = parent
        self.location = EarthLocation(lat=32.7908*u.deg, lon=79.0002*u.deg, height=4507*u.m)
        self.current_time = Time.now()
        self.zoom_level = 1.0
        self.show_constellation_images = True
        self.show_planets = True
        self.show_constellation_lines = True
        self.show_milky_way = True
        self.weather_data = None
        self.setup_ui()
        self.update_timer = QTimer(self)
        self.update_timer.timeout.connect(self.update_sky)
        self.update_timer.start(1000)
        self.fetch_weather_data()

        # Initialize data
        self.stars = self.generate_star_positions()
        self.planets = self.generate_planet_positions()
        self.deep_sky_objects = self.generate_deep_sky_objects()
        self.constellations = self.generate_constellation_data()
        self.selected_object = None
        self.current_ra = 0
        self.current_dec = 0

    def setup_ui(self):
        layout = QVBoxLayout()

        # Time and weather controls
        control_bar = QHBoxLayout()

        # Time controls
        time_controls = QHBoxLayout()
        self.time_label = QLabel(f"Time: {self.current_time.iso[:-7]}")
        self.time_speed = QComboBox()
        self.time_speed.addItems(["Real Time", "10x", "100x", "Paused"])
        self.time_speed.currentIndexChanged.connect(self.change_time_speed)
        time_controls.addWidget(self.time_label)
        time_controls.addWidget(QLabel("Time Speed:"))
        time_controls.addWidget(self.time_speed)

        # Weather display
        self.weather_label = QLabel("Weather: Loading...")
        self.weather_label.setStyleSheet(f"color: {COLOR_TEXT_LIGHT}; font-style: italic;")

        control_bar.addLayout(time_controls)
        control_bar.addStretch()
        control_bar.addWidget(self.weather_label)

        # Search bar
        search_bar = QHBoxLayout()
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search for a planet, star or constellation...")
        self.search_button = QPushButton("Search")
        self.search_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT3};
                color: white;
                padding: 5px;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT1};
            }}
        """)
        search_bar.addWidget(self.search_input)
        search_bar.addWidget(self.search_button)

        # Visualization options
        options_group = QGroupBox("Visualization Options")
        options_layout = QHBoxLayout()

        self.constellation_images_check = QCheckBox("Constellation Images")
        self.constellation_images_check.setChecked(True)
        self.constellation_images_check.stateChanged.connect(self.toggle_constellation_images)

        self.planets_check = QCheckBox("Planets Images")
        self.planets_check.setChecked(True)
        self.planets_check.stateChanged.connect(self.toggle_planets)

        self.constellation_lines_check = QCheckBox("Constellation Lines")
        self.constellation_lines_check.setChecked(True)
        self.constellation_lines_check.stateChanged.connect(self.toggle_constellation_lines)

        self.milky_way_check = QCheckBox("Milky Way")
        self.milky_way_check.setChecked(True)
        self.milky_way_check.stateChanged.connect(self.toggle_milky_way)

        options_layout.addWidget(self.constellation_images_check)
        options_layout.addWidget(self.planets_check)
        options_layout.addWidget(self.constellation_lines_check)
        options_layout.addWidget(self.milky_way_check)
        options_group.setLayout(options_layout)

        # Create figure and canvas with dark theme
        self.figure = Figure(figsize=(10, 8), facecolor=COLOR_DARK)
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111, projection='aitoff')
        self.ax.set_facecolor(COLOR_DARK)

        # Zoom controls
        zoom_controls = QHBoxLayout()
        self.zoom_in_btn = QPushButton("+")
        self.zoom_in_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT1};
                color: white;
                font-weight: bold;
                min-width: 30px;
                max-width: 30px;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT3};
            }}
        """)
        self.zoom_in_btn.clicked.connect(self.zoom_in)

        self.zoom_out_btn = QPushButton("-")
        self.zoom_out_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT1};
                color: white;
                font-weight: bold;
                min-width: 30px;
                max-width: 30px;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT3};
            }}
        """)
        self.zoom_out_btn.clicked.connect(self.zoom_out)

        self.zoom_reset_btn = QPushButton("Reset Zoom")
        self.zoom_reset_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT3};
                color: white;
                padding: 5px;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT1};
            }}
        """)
        self.zoom_reset_btn.clicked.connect(self.reset_zoom)

        zoom_controls.addWidget(self.zoom_out_btn)
        zoom_controls.addWidget(self.zoom_reset_btn)
        zoom_controls.addWidget(self.zoom_in_btn)
        zoom_controls.addStretch()

        # Connect click event
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_hover)

        # Control buttons
        control_layout = QHBoxLayout()

        self.goto_btn = QPushButton("GoTo Selected Object")
        self.goto_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT1};
                color: white;
                padding: 8px;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT3};
            }}
        """)
        self.goto_btn.clicked.connect(self.goto_selected_object)
        self.goto_btn.setEnabled(False)

        self.center_btn = QPushButton("Center View")
        self.center_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT3};
                color: white;
                padding: 8px;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT1};
            }}
        """)
        self.center_btn.clicked.connect(self.center_view)

        control_layout.addWidget(self.goto_btn)
        control_layout.addWidget(self.center_btn)
        control_layout.addStretch()

        # Selected object info
        self.selected_object_info = QLabel("No object selected")
        self.selected_object_info.setStyleSheet(f"color: {COLOR_TEXT_LIGHT}; font-weight: bold;")
        self.selected_object_info.setWordWrap(True)

        # Current position info
        self.current_pos_info = QLabel("Current: RA 00:00:00, DEC +00:00:00")
        self.current_pos_info.setStyleSheet(f"color: {COLOR_TEXT_LIGHT};")

        # Hover info
        self.hover_info = QLabel("Hover over objects for details")
        self.hover_info.setStyleSheet(f"color: {COLOR_TEXT_MEDIUM}; font-style: italic;")

        # Build layout
        layout.addLayout(control_bar)
        layout.addLayout(search_bar)
        layout.addWidget(options_group)
        layout.addWidget(self.canvas)
        layout.addLayout(zoom_controls)
        layout.addWidget(self.current_pos_info)
        layout.addWidget(self.hover_info)
        layout.addWidget(self.selected_object_info)
        layout.addLayout(control_layout)
        self.setLayout(layout)

        # Initialize data
        self.stars = self.generate_star_positions()
        self.planets = self.generate_planet_positions()
        self.deep_sky_objects = self.generate_deep_sky_objects()
        self.constellations = self.generate_constellation_data()
        self.selected_object = None
        self.current_ra = 0
        self.current_dec = 0
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_hover)

    def fetch_weather_data(self):
        try:
            # Using OpenWeatherMap API (you'll need an API key)
            response = requests.get(
                "http://api.openweathermap.org/data/2.5/weather",
                params={
                    'lat': self.location.lat.value,
                    'lon': self.location.lon.value,
                    'appid': 'YOUR_API_KEY',  # Replace with your API key
                    'units': 'metric'
                }
            )
            self.weather_data = response.json()
            self.update_weather_display()
        except Exception as e:
            print(f"Weather data error: {e}")
            self.weather_data = None
            self.weather_label.setText("Weather: Data unavailable")

    def update_weather_display(self):
        if self.weather_data:
            temp = self.weather_data['main']['temp']
            wind_speed = self.weather_data['wind']['speed']
            conditions = self.weather_data['weather'][0]['description']
            self.weather_label.setText(
                f"Weather: {conditions.capitalize()}, {temp:.1f}°C, Wind: {wind_speed:.1f} m/s"
            )
        else:
            self.weather_label.setText("Weather: Data unavailable")

    def toggle_constellation_images(self, state):
        self.show_constellation_images = state == Qt.Checked
        self.plot_sky()

    def toggle_planets(self, state):
        self.show_planets = state == Qt.Checked
        self.plot_sky()

    def toggle_constellation_lines(self, state):
        self.show_constellation_lines = state == Qt.Checked
        self.plot_sky()

    def toggle_milky_way(self, state):
        self.show_milky_way = state == Qt.Checked
        self.plot_sky()

    def zoom_in(self):
        self.zoom_level *= 1.2
        self.plot_sky()

    def zoom_out(self):
        self.zoom_level /= 1.2
        self.plot_sky()

    def reset_zoom(self):
        self.zoom_level = 1.0
        self.plot_sky()

    def change_time_speed(self, index):
        speeds = [1, 10, 100, 0]
        interval = 1000 // speeds[index] if speeds[index] > 0 else 0
        self.update_timer.setInterval(interval)

    def update_sky(self):
        if self.time_speed.currentText() != "Paused":
            time_step = 1  # seconds
            if self.time_speed.currentText() == "10x":
                time_step = 10
            elif self.time_speed.currentText() == "100x":
                time_step = 100

            self.current_time = Time(self.current_time.datetime + timedelta(seconds=time_step))
            self.time_label.setText(f"Time: {self.current_time.iso[:-7]}")
            self.plot_sky()

    def generate_star_positions(self):
        # Generate star positions with more realistic data
        stars = []

        # Bright stars data (name, ra, dec, mag, spectral_type, distance_ly)
        bright_stars = [
            ("Sirius", "06:45:08.9", "-16:42:58", -1.46, "A1V", 8.6),
            ("Canopus", "06:23:57.1", "-52:41:44", -0.74, "F0II", 310),
            ("Arcturus", "14:15:39.7", "+19:10:56", -0.05, "K1.5III", 36.7),
            ("Vega", "18:36:56.3", "+38:47:01", 0.03, "A0V", 25),
            ("Capella", "05:16:41.4", "+45:59:53", 0.08, "G3III", 42.2),
            ("Rigel", "05:14:32.3", "-08:12:06", 0.18, "B8Ia", 860),
            ("Procyon", "07:39:18.1", "+05:13:30", 0.40, "F5IV-V", 11.46),
            ("Betelgeuse", "05:55:10.3", "+07:24:25", 0.45, "M2Iab", 643),
            ("Achernar", "01:37:42.8", "-57:14:12", 0.45, "B6Vep", 139),
            ("Hadar", "14:03:49.4", "-60:22:23", 0.61, "B1III", 390)
        ]

        for name, ra, dec, mag, spectral_type, distance in bright_stars:
            coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            stars.append({
                'name': name,
                'ra': ra,
                'dec': dec,
                'mag': mag,
                'type': 'star',
                'spectral_type': spectral_type,
                'distance': distance,
                'constellation': get_constellation(coord)
            })

        # Add fainter stars (100 stars)
        for i in range(100):
            ra_h = random.uniform(0, 24)
            dec_deg = random.uniform(-90, 90)
            mag = random.uniform(3, 6)
            name = f"Star-{i+1}"

            # Convert to hms/dms format
            ra = f"{int(ra_h)}:{int((ra_h%1)*60)}:{int(((ra_h%1)*60)%1*60):.1f}"
            dec_sign = '+' if dec_deg >= 0 else '-'
            dec_deg_abs = abs(dec_deg)
            dec = f"{dec_sign}{int(dec_deg_abs)}:{int((dec_deg_abs%1)*60)}:{int(((dec_deg_abs%1)*60)%1*60):.1f}"

            coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            stars.append({
                'name': name,
                'ra': ra,
                'dec': dec,
                'mag': mag,
                'type': 'star',
                'spectral_type': random.choice(["G2V", "K5III", "M0V", "B3V", "A0V", "F5IV"]),
                'distance': random.uniform(10, 1000),
                'constellation': get_constellation(coord)
            })

        return stars

    def generate_planet_positions(self):
        # Generate positions for major planets
        planets = [
            {
                'name': 'Mercury',
                'type': 'planet',
                'image': 'assets/mercury.png',
                'mag': -0.4,
                'size': 0.1,
                'description': 'Innermost planet, smallest in solar system'
            },
            {
                'name': 'Venus',
                'type': 'planet',
                'image': 'assets/venus.png',
                'mag': -4.6,
                'size': 0.2,
                'description': 'Second planet, hottest in solar system'
            },
            {
                'name': 'Mars',
                'type': 'planet',
                'image': 'assets/mars.png',
                'mag': -2.0,
                'size': 0.15,
                'description': 'Fourth planet, known as the Red Planet'
            },
            {
                'name': 'Jupiter',
                'type': 'planet',
                'image': 'assets/jupiter.png',
                'mag': -2.7,
                'size': 0.3,
                'description': 'Largest planet, gas giant with Great Red Spot'
            },
            {
                'name': 'Saturn',
                'type': 'planet',
                'image': 'assets/saturn.png',
                'mag': 0.7,
                'size': 0.25,
                'description': 'Second largest, known for its ring system'
            }
        ]

        # Calculate current positions (simplified)
        for planet in planets:
            # Simplified calculation - in a real app you'd use ephemeris data
            base_ra = random.uniform(0, 24)
            base_dec = random.uniform(-30, 30)

            planet['ra'] = f"{int(base_ra)}:{int((base_ra%1)*60)}:{int(((base_ra%1)*60)%1*60):.1f}"
            planet['dec'] = f"{'+' if base_dec >=0 else '-'}{abs(int(base_dec))}:{int((abs(base_dec)%1)*60)}:{int(((abs(base_dec)%1)*60)%1*60):.1f}"

        return planets

    def generate_constellation_data(self):
        # Constellation data with lines and images
        constellations = {
            'Orion': {
                'lines': [
                    [("05:40:45.5", "-02:23:59"), ("05:36:12.8", "-01:12:07")],  # Betelgeuse to Bellatrix
                    [("05:36:12.8", "-01:12:07"), ("05:25:07.9", "-02:22:58")],  # Bellatrix to Alnilam
                    [("05:25:07.9", "-02:22:58"), ("05:16:00.3", "-08:12:06")],  # Alnilam to Rigel
                    [("05:25:07.9", "-02:22:58"), ("05:14:32.3", "-08:12:06")],  # Alnilam to Saiph
                ],
                'image': 'assets/orion.png',
                'center': ("05:35:00", "-05:00:00"),
                'size': 0.3,
                'description': 'The Hunter, one of the most recognizable constellations'
            },
            'Ursa Major': {
                'lines': [
                    [("11:03:43.6", "+61:45:04"), ("11:53:49.8", "+53:41:41")],  # Dubhe to Merak
                    [("11:53:49.8", "+53:41:41"), ("12:15:25.6", "+57:01:57")],  # Merak to Phecda
                ],
                'image': 'assets/ursa_major.png',
                'center': ("12:00:00", "+55:00:00"),
                'size': 0.4,
                'description': 'The Great Bear, contains the Big Dipper asterism'
            },
            # Add more constellations as needed
        }
        return constellations

    def plot_sky(self):
        self.ax.clear()

        # Set the zoom level
        #self.ax.set_xlim(-np.pi/self.zoom_level, np.pi/self.zoom_level)
        #self.ax.set_ylim(-np.pi/2/self.zoom_level, np.pi/2/self.zoom_level)

        # Plot grid
        self.ax.grid(True, color=COLOR_LIGHT, linestyle='--', linewidth=0.5)

        # Plot Milky Way background if enabled
        if self.show_milky_way:
            try:
                milky_way_img = imread('assets/milky_way.png')
                self.ax.imshow(milky_way_img, extent=[-np.pi, np.pi, -np.pi/2, np.pi/2],
                              alpha=0.2, aspect='auto', origin='upper')
            except:
                pass  # Skip if image not available

        # Calculate current Alt/Az for all objects
        all_objects = self.stars + self.deep_sky_objects
        if self.show_planets:
            all_objects += self.planets

        coords = SkyCoord([obj['ra'] for obj in all_objects],
                         [obj['dec'] for obj in all_objects],
                         unit=(u.hourangle, u.deg))

        # Transform to current time and location
        altaz = coords.transform_to(AltAz(obstime=self.current_time, location=self.location))

        # Only plot objects above horizon
        above_horizon = altaz.alt > 0*u.deg
        visible_objects = [obj for obj, visible in zip(all_objects, above_horizon) if visible]
        visible_coords = coords[above_horizon]
        visible_altaz = altaz[above_horizon]

        # Convert to radians for plotting
        ras = visible_coords.ra.wrap_at(180*u.deg).radian
        decs = visible_coords.dec.radian
        mags = [obj['mag'] for obj in visible_objects]
        names = [obj.get('name') or obj.get('id') for obj in visible_objects]
        types = [obj.get('type', 'star') for obj in visible_objects]

        # Plot objects with size based on magnitude and type
        sizes = []
        colors = []
        for obj in visible_objects:
            if obj.get('object_type') == 'deep_sky':
                # Deep sky objects are larger and blueish
                sizes.append(80 * (10 - obj['mag']) / 10)
                colors.append(COLOR_ACCENT2)
            elif obj.get('type') == 'planet':
                # Planets are larger and yellow
                sizes.append(150 * (10 - obj['mag']) / 10)
                colors.append('#FFFF00')
            else:
                # Stars - size based on magnitude, color based on spectral type
                sizes.append(50 * (10 - obj['mag']) / 10)
                # Simple color based on spectral type
                if obj['spectral_type'][0] in ('O', 'B'):
                    colors.append('#9bb4ff')  # Blue
                elif obj['spectral_type'][0] in ('A', 'F'):
                    colors.append('#ffffff')  # White
                elif obj['spectral_type'][0] in ('G', 'K'):
                    colors.append('#ffd700')  # Yellow
                else:  # M
                    colors.append('#ff8c00')  # Orange

        # Plot all objects
        self.ax.scatter(ras, decs, s=sizes, c=colors, alpha=0.8, edgecolors='none')

        # Plot constellation images if enabled
        if self.show_constellation_images:
            for const_name, const_data in self.constellations.items():
                try:
                    coord = SkyCoord(const_data['center'][0], const_data['center'][1],
                                    unit=(u.hourangle, u.deg))
                    altaz = coord.transform_to(AltAz(obstime=self.current_time, location=self.location))

                    if altaz.alt > 0*u.deg:
                        img = Image.open(const_data['image'])
                        img = img.resize((50, 50))  # Resize image
                        img = np.array(img)

                        imagebox = OffsetImage(img, zoom=const_data['size'], alpha=0.5)
                        ab = AnnotationBbox(imagebox,
                                          (coord.ra.wrap_at(180*u.deg).radian,
                                           coord.dec.radian),
                                          frameon=False)
                        self.ax.add_artist(ab)
                except Exception as e:
                    print(f"Error loading constellation image: {e}")

        # Plot constellation lines if enabled
        if self.show_constellation_lines:
            for const_name, const_data in self.constellations.items():
                for line in const_data['lines']:
                    start = SkyCoord(line[0][0], line[0][1], unit=(u.hourangle, u.deg))
                    end = SkyCoord(line[1][0], line[1][1], unit=(u.hourangle, u.deg))

                    # Only plot if both stars are above horizon
                    start_altaz = start.transform_to(AltAz(obstime=self.current_time, location=self.location))
                    end_altaz = end.transform_to(AltAz(obstime=self.current_time, location=self.location))

                    if start_altaz.alt > 0*u.deg and end_altaz.alt > 0*u.deg:
                        x = [start.ra.wrap_at(180*u.deg).radian, end.ra.wrap_at(180*u.deg).radian]
                        y = [start.dec.radian, end.dec.radian]
                        self.ax.plot(x, y, color=COLOR_TEXT_MEDIUM, alpha=0.5, linewidth=1)

        # Label bright objects
        for ra, dec, name, obj_type in zip(ras, decs, names, types):
            if name and (abs(dec) < np.pi/2):  # Only label objects not at poles
                fontsize = 8 if obj_type == 'star' else 9
                color = COLOR_TEXT_LIGHT if obj_type == 'star' else COLOR_ACCENT2
                self.ax.text(ra, dec, name, color=color, fontsize=fontsize,
                            ha='center', va='center')

        # Highlight selected object if one is selected
        if self.selected_object:
            coord = SkyCoord(self.selected_object['ra'], self.selected_object['dec'],
                            unit=(u.hourangle, u.deg))
            self.ax.scatter(
                coord.ra.wrap_at(180*u.deg).radian,
                coord.dec.radian,
                s=200, color='red', alpha=0.7, marker='*')

        # Update current position marker
        if not np.isnan(self.current_ra) and not np.isnan(self.current_dec):
            coord = SkyCoord(self.current_ra, self.current_dec, unit=(u.hourangle, u.deg))
            self.ax.scatter(
                coord.ra.wrap_at(180*u.deg).radian,
                coord.dec.radian,
                s=100, color='yellow', alpha=0.8, marker='x')

        # Set labels
        self.ax.set_xticklabels(['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h'])

        # Add title with current time
        self.ax.set_title(f"Sky View - {self.current_time.iso[:-7]}", color=COLOR_TEXT_LIGHT)

        self.canvas.draw()

    def on_click(self, event):
        if event.inaxes != self.ax:
            return

        # Convert click coordinates to RA/DEC
        x, y = event.xdata, event.ydata
        if x is None or y is None:  # Click was outside the plot area
            return

        ra_rad = x % (2*np.pi)
        dec_rad = y

        # Convert to degrees
        ra_deg = np.degrees(ra_rad) % 360
        dec_deg = np.degrees(dec_rad)

        # Update current position
        self.current_ra = ra_deg / 15  # Convert to hours
        self.current_dec = dec_deg
        self.update_position_display()

        # Find nearest object
        min_dist = float('inf')
        nearest_obj = None

        all_objects = self.stars + self.deep_sky_objects + self.planets
        for obj in all_objects:
            coord = SkyCoord(obj['ra'], obj['dec'], unit=(u.hourangle, u.deg))
            dist = np.sqrt((coord.ra.deg - ra_deg)**2 + (coord.dec.deg - dec_deg)**2)

            if dist < min_dist and dist < 5:  # 5 degree threshold
                min_dist = dist
                nearest_obj = obj

        if nearest_obj:
            self.selected_object = nearest_obj
            self.show_object_details(nearest_obj)
            self.goto_btn.setEnabled(True)
            self.plot_sky()

    def update_position_display(self):
        ra_h = int(self.current_ra)
        ra_m = int((self.current_ra % 1) * 60)
        ra_s = int((((self.current_ra % 1) * 60) % 1 * 60))

        dec_sign = '+' if self.current_dec >= 0 else '-'
        dec_deg = abs(int(self.current_dec))
        dec_m = int((abs(self.current_dec) % 1 * 60))
        dec_s = int(((abs(self.current_dec) % 1 * 60) % 1 * 60))

        self.current_pos_info.setText(
            f"Current: RA {ra_h:02d}:{ra_m:02d}:{ra_s:02d}, DEC {dec_sign}{dec_deg:02d}:{dec_m:02d}:{dec_s:02d}"
        )

    def show_object_details(self, obj):
        if obj.get('type') == 'planet':
            # Planet details
            details = [
                f"Name: {obj.get('name', 'N/A')} (Planet)",
                f"RA: {obj['ra']}",
                f"DEC: {obj['dec']}",
                f"Magnitude: {obj.get('mag', 'N/A')}",
                f"Description: {obj.get('description', 'N/A')}",
                f"Current Time: {self.current_time.iso[:-7]}"
            ]

            # Add weather data if available
            if self.weather_data:
                details.extend([
                    f"Current Temperature: {self.weather_data['main']['temp']:.1f}°C",
                    f"Wind Speed: {self.weather_data['wind']['speed']:.1f} m/s",
                    f"Weather: {self.weather_data['weather'][0]['description'].capitalize()}"
                ])
        elif obj.get('object_type') == 'deep_sky':
            # Deep sky object details
            details = [
                f"Name: {obj.get('name', 'N/A')}",
                f"Catalog ID: {obj.get('id', 'N/A')}",
                f"Type: {obj.get('type', 'N/A')}",
                f"RA: {obj['ra']}",
                f"DEC: {obj['dec']}",
                f"Magnitude: {obj.get('mag', 'N/A')}",
                f"Size: {obj.get('size', 'N/A')}",
                f"Constellation: {obj.get('constellation', 'N/A')}",
                f"Current Time: {self.current_time.iso[:-7]}"
            ]
        else:
            # Star details
            coord = SkyCoord(obj['ra'], obj['dec'], unit=(u.hourangle, u.deg))
            altaz = coord.transform_to(AltAz(obstime=self.current_time, location=self.location))

            details = [
                f"Name: {obj.get('name', 'N/A')}",
                f"RA: {obj['ra']}",
                f"DEC: {obj['dec']}",
                f"Magnitude: {obj.get('mag', 'N/A')}",
                f"Spectral Type: {obj.get('spectral_type', 'N/A')}",
                f"Distance: {obj.get('distance', 'N/A')} light years",
                f"Constellation: {obj.get('constellation', 'N/A')}",
                f"Current Alt/Az: {altaz.alt.deg:.1f}°, {altaz.az.deg:.1f}°",
                f"Current Time: {self.current_time.iso[:-7]}"
            ]

            # Add weather data if available
            if self.weather_data:
                details.extend([
                    f"Current Temperature: {self.weather_data['main']['temp']:.1f}°C",
                    f"Wind Speed: {self.weather_data['wind']['speed']:.1f} m/s",
                    f"Weather: {self.weather_data['weather'][0]['description'].capitalize()}"
                ])

        self.selected_object_info.setText("\n".join(details))

    def on_hover(self, event):
        if event.inaxes != self.ax:
            self.hover_info.setText("Hover over objects for details")
            return

        # Convert hover coordinates to RA/DEC
        x, y = event.xdata, event.ydata
        if x is None or y is None:
            self.hover_info.setText("Hover over objects for details")
            return

        ra_rad = x % (2*np.pi)
        dec_rad = y

        # Convert to degrees
        ra_deg = np.degrees(ra_rad) % 360
        dec_deg = np.degrees(dec_rad)

        # Find nearest object
        min_dist = float('inf')
        nearest_obj = None

        all_objects = self.stars + self.deep_sky_objects
        if self.show_planets:
            all_objects += self.planets

        for obj in all_objects:
            coord = SkyCoord(obj['ra'], obj['dec'], unit=(u.hourangle, u.deg))
            dist = np.sqrt((coord.ra.deg - ra_deg)**2 + (coord.dec.deg - dec_deg)**2)

            if dist < min_dist and dist < 2:  # 2 degree threshold for hover
                min_dist = dist
                nearest_obj = obj

        if nearest_obj:
            if nearest_obj.get('object_type') == 'deep_sky':
                info = f"{nearest_obj.get('id', '')} {nearest_obj.get('name', '')} - {nearest_obj.get('type', '')}"
            elif nearest_obj.get('type') == 'planet':
                info = f"Planet: {nearest_obj.get('name', '')} - Mag: {nearest_obj.get('mag', '')}"
            else:
                info = f"{nearest_obj.get('name', '')} - Mag: {nearest_obj.get('mag', '')}"
            self.hover_info.setText(info)
        else:
            self.hover_info.setText("Hover over objects for details")

    def goto_selected_object(self):
        if not self.selected_object:
            return

        if not self.parent_window or not hasattr(self.parent_window, 'opcua_thread'):
            QMessageBox.warning(self, "Error", "OPC-UA connection not available")
            return

        try:
            coord = SkyCoord(self.selected_object['ra'], self.selected_object['dec'],
                            unit=(u.hourangle, u.deg))
            altaz = coord.transform_to(AltAz(obstime=self.current_time, location=self.parent_window.location))
            az_target = altaz.az.deg
            alt_target = altaz.alt.deg

            # Movement command sequence
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_Acc",
                1.0,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_velocity",
                30.0,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_Acc",
                1.0,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_velocity",
                3.0,
                ua.VariantType.Double
            )

            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_moveAbsolute",
                False,
                ua.VariantType.Boolean
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.EL_moveAbsolute",
                False,
                ua.VariantType.Boolean
            )

            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.Az_MoveAbsolutePosition",
                az_target,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.El_MoveAbsolutePosition",
                alt_target,
                ua.VariantType.Double
            )

            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_moveAbsolute",
                True,
                ua.VariantType.Boolean
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.EL_moveAbsolute",
                True,
                ua.VariantType.Boolean
            )

            QMessageBox.information(self, "GoTo Started",
                                  f"Slewing to {self.selected_object.get('name', self.selected_object.get('id', 'object'))}")

            # Update the mount tab with the new coordinates
            if hasattr(self.parent_window, 'mount_tab'):
                self.parent_window.mount_tab.goto_controls.ra_input.setText(self.selected_object['ra'])
                self.parent_window.mount_tab.goto_controls.dec_input.setText(self.selected_object['dec'])
                self.parent_window.mount_tab.goto_controls.az_display.setText(f"{az_target:.6f}")
                self.parent_window.mount_tab.goto_controls.alt_display.setText(f"{alt_target:.6f}")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to initiate GoTo: {str(e)}")

    def center_view(self):
        self.zoom_level = 1.0
        self.plot_sky()

    def generate_deep_sky_objects(self):
        # Generate some deep sky objects
        deep_sky = [
            {
                'id': 'M1',
                'name': 'Crab Nebula',
                'ra': '05:34:31.94',
                'dec': '+22:00:52.2',
                'mag': 8.4,
                'type': 'Supernova Remnant',
                'object_type': 'deep_sky',
                'constellation': 'Taurus',
                'size': '6x4 arcmin'
            },
            {
                'id': 'M31',
                'name': 'Andromeda Galaxy',
                'ra': '00:42:44.3',
                'dec': '+41:16:09',
                'mag': 3.4,
                'type': 'Galaxy',
                'object_type': 'deep_sky',
                'constellation': 'Andromeda',
                'size': '190x60 arcmin'
            },
            {
                'id': 'M42',
                'name': 'Orion Nebula',
                'ra': '05:35:17.3',
                'dec': '-05:23:28',
                'mag': 4.0,
                'type': 'Diffuse Nebula',
                'object_type': 'deep_sky',
                'constellation': 'Orion',
                'size': '85x60 arcmin'
            }
        ]
        return deep_sky

class StatusIndicator(QLabel):
    def __init__(self, text, parent=None):
        super().__init__(text, parent)
        self.setAlignment(Qt.AlignCenter)
        self.setMinimumWidth(60)
        self.set_status("off")

    def set_status(self, state):
        state_colors = {
            "on": COLOR_ACCENT1,
            "off": "#8B0000",
            "warning": "#FFA500",
            "error": "#d32f2f",
            "active": COLOR_ACCENT3,
            "closed": "#8B0000",
            "open": COLOR_ACCENT1,
            "moving": "#FFA500",
            "idle": COLOR_TEXT_MEDIUM,
            "none": COLOR_TEXT_MEDIUM,
            "exposing": COLOR_ACCENT1
        }

        state_text = state.upper()
        color = state_colors.get(state.lower(), COLOR_TEXT_MEDIUM)

        self.setText(state_text)
        self.setStyleSheet(f"""
            QLabel {{
                background-color: {color};
                color: white;
                padding: 3px;
                border-radius: 4px;
                font-weight: bold;
                min-width: 60px;
                max-width: 60px;
            }}
        """)

class TrackingThread(QThread):
    position_updated = pyqtSignal(float, float)  # az, alt

    def __init__(self, ra_str, dec_str, location):
        super().__init__()
        self.ra_str = ra_str
        self.dec_str = dec_str
        self.location = location
        self.running = False

    def run(self):
        self.running = True
        while self.running:
            try:
                coord = SkyCoord(self.ra_str, self.dec_str, unit=(u.hourangle, u.deg), frame='icrs')
                altaz = coord.transform_to(AltAz(obstime=Time.now(), location=self.location))
                az = round(altaz.az.deg, 6)
                alt = round(altaz.alt.deg, 6)
                self.position_updated.emit(az, alt)
            except Exception as e:
                print(f"Tracking error: {e}")

            self.msleep(10)  # 10ms interval for continuous updates

    def stop(self):
        self.running = False
        self.wait()

class OPCUAClientThread(QThread):
    data_updated = pyqtSignal(dict)
    connection_status = pyqtSignal(bool)
    error_data = pyqtSignal(float)

    def __init__(self, server_url):
        super().__init__()
        self.server_url = server_url
        self.running = False
        self.client = None
        self.location = EarthLocation(lat=32.7908*u.deg, lon=79.0002*u.deg, height=4507*u.m)

    def run(self):
        try:
            self.client = Client(self.server_url)
            self.client.connect()
            self.connection_status.emit(True)
            self.running = True

            error_history = []
            target_az = 0
            target_alt = 0
            current_az = 0
            current_alt = 0

            while self.running:
                try:
                    data = {
                        'az_enable': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.MC_Power_Az.Enable"),
                        'el_enable': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.MC_Power_El.Enable"),
                        'az_pos': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.Az_actpos"),
                        'el_pos': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.El_actpos"),
                        'az_done': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.MC_MoveAbsolute_Az.Done"),
                        'el_done': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.MC_MoveAbsolute_El.Done"),
                        'error': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.FollowingError"),
                        'fwhm': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.MASS_DIMM_FWHM"),
                        'r0': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.MASS_DIMM_r0"),
                        'tau0': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.MASS_DIMM_tau0"),
                        'theta0': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.MASS_DIMM_theta0"),
                        'az_target': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.Az_target"),
                        'el_target': self.read_node("ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.El_target")
                    }

                    self.data_updated.emit(data)

                    if data['error'] is not None:
                        error_history.append(data['error'])
                        if len(error_history) > 100:
                            error_history = error_history[-100:]
                        self.error_data.emit(data['error'])

                    if data['az_target'] is not None:
                        target_az = data['az_target']
                    if data['el_target'] is not None:
                        target_alt = data['el_target']
                    if data['az_pos'] is not None:
                        current_az = data['az_pos']
                    if data['el_pos'] is not None:
                        current_alt = data['el_pos']

                    # Calculate following error as difference between target and current position
                    following_error = np.sqrt((target_az - current_az)**2 + (target_alt - current_alt)**2)
                    self.error_data.emit(following_error)

                except Exception as e:
                    print(f"Error reading OPC-UA data: {e}")

                self.msleep(100)

        except Exception as e:
            print(f"OPC-UA Error: {e}")
            self.connection_status.emit(False)
        finally:
            if self.client:
                try:
                    self.client.disconnect()
                except Exception as e:
                    print(f"Error during disconnect: {e}")

    def read_node(self, node_str):
        try:
            node = self.client.get_node(node_str)
            return node.get_value()
        except Exception:
            return None

    def write_node(self, node_str, value, variant_type):
        try:
            node = self.client.get_node(node_str)
            var = ua.Variant(value, variant_type)
            node.set_value(var)
            return True
        except Exception as e:
            print(f"Write error: {e}")
            return False

    def stop(self):
        self.running = False
        self.wait()

class ErrorPlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure(figsize=(12, 5), dpi=100)
        self.axes = fig.add_subplot(111)
        super().__init__(fig)
        self.setParent(parent)
        self.setStyleSheet(f"border: 1px solid {COLOR_LIGHT};")
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.time_data = np.linspace(0, 10, 100)
        self.error_data = np.zeros(100)
        self.plot_empty()

    def plot_empty(self):
        self.axes.clear()
        self.axes.set_facecolor(COLOR_DARK)
        for spine in self.axes.spines.values():
            spine.set_color(COLOR_TEXT_MEDIUM)
            spine.set_linewidth(0.5)
        self.axes.tick_params(axis='both', colors=COLOR_TEXT_MEDIUM, labelsize=8)
        self.axes.grid(True, color=COLOR_LIGHT, linestyle='--', linewidth=0.5)
        self.axes.set_title("Following Error Plot (Disconnected)", color=COLOR_TEXT_LIGHT, fontsize=10)
        self.axes.set_xlabel("Time (s)", color=COLOR_TEXT_MEDIUM, fontsize=9)
        self.axes.set_ylabel("Error (arcsec)", color=COLOR_TEXT_MEDIUM, fontsize=9)
        self.axes.set_ylim(-100, 100)
        self.draw()

    def update_plot(self, error_value):
        self.error_data = np.roll(self.error_data, -1)
        self.error_data[-1] = error_value

        self.axes.clear()
        self.axes.plot(self.time_data, self.error_data, color=COLOR_ACCENT2, linewidth=1.5)
        self.axes.set_facecolor(COLOR_DARK)

        for spine in self.axes.spines.values():
            spine.set_color(COLOR_TEXT_MEDIUM)
            spine.set_linewidth(0.5)

        self.axes.tick_params(axis='both', colors=COLOR_TEXT_MEDIUM, labelsize=8)
        self.axes.grid(True, color=COLOR_LIGHT, linestyle='--', linewidth=0.5)
        self.axes.set_title("Following Error Plot", color=COLOR_TEXT_LIGHT, fontsize=10)
        self.axes.set_xlabel("Time (s)", color=COLOR_TEXT_MEDIUM, fontsize=9)
        self.axes.set_ylabel("Error (arcsec)", color=COLOR_TEXT_MEDIUM, fontsize=9)

        y_margin = max(100, np.max(np.abs(self.error_data)) * 1.1)
        self.axes.set_ylim(-y_margin, y_margin)
        self.draw()

class GotoControls(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.location = EarthLocation(lat=32.7908*u.deg, lon=79.0002*u.deg, height=4507*u.m)
        self.setup_ui()
        self.update_time()
        conf.auto_download = True
        self.tracking_thread = None
        self.parent_window = parent

    def setup_ui(self):
        layout = QVBoxLayout()

        goto_group = QGroupBox("GOTO Controls")
        goto_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        goto_layout = QVBoxLayout()

        radec_frame = QFrame()
        radec_layout = QVBoxLayout()

        ra_hbox = QHBoxLayout()
        ra_hbox.addWidget(QLabel("RA:"))
        self.ra_input = QLineEdit("00:00:00")
        ra_hbox.addWidget(self.ra_input)
        radec_layout.addLayout(ra_hbox)

        dec_hbox = QHBoxLayout()
        dec_hbox.addWidget(QLabel("DEC:"))
        self.dec_input = QLineEdit("+00:00:00")
        dec_hbox.addWidget(self.dec_input)
        radec_layout.addLayout(dec_hbox)

        # Time display
        time_hbox = QHBoxLayout()
        self.utc_label = QLabel("UTC: -")
        self.sidereal_label = QLabel("LST: -")
        time_hbox.addWidget(self.utc_label)
        time_hbox.addWidget(self.sidereal_label)
        radec_layout.addLayout(time_hbox)

        # Convert button
        convert_btn = QPushButton("Convert")
        convert_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT3};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT1};
            }}
        """)
        convert_btn.clicked.connect(self.convert_coordinates)
        radec_layout.addWidget(convert_btn)

        radec_btn_layout = QHBoxLayout()
        self.radec_go_btn = QPushButton("GO")
        self.radec_go_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT1};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT3};
            }}
        """)
        self.radec_stop_btn = QPushButton("STOP")
        self.radec_stop_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: #8B0000;
                color: white;
            }}
            QPushButton:hover {{
                background-color: #A52A2A;
            }}
        """)
        self.track_btn = QPushButton("Track")
        self.track_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT3};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT1};
            }}
        """)
        radec_btn_layout.addWidget(self.radec_go_btn)
        radec_btn_layout.addWidget(self.radec_stop_btn)
        radec_btn_layout.addWidget(self.track_btn)
        radec_layout.addLayout(radec_btn_layout)

        radec_frame.setLayout(radec_layout)
        goto_layout.addWidget(radec_frame)

        # Add Alt/Az display from coordinate conversion
        altaz_frame = QFrame()
        altaz_layout = QGridLayout()

        altaz_layout.addWidget(QLabel("Azimuth:"), 0, 0)
        self.az_display = QLabel("0.000")
        self.az_display.setAlignment(Qt.AlignRight)
        altaz_layout.addWidget(self.az_display, 0, 1)

        altaz_layout.addWidget(QLabel("Altitude:"), 1, 0)
        self.alt_display = QLabel("0.000")
        self.alt_display.setAlignment(Qt.AlignRight)
        altaz_layout.addWidget(self.alt_display, 1, 1)

        # Tracking status label
        self.tracking_status = QLabel("Tracking: Not Active")
        altaz_layout.addWidget(self.tracking_status, 2, 0, 1, 2)

        altaz_frame.setLayout(altaz_layout)
        goto_layout.addWidget(altaz_frame)

        goto_group.setLayout(goto_layout)
        layout.addWidget(goto_group)

        self.setLayout(layout)

        # Connect signals
        self.radec_go_btn.clicked.connect(self.goto_position)
        self.radec_stop_btn.clicked.connect(self.stop_motion)
        self.track_btn.clicked.connect(self.toggle_tracking)

    def update_time(self):
        current_time = Time.now()
        self.utc_label.setText(f"UTC: {current_time.iso[:-7]}")
        lst = current_time.sidereal_time('mean', longitude=self.location.lon)
        self.sidereal_label.setText(f"LST: {lst.to_string(sep=':', precision=0)}")
        QTimer.singleShot(1000, self.update_time)

    def convert_coordinates(self):
        try:
            ra_str = self.ra_input.text().strip()
            dec_str = self.dec_input.text().strip()

            coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg), frame='icrs')
            altaz = coord.transform_to(AltAz(obstime=Time.now(), location=self.location))

            self.az_display.setText(f"{altaz.az.deg:.6f}")
            self.alt_display.setText(f"{altaz.alt.deg:.6f}")
            return altaz.az.deg, altaz.alt.deg

        except Exception as e:
            QMessageBox.warning(self, "Conversion Error", f"Invalid coordinates or conversion error:\n{str(e)}")
            return None, None

    def goto_position(self):
        az, alt = self.convert_coordinates()
        if az is None or alt is None:
            return

        if not hasattr(self.parent_window, 'opcua_thread') or not self.parent_window.opcua_thread.isRunning():
            QMessageBox.information(self, "Coordinates Converted",
                                  f"Coordinates converted to:\nAzimuth: {az:.6f}°\nAltitude: {alt:.6f}°\n\n"
                                  "Telescope is not connected - cannot move.")
            return

        try:
            # Movement command sequence
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_Acc",
                1.0,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_velocity",
                30.0,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_Acc",
                1.0,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_velocity",
                3.0,
                ua.VariantType.Double
            )

            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_moveAbsolute",
                False,
                ua.VariantType.Boolean
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.EL_moveAbsolute",
                False,
                ua.VariantType.Boolean
            )

            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.Az_MoveAbsolutePosition",
                az,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.El_MoveAbsolutePosition",
                alt,
                ua.VariantType.Double
            )

            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_moveAbsolute",
                True,
                ua.VariantType.Boolean
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.EL_moveAbsolute",
                True,
                ua.VariantType.Boolean
            )

            QMessageBox.information(self, "GoTo Started", f"Slewing to Az={az:.3f}°, Alt={alt:.3f}°")

            # Update mount status if available
            if hasattr(self.parent_window, 'status_indicators'):
                self.parent_window.status_indicators['moving'].set_status("moving")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to initiate GoTo: {str(e)}")

    def stop_motion(self):
        if not hasattr(self.parent_window, 'opcua_thread') or not self.parent_window.opcua_thread.isRunning():
            QMessageBox.information(self, "Not Connected", "Telescope is not connected - cannot stop motion.")
            return

        try:
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_halt2",
                True,
                ua.VariantType.Boolean
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_halt",
                True,
                ua.VariantType.Boolean
            )

            # Update mount status if available
            if hasattr(self.parent_window, 'status_indicators'):
                self.parent_window.status_indicators['moving'].set_status("idle")

            QMessageBox.information(self, "Movement Stopped", "Telescope movement stopped")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to stop movement: {str(e)}")

    def toggle_tracking(self):
        if self.tracking_thread and self.tracking_thread.isRunning():
            self.stop_tracking()
        else:
            self.start_tracking()

    def start_tracking(self):
        ra_str = self.ra_input.text().strip()
        dec_str = self.dec_input.text().strip()
        if not ra_str or not dec_str:
            QMessageBox.warning(self, "Input Error", "Please enter both RA and DEC coordinates")
            return

        try:
            # Validate input
            coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg), frame='icrs')

            if self.tracking_thread and self.tracking_thread.isRunning():
                self.tracking_thread.stop()
                self.tracking_thread.wait()

            self.tracking_thread = TrackingThread(ra_str, dec_str, self.location)
            self.tracking_thread.position_updated.connect(self.update_tracking_position)
            self.tracking_thread.start()

            self.track_btn.setText("Stop Tracking")
            self.tracking_status.setText("Tracking: Active (Computing)")

            # Update mount status if connected
            if hasattr(self.parent_window, 'status_indicators'):
                self.parent_window.status_indicators['tracking'].set_status("on")

            QMessageBox.information(self, "Tracking Started",
                                   f"Tracking RA={ra_str}, DEC={dec_str}\n"
                                   "Telescope will follow when connected.")

        except Exception as e:
            QMessageBox.warning(self, "Invalid Coordinates", f"Cannot start tracking:\n{e}")

    def stop_tracking(self):
        if self.tracking_thread and self.tracking_thread.isRunning():
            self.tracking_thread.stop()
            self.tracking_thread.wait()

        self.track_btn.setText("Track")
        self.tracking_status.setText("Tracking: Not Active")

        # Update mount status if connected
        if hasattr(self.parent_window, 'status_indicators'):
            self.parent_window.status_indicators['tracking'].set_status("off")

        QMessageBox.information(self, "Tracking Stopped", "Tracking stopped")

    def update_tracking_position(self, az, alt):
        """Update Az/Alt display and optionally send to telescope if connected"""
        try:
            # Always update the display with calculated positions
            self.az_display.setText(f"{az:.6f}")
            self.alt_display.setText(f"{alt:.6f}")
            self.tracking_status.setText(f"Tracking: Active (Az={az:.2f}°, Alt={alt:.2f}°)")

            # Only send to OPC UA if connected and tracking is active
            if (self.tracking_thread and self.tracking_thread.isRunning() and
                hasattr(self.parent_window, 'opcua_thread') and
                self.parent_window.opcua_thread and
                self.parent_window.opcua_thread.isRunning() and
                self.parent_window.opcua_thread.client):
                try:
                    self.parent_window.opcua_thread.write_node(
                        "ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.Az_MoveAbsolutePosition",
                        az,
                        ua.VariantType.Double)
                    self.parent_window.opcua_thread.write_node(
                        "ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.El_MoveAbsolutePosition",
                        alt,
                        ua.VariantType.Double)
                except Exception as opc_error:
                    print(f"[OPC-UA Write Error] {opc_error}")

        except Exception as e:
            print(f"[Tracking Update Error] {e}")

class MassDimmWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.connected = False
        self.setup_ui()
        self.setup_timers()

    def setup_ui(self):
        layout = QVBoxLayout()

        status_group = QGroupBox("MASS-DIMM Status")
        status_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        status_layout = QGridLayout()

        self.status_indicators = {
            'power': StatusIndicator("OFF"),
            'cooling': StatusIndicator("OFF"),
            'ao': StatusIndicator("OFF"),
            'shutter': StatusIndicator("CLOSED")
        }

        for i, (label, widget) in enumerate(self.status_indicators.items()):
            status_layout.addWidget(QLabel(label.capitalize() + ":"), i//2, (i%2)*2)
            status_layout.addWidget(widget, i//2, (i%2)*2+1)

        status_group.setLayout(status_layout)
        layout.addWidget(status_group)

        seeing_group = QGroupBox("Atmospheric Parameters")
        seeing_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        seeing_layout = QGridLayout()
        self.seeing_labels = {
            'fwhm': QLabel("0.00 arcsec"),
            'r0': QLabel("0.00 m"),
            'tau0': QLabel("0.00 ms"),
            'theta0': QLabel("0.00 arcsec"),
            'seeing': QLabel("0.00 arcsec"),
            'strehl': QLabel("0.00")
        }

        seeing_layout.addWidget(QLabel("FWHM:"), 0, 0)
        seeing_layout.addWidget(self.seeing_labels['fwhm'], 0, 1)
        seeing_layout.addWidget(QLabel("r0 (Fried):"), 1, 0)
        seeing_layout.addWidget(self.seeing_labels['r0'], 1, 1)
        seeing_layout.addWidget(QLabel("τ0 (Coherence):"), 2, 0)
        seeing_layout.addWidget(self.seeing_labels['tau0'], 2, 1)
        seeing_layout.addWidget(QLabel("θ0 (Isoplanatic):"), 3, 0)
        seeing_layout.addWidget(self.seeing_labels['theta0'], 3, 1)
        seeing_layout.addWidget(QLabel("Seeing:"), 4, 0)
        seeing_layout.addWidget(self.seeing_labels['seeing'], 4, 1)
        seeing_layout.addWidget(QLabel("Strehl Ratio:"), 5, 0)
        seeing_layout.addWidget(self.seeing_labels['strehl'], 5, 1)

        for i in range(6):
            seeing_layout.itemAtPosition(i, 1).widget().setAlignment(Qt.AlignRight)

        seeing_group.setLayout(seeing_layout)
        layout.addWidget(seeing_group)

        button_layout = QHBoxLayout()
        self.start_btn = QPushButton("Start")
        self.start_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT1};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT3};
            }}
        """)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: #8B0000;
                color: white;
            }}
            QPushButton:hover {{
                background-color: #A52A2A;
            }}
        """)
        self.calibrate_btn = QPushButton("Calibrate")
        self.calibrate_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT3};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT1};
            }}
        """)

        button_layout.addWidget(self.start_btn)
        button_layout.addWidget(self.stop_btn)
        button_layout.addWidget(self.calibrate_btn)
        layout.addLayout(button_layout)

        self.plot = FigureCanvas(Figure(figsize=(10, 4)))
        self.plot.setStyleSheet(f"border: 1px solid {COLOR_LIGHT};")
        self.plot.axes = self.plot.figure.add_subplot(111)
        self.plot.axes.set_facecolor(COLOR_DARK)
        for spine in self.plot.axes.spines.values():
            spine.set_color(COLOR_TEXT_MEDIUM)
            spine.set_linewidth(0.5)
        self.plot.axes.tick_params(axis='both', colors=COLOR_TEXT_MEDIUM)
        self.plot.axes.set_title("Seeing Monitor (Disconnected)", color=COLOR_TEXT_LIGHT)
        self.plot.axes.set_xlabel("Time (min)", color=COLOR_TEXT_MEDIUM)
        self.plot.axes.set_ylabel("FWHM (arcsec)", color=COLOR_TEXT_MEDIUM)
        self.plot.axes.grid(True, color=COLOR_LIGHT, linestyle='--')
        self.plot.axes.set_ylim(0, 2.0)

        layout.addWidget(self.plot)
        self.setLayout(layout)

        self.start_btn.clicked.connect(self.start_measurement)
        self.stop_btn.clicked.connect(self.stop_measurement)
        self.calibrate_btn.clicked.connect(self.calibrate)

    def setup_timers(self):
        self.measurement_timer = QTimer(self)
        self.measurement_timer.timeout.connect(self.update_measurement)
        self.running = False
        self.measurement_data = []
        self.measurement_times = []

    def set_connected(self, connected):
        self.connected = connected
        if not connected:
            self.stop_measurement()
            self.plot.axes.set_title("Seeing Monitor (Disconnected)", color=COLOR_TEXT_LIGHT)
            self.plot.draw()

    def start_measurement(self):
        if not self.connected:
            QMessageBox.warning(self, "Connection Error", "Not connected to telescope")
            return

        if not self.running:
            self.running = True
            self.measurement_timer.start(1000)
            self.status_indicators['power'].set_status("on")
            self.status_indicators['cooling'].set_status("active")
            self.start_btn.setEnabled(False)
            self.stop_btn.setEnabled(True)
            self.plot.axes.set_title("Seeing Monitor", color=COLOR_TEXT_LIGHT)
            self.plot.draw()

    def stop_measurement(self):
        if self.running:
            self.running = False
            self.measurement_timer.stop()
            self.status_indicators['power'].set_status("off")
            self.status_indicators['cooling'].set_status("off")
            self.start_btn.setEnabled(True)
            self.stop_btn.setEnabled(False)

    def calibrate(self):
        if not self.connected:
            QMessageBox.warning(self, "Connection Error", "Not connected to telescope")
            return

        QMessageBox.information(self, "Calibration", "Starting MASS-DIMM calibration...")
        self.status_indicators['ao'].set_status("warning")
        QTimer.singleShot(3000, lambda: self.status_indicators['ao'].set_status("on"))

    def update_measurement(self):
        if not self.connected or not hasattr(self, 'last_fwhm'):
            return

        current_fwhm = self.last_fwhm
        r0 = self.last_r0
        tau0 = self.last_tau0
        theta0 = self.last_theta0

        self.seeing_labels['fwhm'].setText(f"{current_fwhm:.2f} arcsec")
        self.seeing_labels['r0'].setText(f"{r0:.2f} m")
        self.seeing_labels['tau0'].setText(f"{tau0:.2f} ms")
        self.seeing_labels['theta0'].setText(f"{theta0:.2f} arcsec")
        self.seeing_labels['seeing'].setText(f"{current_fwhm * 1.2:.2f} arcsec")
        self.seeing_labels['strehl'].setText(f"{0.8 - current_fwhm/10:.2f}")

        self.measurement_times.append(len(self.measurement_times))
        self.measurement_data.append(current_fwhm)

        if len(self.measurement_times) > 60:
            self.measurement_times = self.measurement_times[-60:]
            self.measurement_data = self.measurement_data[-60:]

        self.plot.axes.clear()
        self.plot.axes.plot(self.measurement_times, self.measurement_data, color=COLOR_ACCENT2, linewidth=1.5)
        self.plot.axes.set_facecolor(COLOR_DARK)
        for spine in self.plot.axes.spines.values():
            spine.set_color(COLOR_TEXT_MEDIUM)
            spine.set_linewidth(0.5)
        self.plot.axes.tick_params(axis='both', colors=COLOR_TEXT_MEDIUM)
        self.plot.axes.set_title("Seeing Monitor", color=COLOR_TEXT_LIGHT)
        self.plot.axes.set_xlabel("Time (min)", color=COLOR_TEXT_MEDIUM)
        self.plot.axes.set_ylabel("FWHM (arcsec)", color=COLOR_TEXT_MEDIUM)
        self.plot.axes.grid(True, color=COLOR_LIGHT, linestyle='--')
        self.plot.axes.set_ylim(0, 2.0)
        self.plot.draw()

    def update_data(self, data):
        if data['fwhm'] is not None:
            self.last_fwhm = data['fwhm']
        if data['r0'] is not None:
            self.last_r0 = data['r0']
        if data['tau0'] is not None:
            self.last_tau0 = data['tau0']
        if data['theta0'] is not None:
            self.last_theta0 = data['theta0']

class FilterWheelTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout()

        filter_group = QGroupBox("Filter Selection")
        filter_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        filter_layout = QVBoxLayout()

        self.filter_combo = QComboBox()
        self.filter_combo.addItems(["Clear", "U", "B", "V", "R", "I", "Hα", "OIII", "SII"])

        move_layout = QHBoxLayout()
        move_btn = QPushButton("Move to Filter")
        move_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT1};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT3};
            }}
        """)
        stop_btn = QPushButton("Stop")
        stop_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: #8B0000;
                color: white;
            }}
            QPushButton:hover {{
                background-color: #A52A2A;
            }}
        """)

        move_layout.addWidget(self.filter_combo)
        move_layout.addWidget(move_btn)
        move_layout.addWidget(stop_btn)

        self.current_filter = QLabel("Current Filter: None")
        self.filter_pos = QLabel("Position: 0")

        filter_layout.addLayout(move_layout)
        filter_layout.addWidget(self.current_filter)
        filter_layout.addWidget(self.filter_pos)
        filter_group.setLayout(filter_layout)
        layout.addWidget(filter_group)

        status_group = QGroupBox("Filter Wheel Status")
        status_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        status_layout = QGridLayout()

        self.status_indicators = {
            'power': StatusIndicator("OFF"),
            'moving': StatusIndicator("IDLE"),
            'home': StatusIndicator("UNKNOWN"),
            'error': StatusIndicator("OK")
        }

        for i, (label, widget) in enumerate(self.status_indicators.items()):
            status_layout.addWidget(QLabel(label.capitalize() + ":"), i//2, (i%2)*2)
            status_layout.addWidget(widget, i//2, (i%2)*2+1)

        status_group.setLayout(status_layout)
        layout.addWidget(status_group)

        info_group = QGroupBox("Filter Information")
        info_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        info_layout = QGridLayout()

        self.filter_info = {
            'name': QLabel("Name: -"),
            'type': QLabel("Type: -"),
            'wavelength': QLabel("λ: - nm"),
            'bandwidth': QLabel("Δλ: - nm"),
            'transmission': QLabel("Transmission: -%")
        }

        for i, (key, widget) in enumerate(self.filter_info.items()):
            info_layout.addWidget(widget, i//2, i%2)

        info_group.setLayout(info_layout)
        layout.addWidget(info_group)

        self.setLayout(layout)

class Camera:
    def __init__(self, name):
        self.name = name
        self._is_initialized = False
        self._temperature = 25.0
        self._cooling = False
        self._filter_position = 0
        self._filters = ["Clear", "U", "B", "V", "R", "I", "Hα", "OIII", "SII"]
        self._exposure_progress = 0
        self._is_exposing = False
        self._image_dir = "images"
        self._file_extension = "fits"
        self._serial_number = f"{random.randint(100000, 999999)}"
        os.makedirs(self._image_dir, exist_ok=True)

    @property
    def uid(self):
        return self._serial_number[:6]

    @property
    def is_initialized(self):
        return self._is_initialized

    @property
    def temperature(self):
        return self._temperature

    @property
    def filter(self):
        return self._filters[self._filter_position]

    @property
    def filters(self):
        return self._filters

    def connect(self):
        if self._is_initialized:
            return True

        time.sleep(1)
        self._is_initialized = True
        return True

    def disconnect(self):
        if not self._is_initialized:
            return True

        time.sleep(0.5)
        self._is_initialized = False
        return True

    def set_temperature(self, temp):
        if not self._is_initialized:
            return False

        self._cooling = True
        self._temperature = temp
        return True

    def set_filter(self, position):
        if not self._is_initialized:
            return False

        if position < 0 or position >= len(self._filters):
            return False

        time.sleep(1)
        self._filter_position = position
        return True

    def take_exposure(self, exposure_time, filename=None, dark=False, flat=False):
        if not self._is_initialized:
            return None

        self._is_exposing = True
        self._exposure_progress = 0

        if filename is None:
            filename = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}.{self._file_extension}"

        file_path = os.path.join(self._image_dir, filename)

        exposure_event = Event()
        exposure_thread = Thread(target=self._simulate_exposure, args=(exposure_time, file_path, exposure_event))
        exposure_thread.start()

        return file_path, exposure_event

    def _simulate_exposure(self, exposure_time, file_path, event):
        steps = int(exposure_time * 10)
        for i in range(steps + 1):
            if not self._is_exposing:
                break
            self._exposure_progress = (i / steps) * 100
            time.sleep(0.1)

        if self._is_exposing:
            with open(file_path, 'w') as f:
                f.write("SIMULATED FITS FILE")

        self._is_exposing = False
        self._exposure_progress = 0
        event.set()

    def abort_exposure(self):
        self._is_exposing = False
        return True

    def get_exposure_progress(self):
        return self._exposure_progress

    def get_status(self):
        return {
            'connected': self._is_initialized,
            'temperature': self._temperature,
            'cooling': self._cooling,
            'filter': self.filter,
            'filter_position': self._filter_position,
            'exposing': self._is_exposing,
            'progress': self._exposure_progress,
            'serial': self._serial_number
        }

class CameraTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.camera = None
        self.current_image = None
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout()

        control_group = QGroupBox("Camera Control")
        control_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        control_layout = QGridLayout()

        control_layout.addWidget(QLabel("Camera:"), 0, 0)
        self.camera_combo = QComboBox()
        self.camera_combo.addItems(["Simulator", "Primary", "Guide"])
        control_layout.addWidget(self.camera_combo, 0, 1)

        self.connect_btn = QPushButton("Connect")
        self.connect_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT3};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT1};
            }}
        """)
        self.connect_btn.clicked.connect(self.toggle_connection)
        control_layout.addWidget(self.connect_btn, 0, 2)

        control_layout.addWidget(QLabel("Temperature (°C):"), 1, 0)
        self.temp_input = QLineEdit("25.0")
        self.temp_input.setValidator(QDoubleValidator(-50, 50, 1))
        control_layout.addWidget(self.temp_input, 1, 1)

        self.temp_set_btn = QPushButton("Set")
        self.temp_set_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_MEDIUM};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_LIGHT};
            }}
        """)
        self.temp_set_btn.clicked.connect(self.set_temperature)
        control_layout.addWidget(self.temp_set_btn, 1, 2)

        control_layout.addWidget(QLabel("Filter:"), 2, 0)
        self.filter_combo = QComboBox()
        self.filter_combo.addItems(["Clear", "U", "B", "V", "R", "I", "Hα", "OIII", "SII"])
        control_layout.addWidget(self.filter_combo, 2, 1)

        self.filter_set_btn = QPushButton("Set")
        self.filter_set_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_MEDIUM};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_LIGHT};
            }}
        """)
        self.filter_set_btn.clicked.connect(self.set_filter)
        control_layout.addWidget(self.filter_set_btn, 2, 2)

        control_layout.addWidget(QLabel("Exposure (s):"), 3, 0)
        self.exp_time_input = QLineEdit("5.0")
        self.exp_time_input.setValidator(QDoubleValidator(0.1, 3600, 1))
        control_layout.addWidget(self.exp_time_input, 3, 1)

        control_layout.addWidget(QLabel("Gain:"), 4, 0)
        self.gain_input = QLineEdit("50")
        self.gain_input.setValidator(QDoubleValidator(0, 100, 0))
        control_layout.addWidget(self.gain_input, 4, 1)

        control_layout.addWidget(QLabel("Offset:"), 5, 0)
        self.offset_input = QLineEdit("30")
        self.offset_input.setValidator(QDoubleValidator(0, 100, 0))
        control_layout.addWidget(self.offset_input, 5, 1)

        control_layout.addWidget(QLabel("Filename:"), 6, 0)
        self.filename_input = QLineEdit()
        self.filename_input.setPlaceholderText("auto")
        control_layout.addWidget(self.filename_input, 6, 1, 1, 2)

        self.capture_btn = QPushButton("Capture")
        self.capture_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT1};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT3};
            }}
        """)
        self.capture_btn.clicked.connect(self.capture_image)
        control_layout.addWidget(self.capture_btn, 7, 0)

        self.abort_btn = QPushButton("Abort")
        self.abort_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: #8B0000;
                color: white;
            }}
            QPushButton:hover {{
                background-color: #A52A2A;
            }}
        """)
        self.abort_btn.clicked.connect(self.abort_capture)
        control_layout.addWidget(self.abort_btn, 7, 1)

        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_MEDIUM};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_LIGHT};
            }}
        """)
        self.browse_btn.clicked.connect(self.browse_images)
        control_layout.addWidget(self.browse_btn, 7, 2)

        control_group.setLayout(control_layout)
        layout.addWidget(control_group)

        status_group = QGroupBox("Camera Status")
        status_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        status_layout = QGridLayout()

        self.status_indicators = {
            'power': StatusIndicator("OFF"),
            'cooling': StatusIndicator("OFF"),
            'filter': StatusIndicator("NONE"),
            'exposing': StatusIndicator("IDLE")
        }

        for i, (label, widget) in enumerate(self.status_indicators.items()):
            status_layout.addWidget(QLabel(label.capitalize() + ":"), i//2, (i%2)*2)
            status_layout.addWidget(widget, i//2, (i%2)*2+1)

        self.temp_display = QLabel("Temperature: - °C")
        status_layout.addWidget(self.temp_display, 2, 0, 1, 2)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        status_layout.addWidget(self.progress_bar, 3, 0, 1, 2)

        self.camera_info = QLabel("Camera: Not connected")
        status_layout.addWidget(self.camera_info, 4, 0, 1, 2)

        status_group.setLayout(status_layout)
        layout.addWidget(status_group)

        image_group = QGroupBox("Image Preview")
        image_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        image_layout = QVBoxLayout()

        self.image_display = QLabel("No image available")
        self.image_display.setAlignment(Qt.AlignCenter)
        self.image_display.setStyleSheet(f"""
            QLabel {{
                border: 1px solid {COLOR_LIGHT};
                background-color: {COLOR_DARK};
                min-height: 300px;
            }}
        """)
        image_layout.addWidget(self.image_display)

        self.image_info = QLabel("Image info: None loaded")
        image_layout.addWidget(self.image_info)

        image_group.setLayout(image_layout)
        layout.addWidget(image_group)

        self.setLayout(layout)
        self.set_controls_enabled(False)

        self.update_timer = QTimer(self)
        self.update_timer.timeout.connect(self.update_status)
        self.update_timer.start(500)

    def set_controls_enabled(self, enabled):
        self.temp_input.setEnabled(enabled)
        self.temp_set_btn.setEnabled(enabled)
        self.filter_combo.setEnabled(enabled)
        self.filter_set_btn.setEnabled(enabled)
        self.exp_time_input.setEnabled(enabled)
        self.gain_input.setEnabled(enabled)
        self.offset_input.setEnabled(enabled)
        self.filename_input.setEnabled(enabled)
        self.capture_btn.setEnabled(enabled)
        self.abort_btn.setEnabled(enabled)
        self.browse_btn.setEnabled(enabled)

    def toggle_connection(self):
        if self.camera and self.camera.is_initialized:
            self.disconnect_camera()
        else:
            self.connect_camera()

    def connect_camera(self):
        camera_name = self.camera_combo.currentText()
        try:
            self.camera = Camera(camera_name)
            if self.camera.connect():
                self.connect_btn.setText("Disconnect")
                self.set_controls_enabled(True)
                self.status_indicators['power'].set_status("on")
                self.camera_info.setText(f"Camera: {camera_name} (SN: {self.camera.uid})")
                QMessageBox.information(self, "Connected", f"Successfully connected to {camera_name} camera")
            else:
                self.camera = None
                QMessageBox.warning(self, "Error", "Failed to connect to camera")
        except Exception as e:
            self.camera = None
            QMessageBox.critical(self, "Error", f"Connection error: {str(e)}")

    def disconnect_camera(self):
        try:
            if self.camera and self.camera.disconnect():
                self.connect_btn.setText("Connect")
                self.set_controls_enabled(False)
                self.status_indicators['power'].set_status("off")
                self.status_indicators['cooling'].set_status("off")
                self.status_indicators['filter'].set_status("none")
                self.status_indicators['exposing'].set_status("idle")
                self.temp_display.setText("Temperature: - °C")
                self.progress_bar.setValue(0)
                self.camera_info.setText("Camera: Not connected")
                QMessageBox.information(self, "Disconnected", "Camera disconnected successfully")
            else:
                QMessageBox.warning(self, "Error", "Failed to disconnect camera")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Disconnection error: {str(e)}")
        finally:
            self.camera = None

    def set_temperature(self):
        if not self.camera:
            return

        try:
            temp = float(self.temp_input.text())
            if self.camera.set_temperature(temp):
                self.status_indicators['cooling'].set_status("on")
                QMessageBox.information(self, "Temperature Set", f"Cooling set to {temp}°C")
            else:
                QMessageBox.warning(self, "Error", "Failed to set temperature")
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid temperature")

    def set_filter(self):
        if not self.camera:
            return

        filter_pos = self.filter_combo.currentIndex()
        if self.camera.set_filter(filter_pos):
            self.status_indicators['filter'].set_status(self.camera.filter.lower())
            QMessageBox.information(self, "Filter Changed", f"Filter set to {self.camera.filter}")
        else:
            QMessageBox.warning(self, "Error", "Failed to change filter")

    def capture_image(self):
        if not self.camera:
            return

        try:
            exp_time = float(self.exp_time_input.text())
            filename = self.filename_input.text() or None

            file_path, exposure_event = self.camera.take_exposure(exp_time, filename)

            self.status_indicators['exposing'].set_status("exposing")
            self.exposure_event = exposure_event
            QTimer.singleShot(100, self.check_exposure)

        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter valid exposure settings")

    def check_exposure(self):
        if not self.camera:
            return

        progress = self.camera.get_exposure_progress()
        self.progress_bar.setValue(progress)

        if progress < 100:
            QTimer.singleShot(100, self.check_exposure)
        else:
            self.exposure_complete()

    def exposure_complete(self):
        self.status_indicators['exposing'].set_status("idle")
        self.progress_bar.setValue(0)

        self.current_image = "sample_image.jpg"
        pixmap = QPixmap("assets/camera_icon.png")
        self.image_display.setPixmap(pixmap.scaled(
            self.image_display.width(),
            self.image_display.height(),
            Qt.KeepAspectRatio))

        self.image_info.setText(f"Image captured: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        QMessageBox.information(self, "Capture Complete", "Image capture finished successfully")

    def abort_capture(self):
        if self.camera and self.camera.abort_exposure():
            self.status_indicators['exposing'].set_status("idle")
            self.progress_bar.setValue(0)
            QMessageBox.information(self, "Capture Aborted", "Exposure aborted")
        else:
            QMessageBox.warning(self, "Error", "No exposure to abort or failed to abort")

    def browse_images(self):
        if not self.camera:
            return

        options = QFileDialog.Options()
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open Image", self.camera._image_dir,
            "FITS Files (*.fits *.fit);;All Files (*)", options=options)

        if filename:
            pixmap = QPixmap("assets/camera_icon.png")
            self.image_display.setPixmap(pixmap.scaled(
                self.image_display.width(),
                self.image_display.height(),
                Qt.KeepAspectRatio))

            self.image_info.setText(f"Image loaded: {os.path.basename(filename)}")
            self.current_image = filename

    def update_status(self):
        if not self.camera:
            return

        status = self.camera.get_status()

        self.temp_display.setText(f"Temperature: {status['temperature']:.1f} °C")

        current_filter = self.filter_combo.currentText()
        if current_filter != status['filter']:
            self.filter_combo.setCurrentText(status['filter'])
            self.status_indicators['filter'].set_status(status['filter'].lower())

        cooling_status = "on" if status['cooling'] else "off"
        self.status_indicators['cooling'].set_status(cooling_status)

        exposing_status = "exposing" if status['exposing'] else "idle"
        self.status_indicators['exposing'].set_status(exposing_status)

class CatalogTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.location = parent.location if hasattr(parent, 'location') else None
        self.main_window = self._find_main_window()  # Fixed this line
        self.setup_ui()
        self.load_catalog_data()
    def _find_main_window(self):
        parent = self.parent()
        while parent is not None:
            if isinstance(parent, TelescopeControlUI):
                return parent
            parent = parent.parent()
        return None

    def setup_ui(self):
        layout = QVBoxLayout()

        # Add splitter to show both catalog and sky viewer
        splitter = QSplitter(Qt.Vertical)

        # Top part - catalog search and table
        catalog_widget = QWidget()
        catalog_layout = QVBoxLayout()

        search_group = QGroupBox("Catalog Search")
        search_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        search_layout = QHBoxLayout()

        self.catalog_combo = QComboBox()
        self.catalog_combo.addItems(["Messier", "NGC", "IC", "SAO", "Bright Stars"])
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search by name or ID...")
        search_btn = QPushButton("Search")
        search_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT3};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT1};
            }}
        """)
        goto_btn = QPushButton("GoTo Selected")
        goto_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT1};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT3};
            }}
        """)
        goto_btn.clicked.connect(self.goto_selected_object)

        search_layout.addWidget(self.catalog_combo)
        search_layout.addWidget(self.search_input)
        search_layout.addWidget(search_btn)
        search_layout.addWidget(goto_btn)

        search_group.setLayout(search_layout)
        catalog_layout.addWidget(search_group)

        self.catalog_table = QTableWidget()
        self.catalog_table.setColumnCount(6)
        self.catalog_table.setHorizontalHeaderLabels(["ID", "Name", "Type", "RA", "DEC", "Magnitude"])
        self.catalog_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.catalog_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.catalog_table.setSelectionMode(QTableWidget.SingleSelection)
        self.catalog_table.setStyleSheet(f"""
            QTableWidget {{
                border: 1px solid {COLOR_LIGHT};
                gridline-color: {COLOR_LIGHT};
                color: {COLOR_TEXT_LIGHT};
            }}
            QHeaderView::section {{
                background-color: {COLOR_MEDIUM};
                color: {COLOR_TEXT_LIGHT};
                padding: 5px;
                border: none;
            }}
        """)
        catalog_layout.addWidget(self.catalog_table)

        details_group = QGroupBox("Object Details")
        details_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        details_layout = QGridLayout()
        self.details_labels = {
            'id': QLabel("ID: -"),
            'name': QLabel("Name: -"),
            'type': QLabel("Type: -"),
            'ra': QLabel("RA: -"),
            'dec': QLabel("DEC: -"),
            'mag': QLabel("Magnitude: -"),
            'size': QLabel("Size: -"),
            'const': QLabel("Constellation: -"),
            'alt': QLabel("Altitude: -"),
            'az': QLabel("Azimuth: -")
        }

        for i, (key, widget) in enumerate(self.details_labels.items()):
            details_layout.addWidget(widget, i//3, i%3)

        details_group.setLayout(details_layout)
        catalog_layout.addWidget(details_group)

        catalog_widget.setLayout(catalog_layout)

        # Bottom part - sky viewer
        self.sky_viewer = SkyViewer(self)

        splitter.addWidget(catalog_widget)
        splitter.addWidget(self.sky_viewer)
        splitter.setSizes([400, 400])

        layout.addWidget(splitter)
        self.setLayout(layout)

        search_btn.clicked.connect(self.search_catalog)
        self.catalog_table.itemSelectionChanged.connect(self.show_object_details)

        self.setLayout(layout)

    def load_catalog_data(self):
        self.catalogs = {
            "Messier": [
                               {"id": "M1", "name": "Crab Nebula", "type": "Supernova Remnant",
                                "ra": "05:34:31.94", "dec": "+22:00:52.2", "mag": 8.4, "size": "6x4 arcmin", "const": "Taurus"},
                               {"id": "M2", "name": "", "type": "Globular Cluster",
                                "ra": "21:33:27", "dec": "-00:49:24", "mag": 6.5, "size": "16 arcmin", "const": "Aquarius"},
                               {"id": "M3", "name": "", "type": "Globular Cluster",
                                "ra": "13:42:11", "dec": "+28:22:38", "mag": 6.2, "size": "18 arcmin", "const": "Canes Venatici"},
                               {"id": "M4", "name": "", "type": "Globular Cluster",
                                "ra": "16:23:35", "dec": "-26:31:33", "mag": 5.6, "size": "26 arcmin", "const": "Scorpius"},
                               {"id": "M5", "name": "", "type": "Globular Cluster",
                                "ra": "15:18:33", "dec": "+02:04:58", "mag": 5.6, "size": "17 arcmin", "const": "Serpens"},
                               {"id": "M6", "name": "Butterfly Cluster", "type": "Open Cluster",
                                "ra": "17:40:20", "dec": "-32:15:15", "mag": 4.2, "size": "25 arcmin", "const": "Scorpius"},
                               {"id": "M7", "name": "Ptolemy's Cluster", "type": "Open Cluster",
                                "ra": "17:53:51", "dec": "-34:47:34", "mag": 3.3, "size": "80 arcmin", "const": "Scorpius"},
                               {"id": "M8", "name": "Lagoon Nebula", "type": "Diffuse Nebula",
                                "ra": "18:03:37", "dec": "-24:23:12", "mag": 6.0, "size": "90x40 arcmin", "const": "Sagittarius"},
                               {"id": "M9", "name": "", "type": "Globular Cluster",
                                "ra": "17:19:11", "dec": "-18:30:59", "mag": 7.7, "size": "9 arcmin", "const": "Ophiuchus"},
                               {"id": "M10", "name": "", "type": "Globular Cluster",
                                "ra": "16:57:09", "dec": "-04:05:58", "mag": 6.4, "size": "15 arcmin", "const": "Ophiuchus"},
                               {"id": "M11", "name": "Wild Duck Cluster", "type": "Open Cluster",
                                "ra": "18:51:06", "dec": "-06:16:00", "mag": 5.8, "size": "14 arcmin", "const": "Scutum"},
                               {"id": "M12", "name": "", "type": "Globular Cluster",
                                "ra": "16:47:14", "dec": "-01:56:52", "mag": 6.7, "size": "16 arcmin", "const": "Ophiuchus"},
                               {"id": "M13", "name": "Great Hercules Cluster", "type": "Globular Cluster",
                                "ra": "16:41:41", "dec": "+36:27:35", "mag": 5.8, "size": "20 arcmin", "const": "Hercules"},
                               {"id": "M14", "name": "", "type": "Globular Cluster",
                                "ra": "17:37:36", "dec": "-03:14:45", "mag": 7.6, "size": "11 arcmin", "const": "Ophiuchus"},
                               {"id": "M15", "name": "", "type": "Globular Cluster",
                                "ra": "21:29:58", "dec": "+12:10:01", "mag": 6.2, "size": "18 arcmin", "const": "Pegasus"},
                               {"id": "M16", "name": "Eagle Nebula", "type": "Diffuse Nebula",
                                "ra": "18:18:48", "dec": "-13:49:00", "mag": 6.0, "size": "35x28 arcmin", "const": "Serpens"},
                               {"id": "M17", "name": "Omega Nebula", "type": "Diffuse Nebula",
                                "ra": "18:20:26", "dec": "-16:10:36", "mag": 6.0, "size": "20x15 arcmin", "const": "Sagittarius"},
                               {"id": "M18", "name": "", "type": "Open Cluster",
                                "ra": "18:19:54", "dec": "-17:06:00", "mag": 6.9, "size": "9 arcmin", "const": "Sagittarius"},
                               {"id": "M19", "name": "", "type": "Globular Cluster",
                                "ra": "17:02:37", "dec": "-26:16:05", "mag": 6.8, "size": "13 arcmin", "const": "Ophiuchus"},
                               {"id": "M20", "name": "Trifid Nebula", "type": "Diffuse Nebula",
                                "ra": "18:02:23", "dec": "-23:01:48", "mag": 6.3, "size": "28 arcmin", "const": "Sagittarius"},
                               {"id": "M21", "name": "", "type": "Open Cluster",
                                "ra": "18:04:36", "dec": "-22:29:00", "mag": 5.9, "size": "13 arcmin", "const": "Sagittarius"},
                               {"id": "M22", "name": "", "type": "Globular Cluster",
                                "ra": "18:36:24", "dec": "-23:54:12", "mag": 5.1, "size": "24 arcmin", "const": "Sagittarius"},
                               {"id": "M23", "name": "", "type": "Open Cluster",
                                "ra": "17:56:48", "dec": "-19:01:00", "mag": 5.5, "size": "27 arcmin", "const": "Sagittarius"},
                               {"id": "M24", "name": "Sagittarius Star Cloud", "type": "Star Cloud",
                                "ra": "18:17:00", "dec": "-18:29:00", "mag": 4.6, "size": "90 arcmin", "const": "Sagittarius"},
                               {"id": "M25", "name": "", "type": "Open Cluster",
                                "ra": "18:31:42", "dec": "-19:07:00", "mag": 4.6, "size": "32 arcmin", "const": "Sagittarius"},
                               {"id": "M26", "name": "", "type": "Open Cluster",
                                "ra": "18:45:18", "dec": "-09:23:00", "mag": 8.0, "size": "15 arcmin", "const": "Scutum"},
                               {"id": "M27", "name": "Dumbbell Nebula", "type": "Planetary Nebula",
                                "ra": "19:59:36", "dec": "+22:43:00", "mag": 7.4, "size": "8x5.7 arcmin", "const": "Vulpecula"},
                               {"id": "M28", "name": "", "type": "Globular Cluster",
                                "ra": "18:24:32", "dec": "-24:52:11", "mag": 6.8, "size": "11 arcmin", "const": "Sagittarius"},
                               {"id": "M29", "name": "", "type": "Open Cluster",
                                "ra": "20:23:56", "dec": "+38:31:24", "mag": 6.6, "size": "7 arcmin", "const": "Cygnus"},
                               {"id": "M30", "name": "", "type": "Globular Cluster",
                                "ra": "21:40:22", "dec": "-23:10:48", "mag": 7.2, "size": "12 arcmin", "const": "Capricornus"},
                               {"id": "M31", "name": "Andromeda Galaxy", "type": "Galaxy",
                                "ra": "00:42:44.3", "dec": "+41:16:09", "mag": 3.4, "size": "190x60 arcmin", "const": "Andromeda"},
                               {"id": "M32", "name": "", "type": "Dwarf Galaxy",
                                "ra": "00:42:41.8", "dec": "+40:51:55", "mag": 8.1, "size": "8x6 arcmin", "const": "Andromeda"},
                               {"id": "M33", "name": "Triangulum Galaxy", "type": "Galaxy",
                                "ra": "01:33:50.9", "dec": "+30:39:37", "mag": 5.7, "size": "62x39 arcmin", "const": "Triangulum"},
                               {"id": "M34", "name": "", "type": "Open Cluster",
                                "ra": "02:42:06", "dec": "+42:45:00", "mag": 5.2, "size": "35 arcmin", "const": "Perseus"},
                               {"id": "M35", "name": "", "type": "Open Cluster",
                                "ra": "06:08:54", "dec": "+24:20:00", "mag": 5.1, "size": "28 arcmin", "const": "Gemini"},
                               {"id": "M36", "name": "", "type": "Open Cluster",
                                "ra": "05:36:12", "dec": "+34:08:24", "mag": 6.0, "size": "12 arcmin", "const": "Auriga"},
                               {"id": "M37", "name": "", "type": "Open Cluster",
                                "ra": "05:52:18", "dec": "+32:33:00", "mag": 5.6, "size": "24 arcmin", "const": "Auriga"},
                               {"id": "M38", "name": "", "type": "Open Cluster",
                                "ra": "05:28:42", "dec": "+35:51:00", "mag": 6.4, "size": "21 arcmin", "const": "Auriga"},
                               {"id": "M39", "name": "", "type": "Open Cluster",
                                "ra": "21:32:12", "dec": "+48:26:00", "mag": 4.6, "size": "32 arcmin", "const": "Cygnus"},
                               {"id": "M40", "name": "Winnecke 4", "type": "Double Star",
                                "ra": "12:22:12", "dec": "+58:04:59", "mag": 8.4, "size": "-", "const": "Ursa Major"},
                               {"id": "M41", "name": "", "type": "Open Cluster",
                                "ra": "06:46:00", "dec": "-20:44:00", "mag": 4.5, "size": "38 arcmin", "const": "Canis Major"},
                               {"id": "M42", "name": "Orion Nebula", "type": "Diffuse Nebula",
                                "ra": "05:35:17.3", "dec": "-05:23:28", "mag": 4.0, "size": "85x60 arcmin", "const": "Orion"},
                               {"id": "M43", "name": "De Mairan's Nebula", "type": "Diffuse Nebula",
                                "ra": "05:35:31", "dec": "-05:16:12", "mag": 9.0, "size": "20x15 arcmin", "const": "Orion"},
                               {"id": "M44", "name": "Beehive Cluster", "type": "Open Cluster",
                                "ra": "08:40:24", "dec": "+19:41:00", "mag": 3.1, "size": "95 arcmin", "const": "Cancer"},
                               {"id": "M45", "name": "Pleiades", "type": "Open Cluster",
                                "ra": "03:47:24", "dec": "+24:07:00", "mag": 1.6, "size": "110 arcmin", "const": "Taurus"},
                               {"id": "M46", "name": "", "type": "Open Cluster",
                                "ra": "07:41:42", "dec": "-14:48:36", "mag": 6.1, "size": "27 arcmin", "const": "Puppis"},
                               {"id": "M47", "name": "", "type": "Open Cluster",
                                "ra": "07:36:36", "dec": "-14:28:48", "mag": 4.2, "size": "30 arcmin", "const": "Puppis"},
                               {"id": "M48", "name": "", "type": "Open Cluster",
                                "ra": "08:13:42", "dec": "-05:45:00", "mag": 5.8, "size": "54 arcmin", "const": "Hydra"},
                               {"id": "M49", "name": "", "type": "Galaxy",
                                "ra": "12:29:46.7", "dec": "+08:00:02", "mag": 8.4, "size": "9x7.5 arcmin", "const": "Virgo"},
                               {"id": "M50", "name": "", "type": "Open Cluster",
                                "ra": "07:03:00", "dec": "-08:20:00", "mag": 5.9, "size": "16 arcmin", "const": "Monoceros"},
                               {"id": "M51", "name": "Whirlpool Galaxy", "type": "Galaxy",
                                "ra": "13:29:52.7", "dec": "+47:11:43", "mag": 8.4, "size": "11x7 arcmin", "const": "Canes Venatici"},
                               {"id": "M52", "name": "", "type": "Open Cluster",
                                "ra": "23:24:48", "dec": "+61:35:00", "mag": 6.9, "size": "13 arcmin", "const": "Cassiopeia"},
                               {"id": "M53", "name": "", "type": "Globular Cluster",
                                "ra": "13:12:55", "dec": "+18:10:09", "mag": 7.6, "size": "13 arcmin", "const": "Coma Berenices"},
                               {"id": "M54", "name": "", "type": "Globular Cluster",
                                "ra": "18:55:03", "dec": "-30:28:42", "mag": 7.6, "size": "9 arcmin", "const": "Sagittarius"},
                               {"id": "M55", "name": "", "type": "Globular Cluster",
                                "ra": "19:39:59", "dec": "-30:57:44", "mag": 6.3, "size": "19 arcmin", "const": "Sagittarius"},
                               {"id": "M56", "name": "", "type": "Globular Cluster",
                                "ra": "19:16:35", "dec": "+30:11:04", "mag": 8.3, "size": "7 arcmin", "const": "Lyra"},
                               {"id": "M57", "name": "Ring Nebula", "type": "Planetary Nebula",
                                "ra": "18:53:35.08", "dec": "+33:01:45.0", "mag": 8.8, "size": "2.5x2.0 arcmin", "const": "Lyra"},
                               {"id": "M58", "name": "", "type": "Galaxy",
                                "ra": "12:37:43.5", "dec": "+11:49:05", "mag": 9.7, "size": "5.5x4.5 arcmin", "const": "Virgo"},
                               {"id": "M59", "name": "", "type": "Galaxy",
                                "ra": "12:42:02.3", "dec": "+11:38:49", "mag": 9.6, "size": "5x3.5 arcmin", "const": "Virgo"},
                               {"id": "M60", "name": "", "type": "Galaxy",
                                "ra": "12:43:39.6", "dec": "+11:33:09", "mag": 8.8, "size": "7x6 arcmin", "const": "Virgo"},
                               {"id": "M61", "name": "", "type": "Galaxy",
                                "ra": "12:21:54.9", "dec": "+04:28:25", "mag": 9.7, "size": "6x5.5 arcmin", "const": "Virgo"},
                               {"id": "M62", "name": "", "type": "Globular Cluster",
                                "ra": "17:01:12", "dec": "-30:06:44", "mag": 6.5, "size": "14 arcmin", "const": "Ophiuchus"},
                               {"id": "M63", "name": "Sunflower Galaxy", "type": "Galaxy",
                                "ra": "13:15:49.3", "dec": "+42:01:45", "mag": 8.6, "size": "12x7.5 arcmin", "const": "Canes Venatici"},
                               {"id": "M64", "name": "Black Eye Galaxy", "type": "Galaxy",
                                "ra": "12:56:43.7", "dec": "+21:40:58", "mag": 8.5, "size": "10x5 arcmin", "const": "Coma Berenices"},
                               {"id": "M65", "name": "", "type": "Galaxy",
                                "ra": "11:18:55.9", "dec": "+13:05:32", "mag": 9.3, "size": "8x1.5 arcmin", "const": "Leo"},
                               {"id": "M66", "name": "", "type": "Galaxy",
                                "ra": "11:20:15.0", "dec": "+12:59:30", "mag": 8.9, "size": "9x4 arcmin", "const": "Leo"},
                               {"id": "M67", "name": "", "type": "Open Cluster",
                                "ra": "08:50:24", "dec": "+11:49:00", "mag": 6.1, "size": "30 arcmin", "const": "Cancer"},
                               {"id": "M68", "name": "", "type": "Globular Cluster",
                                "ra": "12:39:27", "dec": "-26:44:38", "mag": 7.8, "size": "11 arcmin", "const": "Hydra"},
                               {"id": "M69", "name": "", "type": "Globular Cluster",
                                "ra": "18:31:23", "dec": "-32:20:53", "mag": 7.6, "size": "7 arcmin", "const": "Sagittarius"},
                               {"id": "M70", "name": "", "type": "Globular Cluster",
                                "ra": "18:43:12", "dec": "-32:17:31", "mag": 7.9, "size": "8 arcmin", "const": "Sagittarius"},
                               {"id": "M71", "name": "", "type": "Globular Cluster",
                                "ra": "19:53:46", "dec": "+18:46:42", "mag": 8.2, "size": "7 arcmin", "const": "Sagitta"},
                               {"id": "M72", "name": "", "type": "Globular Cluster",
                                "ra": "20:53:27", "dec": "-12:32:14", "mag": 9.3, "size": "6 arcmin", "const": "Aquarius"},
                               {"id": "M73", "name": "", "type": "Asterism",
                                "ra": "20:58:54", "dec": "-12:38:00", "mag": 8.9, "size": "2.8 arcmin", "const": "Aquarius"},
                               {"id": "M74", "name": "", "type": "Galaxy",
                                "ra": "01:36:41.8", "dec": "+15:47:01", "mag": 9.4, "size": "10x9 arcmin", "const": "Pisces"},
                               {"id": "M75", "name": "", "type": "Globular Cluster",
                                "ra": "20:06:04", "dec": "-21:55:16", "mag": 8.5, "size": "6 arcmin", "const": "Sagittarius"},
                               {"id": "M76", "name": "Little Dumbbell Nebula", "type": "Planetary Nebula",
                                "ra": "01:42:19.9", "dec": "+51:34:31", "mag": 10.1, "size": "2.7x1.8 arcmin", "const": "Perseus"},
                               {"id": "M77", "name": "Cetus A", "type": "Galaxy",
                                "ra": "02:42:40.7", "dec": "-00:00:48", "mag": 8.9, "size": "7x6 arcmin", "const": "Cetus"},
                               {"id": "M78", "name": "", "type": "Diffuse Nebula",
                                "ra": "05:46:46", "dec": "+00:00:50", "mag": 8.3, "size": "8x6 arcmin", "const": "Orion"},
                               {"id": "M79", "name": "", "type": "Globular Cluster",
                                "ra": "05:24:10", "dec": "-24:31:27", "mag": 7.7, "size": "8 arcmin", "const": "Lepus"},
                               {"id": "M80", "name": "", "type": "Globular Cluster",
                                "ra": "16:17:02", "dec": "-22:58:30", "mag": 7.3, "size": "10 arcmin", "const": "Scorpius"},
                               {"id": "M81", "name": "Bode's Galaxy", "type": "Galaxy",
                                "ra": "09:55:33", "dec": "+69:03:55", "mag": 6.9, "size": "26x14 arcmin", "const": "Ursa Major"},
                               {"id": "M82", "name": "Cigar Galaxy", "type": "Galaxy",
                                "ra": "09:55:52", "dec": "+69:40:47", "mag": 8.4, "size": "11x4 arcmin", "const": "Ursa Major"},
                               {"id": "M83", "name": "Southern Pinwheel", "type": "Galaxy",
                                "ra": "13:37:00", "dec": "-29:51:56", "mag": 7.5, "size": "13x12 arcmin", "const": "Hydra"},
                               {"id": "M84", "name": "", "type": "Galaxy",
                                "ra": "12:25:03.7", "dec": "+12:53:13", "mag": 9.1, "size": "6x5 arcmin", "const": "Virgo"},
                               {"id": "M85", "name": "", "type": "Galaxy",
                                "ra": "12:25:24.0", "dec": "+18:11:28", "mag": 9.1, "size": "7x5 arcmin", "const": "Coma Berenices"},
                               {"id": "M86", "name": "", "type": "Galaxy",
                                "ra": "12:26:11.7", "dec": "+12:56:46", "mag": 8.9, "size": "8x5 arcmin", "const": "Virgo"},
                               {"id": "M87", "name": "Virgo A", "type": "Galaxy",
                                "ra": "12:30:49.4", "dec": "+12:23:28", "mag": 8.6, "size": "7 arcmin", "const": "Virgo"},
                               {"id": "M88", "name": "", "type": "Galaxy",
                                "ra": "12:31:59.2", "dec": "+14:25:14", "mag": 9.6, "size": "7x4 arcmin", "const": "Coma Berenices"},
                               {"id": "M89", "name": "", "type": "Galaxy",
                                "ra": "12:35:39.8", "dec": "+12:33:23", "mag": 9.8, "size": "5x5 arcmin", "const": "Virgo"},
                               {"id": "M90", "name": "", "type": "Galaxy",
                                "ra": "12:36:49.8", "dec": "+13:09:46", "mag": 9.5, "size": "10x4 arcmin", "const": "Virgo"},
                               {"id": "M91", "name": "", "type": "Galaxy",
                                "ra": "12:35:26.4", "dec": "+14:29:47", "mag": 10.2, "size": "5x4 arcmin", "const": "Coma Berenices"},
                               {"id": "M92", "name": "", "type": "Globular Cluster",
                                "ra": "17:17:07", "dec": "+43:08:09", "mag": 6.4, "size": "14 arcmin", "const": "Hercules"},
                               {"id": "M93", "name": "", "type": "Open Cluster",
                                "ra": "07:44:30", "dec": "-23:51:24", "mag": 6.0, "size": "22 arcmin", "const": "Puppis"},
                               {"id": "M94", "name": "Cat's Eye Galaxy", "type": "Galaxy",
                                "ra": "12:50:53.1", "dec": "+41:07:14", "mag": 8.2, "size": "11x9 arcmin", "const": "Canes Venatici"},
                               {"id": "M95", "name": "", "type": "Galaxy",
                                "ra": "10:43:57.7", "dec": "+11:42:14", "mag": 9.7, "size": "7x5 arcmin", "const": "Leo"},
                               {"id": "M96", "name": "", "type": "Galaxy",
                                "ra": "10:46:45.7", "dec": "+11:49:12", "mag": 9.2, "size": "7x5 arcmin", "const": "Leo"},
                               {"id": "M97", "name": "Owl Nebula", "type": "Planetary Nebula",
                                "ra": "11:14:47.7", "dec": "+55:01:08", "mag": 9.9, "size": "3.4x3.3 arcmin", "const": "Ursa Major"},
                               {"id": "M98", "name": "", "type": "Galaxy",
                                "ra": "12:13:48.3", "dec": "+14:54:01", "mag": 10.1, "size": "9x3 arcmin", "const": "Coma Berenices"},
                               {"id": "M99", "name": "", "type": "Galaxy",
                                "ra": "12:18:49.6", "dec": "+14:25:00", "mag": 9.9, "size": "5x5 arcmin", "const": "Coma Berenices"},
                               {"id": "M100", "name": "", "type": "Galaxy",
                                "ra": "12:22:54.9", "dec": "+15:49:21", "mag": 9.3, "size": "7x6 arcmin", "const": "Coma Berenices"},
                               {"id": "M101", "name": "Pinwheel Galaxy", "type": "Galaxy",
                                "ra": "14:03:12.6", "dec": "+54:20:57", "mag": 7.9, "size": "29x27 arcmin", "const": "Ursa Major"},
                               {"id": "M102", "name": "Spindle Galaxy", "type": "Galaxy",
                                "ra": "15:06:29.5", "dec": "+55:45:48", "mag": 9.9, "size": "5x2 arcmin", "const": "Draco"},
                               {"id": "M103", "name": "", "type": "Open Cluster",
                                "ra": "01:33:24", "dec": "+60:39:00", "mag": 7.4, "size": "6 arcmin", "const": "Cassiopeia"},
                               {"id": "M104", "name": "Sombrero Galaxy", "type": "Galaxy",
                                "ra": "12:39:59", "dec": "-11:37:23", "mag": 8.0, "size": "8.7x3.5 arcmin", "const": "Virgo"},
                               {"id": "M105", "name": "", "type": "Galaxy",
                                "ra": "10:47:49.6", "dec": "+12:34:54", "mag": 9.3, "size": "5x4 arcmin", "const": "Leo"},
                               {"id": "M106", "name": "", "type": "Galaxy",
                                "ra": "12:18:57.5", "dec": "+47:18:14", "mag": 8.4, "size": "19x7 arcmin", "const": "Canes Venatici"}
                   ],
            "NGC": [
                                   {"id": "NGC7000", "name": "North America Nebula", "type": "Diffuse Nebula",
                                    "ra": "20:58:47", "dec": "+44:19:48", "mag": 4.0, "size": "120x100 arcmin", "const": "Cygnus"},
                                   {"id": "NGC1976", "name": "Orion Nebula", "type": "Diffuse Nebula",
                                    "ra": "05:35:17.3", "dec": "-05:23:28", "mag": 4.0, "size": "85x60 arcmin", "const": "Orion"},
                                   {"id": "NGC7293", "name": "Helix Nebula", "type": "Planetary Nebula",
                                    "ra": "22:29:38", "dec": "-20:50:13", "mag": 7.3, "size": "16 arcmin", "const": "Aquarius"},
                                   {"id": "NGC253", "name": "Sculptor Galaxy", "type": "Galaxy",
                                    "ra": "00:47:33", "dec": "-25:17:18", "mag": 7.1, "size": "27x7 arcmin", "const": "Sculptor"},
                                   {"id": "NGC5128", "name": "Centaurus A", "type": "Galaxy",
                                    "ra": "13:25:27", "dec": "-43:01:09", "mag": 6.6, "size": "25x20 arcmin", "const": "Centaurus"},
                                   {"id": "NGC3628", "name": "Hamburger Galaxy", "type": "Galaxy",
                                    "ra": "11:20:17", "dec": "+13:35:23", "mag": 9.5, "size": "15x3 arcmin", "const": "Leo"},
                                   {"id": "NGC4565", "name": "Needle Galaxy", "type": "Galaxy",
                                    "ra": "12:36:21", "dec": "+25:59:14", "mag": 9.6, "size": "15x2 arcmin", "const": "Coma Berenices"},
                                   {"id": "NGC6543", "name": "Cat's Eye Nebula", "type": "Planetary Nebula",
                                    "ra": "17:58:33", "dec": "+66:38:00", "mag": 8.1, "size": "20 arcsec", "const": "Draco"},
                                   {"id": "NGC6960", "name": "Western Veil Nebula", "type": "Supernova Remnant",
                                    "ra": "20:45:38", "dec": "+30:42:30", "mag": 7.0, "size": "70x6 arcmin", "const": "Cygnus"},
                                   {"id": "NGC6992", "name": "Eastern Veil Nebula", "type": "Supernova Remnant",
                                    "ra": "20:56:19", "dec": "+31:43:00", "mag": 7.0, "size": "70x6 arcmin", "const": "Cygnus"}
                               ],
            "IC": [
                                   {"id": "IC434", "name": "Horsehead Nebula", "type": "Dark Nebula",
                                    "ra": "05:40:59", "dec": "-02:27:30", "mag": 6.8, "size": "8x6 arcmin", "const": "Orion"},
                                   {"id": "IC2118", "name": "Witch Head Nebula", "type": "Reflection Nebula",
                                    "ra": "05:02:00", "dec": "-07:54:00", "mag": 13.0, "size": "180x60 arcmin", "const": "Eridanus"},
                                   {"id": "IC405", "name": "Flaming Star Nebula", "type": "Emission Nebula",
                                    "ra": "05:16:05", "dec": "+34:27:49", "mag": 6.0, "size": "37x19 arcmin", "const": "Auriga"},
                                   {"id": "IC443", "name": "Jellyfish Nebula", "type": "Supernova Remnant",
                                    "ra": "06:17:00", "dec": "+22:21:00", "mag": 12.0, "size": "50 arcmin", "const": "Gemini"},
                                   {"id": "IC1396", "name": "Elephant Trunk Nebula", "type": "Emission Nebula",
                                    "ra": "21:39:00", "dec": "+57:30:00", "mag": 3.5, "size": "170 arcmin", "const": "Cepheus"}
                               ],
            "SAO": [
                                   {"id": "SAO113271", "name": "Polaris", "type": "Star",
                                    "ra": "02:31:49", "dec": "+89:15:51", "mag": 1.97, "size": "-", "const": "Ursa Minor"},
                                   {"id": "SAO67174", "name": "Vega", "type": "Star",
                                    "ra": "18:36:56", "dec": "+38:47:01", "mag": 0.03, "size": "-", "const": "Lyra"},
                                   {"id": "SAO28737", "name": "Capella", "type": "Star",
                                    "ra": "05:16:41", "dec": "+45:59:53", "mag": 0.08, "size": "-", "const": "Auriga"},
                                   {"id": "SAO79294", "name": "Arcturus", "type": "Star",
                                    "ra": "14:15:40", "dec": "+19:10:56", "mag": -0.05, "size": "-", "const": "Bootes"},
                                   {"id": "SAO159459", "name": "Rigel Kentaurus", "type": "Star",
                                    "ra": "14:39:36", "dec": "-60:50:02", "mag": -0.01, "size": "-", "const": "Centaurus"}
                               ],
            "Bright Stars": [
                                   {"id": "Alpha CMa", "name": "Sirius", "type": "Star",
                                    "ra": "06:45:08.9", "dec": "-16:42:58", "mag": -1.46, "size": "-", "const": "Canis Major"},
                                   {"id": "Alpha Lyr", "name": "Vega", "type": "Star",
                                    "ra": "18:36:56.3", "dec": "+38:47:01", "mag": 0.03, "size": "-", "const": "Lyra"},
                                   {"id": "Alpha Boo", "name": "Arcturus", "type": "Star",
                                    "ra": "14:15:39.7", "dec": "+19:10:56", "mag": -0.05, "size": "-", "const": "Bootes"},
                                   {"id": "Alpha Tau", "name": "Aldebaran", "type": "Star",
                                    "ra": "04:35:55.2", "dec": "+16:30:33", "mag": 0.85, "size": "-", "const": "Taurus"},
                                   {"id": "Beta Ori", "name": "Rigel", "type": "Star",
                                    "ra": "05:14:32.3", "dec": "-08:12:06", "mag": 0.18, "size": "-", "const": "Orion"},
                                   {"id": "Alpha Aur", "name": "Capella", "type": "Star",
                                    "ra": "05:16:41.4", "dec": "+45:59:53", "mag": 0.08, "size": "-", "const": "Auriga"},
                                   {"id": "Beta Gem", "name": "Pollux", "type": "Star",
                                    "ra": "07:45:19.4", "dec": "+28:01:34", "mag": 1.14, "size": "-", "const": "Gemini"},
                                   {"id": "Alpha CMi", "name": "Procyon", "type": "Star",
                                    "ra": "07:39:18.1", "dec": "+05:13:30", "mag": 0.34, "size": "-", "const": "Canis Minor"},
                                   {"id": "Alpha Ori", "name": "Betelgeuse", "type": "Star",
                                    "ra": "05:55:10.3", "dec": "+07:24:25", "mag": 0.42, "size": "-", "const": "Orion"},
                                   {"id": "Alpha Eri", "name": "Achernar", "type": "Star",
                                    "ra": "01:37:42.8", "dec": "-57:14:12", "mag": 0.45, "size": "-", "const": "Eridanus"},
                                   {"id": "Beta Tau", "name": "Elnath", "type": "Star",
                                    "ra": "05:26:17.5", "dec": "+28:36:27", "mag": 1.65, "size": "-", "const": "Taurus"},
                                   {"id": "Alpha Vir", "name": "Spica", "type": "Star",
                                    "ra": "13:25:11.6", "dec": "-11:09:41", "mag": 0.98, "size": "-", "const": "Virgo"},
                                   {"id": "Alpha Sco", "name": "Antares", "type": "Star",
                                    "ra": "16:29:24.4", "dec": "-26:25:55", "mag": 0.96, "size": "-", "const": "Scorpius"},
                                   {"id": "Alpha Cyg", "name": "Deneb", "type": "Star",
                                    "ra": "20:41:25.9", "dec": "+45:16:49", "mag": 1.25, "size": "-", "const": "Cygnus"},
                                   {"id": "Alpha PsA", "name": "Fomalhaut", "type": "Star",
                                    "ra": "22:57:39.0", "dec": "-29:37:20", "mag": 1.16, "size": "-", "const": "Piscis Austrinus"}
                   ]
               }
    def search_catalog(self):
        catalog = self.catalog_combo.currentText()
        search_term = self.search_input.text().strip().lower()

        if catalog not in self.catalogs:
            return

        self.catalog_table.setRowCount(0)

        for obj in self.catalogs[catalog]:
            if (search_term in obj['id'].lower() or
                (obj['name'] and search_term in obj['name'].lower())):
                row = self.catalog_table.rowCount()
                self.catalog_table.insertRow(row)

                self.catalog_table.setItem(row, 0, QTableWidgetItem(obj['id']))
                self.catalog_table.setItem(row, 1, QTableWidgetItem(obj.get('name', '')))
                self.catalog_table.setItem(row, 2, QTableWidgetItem(obj['type']))
                self.catalog_table.setItem(row, 3, QTableWidgetItem(obj['ra']))
                self.catalog_table.setItem(row, 4, QTableWidgetItem(obj['dec']))
                self.catalog_table.setItem(row, 5, QTableWidgetItem(str(obj['mag'])))

    def show_object_details(self):
        selected = self.catalog_table.selectedItems()
        if not selected:
            return

        catalog = self.catalog_combo.currentText()
        obj_id = selected[0].text()

        obj = None
        for candidate in self.catalogs.get(catalog, []):
            if candidate['id'] == obj_id:
                obj = candidate
                break

        if not obj:
            return

        self.details_labels['id'].setText(f"ID: {obj['id']}")
        self.details_labels['name'].setText(f"Name: {obj.get('name', '-')}")
        self.details_labels['type'].setText(f"Type: {obj['type']}")
        self.details_labels['ra'].setText(f"RA: {obj['ra']}")
        self.details_labels['dec'].setText(f"DEC: {obj['dec']}")
        self.details_labels['mag'].setText(f"Magnitude: {obj['mag']}")
        self.details_labels['size'].setText(f"Size: {obj.get('size', '-')}")
        self.details_labels['const'].setText(f"Constellation: {obj.get('const', '-')}")

        try:
            coord = SkyCoord(obj['ra'], obj['dec'], unit=(u.hourangle, u.deg), frame='icrs')
            altaz = coord.transform_to(AltAz(obstime=Time.now(), location=self.location))
            self.details_labels['alt'].setText(f"Altitude: {altaz.alt.deg:.2f}°")
            self.details_labels['az'].setText(f"Azimuth: {altaz.az.deg:.2f}°")
        except Exception as e:
            print(f"Coordinate conversion error: {e}")
            self.details_labels['alt'].setText("Altitude: -")
            self.details_labels['az'].setText("Azimuth: -")

    def goto_selected_object(self):
        if not self.selected_object:
            return

        if not self.parent_window or not hasattr(self.parent_window, 'opcua_thread'):
            QMessageBox.warning(self, "Error", "OPC-UA connection not available")
            return

        try:
            coord = SkyCoord(self.selected_object['ra'], self.selected_object['dec'],
                            unit=(u.hourangle, u.deg))
            altaz = coord.transform_to(AltAz(obstime=self.current_time, location=self.parent_window.location))
            az_target = altaz.az.deg
            alt_target = altaz.alt.deg
            # Movement command sequence
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_Acc",
                1.0,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_velocity",
                30.0,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_Acc",
                1.0,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_velocity",
                3.0,
                ua.VariantType.Double
            )

            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_moveAbsolute",
                False,
                ua.VariantType.Boolean
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.EL_moveAbsolute",
                False,
                ua.VariantType.Boolean
            )

            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.Az_MoveAbsolutePosition",
                az_target,
                ua.VariantType.Double
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.El_MoveAbsolutePosition",
                alt_target,
                ua.VariantType.Double
            )

            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_moveAbsolute",
                True,
                ua.VariantType.Boolean
            )
            self.parent_window.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.EL_moveAbsolute",
                True,
                ua.VariantType.Boolean
            )

            QMessageBox.information(self, "GoTo Started",
                                  f"Slewing to {self.selected_object.get('name', self.selected_object.get('id', 'object'))}")

            # Update the mount tab with the new coordinates
            if hasattr(self.parent_window, 'mount_tab'):
                self.parent_window.mount_tab.goto_controls.ra_input.setText(self.selected_object['ra'])
                self.parent_window.mount_tab.goto_controls.dec_input.setText(self.selected_object['dec'])
                self.parent_window.mount_tab.goto_controls.az_display.setText(f"{az_target:.6f}")
                self.parent_window.mount_tab.goto_controls.alt_display.setText(f"{alt_target:.6f}")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to initiate GoTo: {str(e)}")


class TelescopeControlUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Telescope Control System")
        self.setGeometry(100, 100, 1280, 720)
        self.setup_style()

        self.location = EarthLocation(lat=32.7908*u.deg, lon=79.0002*u.deg, height=4507*u.m)
        self.opcua_thread = OPCUAClientThread("opc.tcp://192.168.10.50:4840")
        self.opcua_thread.data_updated.connect(self.update_opcua_data)
        self.opcua_thread.connection_status.connect(self.update_connection_status)
        self.opcua_thread.error_data.connect(self.update_error_plot)

        self.init_ui()
        self.update_time()
        self.tracking_thread = None

    def setup_style(self):
        self.setStyleSheet(f"""
            QMainWindow {{
                background-color: {COLOR_BACKGROUND};
                color: {COLOR_TEXT_LIGHT};
            }}
            QGroupBox {{
                border: 1px solid {COLOR_LIGHT};
                border-radius: 4px;
                margin-top: 10px;
                padding-top: 15px;
                background-color: {COLOR_DARK};
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 3px;
                color: {COLOR_TEXT_MEDIUM};
            }}
            QPushButton {{
                background-color: {COLOR_MEDIUM};
                color: {COLOR_TEXT_LIGHT};
                border: 1px solid {COLOR_LIGHT};
                padding: 5px;
                min-width: 70px;
            }}
            QPushButton:hover {{
                background-color: {COLOR_LIGHT};
            }}
            QPushButton:pressed {{
                background-color: {COLOR_DARK};
            }}
            QLineEdit {{
                background-color: {COLOR_DARK};
                color: {COLOR_TEXT_LIGHT};
                border: 1px solid {COLOR_LIGHT};
                padding: 5px;
                selection-background-color: {COLOR_MEDIUM};
            }}
            QTabWidget::pane {{
                border: 1px solid {COLOR_LIGHT};
                background: {COLOR_DARK};
            }}
            QTabBar::tab {{
                background: {COLOR_MEDIUM};
                color: {COLOR_TEXT_LIGHT};
                padding: 8px;
                border: 1px solid {COLOR_LIGHT};
                border-bottom: none;
                border-top-left-radius: 4px;
                border-top-right-radius: 4px;
            }}
            QTabBar::tab:selected {{
                background: {COLOR_LIGHT};
            }}
            QLabel {{
                color: {COLOR_TEXT_LIGHT};
            }}
            QCheckBox {{
                color: {COLOR_TEXT_LIGHT};
            }}
            QCheckBox::indicator {{
                width: 16px;
                height: 16px;
            }}
            QScrollArea {{
                border: none;
            }}
            QFrame {{
                background-color: {COLOR_DARK};
            }}
            QComboBox {{
                background-color: {COLOR_MEDIUM};
                color: {COLOR_TEXT_LIGHT};
                border: 1px solid {COLOR_LIGHT};
                padding: 3px;
            }}
            QProgressBar {{
                border: 1px solid {COLOR_LIGHT};
                border-radius: 3px;
                text-align: center;
                color: {COLOR_TEXT_LIGHT};
            }}
            QProgressBar::chunk {{
                background-color: {COLOR_ACCENT1};
                width: 10px;
            }}
            QSlider::groove:horizontal {{
                border: 1px solid {COLOR_LIGHT};
                height: 8px;
                background: {COLOR_DARK};
                margin: 2px 0;
            }}
            QSlider::handle:horizontal {{
                background: {COLOR_TEXT_LIGHT};
                border: 1px solid {COLOR_LIGHT};
                width: 18px;
                margin: -4px 0;
                border-radius: 3px;
            }}
            QTableWidget {{
                border: 1px solid {COLOR_LIGHT};
                gridline-color: {COLOR_LIGHT};
                color: {COLOR_TEXT_LIGHT};
            }}
            QHeaderView::section {{
                background-color: {COLOR_MEDIUM};
                color: {COLOR_TEXT_LIGHT};
                padding: 5px;
                border: none;
            }}
        """)

    def init_ui(self):
        main_widget = QWidget()
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)

        header_layout = QHBoxLayout()

        logo_label = QLabel()
        logo_pixmap = QPixmap("logo.jpeg") if os.path.exists("logo.jpeg") else QPixmap()
        if not logo_pixmap.isNull():
            logo_pixmap = logo_pixmap.scaled(40, 40, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        logo_label.setPixmap(logo_pixmap)
        logo_label.setFixedSize(45, 45)

        self.header_label = QLabel("Telescope Control System")
        self.header_label.setFont(QFont("Arial", 16, QFont.Bold))

        self.time_label = QLabel()
        self.time_label.setFont(QFont("Arial", 10))

        header_layout.addWidget(logo_label)
        header_layout.addWidget(self.header_label)
        header_layout.addStretch()
        header_layout.addWidget(self.time_label)

        main_layout.addLayout(header_layout)
        tabs = QTabWidget()
        tabs.setStyleSheet(f"""
            QTabBar::tab {{
                height: 30px;
                width: 120px;
                font-weight: bold;
            }}
        """)

        self.mount_tab = QWidget()
        mount_tab_layout = QHBoxLayout()
        mount_tab_layout.setContentsMargins(5, 5, 5, 5)
        mount_tab_layout.setSpacing(10)

        left_panel = QVBoxLayout()
        left_panel.setSpacing(10)

        data_group = QGroupBox("Mount Data")
        data_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        data_layout = QGridLayout()

        self.mount_data_labels = {
            'following_error': QLabel("0.0"),
            'altitude': QLabel("0.0000000"),
            'azimuth': QLabel("0.0000000"),
            'status': QLabel("Disconnected")
        }

        data_layout.addWidget(QLabel("Following Error:"), 0, 0)
        data_layout.addWidget(self.mount_data_labels['following_error'], 0, 1)
        data_layout.addWidget(QLabel("Altitude:"), 1, 0)
        data_layout.addWidget(self.mount_data_labels['altitude'], 1, 1)
        data_layout.addWidget(QLabel("Azimuth:"), 2, 0)
        data_layout.addWidget(self.mount_data_labels['azimuth'], 2, 1)
        data_layout.addWidget(QLabel("Status:"), 3, 0)
        data_layout.addWidget(self.mount_data_labels['status'], 3, 1)

        for i in range(4):
            data_layout.itemAtPosition(i, 1).widget().setAlignment(Qt.AlignRight)

        data_group.setLayout(data_layout)
        left_panel.addWidget(data_group)

        self.error_plot = ErrorPlotCanvas(self)
        self.error_plot.setMinimumHeight(300)
        left_panel.addWidget(self.error_plot)

        # GOTO controls
        self.goto_controls = GotoControls(self)
        left_panel.addWidget(self.goto_controls)

        mount_tab_layout.addLayout(left_panel, 2)

        right_panel = QVBoxLayout()
        right_panel.setSpacing(10)

        control_group = QGroupBox("Mount Controls")
        control_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        control_layout = QGridLayout()

        self.connect_btn = QPushButton("Connect")
        self.connect_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT3};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT1};
            }}
        """)
        self.connect_btn.clicked.connect(self.toggle_connection)
        self.home_btn = QPushButton("Home")
        self.home_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_ACCENT1};
                color: white;
            }}
            QPushButton:hover {{
                background-color: {COLOR_ACCENT3};
            }}
        """)
        self.home_btn.clicked.connect(self.home_mount)
        self.park_btn = QPushButton("Park")
        self.park_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: #FFA500;
                color: white;
            }}
            QPushButton:hover {{
                background-color: #d6a01e;
            }}
        """)
        self.park_btn.clicked.connect(self.park_mount)

        control_layout.addWidget(self.connect_btn, 0, 0)
        control_layout.addWidget(self.home_btn, 0, 1)
        control_layout.addWidget(self.park_btn, 1, 0, 1, 2)

        self.direction_btns = {}
        for i, (label, tooltip) in enumerate([
            ("←", "Move West"),
            ("→", "Move East"),
            ("↑", "Move North"),
            ("↓", "Move South")
        ]):
            btn = QPushButton(label)
            btn.setToolTip(tooltip)
            btn.setStyleSheet(f"""
                QPushButton {{
                    background-color: {COLOR_MEDIUM};
                    color: {COLOR_TEXT_LIGHT};
                    font-size: 16px;
                }}
                QPushButton:hover {{
                    background-color: {COLOR_LIGHT};
                }}
            """)
            btn.setFixedSize(60, 40)
            btn.clicked.connect(lambda _, d=label: self.move_direction(d))
            control_layout.addWidget(btn, 2 + i//2, i%2)
            self.direction_btns[label] = btn

        stop_btn_layout = QHBoxLayout()
        self.stop_az_btn = QPushButton("STOP AZ")
        self.stop_az_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_LIGHT};
                color: {COLOR_TEXT_LIGHT};
            }}
            QPushButton:hover {{
                background-color: {COLOR_MEDIUM};
            }}
        """)
        self.stop_az_btn.clicked.connect(self.stop_movement)
        self.stop_alt_btn = QPushButton("STOP ALT")
        self.stop_alt_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLOR_LIGHT};
                color: {COLOR_TEXT_LIGHT};
            }}
            QPushButton:hover {{
                background-color: {COLOR_MEDIUM};
            }}
        """)
        self.stop_alt_btn.clicked.connect(self.stop_movement)

        stop_btn_layout.addWidget(self.stop_az_btn)
        stop_btn_layout.addWidget(self.stop_alt_btn)
        control_layout.addLayout(stop_btn_layout, 6, 0, 1, 2)

        control_group.setLayout(control_layout)
        right_panel.addWidget(control_group)

        status_group = QGroupBox("Mount Status")
        status_group.setStyleSheet(f"QGroupBox {{ color: {COLOR_TEXT_MEDIUM}; }}")
        status_layout = QGridLayout()

        self.status_indicators = {
            'power': StatusIndicator("OFF"),
            'moving': StatusIndicator("IDLE"),
            'limits': StatusIndicator("OK"),
            'tracking': StatusIndicator("OFF")
        }

        for i, (label, widget) in enumerate(self.status_indicators.items()):
            status_layout.addWidget(QLabel(label.capitalize() + ":"), i//2, (i%2)*2)
            status_layout.addWidget(widget, i//2, (i%2)*2+1)

        status_group.setLayout(status_layout)
        right_panel.addWidget(status_group)

        mount_tab_layout.addLayout(right_panel, 1)
        self.mount_tab.setLayout(mount_tab_layout)
        tabs.addTab(self.mount_tab, "Mount")

        self.mass_dimm_tab = MassDimmWidget()
        tabs.addTab(self.mass_dimm_tab, "MASS-DIMM")

        tabs.addTab(FilterWheelTab(), "Filterwheel")
        tabs.addTab(CameraTab(), "Camera")
        tabs.addTab(CatalogTab(self), "Catalog")

        main_layout.addWidget(tabs)
        main_widget.setLayout(main_layout)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(main_widget)
        self.setCentralWidget(scroll_area)

        self.set_controls_enabled(False)

    def update_time(self):
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.time_label.setText(f"Local Time: {current_time}")
        QTimer.singleShot(1000, self.update_time)

    def update_opcua_data(self, data):
        if data['az_pos'] is not None:
            self.mount_data_labels['azimuth'].setText(f"{data['az_pos']:.7f}")
        if data['el_pos'] is not None:
            self.mount_data_labels['altitude'].setText(f"{data['el_pos']:.7f}")

        if data['az_done'] is not None:
            status = "idle" if data['az_done'] else "moving"
            self.status_indicators['moving'].set_status(status)

        if data['az_enable'] is not None:
            self.status_indicators['power'].set_status("on" if data['az_enable'] else "off")

        self.mass_dimm_tab.update_data(data)

    def update_error_plot(self, error_value):
        self.error_plot.update_plot(error_value)
        self.mount_data_labels['following_error'].setText(f"{error_value:.2f}")

    def update_connection_status(self, connected):
        if connected:
            self.mount_data_labels['status'].setText("Connected")
            self.set_controls_enabled(True)
            self.error_plot.axes.set_title("Following Error Plot", color=COLOR_TEXT_LIGHT)
            self.mass_dimm_tab.set_connected(True)
            QMessageBox.information(self, "Connected", "Telescope connection established")
        else:
            self.mount_data_labels['status'].setText("Disconnected")
            self.set_controls_enabled(False)
            self.error_plot.plot_empty()
            self.mass_dimm_tab.set_connected(False)
            QMessageBox.warning(self, "Disconnected", "Telescope connection lost")

    def set_controls_enabled(self, enabled):
        self.home_btn.setEnabled(enabled)
        self.park_btn.setEnabled(enabled)
        self.goto_controls.radec_go_btn.setEnabled(enabled)
        self.goto_controls.radec_stop_btn.setEnabled(enabled)
        self.goto_controls.track_btn.setEnabled(enabled)

        for btn in self.direction_btns.values():
            btn.setEnabled(enabled)

        self.stop_az_btn.setEnabled(enabled)
        self.stop_alt_btn.setEnabled(enabled)

    def toggle_connection(self):
        if self.opcua_thread.isRunning():
            self.opcua_thread.stop()
            self.connect_btn.setText("Connect")
        else:
            self.opcua_thread.start()
            self.connect_btn.setText("Disconnect")

    def home_mount(self):
        if not self.opcua_thread.isRunning():
            QMessageBox.warning(self, "Not Connected", "Please connect to the telescope first")
            return

        reply = QMessageBox.question(
            self, "Confirm Homing",
            "Are you sure you want to home the mount?",
            QMessageBox.Yes | QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            try:
                success_az = self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_starthome",
                    True,
                    ua.VariantType.Boolean
                )

                success_el = self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_starthome",
                    True,
                    ua.VariantType.Boolean
                )

                if success_az and success_el:
                    QMessageBox.information(self, "Homing Started", "Mount homing procedure initiated")
                    self.status_indicators['moving'].set_status("moving")
                else:
                    QMessageBox.warning(self, "Homing Failed", "Failed to send homing command")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to initiate homing: {str(e)}")

    def park_mount(self):
        if not self.opcua_thread.isRunning():
            QMessageBox.warning(self, "Not Connected", "Please connect to the telescope first")
            return

        reply = QMessageBox.question(
            self, "Confirm Parking",
            "Are you sure you want to park the mount?",
            QMessageBox.Yes | QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            try:
                # Movement command sequence
                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_halt2",
                    False,
                    ua.VariantType.Boolean
                )
                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_halt",
                    False,
                    ua.VariantType.Boolean
                )

                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_Acc",
                    1.0,
                    ua.VariantType.Double
                )
                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_velocity",
                    30.0,
                    ua.VariantType.Double
                )
                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_Acc",
                    1.0,
                    ua.VariantType.Double
                )
                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_velocity",
                    3.0,
                    ua.VariantType.Double
                )

                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_moveAbsolute",
                    False,
                    ua.VariantType.Boolean
                )
                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.EL_moveAbsolute",
                    False,
                    ua.VariantType.Boolean
                )

                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.Az_MoveAbsolutePosition",
                    0.0,
                    ua.VariantType.Double
                )
                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.El_MoveAbsolutePosition",
                    45.0,
                    ua.VariantType.Double
                )

                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_moveAbsolute",
                    True,
                    ua.VariantType.Boolean
                )
                self.opcua_thread.write_node(
                    "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.EL_moveAbsolute",
                    True,
                    ua.VariantType.Boolean
                )

                QMessageBox.information(self, "Parking Started", "Mount parking procedure initiated")
                self.status_indicators['moving'].set_status("moving")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to initiate parking: {str(e)}")

    def move_direction(self, direction):
        if not self.opcua_thread.isRunning():
            QMessageBox.warning(self, "Not Connected", "Please connect to the telescope first")
            return

        directions = {
            "←": ("Move West", "Az", -1),
            "→": ("Move East", "Az", 1),
            "↑": ("Move North", "El", 1),
            "↓": ("Move South", "El", -1)
        }

        if direction in directions:
            name, axis, sign = directions[direction]
            try:
                velocity = 1.0 * sign
                success = self.opcua_thread.write_node(
                    f"ns=4;s=|var|PAC320-CWE21-3A.Application.PLC_PRG.MC_MoveVelocity_{axis}.Velocity",
                    velocity,
                    ua.VariantType.Double
                )

                if success:
                    self.status_indicators['moving'].set_status("moving")
                else:
                    QMessageBox.warning(self, "Move Failed", f"Failed to send {name} command")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to initiate movement: {str(e)}")

    def stop_movement(self):
        try:
            self.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.Az_halt2",
                True,
                ua.VariantType.Boolean
            )
            self.opcua_thread.write_node(
                "ns=4;s=|var|PAC320-CWE21-3A.Application.global_variable.El_halt",
                True,
                ua.VariantType.Boolean
            )
            self.status_indicators['moving'].set_status("idle")
            QMessageBox.information(self, "Movement Stopped", "Telescope movement stopped")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to stop movement: {str(e)}")

    def closeEvent(self, event):
        if hasattr(self, 'tracking_thread') and self.tracking_thread and self.tracking_thread.isRunning():
            self.tracking_thread.stop()
            self.tracking_thread.wait()

        if hasattr(self, 'opcua_thread'):
            self.opcua_thread.stop()
        event.accept()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle('Fusion')

    font = QFont()
    font.setFamily("Arial")
    font.setPointSize(10)
    app.setFont(font)

    if not os.path.exists("assets"):
     os.makedirs("assets")
    window = TelescopeControlUI()
    window.showMaximized()
    sys.exit(app.exec_())
