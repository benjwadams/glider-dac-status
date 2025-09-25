#!/usr/bin/env python
# -*- coding: utf-8 -*-
import json
import requests
import sys
from app import app
from shapely.geometry import Point, LineString
from shapely.strtree import STRtree
import shapely.geometry as sgeom
from status.profile_plots import iter_deployments, is_recent_data, is_recent_update
from requests.exceptions import RequestException
import numpy as np
from datetime import datetime

import os
os.environ["CARTOPY_USER_BACKGROUNDS"] = "/tmp/cartopy"
os.environ["CARTOPY_DATA_DIR"] = "/tmp/cartopy"

import cartopy.io.shapereader as shpreader


# Load higher-resolution land polygons for better accuracy
# Load into spatial index for performance
shp_geoms = shpreader.Reader(
               shpreader.natural_earth(resolution='50m', category='physical',
               name='land')).geometries()
land_tree = STRtree([geom for geom in shp_geoms])

def get_trajectory(erddap_url):
    '''
    Reads the trajectory information from ERDDAP and returns a GEOJSON-like
    structure. Filters by min_time from deployment date.
    '''
    # Example URL:
    # https://gliders.ioos.us/erddap/tabledap/ru01-20140104T1621.json?latitude,longitude&time&orderBy(%22time%22)

    # get deployment time (e.g., 20250611T0000)
    min_time = erddap_url.split("/")[-1].replace(".html", "").split("-")[-1]

    # fix url with json extension
    url = erddap_url.replace("html", "json")

    # ERDDAP requires the variable being sorted to be present in the variable
    # list. The time variable will be removed before converting to GeoJSON

    valid_response = False
    for qc_append in ("qartod_location_test_flag,", ""):
        url_append = url + f"?longitude,latitude,{qc_append}time&orderBy(%22time%22)"
        try:
            response = requests.get(url_append, timeout=180, allow_redirects=True)
            response.raise_for_status()
        except RequestException as e:
            print(e)
            continue
        else:
            valid_response = True
            break

    if not valid_response:
        app.logger.error(f"Failed to fetch trajectory: {url_append}")

    data = response.json()

    # Map rows into lon/lat/time/flag
    col_names = data["table"]["columnNames"]
    rows = data["table"]["rows"]

    # Identify column indices dynamically
    lon_idx = col_names.index("longitude")
    lat_idx = col_names.index("latitude")
    time_idx = col_names.index("time")
    flag_idx = col_names.index("qartod_location_test_flag") if "qartod_location_test_flag" in col_names else None

    geo_data = {
        "type": "LineString",
        "coordinates": [(r[lon_idx], r[lat_idx]) for r in rows],
        "time": [r[time_idx] for r in rows],
        "flag": [r[flag_idx] for r in rows] if flag_idx is not None else None,
    }

    # Call your parse function with min_time filter
    geometry = parse_geometry_with_checks(geo_data, geo_data["flag"] is not None, min_time)

    # Simplify trajectory
    coords = LineString(geometry["coordinates"])
    trajectory = coords.simplify(0.02, preserve_topology=False)

    geometry = {
        "type": "LineString",
        "coordinates": list(trajectory.coords),
        "properties": {
            "oceansmap_type": "glider",
        }
    }
    return geometry


def get_path(deployment):
    '''
    Returns the path to the trajectory file

    :param dict deployment: Dictionary containing the deployment metadata
    '''
    trajectory_dir = app.config.get('TRAJECTORY_DIR')
    username = deployment['username']
    name = deployment['name']
    dir_path = os.path.join(trajectory_dir, username)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    file_path = os.path.join(dir_path, name + '.json')
    return file_path


def write_trajectory(deployment, geo_data):
    '''
    Writes a geojson like python structure to the appropriate data file

    :param dict deployment: Dictionary containing the deployment metadata
    :param dict geometry: A GeoJSON Geometry object
    '''
    file_path = get_path(deployment)
    with open(file_path, 'w') as f:
        f.write(json.dumps(geo_data))


def parse_geometry_with_checks(geometry: dict, has_flag: bool, min_time: str = None):
    """
    Filters out bad coordinate pairs based on:
      - minimum time threshold (if provided),
      - flags,
      - land masking,
      - outlier detection.

    Returns geometry with only 'coordinates'.
    """

    coords = []
    times = geometry.get("time")

    # --- Step 0: Time filtering ---
    if min_time and times:
        min_dt = datetime.strptime(min_time, "%Y%m%dT%H%M")
        filtered = [((lon, lat), t) for (lon, lat), t in zip(geometry['coordinates'], times)
                    if datetime.strptime(t, "%Y-%m-%dT%H:%M:%SZ") >= min_dt]
        filtered_coords = [lonlat for lonlat, _ in filtered]
    else:
        filtered_coords = geometry['coordinates']

    # --- Step 1: Filter by flags and missing values ---
    if has_flag:
        filtered_coords = [
            (lon, lat) for (lon, lat), flag in zip(filtered_coords, geometry['flag'])
            if (flag is None or flag == 1) and lon is not None and lat is not None
        ]
    else:
        filtered_coords = [(lon, lat) for lon, lat in filtered_coords
                           if lon is not None and lat is not None]

    # --- Step 2: Remove points that fall on land ---
    sea_coords = exclude_trajectory_path_from_landmass(filtered_coords)
    if sea_coords.size < len(filtered_coords):
        logger.warning(f"Detected and removed possible land overlap for deployment")

    # --- Step 3: Remove statistical outliers ---
    if sea_coords.size > 0:
        lons, lats = sea_coords.T

        z_lon = np.abs((lons - lons.mean()) / lons.std())
        z_lat = np.abs((lats - lats.mean()) / lats.std())

        threshold = 3
        outlier_filtered_coords = sea_coords[(z_lon <= threshold) &
                                             (z_lat <= threshold)]
    # empty coordinates should be rare, but will avoid throwing an exception
    else:
        outlier_filtered_coords = sea_coords
    return {'coordinates': outlier_filtered_coords.tolist()}


def exclude_trajectory_path_from_landmass(
    trajectory_coords: List[Tuple[float, float]]
) -> np.ndarray:
    """
    Remove trajectory points that fall on land.

    Parameters
    ----------
    trajectory_coords : List[Tuple[float, float]]
        Sequence of (longitude, latitude) pairs.

    Returns
    -------
    np.ndarray
        Array of (longitude, latitude) pairs that lie over water.
    """
    points = [Point(lon, lat) for lon, lat in trajectory_coords]
    land_vertex_indices, _ = land_tree.query(points, predicate="within")
    coords_array = np.array(trajectory_coords, dtype=float)
    return np.delete(coords_array, land_vertex_indices, axis=0)


def trajectory_exists(deployment):
    '''
    Returns True if the data is within the last week

    :param dict deployment: Dictionary containing the deployment metadata
    '''

    file_path = get_path(deployment)
    return os.path.exists(file_path)


def generate_trajectories(deployments=None):
    '''
    Determine which trajectories need to be built, and write geojson to file
    '''
    # TODO: Use a less brute force approach to filtering
    for deployment in iter_deployments():
        if deployments is not None and deployment["name"] not in deployments:
            continue
        try:
            # Only add if the deployment has been recently updated or the data is recent
            recent_update = is_recent_update(deployment['updated'])
            recent_data = is_recent_data(deployment)
            existing_trajectory = trajectory_exists(deployment)
            if (not deployment["name"].endswith("-delayed") and
                (recent_update or recent_data or not existing_trajectory
                or not deployment["completed"])):
                app.logger.info(f"Fetching trajectory for {deployment['name']}")
                geo_data = get_trajectory(deployment['erddap'])
                write_trajectory(deployment, geo_data)
        except Exception:
            from traceback import print_exc
            print_exc()
    return 0


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description=generate_trajectories.__doc__)
    parser.add_argument(
        '-d', '--deployment',
        action='append',
        help='Which deployment to build'
    )
    args = parser.parse_args()
    sys.exit(generate_trajectories(args.deployment))
