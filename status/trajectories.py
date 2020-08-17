#!/usr/bin/env python
# -*- coding: utf-8 -*-
import json
import requests
import os
import sys
from app import app
from shapely.geometry import LineString
from status.profile_plots import iter_deployments, is_recent_data, is_recent_update


def get_trajectory(erddap_url):
    '''
    Reads the trajectory information from ERDDAP and returns a GEOJSON like
    structure.
    '''
    # https://gliders.ioos.us/erddap/tabledap/ru01-20140104T1621.json?latitude,longitude&time&orderBy(%22time%22)
    url = erddap_url.replace('html', 'json')
    # ERDDAP requires the variable being sorted to be present in the variable
    # list.  The time variable will be removed before converting to GeoJSON
    url += '?longitude,latitude,time&orderBy(%22time%22)'
    response = requests.get(url, timeout=180)
    if response.status_code != 200:
        raise IOError("Failed to fetch trajectories: {}".format(erddap_url))
    data = response.json()
    geo_data = {
        'type': 'LineString',
        'coordinates': [c[0:2] for c in data['table']['rows']]
    }

    geometry = parse_geometry(geo_data)
    coords = LineString(geometry['coordinates'])
    trajectory = coords.simplify(0.02, preserve_topology=False)
    geometry = {
        'type': 'LineString',
        'coordinates': list(trajectory.coords),
        'properties': {
            'oceansmap_type': 'glider'
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


def parse_geometry(geometry):
    '''
    Filters out potentially bad coordinate pairs as returned from
    GliderDAC. Returns a safe geometry object.

    :param dict geometry: A GeoJSON Geometry object
    '''
    coords = []
    for lon, lat in geometry['coordinates']:
        if lon is None or lat is None:
            continue
        coords.append([lon, lat])
    return {'coordinates': coords}


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
    for deployment in iter_deployments():
        try:
            # Only add if the deployment has been recently updated or the data is recent
            recent_update = is_recent_update(deployment['updated'])
            recent_data = is_recent_data(deployment)
            existing_trajectory = trajectory_exists(deployment)
            if (recent_update or recent_data or not existing_trajectory):
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
