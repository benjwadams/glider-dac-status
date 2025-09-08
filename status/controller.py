#!/usr/bin/env python
'''
status.controller

The controller definition for the status application. This mostly contains the
routes and logic API
'''

import requests
from flask import jsonify, current_app
from status import api
from status.trajectories import get_trajectory


@api.route('/test')
def test():
    return jsonify(message="Running")


# --------------------------------------------------------------------------------
# Proxies - Use with caution
# --------------------------------------------------------------------------------
@api.route('/deployment', methods=['GET'])
def get_deployments():
    url = current_app.config.get('DAC_API')
    response = requests.get(url)
    return response.content, response.status_code, dict(response.headers)


@api.route('/deployment/<string:username>/<string:deployment_name>', methods=['GET'])
def get_deployment(username, deployment_name):
    url = current_app.config.get('DAC_API')
    url += '/%s/%s' % (username, deployment_name)
    response = requests.get(url)
    return response.content, response.status_code, dict(response.headers)


@api.route('/track/<string:username>/<string:deployment_name>')
def track(username, deployment_name):
    url = current_app.config.get('DAC_API')
    url += '/%s/%s' % (username, deployment_name)
    response = requests.get(url)
    if response.status_code != 200:
        return jsonify(error="Unable to read from DAC API"), 500
    deployment = response.json()
    erddap_url = deployment['erddap']
    geo_data = get_trajectory(erddap_url)
    return jsonify(**geo_data)
