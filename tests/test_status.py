import pytest
from app import app
from flask import jsonify, Flask
from datetime import datetime, timezone, timedelta
import json
from contextlib import contextmanager
from status.trajectories import exclude_trajectory_path_from_landmass
import numpy as np

@app.route('/api/deployments')
def deployments_mock():
    return jsonify(_generate_status_json())

@pytest.fixture
def client():
    app.config["TESTING"] = True
    app.config["DEBUG"] = False
    with app.test_client() as client:
        yield client

def test_index(client):
    resp = client.get("/")
    assert resp.status == "200 OK"

def test_deployments_api(client):
    with app.app_context():
        with deployments_mock():
            resp = client.get('/api/deployments')
            json_data = json.loads(resp.data)
            assert (json_data["results"][0]["name"] ==
                            "test-20200101T0000Z")

def _generate_status_json(data_type="api_return"):
    dt_now = datetime.now(tz=timezone.utc)
    dt_yesterday = dt_now - timedelta(1)
    status_dict = {
                    "results": [
                        {"name": "test-20200101T0000Z",
                        "datasetID": "test-20200101T0000Z",
                        "id": "5ddddddddddddddddddddddd",
                        "attribution": "None, test data only.",
                        "attribution": "Test Dataset",
                        "username": "test_user",
                        "archive_safe": True,
                        "checksum": "2385f68e7603c9d9aebb6c327353288a",
                        "completed": True,
                        "deployment_date": dt_now.timestamp(),
                        "deployment_dir": "test/test-20200101T0000Z",
                        "dac_url": "https://gliders.ioos.us/providers/deployment/5ddddddddddddddddddddddd",
                        "glider_name": "test",
                        "east": 180.0,
                        "west": -180.0,
                        "north": 90.0,
                        "south": -90.0,
                        "start": (dt_now - timedelta(1)).timestamp(),
                        "end": dt_now.timestamp(),
                        "delayed_mode": False,
                        "estimated_deploy_date": None,
                        "estimated_deploy_location": None,
                        "dac_url": "https://gliders.ioos.us/providers/deployment/5ddddddddddddddddddddddd",
                        "graph": "https://gliders.ioos.us/erddap/tabledap/test-20200101T0000Z.graph",
                        "fgdc": "https://gliders.ioos.us/erddap/metadata/fgdc/xml/test-20200101T0000Z_fgdc.xml",
                        "institution": "unit_test",
                        "iso": "https://gliders.ioos.us/erdap/tabledap/test-20200101T0000Z.iso19115",
                        "latest_nc_file": "test_20200105T0000Z.nc",
                        "latest_nc_file": "test_20200105T0000Z.nc",
                        "meta": "https://gliders.ioos.us/erddap/info/test-20200101T0000Z/index",
                        "nc_file_last_update": (dt_now - timedelta(1)).timestamp(),
                        "nc_files_count": 30,
                        "operator": "test",
                        "potential_invalid_files": [],
                        "rss": "https://gliders.ioos.us/erdap/rss/test.rss",
                        "status": None,
                        "subset": "https://gliders.ioos.us/erdap/tabledap/test-20200101T0000Z.subset",
                        "summary": "Testing only dataset",
                        "tds": "https://gliders.ioos.us/thredds/dodsC/catalog/deployments/test/test-20200101T0000Z.subset",
                        "ts0": dt_yesterday.isoformat(),
                        "ts1": dt_now.isoformat(),
                        "updated": dt_now.timestamp(),
                        "wmo_id": None
                        }
                    ]
                    }
    if data_type == "api_return":
        status_dict["num_results"] = len(status_dict["results"])
    elif data_type == "written_json":
        status_dict["meta"] = {"fetch_time":
                                dt_now.strftime("%b %d, %Y %H:%M Z")}


    return status_dict

def test_exclude_trajectory_path_from_landmass():
    """Test that glider track vertices over landmasses get properly excluded"""
    los_angeles_point = (-118.2437, 34.0522)  # Los Angeles, CA
    fake_traj = [
        (-160.0, 25.0),        # Somewhere in the central Pacific Ocean
        (-150.0, 28.0),        # A waypoint farther east
        (-140.0, 30.0),        # Approaching the eastern Pacific
        (-130.0, 32.0),        # Off the California coast
        los_angeles_point,
    ]
    traj_filtered = exclude_trajectory_path_from_landmass(fake_traj)
    assert (np.size(traj_filtered, 0) < len(fake_traj) and
            los_angeles_point not in traj_filtered.tolist())
