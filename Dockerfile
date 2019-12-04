FROM python:3.8.0-buster

LABEL description="Glider DAC Status"

RUN mkdir /usr/local/glider-dac-status
COPY . /usr/local/glider-dac-status

WORKDIR /usr/local/glider-dac-status

RUN apt-get update && \
    apt-get install -y python3-netcdf4 && \
    pip install --no-cache cython gunicorn && \
    pip install --no-cache -r requirements/requirements.txt && \
    useradd glider

USER glider

# volume will need to be mounted for config
ENV SETTINGS_FILE_FOR_DYNACONF="/usr/local/glider-dac-status/config/config.production.yml" \
    ENV_FOR_DYNACONF=production \
    FLASK_ENV=production

EXPOSE 5000

# TODO: Find a way to use gunicorn while preserving config settings such as
# port, etc
#CMD [ "gunicorn", "app:app" ]
CMD [ "python", "app.py" ]
