FROM node:13.0.1-alpine AS buildstep
LABEL maintainer "RPS <devops@rpsgroup.com>"

RUN mkdir -p /web
WORKDIR /web

COPY ./web/ /web/
RUN yarn global add grunt-cli && \
    yarn install && \
    grunt

FROM python:3.12
ARG glider_gid_uid=1000

RUN mkdir -p /glider-dac-status /mpl_config
VOLUME /mpl_config
RUN mkdir /glider-dac-status/logs
COPY app.py config.yml flask_environments.py manage.py /glider-dac-status/
COPY status /glider-dac-status/status
COPY navo /glider-dac-status/navo
COPY requirements/requirements.txt /requirements.txt

WORKDIR /glider-dac-status

RUN apt-get update && \
    apt-get install -y python3-netcdf4 libnetcdf-dev libhdf5-dev && \
    pip install --no-cache 'Cython<3.0' gunicorn && \
    pip install --no-cache -r /requirements.txt && \
    rm -rf /var/lib/apt/lists/* && \
    groupadd -g $glider_gid_uid glider && \
    useradd -u $glider_gid_uid -g $glider_gid_uid glider


ENV FLASK_ENV="PRODUCTION"
COPY --from=buildstep /web/ /glider-dac-status/web
RUN chown -R glider:glider /glider-dac-status/ /mpl_config

USER glider
EXPOSE 5000
