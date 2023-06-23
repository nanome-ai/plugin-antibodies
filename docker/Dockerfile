FROM continuumio/miniconda3

ENV ARGS=''
ENV SETTINGS_DIR=/data/user-settings
WORKDIR /app

ARG CACHEBUST
COPY environment.yml .
RUN conda env update

RUN mkdir -p ${SETTINGS_DIR}
COPY user-settings/default.json ${SETTINGS_DIR}

COPY . .

CMD python run.py ${ARGS}
