FROM continuumio/miniconda3

ENV ARGS=''
WORKDIR /app

# TODO: install non-pip dependencies here

ARG CACHEBUST
COPY environment.yml .
RUN conda env update

COPY . .

CMD python run.py ${ARGS}
